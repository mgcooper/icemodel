function varargout = postprocess(ice1, ice2, opts, varargin)
   %POSTPROCESS Post-process simulation output.
   %
   % Syntax:
   %
   % [ice1, ice2] = icemodel.postprocess(ice1, ice2, opts, simyears)
   % [ice1, ice2] = icemodel.postprocess(ice1, ice2, opts, swd, lwd, ...
   %    albedo, time)
   %
   % Description:
   %
   % [ice1, ice2] = icemodel.postprocess(ice1, ice2, opts, simyears)
   % Loads the met data resolved in opts. SIMYEARS may be a scalar year or a
   % year vector; a scalar year subsets the met data to that year, while a
   % vector returns the full concatenated met series for those years.
   %
   % [ice1, ice2] = icemodel.postprocess(ice1, ice2, opts, swd, lwd, ...
   %    albedo, time)
   %
   % See also:
   %
   %#codegen

   % Process inputs
   if nargin == 7
      [swd, lwd, albedo, time] = deal(varargin{:});
   elseif nargin == 4
      simyears = varargin{1};
      met = icemodel.loadmet(opts);
      if isscalar(simyears)
         ii = year(met.Time) == simyears;
         if countTimeSteps(ice1) == height(met)
            [ice1, ice2] = subsetOutput(ice1, ice2, ii);
         end
         met = met(ii, :);
      end
      swd = met.swd;
      lwd = met.lwd;
      albedo = met.albedo;
      time = met.Time;
   else
      error('unrecognized number of inputs')
   end

   % Load physical constants
   Tf = icemodel.physicalConstant('Tf');

   % Calculate runoff
   if strcmp('skinmodel', opts.smbmodel)
      ice1 = icemodel.surface.diagnose_surface_runoff(ice1, opts.dt);
   elseif strcmp('icemodel', opts.smbmodel)
      ice1 = icemodel.column.diagnose_column_runoff(ice1, ice2, opts);
   end

   % Compute a full state and energy balance
   if ~strcmp(opts.output_profile, 'minimal')
      [ice1, ice2] = computeState(ice1, ice2, opts, swd, lwd, albedo, Tf);
   end

   % Convert surface and subsurface ice temperature from Kelvin to Celsius.
   ice1.Tsfc = min(ice1.Tsfc - Tf, 0);
   ice2.Tice = min(ice2.Tice - Tf, 0);

   % Convert logical flags to single
   if isfield(ice1, 'Tice_converged')
      ice1.Tice_converged = single(ice1.Tice_converged);
   end
   if isfield(ice1, 'Tsfc_converged')
      ice1.Tsfc_converged = single(ice1.Tsfc_converged);
   end

   % Convert ice1 to timetable
   time.TimeZone = 'UTC';
   ice1 = struct2table(ice1);
   ice1 = table2timetable(ice1, 'RowTimes', time);

   % Retime 15-minute data to hourly values with per-channel aggregation.
   % Native timetable bins preserve variable classes and partial windows.
   if opts.dt == 900
      [ice1, bin_start, bin_end] = ...
         icemodel.retimeHourlyFixedStep(ice1);
      [ice1, ice2] = retimeLogical(ice1, ice2, bin_start, bin_end);
      ice2 = retimeIce2(ice2, bin_start, bin_end);
   end

   % Round the data to save disk space, retaining necessary precision
   [ice1, ice2] = roundData(ice1, ice2);

   if ~strcmp(opts.output_profile, 'minimal')
      ice2.Time = ice1.Time; % not added in legacy grid saves, maybe remove.

      % Rename ice1 vars to match the naming conventions i use everywhere else
      oldvars = {'Qsi','Qsr','Qsn','Qli','Qle','Qln','Qh','Qe','Qc','Qn','Tsfc'};
      newvars = {'swd','swu','swn','lwd','lwu','lwn','shf','lhf','chf','netr','tsfc'};
      ice1 = renameTableVars(ice1, ...
         oldvars(ismember(oldvars, ice1.Properties.VariableNames)), ...
         newvars(ismember(oldvars, ice1.Properties.VariableNames)));
   end

   switch nargout
      case 2
         varargout{1} = ice1;
         varargout{2} = ice2;
      case 3
         varargout{1} = ice1;
         varargout{2} = ice2;

         if nargin == 7
            met = icemodel.loadmet(opts);
            met = met(isbetween(met.Time, time(1), time(end)), :);
         end
         varargout{3} = icemodel.processmet(met, newTimeStep="hourly");
   end
end

%%
function n = countTimeSteps(ice1)

   fields = fieldnames(ice1);
   n = size(ice1.(fields{1}), 1);
end

%%
function [ice1, ice2] = subsetOutput(ice1, ice2, ii)

   fields = fieldnames(ice1);
   for n = 1:numel(fields)
      ice1.(fields{n}) = ice1.(fields{n})(ii, :);
   end

   fields = fieldnames(ice2);
   for n = 1:numel(fields)
      thisfield = fields{n};
      if strcmp(thisfield, 'Z')
         continue
      elseif strcmp(thisfield, 'Time')
         ice2.Time = ice2.Time(ii, :);
      else
         ice2.(thisfield) = ice2.(thisfield)(:, ii);
      end
   end
end

%%
function ice2 = retimeIce2(ice2, bin_start, bin_end)
   %RETIMEICE2 Aggregate subsurface outputs over the surface hourly bins.

   % Get the field names of ice2
   fields = fieldnames(ice2);

   % Allocate one subsurface column per surface bin so every channel shares
   % the retained output boundaries, including partial and empty native bins.
   n_bins = numel(bin_start);
   tmp = struct();
   for n = 1:numel(fields)
      tmp.(fields{n}) = nan(size(ice2.Tice, 1), n_bins);
   end
   % Replace Z, if this is not a legacy grid run
   if isfield(ice2, 'Z')
      tmp.Z = ice2.Z;
   end

   % Decide each field's rule once. Every df_ channel is a per-step increment,
   % so it sums; errH is a residual and also sums; everything else averages.
   % Z is the depth grid and is copied, not aggregated.
   do_sum = icemodel.isIncrementChannel(fields) | strcmp(fields, 'errH');
   skip = strcmp(fields, 'Z');

   % Aggregate the exact raw samples assigned to each native hourly bin.
   for n = 1:n_bins
      if bin_start(n) == 0
         continue
      end
      ii = bin_start(n):bin_end(n);

      for m = 1:numel(fields)
         if skip(m)
            continue
         end
         thisfield = fields{m};
         if do_sum(m)
            tmp.(thisfield)(:, n) = sum(ice2.(thisfield)(:, ii), 2);
         else
            tmp.(thisfield)(:, n) = mean(ice2.(thisfield)(:, ii), 2);
         end
      end
   end
   ice2 = tmp;
end

%%
function [ice1, ice2] = retimeLogical( ...
      ice1, ice2, bin_start, bin_end)
   %RETIMELOGICAL Aggregate logical flags over the surface hourly bins.

   % Collapse each two-dimensional flag with the same native bin membership
   % used for numeric subsurface outputs, then retain its surface-layer flag.
   fields = fieldnames(ice2);
   for n = 1:numel(fields)
      thisfield = fields{n};
      if islogical(ice2.(thisfield)(1, 1))
         flag = false(size(ice2.(thisfield), 1), numel(bin_start));
         for mm = 1:numel(bin_start)
            if bin_start(mm) == 0
               continue
            end
            ii = bin_start(mm):bin_end(mm);
            flag(:, mm) = any(ice2.(thisfield)(:, ii), 2);
         end
         ice1.(thisfield) = transpose(flag(1, :));
         ice2 = rmfield(ice2, thisfield);
      end
   end

   % 1-d logical

end

%%
function [ice1, ice2] = roundData(ice1, ice2)

   % Round legacy ice1 channels to five digits. Preserve diagnostic mass-budget
   % ledgers at solver precision so signed closure identities remain testable.
   % Round one variable at a time. Brace EXTRACTION concatenates the selected
   % columns into one array first, and the single-precision columns
   % (Tsfc_converged, Tice_converged) promote the whole block to single, so
   % every double channel would round at single precision. Brace ASSIGNMENT
   % restores each variable's original class, which makes that invisible in
   % the stored types.
   vars1 = ice1.Properties.VariableNames;
   is_budget = ismember(vars1, icemodel.namelists.budgetoutputs());

   % df_rof is a per-step increment like ice2's df_liq/df_evp/df_lyr, so it
   % keeps their 8-digit precision. Rounding it to five digits would zero the
   % small overflow values the closure identities are checked against.
   is_increment = icemodel.isIncrementChannel(vars1);
   keep_precision = is_budget | is_increment;
   round_names = vars1(~keep_precision);
   increment_names = vars1(is_increment & ~is_budget);

   % The isnumeric guards are needed only by this per-variable form, because
   % round rejects a logical column. retimeLogical can
   % move a logical ice2 flag channel into ice1, though no shipped vars2 list
   % currently names one.
   for k = 1:numel(round_names)
      if isnumeric(ice1.(round_names{k}))
         ice1.(round_names{k}) = round(ice1.(round_names{k}), 5);
      end
   end
   for k = 1:numel(increment_names)
      if isnumeric(ice1.(increment_names{k}))
         ice1.(increment_names{k}) = round(ice1.(increment_names{k}), 8);
      end
   end

   % Round the ice2 data
   fields = fieldnames(ice2);
   for mm = 1:numel(fields)
      thisfield = fields{mm};
      switch thisfield
         case {'f_ice','f_liq','k_vap','k_eff'}
            ice2.(thisfield) = round(ice2.(thisfield), 5);
         case {'Tice', 'h_melt','h_freeze'}
            ice2.(thisfield) = round(ice2.(thisfield), 3);
         case {'cp_sno','ro_sno'}
            ice2.(thisfield) = round(ice2.(thisfield), 1);
         case {'df_liq','df_lyr','df_evp','Qsub','Sc','errT','errH'}
            ice2.(thisfield) = round(ice2.(thisfield), 8);
      end
   end

   % For reference, another way to do it
   % persistent lookup
   % if isempty(lookup)
   %    lookup = {
   %       'f_ice', 5; 'f_liq', 5; 'k_vap', 5; 'k_eff', 5;
   %       'Tice', 3; 'h_melt', 3; 'h_freeze', 3;
   %       'cp_sno', 1; 'ro_sno', 1;
   %       'df_liq', 8; 'df_lyr', 8; 'Qsub', 8; 'Sc', 8; 'errT', 8; 'errH', 8
   %       };
   % end
   % [fields, idx] = intersect(lookup(:, 1), fieldnames(ice2));
   % for n = 1:numel(fields)
   %    ice2.(fields{n}) = round(ice2.(fields{n}), lookup{idx(n), 2});
   % end
end

%%
function [ice1, ice2] = computeState(ice1, ice2, opts, swd, lwd, albedo, Tf)

   % Compute bulk density, heat capacity, thermal conductivity, and a full
   % surface and subsurface energy balance. Don't do this for large simulations
   % if time or disk space is limited, instead compute them after the
   % simulation.

   % Compute bulk density (kg/m3), heat capacity (J/kg/K), thermal K (W/m/K)
   T_ice = ice2.Tice;
   f_liq = ice2.f_liq;
   f_ice = ice2.f_ice;

   % Phase-aware effective conductivity for diagnostics
   [k_eff, k_vap] = icemodel.column.bulk_thermal_conductivity(...
      T_ice, f_ice, f_liq);
   ro_sno = icemodel.column.bulk_density(f_ice, f_liq);
   cp_sno = icemodel.column.bulk_specific_heat_capacity( ...
      f_ice, f_liq, ro_sno);

   % Compute a mesh for plotting
   Z = opts.z0_thermal;
   dz = opts.dz_thermal;

   % Assign values to ice2
   ice2.Z      = (dz/2:dz:Z-dz/2)';
   ice2.k_eff  = k_eff;                  % eff. thermal k
   ice2.k_vap  = k_vap;                  % water vapor diffusivity
   ice2.cp_sno = cp_sno;                 % sp. heat cap.
   ice2.ro_sno = ro_sno;                 % ice density

   % Compute the radiative heat fluxes
   lwu = -icemodel.surface.outgoing_longwave_radiation(min(ice1.Tsfc, Tf));

   ice1.albedo = albedo;                        % albedo
   ice1.swd    = swd;                           % shortwave down
   ice1.swu    = swd.*albedo;                   % shortwave up
   ice1.lwd    = lwd;                           % longwave down
   ice1.lwu    = lwu;                           % longwave up
   ice1.swn    = (1-albedo).*swd;               % net shortwave
   ice1.lwn    = ice1.lwd-ice1.lwu;             % net longwave
   ice1.netr   = ice1.swn+ice1.lwn;             % net radiation
   ice1.Qsip   = (1-ice1.chi).*ice1.swd;        % penetrated shortwave
   ice1.Qabs   = (1-albedo).*swd;               % absorbed shortwave
   ice1.Qsrf   = ice1.chi.*ice1.Qabs;           % skin shortwave
   ice1.Qsub   = (1-ice1.chi).*ice1.Qabs;       % subsurf shortwave
   ice1.Qbal   = ice1.Qsub+ice1.Qsrf-ice1.Qabs; % balance
end

% for reference, Qsub can also be computed this way
% ice1.Qsub = sum(ice2.Sc.*dz)';             % subsurf shortwave


% The surface energy balance: Qm = chi*Qsi*(1-albedo) + Qln + Qh + Qe + Qc
% The subsurface shortwave balance: Qsub = (1-chi)*Qsi*(1-albedo)
%
%                      Qsi  Qsr
%                       |    ^
%  skin (wall)          v    |     Qabs = Qsi*(1-albedo) = Qsrf + Qsub
% ----------------------|---/---   Qsrf = chi*Qsi*(1-albedo) = chi*Qabs
%  surface layer        v  /    \
% -------------------- Qsip----- > Qsub = Qsi*(1-chi)*(1-albedo)
%  subsurface layers    |       /       = Qsip*(1-albedo)
%                       v      /        = Qabs*(1-chi)
% ------------------------------
%
% The total absorbed shortwave radiation equals the sum of the 'skin'
% absorbed radiation (Qsrf) and the subsurface absorbed radiation (Qsub)
% The subsurface absorbed radiation is the integral from z=0 to z=infty of
% dQ/dz*dz i.e. sum(Sc.*dz) with Sc the source term in units W/m3. chi is
% just the ratio of Qsrf to Qabs i.e. how much of the total absorbed energy
% is absorbed by the wall. Qsip is the penetrating radiation, some of which
% is absorbed and some of which is reflected (contributes to Qsr)
% in short, if we call the net solar downflux in the ice Q, then:
% Qsub = dQ is the absorbed solar downflux in the ice in each layer
% Sc = dQ/dz is the source term (divergence of the net solar downflux)


% The original method when I saved the grid output
%
% % Retime to hourly
% if opts.dt == 900
%
%    ice1  = retime(ice1,'hourly','mean');
%    feb29 = month(ice1.Time) == 2 & day(ice1.Time) == 29;
%    ice1  = ice1(~feb29,:);
%
%    % init tmp arrays to retime ice2
%    tmp.Tice = nan(size(ice2.Tice,1),numel(ice1.Time));
%
%    if strcmp('icemodel', opts.smbmodel)
%       tmp.df_liq = nan(size(ice2.df_liq,1),numel(ice1.Time));
%       tmp.f_ice = nan(size(ice2.f_ice,1),numel(ice1.Time));
%       tmp.f_liq = nan(size(ice2.f_liq,1),numel(ice1.Time));
%    end
%
%    for n = 1:numel(ice1.Time)
%       % this works b/c we know it's fifteen minute data
%       i1 = n*4-3;
%       i2 = n*4;
%
%       tmp.Tice(:,n) = mean(ice2.Tice(:,i1:i2),2);
%
%       if strcmp('icemodel', opts.smbmodel)
%          tmp.df_liq(:,n) = sum(ice2.df_liq(:,i1:i2),2);
%          tmp.f_liq(:,n) = mean(ice2.f_liq(:,i1:i2),2);
%          tmp.f_ice(:,n) = mean(ice2.f_ice(:,i1:i2),2);
%       end
%    end
%    ice2 = tmp;
% end
