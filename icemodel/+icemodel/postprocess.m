function varargout = postprocess(ice1, ice2, opts, varargin)
   %POSTPROCESS Calculate diagnostic outputs and format simulation output.
   %
   % [ice1, ice2] = icemodel.postprocess(ice1, ice2, opts, simyears)
   % [ice1, ice2] = icemodel.postprocess(ice1, ice2, opts, swd, lwd, albedo, time)
   %
   % [ice1, ice2] = icemodel.postprocess(ice1, ice2, opts, simyears) loads
   % shortwave radiation, longwave radiation, albedo, and time from the
   % meteorological data selected by opts. simyears is a scalar year or a year
   % vector. A scalar year also subsets matching output when it has the full
   % meteorological time axis.
   %
   % [ice1, ice2] = icemodel.postprocess(ice1, ice2, opts, swd, lwd, ...
   %    albedo, time) uses the supplied shortwave radiation, longwave radiation,
   % albedo, and time instead of loading meteorological data. ice1 contains
   % surface output and ice2 contains subsurface output. opts controls the model
   % type, output profile, timestep, and output conventions.
   %
   % See also: icemodel.loadmet, icemodel.retimeHourlyFixedStep,
   % icemodel.processmet
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

   % Load physical constants.
   Tf = icemodel.physicalConstant('Tf');

   % Calculate runoff.
   if strcmp('skinmodel', opts.smbmodel)
      ice1 = icemodel.surface.diagnose_surface_runoff(ice1, opts.dt);
   elseif strcmp('icemodel', opts.smbmodel)
      ice1 = icemodel.column.diagnose_column_runoff(ice1, ice2, opts);
   end

   % Compute a full state and energy balance.
   if ~strcmp(opts.output_profile, 'minimal')
      [ice1, ice2] = computeState(ice1, ice2, opts, swd, lwd, albedo, Tf);
   end

   % Convert surface and subsurface ice temperature from Kelvin to Celsius.
   ice1.Tsfc = min(ice1.Tsfc - Tf, 0);
   ice2.Tice = min(ice2.Tice - Tf, 0);

   % Convert logical flags to single.
   if isfield(ice1, 'Tice_converged')
      ice1.Tice_converged = single(ice1.Tice_converged);
   end
   if isfield(ice1, 'Tsfc_converged')
      ice1.Tsfc_converged = single(ice1.Tsfc_converged);
   end

   % Convert ice1 to timetable.
   time.TimeZone = 'UTC';
   ice1 = struct2table(ice1);
   ice1 = table2timetable(ice1, 'RowTimes', time);

   % Retime 15-minute data into the hourly bins returned by
   % retimeHourlyFixedStep. The bins retain partial and empty windows and each
   % field keeps its aggregation rule and numeric class.
   if opts.dt == 900
      [ice1, bin_start, bin_end] = ...
         icemodel.retimeHourlyFixedStep(ice1);
      [ice1, ice2] = retimeLogical(ice1, ice2, bin_start, bin_end);
      ice2 = retimeIce2(ice2, bin_start, bin_end);
   end

   % Round the data to save disk space, retaining necessary precision.
   [ice1, ice2] = roundData(ice1, ice2);

   if ~strcmp(opts.output_profile, 'minimal')
      ice2.Time = ice1.Time; % not added in legacy grid saves, maybe remove.

      % Rename ice1 vars to match the naming conventions I use everywhere else.
      oldvars = ...
         {'Qsi','Qsr','Qsn','Qli','Qle','Qln','Qh','Qe','Qc','Qn','Tsfc'};
      newvars = ...
         {'swd','swu','swn','lwd','lwu','lwn','shf','lhf','chf','netr','tsfc'};
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

   % Get the field names of ice2.
   fields = fieldnames(ice2);

   % Allocate one subsurface column for each retained surface bin. This gives
   % every field the same hourly, partial, and empty output windows.
   n_bins = numel(bin_start);
   tmp = struct();
   for n = 1:numel(fields)
      tmp.(fields{n}) = nan(size(ice2.Tice, 1), n_bins);
   end
   % Replace Z, if this is not a legacy grid run.
   if isfield(ice2, 'Z')
      tmp.Z = ice2.Z;
   end

   % Decide each field's rule. Every df_ field is a per-step increment, so it
   % sums; errH is a residual and also sums; everything else averages. Z is the
   % depth grid and is copied, not aggregated.
   do_sum = icemodel.isIncrementChannel(fields) | strcmp(fields, 'errH');
   skip = strcmp(fields, 'Z');

   % Aggregate the input samples assigned to each hourly bin.
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
function [ice1, ice2] = retimeLogical(ice1, ice2, bin_start, bin_end)
   %RETIMELOGICAL Aggregate logical flags over the surface hourly bins.

   % Aggregate each two-dimensional flag with the same hourly bin membership
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

   % TODO: 1-d logical
end

%%
function [ice1, ice2] = roundData(ice1, ice2)

   % Round ice1 channels to five digits (legacy method). Keep the diagnostic
   % mass-budget at double precision so the signed budget closures are testable.
   % Round one variable at a time. Table brace extraction concatenates selected
   % columns into one array. If that array includes the single-precision
   % Tsfc_converged or Tice_converged column, MATLAB converts double columns to
   % single before rounding. Per-variable brace assignment preserves each
   % column's original numeric class.
   vars1 = ice1.Properties.VariableNames;
   is_budget = ismember(vars1, icemodel.namelists.budgetoutputs());

   % df_rof is a per-step increment like ice2's df_liq/df_evp/df_lyr, so it
   % keeps their 8-digit precision. Rounding it to five digits would zero the
   % small overflow values the budget closures are checked against.
   is_increment = icemodel.isIncrementChannel(vars1);
   keep_precision = is_budget | is_increment;
   round_names = vars1(~keep_precision);
   increment_names = vars1(is_increment & ~is_budget);

   % The isnumeric guards are needed by this per-variable loop because round
   % rejects a logical column. RETIMELOGICAL can move a logical ice2 field
   % into ice1, so test each field before rounding it.
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
         case {'df_liq','df_lyr','df_evp','df_vap_liq','df_vap_ice', ...
               'Qsub','Sc','errT','errH'}
            ice2.(thisfield) = round(ice2.(thisfield), 8);
      end
   end

   % % For reference, another way to do it:
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
   %computeState compute post-processed state and energy balance.
   %
   % Compute bulk density (kg/m3), heat capacity (J/kg/K), thermal conductivity
   % (W/m/K), and a full surface and subsurface energy balance. Skip this for a
   % large simulation when time or disk space is limited. Compute them after the
   % simulation instead.

   % Pull out the state vectors.
   T_ice = ice2.Tice;
   f_liq = ice2.f_liq;
   f_ice = ice2.f_ice;

   % Compute bulk density, heat capacity, and thermal conductivity. The model
   % uses vapor-free thermal conductivity for the surface conductive heat flux.
   % Subsurface vapor energy transfer is computed at cell faces. For a nodal
   % measurement such as a needle probe, combine the stored diagnostics as
   % k_eff + (1 - f_ice - f_liq) .* k_vap.
   k_eff = icemodel.column.bulk_thermal_conductivity(T_ice, f_ice, f_liq, 0);
   k_vap = icemodel.vapor.vapor_thermal_conductivity(T_ice, f_liq);

   % Bulk density (kg/m3) and heat capacity (J/kg/K)
   ro_sno = icemodel.column.bulk_density(f_ice, f_liq);
   cp_sno = icemodel.column.bulk_specific_heat_capacity(f_ice, f_liq, ro_sno);

   % Compute a mesh for plotting
   Z = opts.z0_thermal;
   dz = opts.dz_thermal;

   % Assign values to ice2
   ice2.Z      = (dz/2:dz:Z-dz/2)';
   ice2.k_eff  = k_eff;                  % eff. thermal k
   ice2.k_vap  = k_vap;                  % raw vapor thermal conductivity
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
%
% The surface energy balance:
% Qm = chi*Qsi*(1-albedo) + Qln + Qh + Qe + Qc
%
% The subsurface shortwave balance:
% Qsub = (1-chi)*Qsi*(1-albedo)
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
% The total absorbed shortwave radiation equals the sum of the 'skin' absorbed
% radiation (Qsrf) and the subsurface absorbed radiation (Qsub) The subsurface
% absorbed radiation is the integral from z=0 to z=infty of dQ/dz*dz i.e.
% sum(Sc.*dz) with Sc the source term in units W/m3. chi is just the ratio of
% Qsrf to Qabs i.e. how much of the total absorbed energy is absorbed by the
% wall. Qsip is the penetrating radiation, some of which is absorbed and some of
% which is reflected (contributes to Qsr) in short, if we call the net solar
% downflux in the ice Q, then: Qsub = dQ is the absorbed solar downflux in the
% ice in each layer Sc = dQ/dz is the source term (divergence of the net solar
% downflux)


% %% The original method when I saved the grid output
%
% % Retime to hourly
% if opts.dt == 900
%    ice1  = retime(ice1,'hourly','mean');
%    feb29 = month(ice1.Time) == 2 & day(ice1.Time) == 29;
%    ice1  = ice1(~feb29,:);
%
%    % init tmp arrays to retime ice2
%    tmp.Tice = nan(size(ice2.Tice,1), numel(ice1.Time));
%
%    if strcmp('icemodel', opts.smbmodel)
%       tmp.df_liq = nan(size(ice2.df_liq,1), numel(ice1.Time));
%       tmp.f_ice = nan(size(ice2.f_ice,1),numel(ice1.Time));
%       tmp.f_liq = nan(size(ice2.f_liq,1),numel(ice1.Time));
%    end
%
%    for n = 1:numel(ice1.Time)
%       % this works b/c we know it's fifteen minute data.
%       i1 = n*4-3; i2 = n*4;
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
