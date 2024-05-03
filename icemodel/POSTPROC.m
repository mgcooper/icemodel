function varargout = POSTPROC(ice1, ice2, opts, varargin)
   %POSTPROC post process the simulation data
   %
   % Syntax:
   %
   % [ice1, ice2] = POSTPROC(ice1, ice2, opts, simyear)
   % [ice1, ice2] = POSTPROC(ice1, ice2, opts, swd, lwd, albedo, time)
   %
   % Description:
   %
   % [ice1, ice2] = POSTPROC(ice1, ice2, opts, simyear) Loads the met file
   % based on the information in opts and the simulation year. The year is
   % needed because the opts struct can be setup for multi-year simulations.
   %
   % [ice1, ice2] = POSTPROC(ice1, ice2, opts, swd, lwd, albedo, time)
   %
   % See also:

   % Process inputs
   if nargin == 7
      [swd, lwd, albedo, time] = deal(varargin{:});
   elseif nargin == 4
      simyear = varargin{1};
      met = icemodel.loadmet(opts, find(opts.simyears == simyear));
      swd = met.swd;
      lwd = met.lwd;
      albedo = met.albedo;
      time = met.Time;
   else
      error('unrecognized number of inputs')
   end

   % Load physical constants
   [Tf, ro_liq, Ls, Lf] = icemodel.physicalConstant('Tf', 'ro_liq', 'Ls', 'Lf');

   % Calculate runoff
   if strcmp('skinmodel', opts.smbmodel)
      ice1 = SRFRUNOFF(ice1, ro_liq, Ls, Lf, opts.dt);
   elseif strcmp('icemodel', opts.smbmodel)
      ice1 = ICERUNOFF(ice1, ice2, opts);
   end
   % ice1.runoff2 = ice1.melt - ice1.freeze;

   % Compute a full state and energy balance
   if ~strcmp(opts.sitename, 'sector')
      [ice1, ice2] = computeState(ice1, ice2, opts, swd, lwd, albedo);
   end

   % Convert surface and subsurface ice temperature from Kelvin to Celsius
   % TODO: This is not right - this adds the triple point temperature. Besides
   % there's no reason to convert to Celsius. But this must wait for the next
   % release.
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

   % Retime 15 min data to hourly
   if opts.dt == 900
      ice1 = retime(ice1, 'hourly', 'mean');
      [ice1, ice2] = retimeLogical(ice1, ice2);
      ice1 = ice1(~(month(ice1.Time) == 2 & day(ice1.Time) == 29), :);
      ice2 = retimeIce2(ice2, ice1.Time);
   end

   % Round the data to save disk space, retaining necessary precision
   [ice1, ice2] = roundData(ice1, ice2);

   if ~strcmp(opts.sitename, 'sector')
      ice2.Time = ice1.Time; % not added in the sector case, maybe remove.

      % Rename ice1 vars to match the naming conventions i use everywhere else
      oldvars = {'Qsi','Qsr','Qsn','Qli','Qle','Qln','Qh','Qe','Qc','Qn','Tsfc'};
      newvars = {'swd','swu','swn','lwd','lwu','lwn','shf','lhf','chf','netr','tsfc'};
      ice1 = renamevars(ice1, ...
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
            met = icemodel.loadmet(opts, find(opts.simyears == year(time(1))));
         end
         varargout{3} = icemodel.processmet(met);
   end
end

%%
function ice2 = retimeIce2(ice2, Time)

   % Get the field names of ice2
   fields = fieldnames(ice2);

   % Initialize temporary arrays to retime ice2
   tmp = struct();
   for n = 1:numel(fields)
      tmp.(fields{n}) = nan(size(ice2.Tice, 1), numel(Time));
   end
   % Replace Z, if this is not a sector run
   if isfield(ice2, 'Z')
      tmp.Z = ice2.Z;
   end

   % Loop over each hour
   for n = 1:numel(Time)
      % this works b/c we know it's fifteen minute data
      i1 = n*4-3;
      i2 = n*4;

      % Iterate over fields to apply appropriate retime method
      for m = 1:numel(fields)
         thisfield = fields{m};
         switch thisfield
            case {'Z'} % might need {'Z', 'LCflag'} and other logicals/etc
               continue
            case {'Tice', 'f_ice', 'f_liq', 'k_eff', 'k_vap', 'ro_sno', 'cp_sno'}
               tmp.(thisfield)(:, n) = mean(ice2.(thisfield)(:, i1:i2), 2);
            case {'df_liq', 'df_evp', 'df_lyr', 'errH'}
               tmp.(thisfield)(:, n) = sum(ice2.(thisfield)(:, i1:i2), 2);
            otherwise % assume averaging is correct
               tmp.(thisfield)(:, n) = mean(ice2.(thisfield)(:, i1:i2), 2);
         end
      end
   end
   ice2 = tmp;
end

%%
function [ice1, ice2] = retimeLogical(ice1, ice2)
   % Retime logical flags

   % 2-d logical
   fields = fieldnames(ice2);
   for n = 1:numel(fields)
      thisfield = fields{n};
      if islogical(ice2.(thisfield)(1, 1))
         flag = false(size(ice2.f_ice));
         idx = 0;
         for mm = 1:4:size(ice2.(thisfield), 2) - 3
            idx = idx+1;
            flag(:, idx) = sum(ice2.(thisfield)(:, mm:mm+3), 2) > 0;
         end
         ice1.(thisfield) = transpose(flag(1, :));
         ice2 = rmfield(ice2, thisfield);
      end
   end

   % 1-d logical

end

%%
function [ice1, ice2] = roundData(ice1, ice2)

   % Round the ice1 data to five digits
   ice1{:, ice1.Properties.VariableNames} = ...
      round(ice1{:, ice1.Properties.VariableNames}, 5);

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
function [ice1, ice2] = computeState(ice1, ice2, opts, swd, lwd, albedo)

   % Compute bulk density, heat capacity, thermal conductivity, and a full
   % surface and subsurface energy balance. Don't do this for large simulations
   % if time or disk space is limited, instead compute them after the
   % simulation.

   % Load physical constants
   [Tf,emissSB,cv_ice,cv_liq,ro_ice,ro_liq,ro_air,k_liq,Ls,Rv,emiss] = ...
      icemodel.physicalConstant('Tf','emissSB','cv_ice','cv_liq','ro_ice', ...
      'ro_liq','ro_air','k_liq','Ls','Rv','emiss');

   % Compute bulk density (kg/m3), heat capacity (J/kg/K), thermal K (W/m/K)
   T_ice = ice2.Tice;
   f_liq = ice2.f_liq;
   f_ice = ice2.f_ice;

   [k_eff, k_vap] = GETGAMMA(T_ice, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);
   ro_sno = ro_ice * f_ice + ro_liq * f_liq + ro_air * (1.0 - f_liq - f_ice);
   cp_sno = (cv_ice * f_ice + cv_liq * f_liq) ./ ro_sno;

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
   ice1.albedo = albedo;                        % albedo
   ice1.swd    = swd;                           % shortwave down
   ice1.swu    = swd.*albedo;                   % shortwave up
   ice1.lwd    = emiss.*lwd;                    % longwave down
   ice1.lwu    = emissSB.*min(ice1.Tsfc,Tf).^4; % longwave up
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


% The original method when I saved the sector output
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
