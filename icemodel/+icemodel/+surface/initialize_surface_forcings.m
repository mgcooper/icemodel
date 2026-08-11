function [tair, swd, lwd, albedo, wspd, rh, psfc, rain, tppt, time, ...
      snow_depth, opts, rainf, snowf] ...
      = initialize_surface_forcings(opts, fileiter)
   %initialize_surface_forcings Load the meteorological forcing vectors.
   %
   %  [tair, swd, lwd, albedo, wspd, rh, psfc, rain, tppt, time] ...
   %     = icemodel.surface.initialize_surface_forcings(opts)
   %  ... = icemodel.surface.initialize_surface_forcings(opts, fileiter)
   %
   % Outputs:
   %  tair   - air temperature [K]
   %  swd    - downwelling shortwave radiation [W m^-2]
   %  lwd    - downwelling longwave radiation [W m^-2]
   %  albedo - surface albedo [1]
   %  wspd   - wind speed [m s^-1]
   %  rh     - relative humidity [%]
   %  psfc   - surface pressure [Pa]
   %  rain   - rainfall-rate placeholder [kg m^-2 s^-1]
   %  tppt   - precipitation wet-bulb temperature [K]
   %  time   - forcing timestamps [datetime]
   %  snow_depth - optional forcing snow depth [m]; NaN when unavailable
   %  opts    - runtime options with time-varying PROMICE observation
   %            heights resolved by icemodel.loadmet
   %  rainf  - phase-source-selected liquid precipitation rate [m s^-1]
   %  snowf  - phase-source-selected solid precipitation rate [m s^-1]
   %
   % The trailing rainf/snowf outputs carry the runtime precipitation phase
   % selection (opts.precip_phase_source, POLICY A10/D-18). 'source' returns
   % the met product's own split as shipped. 'threshold' repartitions the
   % canonical total ppt by air temperature. They come after OPTS, so the
   % positions of the earlier outputs do not change. Under POLICY D-0b they
   % feed no existing physics: the RAIN output stays zero until the
   % advective-rain physics is finished.
   %
   %#codegen

   % The 2nd input is the index into the metfile name list resolved in
   % icemodel.configureRun / icemodel.setopts. If omitted, load and
   % concatenate all files listed in opts.metfname.
   if nargin < 2
      [met, opts] = icemodel.loadmet(opts);
   else
      [met, opts] = icemodel.loadmet(opts, fileiter);
   end

   % Transfer the met data to vectors
   rh = met.rh;
   swd = met.swd;
   lwd = met.lwd;
   tair = met.tair;
   wspd = met.wspd;
   psfc = met.psfc;
   time = met.Time;
   albedo = met.albedo;
   if ismember('snow_depth', met.Properties.VariableNames)
      snow_depth = met.snow_depth;
   else
      snow_depth = nan(height(met), 1);
   end

   % Rainfall forcing is ignored in the core time integration. Keep the
   % zero-rain behavior explicit until rain/snow/ppt forcing support is
   % implemented consistently and enabled as a separate physics change
   % (POLICY D-0b; the rainf/snowf outputs below carry the selected data
   % split without entering the solver).
   rain = 0 * tair;

   % Runtime precipitation phase selection (POLICY A10 / D-18). The option
   % opts.precip_phase_source picks the rainf/snowf split that
   % snowfall-consuming callers receive. The resolution helper enforces the
   % A10 validity contract: nonnegative components that sum to the total.
   % The helper parses options and defaults outside the code-generation
   % subset, so it stays behind the MATLAB-target boundary, like the
   % readiness verifier in icemodel.loadmet. Generated targets return the
   % unresolved NaN sentinel until a generated consumer exists.
   if coder.target('MATLAB')
      % A cached options struct can lack this field. Default to the
      % product's own split (the icemodel.setopts default), so an old
      % struct keeps its behavior.
      phase_source = 'source';
      if isfield(opts, 'precip_phase_source')
         phase_source = opts.precip_phase_source;
      end
      [rainf, snowf] = icemodel.resolvePrecipPhase( ...
         optionalMetColumn(met, 'ppt'), tair, ...
         optionalMetColumn(met, 'rainf'), ...
         optionalMetColumn(met, 'snowf'), phase_source);
   else
      rainf = nan(size(tair));
      snowf = nan(size(tair));
   end

   % TODO: support rainfall and snowfall mass/state evolution. The optional
   % snow-depth input that the THF roughness selector uses is standardized
   % separately as `snow_depth`. It does not imply a full snow-model mass and
   % energy treatment, and it can stay NaN in existing station datasets.
   %
   % Legacy forcing-derivation fallbacks for station datasets that omit
   % `lwd`, `swd`, or `psfc` live under `icemodel.surface`:
   %   empirical_incoming_longwave_radiation
   %   incoming_shortwave_radiation
   %   terrain_adjusted_shortwave_radiation
   %   atmospheric_pressure_from_elevation
   % Wire them in here explicitly if a future forcing workflow needs those
   % derived series instead of direct measured inputs.

   % Solve for wet bulb for use in the advective heat flux calculation.
   tppt = nan(size(rh));
   for n = 1:numel(rh)
      tppt(n) = icemodel.vapor.wet_bulb_temperature(tair(n), rh(n), psfc(n));
   end

end

%%
function values = optionalMetColumn(met, name)
   %OPTIONALMETCOLUMN Return a met column or an all-NaN placeholder.
   %
   % An absent channel stays missing instead of zero-filled, so the phase
   % resolution never creates precipitation from a source that has none.
   if ismember(name, met.Properties.VariableNames)
      values = met.(name);
   else
      values = nan(height(met), 1);
   end
end
