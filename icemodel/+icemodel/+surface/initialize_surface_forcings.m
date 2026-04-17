function [tair, swd, lwd, albedo, wspd, rh, psfc, rain, tppt, time, ...
      snow_depth] = initialize_surface_forcings(opts, fileiter)
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
   %
   %#codegen

   % The 2nd input is the index into the metfile name list resolved in
   % icemodel.configureRun / icemodel.setopts. If omitted, load and
   % concatenate all files listed in opts.metfname.
   if nargin < 2
      met = icemodel.loadmet(opts);
   else
      met = icemodel.loadmet(opts, fileiter);
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

   % Rainfall forcing is ignored in the core time integration.
   % Keep the zero-rain behavior explicit here until rain/snow/ppt forcing
   % support is implemented consistently (Jordan 1991 requires special
   % timestep shortening during accumulation events).
   % When rain is eventually set from forcing data (typically in m/timestep),
   % it must be converted to a rate in m s^-1 before being passed to
   % advective_heat_flux, which expects ppt in m s^-1.
   rain = 0 * tair;

   % TODO: support snowfall, and confirm if forcing files are consistent wrt
   % rain/snow/ppt/prec variable names. The optional forcing-snow-depth hook
   % used by the THF roughness selector is standardized separately as
   % `snow_depth`, but it does not imply a full snow-model mass/energy
   % treatment and may remain NaN in existing station datasets.
   %
   % Legacy forcing-derivation fallbacks for station datasets that omit
   % `lwd`, `swd`, or `psfc` now live under `icemodel.surface`:
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
