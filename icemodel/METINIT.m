function [tair, swd, lwd, albedo, wspd, rh, psfc, rain, tppt, time, ...
      snow_depth] = METINIT(opts, fileiter)
   %METINIT Load and expand the meteorological forcing vectors.
   %
   %  [tair, swd, lwd, albedo, wspd, rh, psfc, rain, tppt, De, br_coefs, time] ...
   %     = METINIT(opts)
   %  ... = METINIT(opts, fileiter)
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
   if isvariable('snow_depth', met)
      snow_depth = met.snow_depth;
   else
      snow_depth = nan(size(tair));
   end

   % Rainfall forcing is still ignored in the core time integration.
   % Keep the legacy zero-rain behavior explicit here until rain/snow/ppt
   % forcing support is implemented consistently.
   rain = 0 * tair;

   % TODO: support snowfall, and confirm if forcing files are consistent wrt
   % rain/snow/ppt/prec variable names. The optional forcing-snow-depth hook
   % used by the THF roughness selector is standardized separately as
   % `snow_depth`, but it does not imply a full snow-model mass/energy
   % treatment and may remain NaN in existing station datasets.

   % Solve for wet bulb
   tppt = nan(size(rh));
   for n = 1:numel(rh)
      tppt(n) = icemodel.vapor.wet_bulb_temperature(tair(n), rh(n), psfc(n));
   end
end
