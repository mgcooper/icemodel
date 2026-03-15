function [tair, swd, lwd, albedo, wspd, rh, psfc, rain, tppt, ...
      De, S, time] = METINIT(opts, fileiter)
   %METINIT initialize the met file
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

   if isvariable('rain', met)
      rain = met.rain / 3600; % opts.dt - the mar rain/snow is mWE / hr not / dt
      rain(isnan(rain)) = 0.0;
      rain = 0 * tair;
   else
      rain = 0 * tair;
   end

   % Solve for wet bulb
   [Ls, cp_air] = icemodel.physicalConstant('Ls', 'cp_air');
   tppt = nan(size(rh));
   for n = 1:numel(rh)
      tppt(n) = SOLVEWB(tair(n), rh(n), Ls, cp_air, psfc(n));
   end

   % Compute the wind transfer and stability coefficients
   [De, S] = WINDCOEF(wspd, opts.z_0, opts.z_tair, opts.z_wind);
end
