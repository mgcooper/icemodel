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

   % Rainfall forcing is still ignored in the core time integration.
   % Keep the legacy zero-rain behavior explicit here until rain/snow/ppt
   % forcing support is implemented consistently.
   rain = 0 * tair;

   % TODO: support snowfall, and confirm if forcing files are consistent wrt
   % rain/snow/ppt/prec variable names

   % Solve for wet bulb
   tppt = nan(size(rh));
   for n = 1:numel(rh)
      tppt(n) = TWETBULB(tair(n), rh(n), psfc(n));
   end

   % The canonical met loader already computes De after all swaps/subsetting.
   De = met.De;
   [~, S] = WINDCOEF(wspd, opts.z_0, opts.z_tair, opts.z_wind);
end
