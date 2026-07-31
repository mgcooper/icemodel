function [x, filled, audit] = fillTwilightClimatology( ...
      times, x, native, latitude, longitude, kwargs)
   %FILLTWILIGHTCLIMATOLOGY Fill one-posting SWD gaps beside known night.
   %
   %  [x, filled, audit] = ...
   %     icemodel.forcing.reconstruct.fillTwilightClimatology( ...
   %     times, x, native, latitude, longitude)
   %
   % A single still-missing civil-twilight posting with exactly one
   % adjacent all-interval darkness posting uses the existing station
   % day-of-year/posting climatology. The untouched native series is the
   % only support pool; finite input samples are never modified.

   arguments
      times datetime
      x (:, 1) double
      native (:, 1) double
      latitude (1, 1) double {mustBeFinite}
      longitude (1, 1) double {mustBeFinite}
      kwargs.min_support (1, 1) double ...
         {mustBeInteger, mustBePositive} = ...
         icemodel.forcing.reconstruct.solarElevationBands(). ...
         twilight_climatology_min_support
   end

   if numel(times) ~= numel(x) || numel(native) ~= numel(x)
      error('icemodel:reconstruct:fillTwilightClimatology:sizeMismatch', ...
         'times, x, and native must have equal lengths');
   end
   if numel(times) < 2
      error('icemodel:reconstruct:fillTwilightClimatology:missingCadence', ...
         'at least two postings are required');
   end

   interval = median(diff(times));
   maximum_elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, latitude, longitude, interval);
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   dark = maximum_elevation <= bands.civil_twilight_deg;
   twilight = maximum_elevation > bands.civil_twilight_deg ...
      & maximum_elevation <= bands.horizon_deg;
   target = find(~isfinite(x) & twilight);
   target = target(target > 1 & target < numel(x));
   if isempty(target)
      filled = false(size(x));
      audit = cell(0, 1);
      return
   end

   one_dark_anchor = xor(dark(target - 1), dark(target + 1));
   finite_light_anchor = ...
      (dark(target - 1) & isfinite(x(target + 1))) ...
      | (dark(target + 1) & isfinite(x(target - 1)));
   target = target(one_dark_anchor & finite_light_anchor);
   if isempty(target)
      filled = false(size(x));
      audit = cell(0, 1);
      return
   end

   [candidate, n_support] = ...
      icemodel.forcing.reconstruct.climatologyFill( ...
      times, native, times(target), min_support=kwargs.min_support);
   valid = icemodel.forcing.reconstruct.scalarValidity("swd", candidate) ...
      & candidate <= bands.twilight_ceiling_wm2;
   target = target(valid);
   candidate = candidate(valid);
   n_support = n_support(valid);
   filled = false(size(x));
   audit = cell(0, 1);
   for k = 1:numel(target)
      x(target(k)) = candidate(k);
      filled(target(k)) = true;
      segment = false(size(x));
      segment(target(k)) = true;
      rows = icemodel.forcing.reconstruct.auditSegments( ...
         times, segment, "swd", "twilight_climatology", sprintf( ...
         ['day-of-year/posting median; support %d; ' ...
         'interval maximum %.3g deg'], n_support(k), ...
         maximum_elevation(target(k))));
      audit = [audit; rows]; %#ok<AGROW>
   end
end
