function [valid, lower_limit, upper_limit] = physicalValidity( ...
      channel, values, times, kwargs)
   %PHYSICALVALIDITY Enforce scalar and relational reconstruction bounds.
   %
   %  valid = icemodel.forcing.reconstruct.physicalValidity( ...
   %     channel, values, times, latitude=67, longitude=-49)
   %  valid = icemodel.forcing.reconstruct.physicalValidity( ...
   %     "swu", values, times, swd=swd)
   %
   % Static channel limits come from physicalBounds. Downwelling shortwave
   % additionally cannot exceed the configured multiple of
   % top-of-atmosphere irradiance
   % (with a small absolute night-noise floor on the ceiling and a civil
   % twilight allowance — POLICY D-28), and upward shortwave cannot exceed
   % the paired downwelling value.
   %
   % Name-value
   %  latitude, longitude : site point, required for the swd solar ceiling.
   %  swd : paired downwelling reference for the swu relational bound.
   %  toa, elevation : optional PRECOMPUTED solar geometry, one value per
   %     sample of `values` — toaIrradiance at `times` and the maximum solar
   %     elevation over each posting interval. Hot-loop callers (the
   %     run x method walk in reconstructSeries, the per-fragment tier-1
   %     checks in fillShortGaps, the per-outage probes in
   %     lastResortProxies) evaluate the geometry once on their full time
   %     axis and pass slices here.
   %  interval : support duration used to compute maximum elevation when
   %     `elevation` is omitted. Zero preserves point-sample behavior for
   %     direct callers; production reconstruction supplies its cadence.
   %
   % Returns
   %  valid : logical mask, one entry per sample.
   %  lower_limit, upper_limit : per-sample effective bounds actually
   %     enforced (scalar registry bounds tightened by the relational
   %     ceiling where one applies). Callers that must CAP an estimate at
   %     the validity limit instead of rejecting it — the seam-blend
   %     fix-up in lastResortProxies — read these so the limit stays
   %     single-sourced here.

   arguments
      channel (1, 1) string
      values (:, 1) double
      times datetime
      kwargs.latitude (1, 1) double = NaN
      kwargs.longitude (1, 1) double = NaN
      kwargs.swd (:, 1) double = nan(0, 1)
      kwargs.toa (:, 1) double = nan(0, 1)
      kwargs.elevation (:, 1) double = nan(0, 1)
      kwargs.interval (1, 1) duration = hours(0)
   end

   if numel(values) ~= numel(times)
      error('icemodel:reconstruct:physicalValidity:sizeMismatch', ...
         'values and times must share one sample axis');
   end
   if kwargs.interval < hours(0)
      error('icemodel:reconstruct:physicalValidity:negativeInterval', ...
         'interval must be nonnegative');
   end

   % Every channel first passes its dedicated scalar registry.
   bounds = icemodel.forcing.reconstruct.physicalBounds(channel);
   valid = isfinite(values) & values >= bounds(1) & values <= bounds(2);
   % Per-sample limits default to the scalar registry; relational branches
   % tighten the ceiling below.
   lower_limit = repmat(bounds(1), numel(values), 1);
   upper_limit = repmat(bounds(2), numel(values), 1);

   if channel == "swd"
      % Solar geometry is required to prove the policy's daylight ceiling;
      % refusing unknown geometry is safer than admitting unbounded flux.
      if ~(isfinite(kwargs.latitude) && isfinite(kwargs.longitude))
         error('icemodel:reconstruct:physicalValidity:missingSolarGeometry', ...
            'swd validity requires latitude and longitude');
      end
      % Precomputed geometry (see the docstring) short-circuits the NOAA
      % evaluation; a wrong-length vector is refused because silently
      % recomputing would hide the caller's slicing bug.
      toa = kwargs.toa;
      if isempty(toa)
         toa = icemodel.forcing.reconstruct.toaIrradiance(times, ...
            kwargs.latitude, kwargs.longitude);
      elseif numel(toa) ~= numel(values)
         error(['icemodel:reconstruct:physicalValidity:' ...
            'precomputedGeometrySize'], ...
            'precomputed toa must contain one value per sample');
      end
      % Pyranometers report small positive thermal-offset noise in
      % darkness, so a bare 1.05x ceiling brands every polar-night
      % sample invalid (EGP winters carry 0.4-1.4 W/m2 at TOA = 0). A
      % noise floor on the ceiling admits instrument reality while still
      % rejecting garbage; like the lwd floor, the bound VALUE is a
      % parameter-level choice under the A15/D-25 principle (ratified as
      % D-28: the 5 W/m2 term is a minimum CEILING in darkness, never a
      % floor on data).
      bands = icemodel.forcing.reconstruct.solarElevationBands();
      ceiling = max(bands.toa_ceiling_multiplier * toa, ...
         bands.toa_ceiling_floor_wm2);
       % Civil twilight (the complete posting reaches between 0 deg and
       % the civil-twilight boundary) scatters real diffuse irradiance of
       % order 8-28 W/m2
      % onto the surface while the geometric TOA model reads exactly
      % zero: that light crossed the atmosphere and IS incident energy,
      % so the ceiling admits up to 50 W/m2 there instead of branding
      % genuine dusk/dawn samples invalid (POLICY D-28; the 50 W/m2
      % allowance is a Section-C-style parameter sized above the
      % observed twilight range). Elevation comes from the same NOAA
      % solar geometry toaIrradiance itself uses, and the boundary comes
      % from the solarElevationBands single source, so the darkness
      % pre-pass, the binned calibration, and this ceiling can never
      % disagree about where twilight is.
       % Use the maximum over the posting support, matching PROMICE staging
       % and darkness reconstruction. A posting that begins in deep night
       % but enters twilight must not be judged from its start instant.
       elevation = kwargs.elevation;
       if isempty(elevation)
          elevation = ...
             icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
             times, kwargs.latitude, kwargs.longitude, kwargs.interval);
      elseif numel(elevation) ~= numel(values)
         error(['icemodel:reconstruct:physicalValidity:' ...
            'precomputedGeometrySize'], ...
            'precomputed elevation must contain one value per sample');
      end
      twilight = elevation(:) < 0 ...
         & elevation(:) >= bands.civil_twilight_deg;
      ceiling(twilight) = max(ceiling(twilight), ...
         bands.twilight_ceiling_wm2);
      valid = valid & values <= ceiling;
      upper_limit = min(upper_limit, ceiling);
   elseif channel == "swu"
      % Upward shortwave is meaningful only beside a finite downwelling
      % reference on the same axis.
      if numel(kwargs.swd) ~= numel(values)
         error('icemodel:reconstruct:physicalValidity:missingShortwaveReference', ...
            'swu validity requires one swd value per sample');
      end
      valid = valid & isfinite(kwargs.swd) & values <= kwargs.swd;
      % min ignores NaN, so a missing swd reference leaves the scalar
      % ceiling in place; validity above already rejects those samples.
      upper_limit = min(upper_limit, kwargs.swd);
   end
end
