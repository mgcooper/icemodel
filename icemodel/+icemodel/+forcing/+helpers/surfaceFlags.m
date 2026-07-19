function flags = surfaceFlags(z, sensors, t, kwargs)
   %SURFACEFLAGS Per-sample quality flags for a PROMICE surface-height series.
   %
   %  flags = icemodel.forcing.helpers.surfaceFlags(z, sensors, t)
   %  flags = ... surfaceFlags(z, sensors, t, transition_times=...)
   %
   % Derives the per-sample flag channels staged alongside a PROMICE/GC-Net
   % surface-height series. We MODIFY no GEUS data here: these are faithful
   % per-sample masks read/derived from the L3 channels, never edits to the
   % series. buildPromiceData attaches them; consumers act on them.
   %
   % GAP FLAG (gap-bridged surface height). The L3 surface height (z_ice_surf or
   % z_surf_combined) is bridged across data gaps by a manual slope when ALL
   % surface-ranging sensors fail (readme, "Surface height estimation"): the
   % trend is preserved but a per-timestep RATE through the gap is not a direct
   % observation. The OLD heuristic flagged only samples where z itself is NaN,
   % which misses these slope-bridged segments (z is finite there, manufactured
   % by interpolation). This flag instead marks a sample gap-bridged when ALL of
   % the underlying contributing sensors (transducer/boom/stake ranges) are NaN
   % yet z is finite -> the surface value at that sample is interpolated, not
   % measured. Leading/trailing samples before the first / after the last finite
   % z are also flagged (no observation to anchor them). Where no sensor channel
   % is available the flag falls back to ~isfinite(z).
   %
   % STATION-TRANSITION FLAG (handover window). A PROMICE "site" can merge
   % several AWS ("stations") over time; at a handover the surface or subsurface
   % series can carry an expected discrete offset (readme Table 1 / Station vs
   % site). When known handover times are supplied (transition_times) this flag
   % marks samples within tol_days of any handover. This is distinct from the gap
   % flag: a transition is a step, not a NaN, so the gap flag never sees it (e.g.
   % a site with a known merge but no surface NaN at the handover).
   %
   % Inputs
   %  z       - double column, the surface-height series (may hold NaN)
   %  sensors - double matrix (nSamples x nSensors) of the underlying L3
   %            surface-ranging channels (transducer_depth, boom_height,
   %            stake_height as available); pass [] when none are present
   %  t       - datetime column, the series time axis
   %
   % Name-value
   %  transition_times : datetime array of known station-handover times ([])
   %  tol_days         : transition-window half-width [days] (default 14)
   %
   % Outputs
   %  flags - struct of per-sample double (0/1) masks aligned to z:
   %          .gap                gap-bridged (interpolated, non-observational)
   %          .station_transition within tol_days of a station handover
   %
   % See also: icemodel.forcing.buildPromiceData,
   %  icemodel.forcing.destepSurface

   arguments
      z (:, 1) double
      sensors double = []
      t (:, 1) datetime = NaT(numel(z), 1)
      kwargs.transition_times (:, 1) datetime = datetime.empty(0, 1)
      kwargs.tol_days (1, 1) double {mustBeNonnegative} = 14
   end

   n = numel(z);
   z_finite = isfinite(z);

   % Gap flag from the underlying sensors: a sample is gap-bridged when z is
   % finite (a value exists) but every contributing sensor is NaN there (no
   % direct measurement, so the value is slope-bridged). Falls back to the z-NaN
   % heuristic only when no sensor channel is available.
   if isempty(sensors)
      gap = ~z_finite;
   else
      all_sensors_nan = all(~isfinite(sensors), 2);
      gap = (z_finite & all_sensors_nan);
      % Leading/trailing samples (before first / after last finite z) carry no
      % observation to anchor them, so they are gap-bridged too.
      first = find(z_finite, 1, 'first');
      last = find(z_finite, 1, 'last');
      if isempty(first)
         gap = true(n, 1);   % no finite z at all: everything is non-observational
      else
         gap(1:first - 1) = true;
         gap(last + 1:end) = true;
         % A sample whose z is NaN is non-observational regardless of sensors.
         gap(~z_finite) = true;
      end
   end

   % Station-transition flag: mark samples within tol_days of any known handover.
   station_transition = false(n, 1);
   if ~isempty(kwargs.transition_times) && ~all(isnat(t))
      for k = 1:numel(kwargs.transition_times)
         station_transition = station_transition | ...
            abs(days(t - kwargs.transition_times(k))) <= kwargs.tol_days;
      end
   end

   flags = struct('gap', double(gap), ...
      'station_transition', double(station_transition));
end
