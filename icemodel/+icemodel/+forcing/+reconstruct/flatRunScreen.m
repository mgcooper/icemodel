function [flagged, findings] = flatRunScreen(met, latitude, longitude, kwargs)
   %FLATRUNSCREEN Flag multi-day buried/rime-encased sensor runs in met data.
   %
   %  [flagged, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
   %     met, 62.3, -49.5)
   %  [flagged, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
   %     met, latitude, longitude, max_tair_range_k=0.3, min_run_days=5)
   %
   % Role
   %  Pre-reconstruction QA screen for the buried / rime-encased station
   %  signature (bead icemodel-g1n.45): a station whose sensors are under
   %  snow or encased in rime keeps reporting, and upstream PROMICE QC can
   %  certify those samples as valid observations. The verified reference
   %  case is FRE 2018-03-25 to 2018-04-15 — 22 consecutive days with
   %  tair daily range < 0.15 K, rh pinned at 85 %, swd < 1 W/m2 under
   %  April daylight, and lwd equal to sigma*tair^4 (the radiometer sees
   %  the rime shell at instrument temperature, a blackbody). The screen
   %  flags multi-day runs where
   %     (a) the daily tair range is near zero (core condition), and
   %     (b) swd stays near zero while top-of-atmosphere irradiance is
   %         high (dead pyranometer under daylight), or
   %     (c) lwd tracks the blackbody emission sigma*tair^4 within a
   %         tolerance (radiometer viewing its own encasement).
   %  Under POLICY A1 the native record is immutable: this screen never
   %  modifies data. Callers must exclude flagged samples from training
   %  pools (climatology, donor fits, step scales) and report the runs;
   %  the natives themselves ship unmodified with the flag as audit.
   %
   % Inputs
   %  met : timetable with a tair channel [K]; swd [W/m2], lwd [W/m2],
   %     and rh [percent] are used as corroborating evidence when
   %     present. Any other channels are ignored.
   %  latitude, longitude : station point (WGS84 degrees, east positive)
   %     for the top-of-atmosphere daylight reference.
   %
   % Name-value (defaults are the single source of these thresholds,
   % following the physicalBounds pattern; retuning is a policy change)
   %  min_run_days : consecutive qualifying days required before a run is
   %     flagged. Default 3 — long enough that persistent synoptic
   %     overcast (1-2 days of small tair range) does not trip the
   %     screen; the FRE reference run lasted 22 days.
   %  max_tair_range_k : daily (max - min) tair below which a day counts
   %     as flat. Default 0.5 K — a free-standing sensor sees several K
   %     of diurnal plus synoptic range; FRE ran below 0.15 K.
   %  max_swd_wm2 : daily-max swd below which the pyranometer counts as
   %     dark. Default 5 W/m2 (matches the physicalValidity night-noise
   %     floor); FRE topped out at 0.9 W/m2 under April daylight.
   %  min_toa_wm2 : daily-max top-of-atmosphere irradiance required
   %     before darkness is evidence — polar night must not corroborate.
   %     Default 100 W/m2.
   %  lwd_blackbody_tol_wm2 : daily mean |lwd - sigma*tair^4| below which
   %     lwd counts as blackbody-locked. Default 15 W/m2 — clear skies
   %     sit 20-90 W/m2 below blackbody, overcast approaches it but the
   %     conjunction with a flat tair separates the cases; FRE ran ~5.
   %  max_rh_range_pct : run-level rh (max - min) below which rh is
   %     reported as pinned. Evidence annotation only, never a flag
   %     condition. Default 1 percent; FRE pinned exactly at 85.
   %  min_daily_coverage : fraction of a day's expected samples that must
   %     be finite before the day is evaluable. Default 0.5 so sparse
   %     days neither trip nor corroborate the screen.
   %  require_corroboration : when true (default, the policy setting) a
   %     flat-tair day must also satisfy (b) or (c) to qualify; set false
   %     to screen on tair alone (diagnostic use).
   %
   % Outputs
   %  flagged : height(met) x 1 logical — true for every sample inside a
   %     flagged run (whole calendar days).
   %  findings : one row per flagged run — start_time, end_time (first
   %     and last flagged sample times), n_days, channels (comma-joined
   %     implicated channels), and the evidence statistics
   %     tair_range_max_k (largest daily range in the run), swd_max_wm2
   %     (largest daily-max swd), toa_max_wm2 (largest daily-max TOA),
   %     lwd_dev_wm2 (run mean of the daily blackbody deviation), and
   %     rh_range_pct (run-level rh range). Stats are NaN when the
   %     channel is absent.
   %
   % See also: icemodel.forcing.reconstruct.physicalValidity,
   %  icemodel.forcing.reconstruct.toaIrradiance,
   %  icemodel.forcing.reconstruct.fillPromiceStation

   arguments
      met timetable
      latitude (1, 1) double {mustBeFinite}
      longitude (1, 1) double {mustBeFinite}
      kwargs.min_run_days (1, 1) double {mustBePositive} = 3
      kwargs.max_tair_range_k (1, 1) double {mustBePositive} = 0.5
      kwargs.max_swd_wm2 (1, 1) double {mustBePositive} = 5
      kwargs.min_toa_wm2 (1, 1) double {mustBePositive} = 100
      kwargs.lwd_blackbody_tol_wm2 (1, 1) double {mustBePositive} = 15
      kwargs.max_rh_range_pct (1, 1) double {mustBePositive} = 1
      kwargs.min_daily_coverage (1, 1) double ...
         {mustBeInRange(kwargs.min_daily_coverage, 0, 1)} = 0.5
      kwargs.require_corroboration (1, 1) logical = true
   end

   % The screen is meaningless without air temperature: the flat-tair
   % condition is the core signature every corroboration hangs off.
   varnames = string(met.Properties.VariableNames);
   if ~any(varnames == "tair")
      error('icemodel:reconstruct:flatRunScreen:missingTair', ...
         'flatRunScreen requires a tair channel in met');
   end

   % Empty findings carry the documented schema so callers can vertcat
   % per-station results without special-casing a no-findings station.
   findings = table('Size', [0, 9], 'VariableTypes', ...
      {'datetime', 'datetime', 'double', 'string', 'double', ...
      'double', 'double', 'double', 'double'}, 'VariableNames', ...
      {'start_time', 'end_time', 'n_days', 'channels', ...
      'tair_range_max_k', 'swd_max_wm2', 'toa_max_wm2', ...
      'lwd_dev_wm2', 'rh_range_pct'});
   times = met.Properties.RowTimes;
   flagged = false(height(met), 1);
   if isempty(times)
      return
   end

   % Group samples into calendar days on the record's own clock; daily
   % aggregation is what turns the diurnal cycle into the range statistic
   % the buried signature suppresses.
   day_start = dateshift(times, 'start', 'day');
   [day_list, ~, gid] = unique(day_start);
   n_days_total = numel(day_list);

   % A day only participates when enough of its expected samples are
   % finite: a mostly-missing day can neither trip nor corroborate. The
   % expected count comes from the record's median cadence.
   if isscalar(times)
      expected_per_day = 1;
   else
      expected_per_day = max(1, round(days(1) / median(diff(times))));
   end
   tair = met.tair;
   tair_count = accumarray(gid(isfinite(tair)), 1, [n_days_total, 1]);
   evaluable = tair_count ./ expected_per_day >= kwargs.min_daily_coverage;

   % Core condition (a): the daily tair range collapses toward zero when
   % the sensor equilibrates inside snow or rime instead of the air.
   tair_range = dailyStat(gid, n_days_total, tair, @max) ...
      - dailyStat(gid, n_days_total, tair, @min);
   cond_tair = evaluable & tair_range < kwargs.max_tair_range_k;

   % Corroboration (b): swd pinned at zero while the top of the
   % atmosphere is bright — a working pyranometer cannot be dark under
   % daylight, but polar night must not count as evidence.
   toa = icemodel.forcing.reconstruct.toaIrradiance( ...
      times, latitude, longitude);
   toa_max = dailyStat(gid, n_days_total, toa, @max);
   swd_max = nan(n_days_total, 1);
   cond_swd = false(n_days_total, 1);
   if any(varnames == "swd")
      swd = met.swd;
      swd_count = accumarray(gid(isfinite(swd)), 1, [n_days_total, 1]);
      swd_max = dailyStat(gid, n_days_total, swd, @max);
      cond_swd = swd_count ./ expected_per_day ...
         >= kwargs.min_daily_coverage ...
         & toa_max >= kwargs.min_toa_wm2 ...
         & swd_max < kwargs.max_swd_wm2;
   end

   % Corroboration (c): lwd locked to the blackbody emission at air
   % temperature — a radiometer under rime radiatively views its own
   % encasement at instrument temperature instead of the sky.
   sigma = icemodel.physicalConstant('SB');
   lwd_dev = nan(n_days_total, 1);
   cond_lwd = false(n_days_total, 1);
   if any(varnames == "lwd")
      % The deviation needs both channels finite on the same sample.
      dev = abs(met.lwd - sigma .* tair .^ 4);
      dev_count = accumarray(gid(isfinite(dev)), 1, [n_days_total, 1]);
      lwd_dev = dailyStat(gid, n_days_total, dev, @mean);
      cond_lwd = dev_count ./ expected_per_day ...
         >= kwargs.min_daily_coverage ...
         & lwd_dev < kwargs.lwd_blackbody_tol_wm2;
   end

   % A qualifying day is flat tair plus at least one radiation
   % corroboration; the diagnostic escape hatch screens on tair alone.
   if kwargs.require_corroboration
      buried_day = cond_tair & (cond_swd | cond_lwd);
   else
      buried_day = cond_tair;
   end

   % Runs must be calendar-consecutive: a new run starts on every
   % qualifying day whose predecessor is not a qualifying immediately
   % preceding day (an axis hole splits a run like a clean day does).
   prev_consecutive = [false; days(diff(day_list)) == 1];
   run_start = buried_day ...
      & ~(prev_consecutive & [false; buried_day(1:end - 1)]);
   run_id = cumsum(run_start);
   run_id(~buried_day) = 0;

   % Keep only runs long enough to outlast synoptic persistence: count
   % the qualifying days under each run id, then drop the short runs so
   % the findings arrays can be preallocated exactly.
   run_ids = unique(run_id(run_id > 0));
   run_lengths = arrayfun(@(id) nnz(run_id == id), run_ids);
   run_ids = run_ids(run_lengths >= kwargs.min_run_days);
   n_runs = numel(run_ids);
   if n_runs == 0
      return
   end

   % One findings row per surviving run, assembled into preallocated
   % arrays and converted to the schema table in a single construction.
   start_time = NaT(n_runs, 1, 'TimeZone', times.TimeZone);
   end_time = start_time;
   n_days = zeros(n_runs, 1);
   channels = strings(n_runs, 1);
   tair_range_max_k = zeros(n_runs, 1);
   swd_max_wm2 = zeros(n_runs, 1);
   toa_max_wm2 = zeros(n_runs, 1);
   lwd_dev_wm2 = zeros(n_runs, 1);
   rh_range_pct = zeros(n_runs, 1);
   for k = 1:n_runs
      in_run = run_id == run_ids(k);
      % Whole calendar days are flagged so partial-day training windows
      % cannot slice into a compromised run.
      in_samples = ismember(gid, find(in_run));
      flagged(in_samples) = true;

      % Implicated channels: tair always (the core condition), swd/lwd
      % when their corroboration fired on any run day, rh when the run
      % holds it pinned (evidence annotation, never a flag condition).
      implicated = "tair";
      if any(cond_swd(in_run))
         implicated = implicated + ",swd";
      end
      if any(cond_lwd(in_run))
         implicated = implicated + ",lwd";
      end
      rh_range = NaN;
      if any(varnames == "rh")
         rh_run = met.rh(in_samples);
         rh_range = max(rh_run, [], 'omitnan') ...
            - min(rh_run, [], 'omitnan');
         if rh_range < kwargs.max_rh_range_pct
            implicated = implicated + ",rh";
         end
      end

      % Evidence statistics summarize the run for the report and the
      % upstream QC escalation without shipping raw samples.
      run_times = times(in_samples);
      start_time(k) = run_times(1);
      end_time(k) = run_times(end);
      n_days(k) = nnz(in_run);
      channels(k) = implicated;
      tair_range_max_k(k) = max(tair_range(in_run));
      swd_max_wm2(k) = max(swd_max(in_run));
      toa_max_wm2(k) = max(toa_max(in_run));
      % NaN daily deviations (low-coverage days admitted via swd) drop
      % out of the run mean; an absent lwd channel stays all-NaN.
      lwd_dev_wm2(k) = mean(lwd_dev(in_run), 'omitnan');
      rh_range_pct(k) = rh_range;
   end
   findings = table(start_time, end_time, n_days, channels, ...
      tair_range_max_k, swd_max_wm2, toa_max_wm2, lwd_dev_wm2, ...
      rh_range_pct);
end

%% Local functions
function stat = dailyStat(gid, n_groups, x, fun)
   %DAILYSTAT Aggregate one channel per day, NaN where no finite samples.
   valid = isfinite(x);
   stat = accumarray(gid(valid), x(valid), [n_groups, 1], fun, NaN);
end
