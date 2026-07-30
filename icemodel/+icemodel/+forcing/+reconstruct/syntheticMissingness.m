function draws = syntheticMissingness(series, channel, runs, kwargs)
   %SYNTHETICMISSINGNESS Draw blocked synthetic gaps into observed segments.
   %
   %  draws = icemodel.forcing.reconstruct.syntheticMissingness( ...
   %     series, "lwd", census.runs, years=split.years_selection, seed=7)
   %
   % Role
   %  Blocked synthetic-missingness sampler for the held-out validation
   %  protocol (gap-fill policy). Gap DURATIONS are
   %  resampled from the real census run-length distribution of the target
   %  channel (optionally per stratum) so the synthetic test matches the
   %  outage structure the engine will face; whole gaps are inserted only
   %  into fully observed spans of the allowed years, never pointwise, so
   %  autocorrelation is honest. Inserted gaps do not overlap each other and
   %  keep an observed context margin on both sides, which the boundary-jump
   %  metric needs.
   %
   % Name-value
   %  years : double vector. Calendar years eligible for insertion — the
   %     split's selection years during method admission, evaluation years
   %     for final grading. Required.
   %  seed : nonnegative integer. Required; the draw is deterministic given
   %     (runs, years, seed).
   %  n_gaps : number of synthetic gaps to attempt (default 25). Fewer are
   %     returned when the observed spans cannot host more without overlap;
   %     the shortfall is reported rather than silently absorbed.
   %  bucket : optional census bucket index filter for the duration pool
   %     (stratified draws pass one bucket at a time).
   %  season : optional season filter ("DJF"/"MAM"/"JJA"/"SON") for the
   %     duration pool.
   %  context_hours : observed margin required on each side of an inserted
   %     gap (default 24).
   %  latitude, longitude, toa_dark_wm2 : solar geometry and the central
   %     meaningful-sun threshold required for SWD draws.
   %
   % Returns
   %  draws : struct with fields
   %     mask : logical column, true where the channel is synthetically
   %        masked (always a subset of currently finite samples).
   %     gaps : table — start_time, end_time, duration_hours, bucket,
   %        season — one row per inserted gap.
   %     requested, inserted : attempted and achieved gap counts.
   %
   % See also: icemodel.forcing.reconstruct.gapCensus,
   %  icemodel.forcing.reconstruct.validationSplit

   arguments
      series timetable
      channel (1, 1) string
      runs table
      kwargs.years (1, :) double {mustBeInteger, mustBeNonempty}
      kwargs.seed (1, 1) double {mustBeInteger, mustBeNonnegative}
      kwargs.n_gaps (1, 1) double {mustBeInteger, mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().synthetic_n_gaps
      kwargs.bucket (1, :) double = []
      kwargs.season (1, :) string = strings(1, 0)
      kwargs.context_hours (1, 1) double {mustBeNonnegative} = ...
         icemodel.forcing.reconstruct.setopts().synthetic_context_hours
      kwargs.latitude (1, 1) double = NaN
      kwargs.longitude (1, 1) double = NaN
      kwargs.toa_dark_wm2 (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().toa_dark_wm2
   end

   times = series.Properties.RowTimes;
   dt_hours = hours(median(diff(times)));

   % The duration pool is the real run-length distribution of this channel,
   % optionally narrowed to one stratum. An empty pool is an error: drawing
   % durations from nothing would silently test the wrong regime.
   pool = runs(runs.channel == channel, :);
   if ~isempty(kwargs.bucket)
      pool = pool(ismember(pool.bucket, kwargs.bucket), :);
   end
   if ~isempty(kwargs.season)
      pool = pool(ismember(pool.season, kwargs.season), :);
   end
   if isempty(pool)
      error('icemodel:reconstruct:syntheticMissingness:emptyDurationPool', ...
         'no census runs for channel %s in the requested stratum', channel);
   end

   % Insertions live only in the requested years and season, on finite
   % samples, away from the series edges by the context margin.
    x = series.(channel);
    eligible = isfinite(x) & ismember(year(times), kwargs.years);
    if channel == "swd"
       if ~(isfinite(kwargs.latitude) && isfinite(kwargs.longitude))
          error(['icemodel:reconstruct:syntheticMissingness:' ...
             'missingSolarGeometry'], ...
             'SWD validation draws require latitude and longitude');
       end
       toa = icemodel.forcing.reconstruct.toaIrradiance( ...
          times, kwargs.latitude, kwargs.longitude);
       eligible = eligible & toa >= kwargs.toa_dark_wm2;
    end
   if ~isempty(kwargs.season)
      eligible = eligible & ismember( ...
         icemodel.forcing.reconstruct.seasonOf(times), kwargs.season);
   end
   context = max(1, round(kwargs.context_hours / dt_hours));

   stream = RandStream('mt19937ar', 'Seed', kwargs.seed);
   duration_hours = pool.duration_hours( ...
      randi(stream, height(pool), kwargs.n_gaps, 1));

   mask = false(height(series), 1);
   gap_rows = cell(kwargs.n_gaps, 1);
   inserted = 0;
   for g = 1:kwargs.n_gaps
      run_len = max(1, round(duration_hours(g) / dt_hours));
      % A viable start needs the whole gap plus both context margins to be
      % eligible, currently unmasked, and contiguous. movsum's shrinking
      % endpoint windows can never reach the full-window sum, so series-edge
      % starts are excluded without any explicit (and crashable) index math.
      window = run_len + 2 * context;
      viable = movsum(eligible & ~mask, [0, window - 1]) == window;
      candidates = find(viable);
      if isempty(candidates)
         continue
      end
      start = candidates(randi(stream, numel(candidates))) + context;
      stop = start + run_len - 1;
      mask(start:stop) = true;
      inserted = inserted + 1;
      gap_rows{inserted} = table(times(start), times(stop), ...
         run_len * dt_hours, ...
         icemodel.forcing.reconstruct.gapDurationBucket( ...
         run_len * dt_hours), ...
         icemodel.forcing.reconstruct.seasonOf(times(start)), ...
         'VariableNames', ...
         {'start_time', 'end_time', 'duration_hours', 'bucket', 'season'});
   end

   gaps = vertcat(gap_rows{1:inserted});
   if inserted == 0
      gaps = table('Size', [0 5], 'VariableTypes', ...
         {'datetime', 'datetime', 'double', 'double', 'string'}, ...
         'VariableNames', {'start_time', 'end_time', 'duration_hours', ...
         'bucket', 'season'});
   else
      % Return gaps in TIME order: consumers extract masked samples with
      % series.(channel)(mask), which is time-ordered, and validationMetrics
      % maps samples to gaps by walking the table sequentially — insertion
      % order would silently misalign that mapping.
      gaps = sortrows(gaps, 'start_time');
   end
   draws = struct('mask', mask, 'gaps', gaps, ...
      'requested', kwargs.n_gaps, 'inserted', inserted);
end
