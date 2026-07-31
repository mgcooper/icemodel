function census = gapCensus(series, kwargs)
   %GAPCENSUS Census contiguous missing runs in one role-contract series.
   %
   %  census = icemodel.forcing.reconstruct.gapCensus(series)
   %  census = icemodel.forcing.reconstruct.gapCensus(series, ...
   %     channels=["tair", "lwd"], latitude=67.0, longitude=-48.8)
   %
   % Role
   %  Family-generic gap census for the reconstruction harness (DesignSpec
   %  2026-07-23-promice-gap-filling-and-ktransect, bead icemodel-g1n.6).
   %  The input is any regular timetable under the engine's role contract —
   %  a PROMICE target, a K-transect donor, a GC-Net record — never a
   %  family-specific type. The census bounds statistics to the observed
   %  record (first..last finite core-channel sample), optionally restricts
   %  shortwave to daylight so diurnal screening does not masquerade as
   %  outage, and returns both the per-run table the synthetic-missingness
   %  sampler resamples and the per-channel summary the policy and report
   %  consume.
   %
   % Name-value
   %  channels : string vector. Channels to census; default every numeric
   %     variable in the series.
   %  edges_hours : numeric vector. Ascending run-duration bucket edges in
   %     wall-clock hours; the policy census uses [0 6 24 72 168 Inf].
   %  cap_hours : wall-clock interpolation cap used to separate tier-1
   %     fixable samples from samples needing a donor or proxy.
   %  cap_hours_by_channel : evidenced per-channel overrides from setopts.
   %  core_channels : string vector. Channels whose joint finite span
   %     defines the observed record for in-record bounding; an absent
   %     configured channel means no valid core record, and an empty
   %     selection disables bounding. Default ["tair", "rh", "psfc"].
   %  daylight_channels : string vector. Channels censused only where TOA
   %     irradiance reaches toa_dark_wm2 (default "swd").
   %  toa_dark_wm2 : central meaningful-sun threshold.
   %  latitude, longitude : scalar site coordinates for the daylight cut.
   %
   % Returns
   %  census : struct with fields
   %     summary : table, one row per channel — n_in_record, n_missing,
   %        pct_missing, n_runs, per-bucket run counts, longest_run_hours,
   %        samples_fixable_within_cap, samples_needing_donor_or_proxy.
   %     runs : table, one row per missing run — channel, start_time,
   %        end_time, duration_hours, bucket, season (the stratum axes the
   %        sampler draws from).
   %     record_start, record_end : datetime bounds of the observed record.
   %     edges_hours : the bucket edges actually used.
   %     cap_hours : the interpolation cap actually used.
   %
   % See also: icemodel.forcing.reconstruct.syntheticMissingness,
   %  icemodel.forcing.reconstruct.validationSplit

   arguments
      series timetable
      kwargs.channels (1, :) string = strings(1, 0)
      kwargs.edges_hours (1, :) double = ...
         icemodel.forcing.reconstruct.bucketEdges()
      kwargs.cap_hours (1, 1) double {mustBePositive, ...
         icemodel.forcing.reconstruct.mustBeCapHours(kwargs.cap_hours)} = ...
         icemodel.forcing.reconstruct.setopts().cap_hours
      kwargs.cap_hours_by_channel (1, 1) struct = ...
         icemodel.forcing.reconstruct.setopts().cap_hours_by_channel
      kwargs.core_channels (1, :) string = ...
         icemodel.forcing.reconstruct.setopts().core_channels
      kwargs.daylight_channels (1, :) string = "swd"
      kwargs.latitude (1, 1) double = NaN
      kwargs.longitude (1, 1) double = NaN
      kwargs.toa_dark_wm2 (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().toa_dark_wm2
   end

   times = series.Properties.RowTimes;
   validated_caps = icemodel.forcing.reconstruct.setopts( ...
      cap_hours=kwargs.cap_hours, ...
      cap_hours_by_channel=kwargs.cap_hours_by_channel);
   cap_hours_by_channel = validated_caps.cap_hours_by_channel;
   names = string(series.Properties.VariableNames);
   channels = kwargs.channels;
   if isempty(channels)
      % Default to every numeric channel so donor and target series census
      % identically without per-family lists.
      numeric = varfun(@isnumeric, series, 'OutputFormat', 'uniform');
      channels = names(numeric);
   end
   missing_channels = setdiff(channels, names);
   if ~isempty(missing_channels)
      error('icemodel:reconstruct:gapCensus:unknownChannel', ...
         'channel(s) not in series: %s', strjoin(missing_channels, ", "));
   end
   dt_hours = hours(median(diff(times)));

   % Bound every statistic to the observed record so pre-deployment and
   % post-retrieval spans never count as outage.
   core = unique(kwargs.core_channels, 'stable');
   in_record = true(height(series), 1);
   if ~isempty(core)
      missing_core = setdiff(core, names);
      if ~isempty(missing_core)
         error('icemodel:reconstruct:gapCensus:noFiniteCore', ...
            'configured core channel(s) absent: %s', ...
            strjoin(missing_core, ", "));
      end
      finite_core = true(height(series), 1);
      for channel = reshape(core, 1, [])
         finite_core = finite_core & isfinite(series.(channel));
      end
      first = find(finite_core, 1, 'first');
      last = find(finite_core, 1, 'last');
      if isempty(first)
         error('icemodel:reconstruct:gapCensus:noFiniteCore', ...
            'no finite samples in core channels %s', strjoin(core, ", "));
      end
      in_record = false(height(series), 1);
      in_record(first:last) = true;
   end
   record_start = times(find(in_record, 1, 'first'));
   record_end = times(find(in_record, 1, 'last'));

   % The meaningful-sun cut applies only when the caller supplies the site
   % point; dark shortwave carries no information and would double-count the
   % diurnal screening pattern as missing data.
   day = true(height(series), 1);
   has_point = isfinite(kwargs.latitude) && isfinite(kwargs.longitude);
   if has_point
      toa = icemodel.forcing.reconstruct.toaIrradiance( ...
         times, kwargs.latitude, kwargs.longitude);
      day = toa >= kwargs.toa_dark_wm2;
   end

   summary_rows = cell(numel(channels), 1);
   run_rows = cell(numel(channels), 1);
   for c = 1:numel(channels)
      channel = channels(c);
      scope = in_record;
      if ismember(channel, kwargs.daylight_channels) && has_point
         scope = in_record & day;
      end
      miss = ~isfinite(series.(channel)) & scope;

      % Contiguous missing runs from the padded difference of the mask.
      d = diff([0; miss(:); 0]);
      starts = find(d == 1);
      stops = find(d == -1) - 1;
      run_hours = (stops - starts + 1) * dt_hours;
       buckets = icemodel.forcing.reconstruct.gapDurationBucket( ...
          run_hours, kwargs.edges_hours);
       % Count the assigned right-closed bucket IDs themselves; histcounts
       % uses the opposite boundary convention at exact policy edges.
       valid_buckets = isfinite(buckets);
       counts = accumarray(buckets(valid_buckets), ...
          ones(nnz(valid_buckets), 1), ...
          [numel(kwargs.edges_hours) - 1, 1]).';

      % Samples the configured interior cap can fix versus samples needing
      % a donor or proxy tier (sized by the POLICY B3 tier-1 cap).
      channel_cap = kwargs.cap_hours;
      if isfield(cap_hours_by_channel, channel)
         channel_cap = cap_hours_by_channel.(channel);
      end
      fixable = sum(run_hours(run_hours <= channel_cap)) / dt_hours;
      summary_rows{c} = [{char(channel), nnz(scope), nnz(miss), ...
         100 * nnz(miss) / max(1, nnz(scope)), numel(run_hours)}, ...
         num2cell(counts), {max([run_hours; 0]), fixable, ...
         nnz(miss) - fixable}];
      run_rows{c} = table( ...
         repmat(channel, numel(starts), 1), ...
         times(starts), times(stops), run_hours, buckets, ...
         icemodel.forcing.reconstruct.seasonOf(times(starts)), ...
         'VariableNames', {'channel', 'start_time', 'end_time', ...
         'duration_hours', 'bucket', 'season'});
   end

   bucket_names = compose("runs_bucket%d", 1:(numel(kwargs.edges_hours) - 1));
   summary = cell2table(vertcat(summary_rows{:}), 'VariableNames', ...
      [{'channel', 'n_in_record', 'n_missing', 'pct_missing', 'n_runs'}, ...
      cellstr(bucket_names), {'longest_run_hours', ...
      'samples_fixable_within_cap', 'samples_needing_donor_or_proxy'}]);
   census = struct('summary', summary, 'runs', vertcat(run_rows{:}), ...
      'record_start', record_start, 'record_end', record_end, ...
      'edges_hours', kwargs.edges_hours, 'cap_hours', kwargs.cap_hours, ...
      'cap_hours_by_channel', cap_hours_by_channel);
end
