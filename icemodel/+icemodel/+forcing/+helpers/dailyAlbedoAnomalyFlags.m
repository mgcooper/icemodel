function [row_flags, report] = dailyAlbedoAnomalyFlags(Time, swd, swu)
   %DAILYALBEDOANOMALYFLAGS Flag transient reflected-shortwave collapses.
   %
   %  [row_flags, report] = ...
   %     icemodel.forcing.helpers.dailyAlbedoAnomalyFlags(Time, swd, swu)
   %
   % The detector compares the daily shortwave-energy ratio
   % `sum(swu)/sum(swd)` with robust baselines before and after each day.
   % A conservative seed must have complete temporal support, at least one
   % third of the expected irradiated samples, at least 2 kWh m-2 incident
   % shortwave, and at least 80% valid reflected-shortwave energy coverage.
   % Seed-containing episodes expand only to nearby, weaker anomalies. This
   % recovery requirement rejects persistent dark-ice seasons while retaining
   % short sensor-collapse episodes on otherwise bright snow and firn.
   %
   % `row_flags` preserves the input vector shape and identifies every timestamp on an
   % anomalous day. `report` contains compact counts/dates plus a daily
   % diagnostic timetable for audit and provenance use. The helper never edits
   % source radiation; callers decide which derived channels to mask.

   arguments
      Time datetime {mustBeVector}
      swd double {mustBeVector}
      swu double {mustBeVector}
   end

   % Reject mismatched inputs before any grouping can silently truncate a
   % source channel.
   input_size = size(Time);
   n_rows = numel(Time);
   if numel(swd) ~= n_rows || numel(swu) ~= n_rows
      error('icemodel:forcing:helpers:dailyAlbedoAnomalyFlags:sizeMismatch', ...
         'Time, swd, and swu must contain the same number of rows.')
   end

   row_flags = false(input_size);
   report = emptyReport();
   if n_rows < 2
      return
   end

   % A naive source axis is declared to be UTC; an already zoned axis is
   % converted to UTC. Reversed axes fail closed because sorting would hide a
   % broken native-grid contract.
   Time.TimeZone = 'UTC';
   Time = Time(:);
   swd = swd(:);
   swu = swu(:);
   original_steps_s = seconds(diff(Time));
   nonzero_steps_s = abs(original_steps_s( ...
      isfinite(original_steps_s) & original_steps_s ~= 0));
   if isempty(nonzero_steps_s)
      return
   end
   cadence_s = median(nonzero_steps_s, 'omitnan');
   tolerance_s = max(1e-6, cadence_s * 1e-9);
   if any(original_steps_s < -tolerance_s)
      return
   end
   [Time, order] = sort(Time);
   swd = swd(order);
   swu = swu(order);

   % Infer a day-dividing dense cadence and its repeated within-day phase.
   % Duplicate, missing, and off-grid slots are rejected day by day below.
   expected_samples = 86400 / cadence_s;
   if ~isfinite(cadence_s) || cadence_s <= 0 || cadence_s >= 86400 ...
         || abs(expected_samples - round(expected_samples)) > tolerance_s
      return
   end
   expected_samples = round(expected_samples);
   if n_rows < expected_samples
      return
   end
   cadence_hours = cadence_s / 3600;
   phase_s = dominantDailyPhase(Time, cadence_s, tolerance_s);

   % Reduce raw shortwave pairs by UTC day. The measured downwelling energy
   % supplies the physically relevant solar-angle weighting.
   day_key = dateshift(Time, 'start', 'day');
   [group, days_utc] = findgroups(day_key);
   ratio = swu ./ swd;
   irradiated = isfinite(swd) & swd >= 10;
   valid_pair = irradiated & isfinite(swu) & swu > 0 ...
      & ratio > 0 & ratio <= 1.2;
   row_count = splitapply(@numel, swd, group);
   pair_count = splitapply(@(keep) nnz(keep), valid_pair, group);
   total_swd = splitapply(@selectedSum, swd, irradiated, group);
   valid_swd = splitapply(@selectedSum, swd, valid_pair, group);
   valid_swu = splitapply(@selectedSum, swu, valid_pair, group);
   alpha = valid_swu ./ valid_swd;
   incident_Whm2 = valid_swd .* cadence_hours;
   energy_coverage = valid_swd ./ total_swd;
   complete_grid = exactDailyGrid(Time, group, numel(days_utc), ...
      expected_samples, cadence_s, phase_s, tolerance_s, row_count);
   alpha(~complete_grid | pair_count == 0 | valid_swd <= 0) = NaN;

   % Baselines use well-supported days but a lower energy floor than seeds so
   % shoulder-season context remains available around a bright anomaly.
   minimum_pairs = max(8, ceil(expected_samples / 3));
   baseline_support = complete_grid & pair_count >= minimum_pairs ...
      & incident_Whm2 >= 500 & energy_coverage >= 0.8 ...
      & isfinite(alpha);
   [pre, post, robust_scale] = rollingBaselines( ...
      days_utc, alpha, baseline_support);

   % Conservative seeds provide specificity. Sensitive candidates are used
   % only to expand a component that already contains a conservative seed.
   seed_support = complete_grid & pair_count >= minimum_pairs ...
      & incident_Whm2 >= 2000 & energy_coverage >= 0.8;
   sensitive_support = complete_grid & pair_count >= minimum_pairs ...
      & incident_Whm2 >= 500 & energy_coverage >= 0.8;
   seed = anomalyDays(alpha, pre, post, robust_scale, seed_support, ...
      0.70, 0.20, 6);
   sensitive = anomalyDays(alpha, pre, post, robust_scale, ...
      sensitive_support, 0.80, 0.10, 4);
   episode = seedContainingEpisodes(days_utc, seed, sensitive, 2);

   % Map daily episode flags back through the sort permutation without changing
   % any raw shortwave sample.
   sorted_flags = ismember(day_key, days_utc(episode));
   restored_flags = false(n_rows, 1);
   restored_flags(order) = sorted_flags;
   row_flags = reshape(restored_flags, input_size);

   % Return compact provenance fields and the complete daily evidence table so
   % builders and artifact audits share one implementation.
   diagnostics = timetable(alpha, pair_count, incident_Whm2, ...
      energy_coverage, complete_grid, baseline_support, pre, post, ...
      robust_scale, seed, sensitive, episode, 'RowTimes', days_utc, ...
      'VariableNames', {'alpha', 'pair_count', 'incident_Whm2', ...
      'energy_coverage', 'complete_grid', 'baseline_support', 'pre', ...
      'post', 'robust_scale', 'seed', 'sensitive', 'episode'});
   diagnostics.Properties.DimensionNames{1} = 'Time';
   report = struct( ...
      'policy', policyText(), ...
      'seed_day_count', nnz(seed), ...
      'episode_day_count', nnz(episode), ...
      'flagged_row_count', nnz(row_flags), ...
      'seed_dates', string(days_utc(seed), 'yyyy-MM-dd'), ...
      'episode_dates', string(days_utc(episode), 'yyyy-MM-dd'), ...
      'diagnostics', diagnostics);
end

function value = selectedSum(values, keep)
   %SELECTEDSUM Sum values only where the paired logical selector is true.
   value = sum(values(keep), 'omitmissing');
end

function phase_s = dominantDailyPhase(Time, cadence_s, tolerance_s)
   %DOMINANTDAILYPHASE Infer the repeated within-day native-grid offset.
   day_start = dateshift(Time, 'start', 'day');
   phase = mod(seconds(Time - day_start), cadence_s);
   phase(abs(phase - cadence_s) <= tolerance_s) = 0;

   % Quantize only at the numerical comparison tolerance so one off-grid sample
   % cannot move the dominant native phase.
   phase = round(phase ./ tolerance_s) .* tolerance_s;
   phase_s = mode(phase);
end

function complete = exactDailyGrid( ...
      Time, day_index, n_days, n_slots, cadence_s, phase_s, tolerance_s, ...
      row_count)
   %EXACTDAILYGRID Require each expected phase-aware native slot exactly once.
   day_start = dateshift(Time, 'start', 'day');
   offset_s = seconds(Time - day_start);
   slot = round((offset_s - phase_s) ./ cadence_s);
   expected_offset_s = phase_s + slot .* cadence_s;
   on_grid = slot >= 0 & slot < n_slots ...
      & abs(offset_s - expected_offset_s) <= tolerance_s;

   % Count every accepted day/slot pair. A duplicate plus a missing timestamp
   % cannot masquerade as a complete day merely because its row count matches.
   slot_index = (day_index(on_grid) - 1) .* n_slots + slot(on_grid) + 1;
   slot_count = zeros(n_days * n_slots, 1);
   if ~isempty(slot_index)
      slot_count = accumarray(slot_index, ones(size(slot_index)), ...
         [n_days * n_slots, 1]);
   end
   unique_grid = all(reshape(slot_count, n_slots, n_days) == 1, 1)';

   % Duplicate or zero steps invalidate their UTC day even when all other rows
   % happen to occupy valid native slots.
   bad_step = find(diff(Time) <= seconds(0)) + 1;
   bad_order = false(n_days, 1);
   if ~isempty(bad_step)
      bad_order = accumarray(day_index(bad_step), ...
         ones(numel(bad_step), 1), [n_days, 1]) > 0;
   end
   complete = row_count == n_slots & unique_grid & ~bad_order;
end

function [pre, post, scale] = rollingBaselines(Time, alpha, eligible)
   %ROLLINGBASELINES Compute separated robust context before and after each day.
   n_days = numel(Time);
   pre = nan(n_days, 1);
   post = nan(n_days, 1);
   scale = nan(n_days, 1);
   for k = 1:n_days
      before = Time >= Time(k) - caldays(30) ...
         & Time <= Time(k) - caldays(4) & eligible;
      after = Time >= Time(k) + caldays(4) ...
         & Time <= Time(k) + caldays(30) & eligible;
      if nnz(before) < 5 || nnz(after) < 5
         continue
      end

      % Use the larger side-specific MAD so one unusually quiet side cannot
      % make the detector hypersensitive.
      before_values = alpha(before);
      after_values = alpha(after);
      pre(k) = median(before_values, 'omitmissing');
      post(k) = median(after_values, 'omitmissing');
      pre_scale = 1.4826 ...
         * median(abs(before_values - pre(k)), 'omitmissing');
      post_scale = 1.4826 ...
         * median(abs(after_values - post(k)), 'omitmissing');
      scale(k) = max(pre_scale, post_scale);
   end
end

function flags = anomalyDays( ...
      alpha, pre, post, scale, support, cap, minimum_drop, sigma)
   %ANOMALYDAYS Apply one two-sided robust transient threshold.
   threshold = max(minimum_drop, sigma .* scale);
   flags = support & alpha < cap ...
      & pre - alpha >= threshold ...
      & post - alpha >= threshold;
end

function episode = seedContainingEpisodes(Time, seed, sensitive, max_gap_days)
   %SEEDCONTAININGEPISODES Retain sensitive components containing a seed.
   candidate = seed | sensitive;
   episode = false(size(candidate));
   indices = find(candidate);
   if isempty(indices)
      return
   end

   % Candidate days separated by at most the settled gap belong to one episode;
   % only a component with a conservative seed survives.
   first = 1;
   for k = 2:numel(indices) + 1
      connected = k <= numel(indices) ...
         && days(Time(indices(k)) - Time(indices(k - 1))) <= max_gap_days;
      if connected
         continue
      end
      component = indices(first:k - 1);
      if any(seed(component))
         episode(component) = true;
      end
      first = k;
   end
end

function report = emptyReport()
   %EMPTYREPORT Return a stable no-evidence report for short or invalid axes.
   diagnostics = timetable();
   diagnostics.Properties.DimensionNames{1} = 'Time';
   report = struct( ...
      'policy', policyText(), ...
      'seed_day_count', 0, ...
      'episode_day_count', 0, ...
      'flagged_row_count', 0, ...
      'seed_dates', strings(0, 1), ...
      'episode_dates', strings(0, 1), ...
      'diagnostics', diagnostics);
end

function text = policyText()
   %POLICYTEXT Describe the settled detector in artifact-safe plain text.
   text = "daily sum(swu)/sum(swd) over swd >= 10 W m-2 valid pairs; " ...
      + "complete UTC grid, >= max(8, one-third daily samples), and >= 80% " ...
      + "valid incident-energy coverage; conservative seeds require >= 2 " ...
      + "kWh m-2, alpha < 0.70, and two-sided drop >= max(0.20, 6 MAD); " ...
      + "seed episodes expand across <= 2-day gaps only to >= 0.5 kWh m-2, " ...
      + "alpha < 0.80, two-sided drop >= max(0.10, 4 MAD) days";
end
