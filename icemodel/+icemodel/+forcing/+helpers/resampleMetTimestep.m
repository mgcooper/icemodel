function met = resampleMetTimestep(met, dt_out)
   %RESAMPLEMETTIMESTEP Resample model met without crossing source outages.
   %
   %  met = icemodel.forcing.helpers.resampleMetTimestep(met, "15m")
   %  met = icemodel.forcing.helpers.resampleMetTimestep(met, "")
   %
   % The 15-minute path treats every source row as the interval mean valid over
   % [t,t+source_dt), matching the model's forward integration contract. It
   % therefore holds each numeric value constant over its declared support;
   % interpolation would change the interval integral. Explicit non-finite
   % values and omitted timestamps remain unavailable, and the final source
   % interval is represented through its exclusive end. Timetable variable
   % attributes and UserData are preserved; UserData gains source-derived
   % support and missing-count provenance for artifact QA. Larger timestamp
   % steps must omit whole native intervals. A compact two-row source may prove
   % a common cadence up to one hour; a longer lone step is rejected because it
   % cannot distinguish native cadence from a gap. dt_out="" is an exact no-op.
   % Compact per-calendar-year summaries preserve exact source/gap/support facts
   % when writemet later divides the fully derived output into yearly artifacts.

   arguments
      met timetable
      dt_out (1, 1) string = "15m"
   end

   if dt_out == ""
      return
   end
   if dt_out ~= "15m"
      error('icemodel:forcing:resampleMetTimestep:unsupportedDt', ...
         'unsupported model-met timestep: %s', dt_out)
   end

   % A prior guarded builder resample remains authoritative when writemet sees
   % its already-15-minute output. Do not replace native-source provenance with
   % counts from the derived timetable.
   metadata = met.Properties.UserData;
   already_guarded = isstruct(metadata) ...
      && isfield(metadata, 'met_resample_policy');
   uniform_15m = height(met) <= 1 ...
      || all(diff(met.Time) == minutes(15));
   if already_guarded
      if uniform_15m
         return
      end
      error('icemodel:forcing:resampleMetTimestep:guardedInputInvalid', ...
         'guarded met has a missing or irregular derived 15-minute row')
   end
   if height(met) > 1 && all(diff(met.Time) == minutes(15))
      met = recordResampleProvenance(met, met, ...
         "native_15m_unchanged", 900, 0, missingCounts(met), ...
         met.Time(end) + minutes(15));
      return
   end

   % One row has no inferable source duration. Preserve that exact sample rather
   % than inventing support; multi-row writers enforce their own requirements.
   if height(met) < 2
      met = recordResampleProvenance(met, met, ...
         "single_sample_unchanged", NaN, 0, missingCounts(met), NaT);
      return
   end

   % A single compact step up to one hour is a supported application cadence.
   % Longer two-row spans cannot distinguish native cadence from omitted rows.
   if height(met) == 2 && seconds(met.Time(2) - met.Time(1)) > 3600
      error('icemodel:forcing:resampleMetTimestep:ambiguousSourceCadence', ...
         ['cannot infer source cadence from a two-row step longer than one ' ...
         'hour; use native output or provide at least three timestamps'])
   end

   % Preserve every source interval while creating the uniform model-met
   % cadence. The support mask independently records the minimum missing output
   % count expected from native values and omitted timestamps.
   source = met;
   source_steps = seconds(diff(source.Time));
   % Omitted source rows appear only as larger steps, so the shortest positive
   % step is the native cadence when the full source is available.
   source_cadence_s = min(source_steps(source_steps > 0), [], 'omitnan');
   validateSourceCadence(source_steps, source_cadence_s);
   support_end = source.Time(end) + seconds(source_cadence_s);
   new_time = (source.Time(1):minutes(15):support_end - minutes(15))';
   [met, source_cadence_s, source_gap_count, expected_missing] = ...
      holdFiniteIntervals(source, new_time, source_cadence_s);
   met = recordResampleProvenance(met, source, ...
      "interval_start_zero_order_hold", source_cadence_s, ...
      source_gap_count, expected_missing, support_end);
end

function validateSourceCadence(source_steps, cadence_s)
   %VALIDATESOURCECADENCE Require 15-minute-aligned forward interval support.
   if ~isfinite(cadence_s) || cadence_s < 900 ...
         || abs(mod(cadence_s, 900)) > 1e-6
      error('icemodel:forcing:resampleMetTimestep:badSourceCadence', ...
         'source cadence must be a finite multiple of 15 minutes')
   end

   % Short or non-grid-aligned steps are ambiguous; larger steps are omitted
   % source intervals only when they span a whole number of native intervals.
   regular = abs(source_steps - cadence_s) <= 1e-6;
   native_multiples = source_steps ./ cadence_s;
   omitted = source_steps > cadence_s + 1e-6 ...
      & abs(native_multiples - round(native_multiples)) <= 1e-9;
   if any(~regular & ~omitted)
      error('icemodel:forcing:resampleMetTimestep:irregularSourceTime', ...
         'source times must follow or omit whole multiples of the native cadence')
   end
end

function [met, cadence_s, gap_count, expected_missing] = ...
      holdFiniteIntervals(source, new_time, cadence_s)
   %HOLDFINITEINTERVALS Repeat each row only over its forward source support.
   source_steps = seconds(diff(source.Time));
   gap_intervals = find(source_steps > 1.5 * cadence_s);
   gap_count = numel(gap_intervals);

   % Retime establishes the output shape and preserves table variable
   % attributes. Numeric channels are replaced below by support-aware holds.
   met = retime(source, new_time, 'fillwithmissing');
   rows_per_interval = round(cadence_s / 900);
   source_rows = round(seconds(source.Time - new_time(1)) / 900) + 1;

   expected_missing = struct();
   names = string(source.Properties.VariableNames);
   for name = reshape(names, 1, [])
      values = source.(char(name));
      if ~isnumeric(values)
         continue
      end

      % Initialize missing so omitted timestamp intervals stay unavailable. A
      % source row then fills exactly source_dt/15m rows, including NaN values.
      held = nan([numel(new_time), size(values, 2)]);
      for k = 1:height(source)
         rows = source_rows(k) + (0:rows_per_interval - 1);
         rows = rows(rows >= 1 & rows <= numel(new_time));
         if isempty(rows)
            continue
         end
         held(rows, :) = repmat(values(k, :), numel(rows), 1);
      end
      met.(char(name)) = held;
      expected_missing.(char(name)) = nnz(~isfinite(held));
   end
end

function counts = missingCounts(T)
   %MISSINGCOUNTS Count unavailable values in each numeric timetable variable.
   counts = struct();
   names = string(T.Properties.VariableNames);
   for name = reshape(names, 1, [])
      values = T.(char(name));
      if isnumeric(values)
         counts.(char(name)) = nnz(~isfinite(values));
      end
   end
end

function met = recordResampleProvenance(met, source, policy, cadence_s, ...
      gap_count, expected_missing, support_end)
   %RECORDRESAMPLEPROVENANCE Persist source facts for read-only artifact QA.
   metadata = source.Properties.UserData;
   if isempty(metadata)
      metadata = struct();
   elseif ~isstruct(metadata)
      error('icemodel:forcing:resampleMetTimestep:badMetadata', ...
         'met.Properties.UserData must be empty or a metadata struct')
   end

   metadata.met_resample_policy = policy;
   metadata.met_resample_source_row_count = height(source);
   metadata.met_resample_source_cadence_seconds = cadence_s;
   metadata.met_resample_source_time_gap_count = gap_count;
   metadata.met_resample_source_missing_counts = missingCounts(source);
   metadata.met_resample_expected_missing_counts = expected_missing;
   metadata.met_resample_time_semantics = "interval_start";
   if height(met) > 0
      metadata.met_resample_support_start_inclusive = met.Time(1);
   else
      metadata.met_resample_support_start_inclusive = NaT;
   end
   metadata.met_resample_support_end_exclusive = support_end;
   gap_intervals = sourceGapIntervals(source, cadence_s);
   metadata.met_resample_source_gap_intervals = gap_intervals;
   metadata.met_resample_yearly_summaries = yearlyResampleSummaries( ...
      met, source, cadence_s, gap_intervals);
   met.Properties.UserData = metadata;
end

function gaps = sourceGapIntervals(source, cadence_s)
   %SOURCEGAPINTERVALS Record omitted native support without inventing rows.
   template = struct('start', NaT, 'end', NaT);
   if height(source) < 2 || ~isfinite(cadence_s)
      gaps = repmat(template, 0, 1);
      return
   end
   source_steps = seconds(diff(source.Time));
   gap_index = find(source_steps > cadence_s + 1e-6);
   gaps = repmat(template, numel(gap_index), 1);
   for n = 1:numel(gap_index)
      k = gap_index(n);
      gaps(n).start = source.Time(k) + seconds(cadence_s);
      gaps(n).end = source.Time(k + 1);
   end
end

function summaries = yearlyResampleSummaries(met, source, cadence_s, gaps)
   %YEARLYRESAMPLESUMMARIES Record exact lineage for later calendar slicing.
   template = struct('year', NaN, 'source_row_count', 0, ...
      'source_time_gap_count', 0, 'source_missing_counts', struct(), ...
      'source_gap_intervals', repmat(struct('start', NaT, 'end', NaT), 0, 1), ...
      'expected_missing_counts', struct(), ...
      'support_start_inclusive', NaT, 'support_end_exclusive', NaT);
   if height(met) == 0 || height(source) == 0 || ~isfinite(cadence_s)
      summaries = repmat(template, 0, 1);
      return
   end

   years_present = unique(year(met.Time));
   summaries = repmat(template, numel(years_present), 1);
   source_support_end = source.Time + seconds(cadence_s);

   for n = 1:numel(years_present)
      yyyy = years_present(n);
      output = met(year(met.Time) == yyyy, :);
      support_start = output.Time(1);
      support_end = output.Time(end) + minutes(15);

      % A source row belongs to every yearly slice touched by its forward
      % interval; a cross-midnight interval can therefore contribute to two.
      source_mask = source.Time < support_end ...
         & source_support_end > support_start;
      gap_overlaps = false(numel(gaps), 1);
      for k = 1:numel(gaps)
         gap_overlaps(k) = gaps(k).start < support_end ...
            && gaps(k).end > support_start;
      end

      summaries(n).year = yyyy;
      summaries(n).source_row_count = nnz(source_mask);
      summaries(n).source_time_gap_count = nnz(gap_overlaps);
      summaries(n).source_missing_counts = missingCounts(source(source_mask, :));
      summaries(n).source_gap_intervals = gaps(gap_overlaps);
      summaries(n).expected_missing_counts = missingCounts(output);
      summaries(n).support_start_inclusive = support_start;
      summaries(n).support_end_exclusive = support_end;
   end
end
