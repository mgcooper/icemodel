function [observed, modeled, aligned] = alignObservationSeries( ...
      observation, candidate)
   %ALIGNOBSERVATIONSERIES Align one observation/model series by its support.
   %
   %  [observed, modeled, aligned] = ...
   %     icemodel.verification.helpers.alignObservationSeries( ...
   %     observation, candidate)
   %
   % OBSERVATION and CANDIDATE are normalized series structs with `values`,
   % `axis`, `axis_kind`, `variable`, and `units` fields. Supported axes are
   % exact timestamps, profile depth, and observation intervals paired with
   % timestamped model samples.

   arguments
      observation (1, 1) struct
      candidate (1, 1) struct
   end

   % Dispatch once on the declared support so every verification consumer uses
   % the same interval, depth, and timestamp semantics.
   if observation.axis_kind == "time" && candidate.axis_kind == "time"
      [observed, modeled, aligned] = alignByTime(observation, candidate);
   elseif observation.axis_kind == "depth" ...
         && candidate.axis_kind == "depth"
      [observed, modeled, aligned] = alignByDepth(observation, candidate);
   elseif observation.axis_kind == "interval" ...
         && candidate.axis_kind == "time"
      [observed, modeled, aligned] = alignIntervalsToTimes( ...
         observation, candidate);
   else
      [observed, modeled, aligned] = deal([], [], table());
   end
end

function [observed, modeled, aligned] = alignIntervalsToTimes( ...
      observation, candidate)
   %ALIGNINTERVALSTOTIMES Integrate complete model support over observations.
   intervals = observation.axis;
   observed_values = double(observation.values(:));
   sample_times = candidate.axis(:);
   sample_values = double(candidate.values(:));
   accumulate_rate = shouldAccumulateInterval(observation, candidate);
   sample_hours = representativeStepHours(sample_times);
   if accumulate_rate
      % A rate sample cannot represent more than one day in the current
      % verification contract, even when sparse source postings are supplied.
      sample_hours = min(sample_hours, 24);
   end
   finite_samples = ~isnat(sample_times) & isfinite(sample_values);

   % Preallocate at the observation count, then compact after rejecting missing
   % values or incomplete interval coverage.
   n_intervals = numel(observed_values);
   keep_rows = false(n_intervals, 1);
   modeled_values = NaN(n_intervals, 1);
   for n = 1:n_intervals
      if accumulate_rate
         keep = sample_times >= intervals(n, 1) ...
            & sample_times < intervals(n, 2) ...
            & finite_samples;
      else
         keep = sample_times >= intervals(n, 1) ...
            & sample_times <= intervals(n, 2) ...
            & finite_samples;
      end
      if ~isfinite(observed_values(n)) || ~any(keep)
         continue
      end
      if numel(unique(sample_times(keep))) ~= nnz(keep)
         % Duplicate model postings would double-count an accumulated flux.
         % Reject the interval rather than choosing or averaging a state that
         % the candidate contract requires to be unique.
         continue
      end
      if accumulate_rate && ~coversInterval(sample_times(keep), ...
            intervals(n, 1), intervals(n, 2), sample_hours)
         continue
      end
      if accumulate_rate
         modeled_values(n) = sum(sample_values(keep)) * sample_hours;
      else
         modeled_values(n) = mean(sample_values(keep));
      end
      keep_rows(n) = true;
   end

   observed = observed_values(keep_rows);
   modeled = modeled_values(keep_rows);
   source_row = find(keep_rows);
   aligned = table(intervals(keep_rows, 1), intervals(keep_rows, 2), ...
      source_row, observed, modeled, 'VariableNames', ...
      {'start_date', 'end_date', 'source_row', 'target', 'candidate'});
end

function tf = coversInterval(sample_times, start_date, end_date, sample_hours)
   %COVERSINTERVAL True when finite samples continuously cover one interval.
   if isempty(sample_times)
      tf = false;
      return
   end

   % The caller has already rejected duplicate postings. Sorting here preserves
   % the ability to detect omitted or non-finite internal samples as gaps.
   sample_times = unique(sort(sample_times(:)));
   coverage_end = end_date - hours(sample_hours);
   if coverage_end < start_date
      coverage_end = start_date;
   end
   if min(sample_times) > start_date || max(sample_times) < coverage_end
      tf = false;
      return
   end

   % A single rate sample covers only a sub-step observation; longer intervals
   % require a continuous sequence at the representative cadence.
   if isscalar(sample_times)
      tf = hours(end_date - start_date) <= sample_hours / 2;
   else
      gaps = hours(diff(sample_times));
      tf = all(gaps <= sample_hours * 1.5);
   end
end

function tf = shouldAccumulateInterval(observation, candidate)
   %SHOULDACCUMULATEINTERVAL True for accumulated-SMB observation semantics.
   variables = [string(observation.variable), string(candidate.variable)];
   units = lower([string(observation.units), string(candidate.units)]);
   tf = any(variables == "smb") || any(contains(units, "/h")) ...
      || any(contains(units, "per hour"));
end

function sample_hours = representativeStepHours(sample_times)
   %REPRESENTATIVESTEPHOURS Return the positive median model timestep.
   valid_times = sort(sample_times(~isnat(sample_times)));
   if numel(valid_times) < 2
      sample_hours = 1;
      return
   end
   steps = hours(diff(valid_times));
   steps = steps(isfinite(steps) & steps > 0);
   if isempty(steps)
      sample_hours = 1;
   else
      sample_hours = median(steps);
   end
end

function [observed, modeled, aligned] = alignByTime(observation, candidate)
   %ALIGNBYTIME Synchronize comparable series on exact common timestamps.
   observation_tt = timetable(observation.axis(:), observation.values(:), ...
      'VariableNames', {'target'});
   candidate_tt = timetable(candidate.axis(:), candidate.values(:), ...
      'VariableNames', {'candidate'});
   sync_tt = synchronize(observation_tt, candidate_tt, ...
      'intersection', 'first');
   valid = isfinite(sync_tt.target) & isfinite(sync_tt.candidate);
   sync_tt = sync_tt(valid, :);
   observed = sync_tt.target;
   modeled = sync_tt.candidate;
   aligned = timetable2table(sync_tt);
end

function [observed, modeled, aligned] = alignByDepth(observation, candidate)
   %ALIGNBYDEPTH Interpolate candidate values only inside common depth support.
   observed_depth = double(observation.axis(:));
   candidate_depth = double(candidate.axis(:));
   observed_values = double(observation.values(:));
   candidate_values = double(candidate.values(:));
   observed_valid = isfinite(observed_depth) & isfinite(observed_values);
   candidate_valid = isfinite(candidate_depth) & isfinite(candidate_values);
   observed_depth = observed_depth(observed_valid);
   observed_values = observed_values(observed_valid);
   candidate_depth = candidate_depth(candidate_valid);
   candidate_values = candidate_values(candidate_valid);
   if isempty(observed_depth) || isempty(candidate_depth)
      [observed, modeled, aligned] = deal([], [], table());
      return
   end

   % Candidate depths must be unique before interpolation. Stable selection
   % preserves the source order for repeated legacy coordinates.
   [candidate_depth, ia] = unique(candidate_depth, 'stable');
   candidate_values = candidate_values(ia);
   keep = observed_depth >= min(candidate_depth) ...
      & observed_depth <= max(candidate_depth);
   observed_depth = observed_depth(keep);
   observed = observed_values(keep);
   if isscalar(candidate_depth)
      modeled = repmat(candidate_values, size(observed));
   else
      modeled = interp1(candidate_depth, candidate_values, ...
         observed_depth, 'linear');
   end

   % Non-finite interpolation results are not valid matched observations.
   valid = isfinite(observed) & isfinite(modeled);
   observed_depth = observed_depth(valid);
   observed = observed(valid);
   modeled = modeled(valid);
   aligned = table(observed_depth(:), observed(:), modeled(:), ...
      'VariableNames', {'depth', 'target', 'candidate'});
end
