function [hourly, bin_start, bin_end] = retimeHourlyFixedStep(TT)
   %RETIMEHOURLYFIXEDSTEP Aggregate fixed-step data to hourly values.
   %
   %  TT = icemodel.retimeHourlyFixedStep(TT)
   %  [TT, bin_start, bin_end] = icemodel.retimeHourlyFixedStep(TT)
   %
   % This helper is used by icemodel.postprocess where model output is
   % known to be 15-minute data (opts.dt == 900). Aligned complete hours use a
   % fixed-array path for compatibility with code generation. Interpreted
   % MATLAB falls back to native timetable retiming for partial, unaligned, or
   % irregular bins. BIN_START and BIN_END identify the inclusive raw-sample
   % bounds for each retained output row; zero bounds denote an empty native
   % timetable bin.
   %
   % The fixed path is written for code generation. Every construct that Coder
   % rejects is inside a coder.target('MATLAB') branch, which the code
   % generator removes. Nothing in this repository compiles this function. The
   % directive below records that intent and lets the analyzer check it.
   %
   % See also: icemodel.postprocess, icemodel.isIncrementChannel
   %
   %#codegen

   % Return early on empty inputs so callers can stay simple.
   bin_start = zeros(0, 1);
   bin_end = zeros(0, 1);
   if isempty(TT)
      hourly = TT;
      return
   end

   % Use the fixed four-sample branch whenever the full 15-minute spacing is
   % known to hold. Generated code rejects all other grids to prevent
   % unsupported native timetable retiming from entering the compiled graph.
   if isFixedStepHourlyCompatible(TT.Properties.RowTimes)
      % Seed the output from the first sample in every block so timetable
      % properties are preserved per variable. This branch preserves the storage
      % class for floating-point variables only. The block mean widens an
      % integer or logical variable to double, and the fallback branch below
      % rejects those instead of widening them. Callers with logical channels
      % must cast before retiming, as icemodel.postprocess does for
      % Tsfc_converged and Tice_converged.
      n_samples = height(TT);
      n_hours = n_samples / 4;
      use_codegen_leap = ~coder.target('MATLAB') && n_hours == 8784;
      if use_codegen_leap
         % Generated WRITEOUTPUT calls contain one complete year. Exclude the
         % 96 raw leap-day samples before aggregation so the timetable output
         % is created once with its final, codegen-constant 8760-row size.
         bin_start = transpose([1:4:5661, 5761:4:n_samples]);
         bin_end = bin_start + 3;
         hourly = TT([1:4:5661, 5761:4:end], :);
         n_output_hours = 8760;
      else
         bin_start = transpose(1:4:n_samples);
         bin_end = bin_start + 3;
         hourly = TT(1:4:end, :);
         n_output_hours = n_hours;
      end

      % Apply the aggregation class one variable at a time. This avoids
      % concatenating mixed single and double variables into one array.
      vars = TT.Properties.VariableNames;
      for n = 1:numel(vars)
         name = vars{n};
         method = aggregationMethod(name);
         if use_codegen_leap
            % The raw values must omit the same leap-day interval as the row
            % labels before they enter four-sample aggregation.
            values = TT.(name)([1:5664, 5761:end], :);
         else
            values = TT.(name);
         end
         hourly.(name) = ...
            aggregateFixedValues(values, n_output_hours, method);
      end
   elseif coder.target('MATLAB')
      [hourly, bin_start, bin_end] = aggregateNativeBins(TT);
   else
      error('icemodel:retimeHourlyFixedStep:unsupportedGrid', ...
         ['generated hourly retiming requires aligned complete ', ...
         '15-minute bins'])
   end

   % Interpreted MATLAB drops the synthetic Feb 29 rows here. The code
   % generation path above already excludes that interval.
   if coder.target('MATLAB')
      keep = ~(month(hourly.Properties.RowTimes) == 2 ...
         & day(hourly.Properties.RowTimes) == 29);
      hourly = hourly(keep, :);
      bin_start = bin_start(keep);
      bin_end = bin_end(keep);
   end
end

function tf = isFixedStepHourlyCompatible(time)
   %ISFIXEDSTEPHOURLYCOMPATIBLE Check the timetable for the fixed retime.

   % Compatibility check: four samples per hour, an hourly first label, and
   % exact quarter-hour spacing.
   tf = numel(time) >= 4 ...
      && mod(numel(time), 4) == 0 ...
      && minute(time(1)) == 0 ...
      && all(diff(time) == minutes(15));
   if coder.target('MATLAB')
      % Interpreted MATLAB can also reject a fractional-minute first sample;
      % SECOND does not support datetime inputs during code generation.
      tf = tf && second(time(1)) == 0;
   end
end

function returned = aggregateFixedValues(values, n_hours, method)
   %AGGREGATEFIXEDVALUES Reduce four-sample blocks using one method code.

   % Preserve trailing variable dimensions while making the four native rows
   % of each hour the leading reduction dimension.
   n_columns = size(values, 2);
   blocks = reshape(values, 4, n_hours, n_columns);

   % RETIME excludes missing numeric samples from reductions. Zero replacement
   % gives SUM its native all-missing result and leaves MEAN as zero divided by
   % zero when no finite samples exist.
   missing = isnan(blocks);
   finite_blocks = blocks;
   finite_blocks(missing) = 0;
   switch method
      case 1
         reduced = sum(finite_blocks, 1) ./ sum(~missing, 1);
      case 2
         reduced = sum(finite_blocks, 1);
      case 3
         % FIRSTVALUE selects the first nonmissing sample in each bin.
         reduced = blocks(1, :, :);
         for n = 2:4
            use_value = isnan(reduced) & ~isnan(blocks(n, :, :));
            candidate = blocks(n, :, :);
            reduced(use_value) = candidate(use_value);
         end
      case 4
         % LASTVALUE selects the final nonmissing sample in each bin.
         reduced = blocks(4, :, :);
         for n = 3:-1:1
            use_value = isnan(reduced) & ~isnan(blocks(n, :, :));
            candidate = blocks(n, :, :);
            reduced(use_value) = candidate(use_value);
         end
      otherwise
         error('icemodel:retimeHourlyFixedStep:aggregationMethod', ...
            'unrecognized fixed-step aggregation method')
   end
   returned = reshape(reduced, n_hours, n_columns);
end

function method = aggregationMethod(name)
   %AGGREGATIONMETHOD Map one channel to mean, sum, first, or last.

   % Numeric method codes keep the generated fixed-array switch constant-sized:
   % 1 means mean, 2 sum, 3 first sample, and 4 last sample.
   method = 1;
   if icemodel.isIncrementChannel(name) ...
         || matchesField(name, icemodel.namelists.budgetoutputs('sum'))
      method = 2;
   elseif matchesField(name, icemodel.namelists.budgetoutputs('first'))
      method = 3;
   elseif matchesField(name, icemodel.namelists.budgetoutputs('last')) ...
         || matchesField(name, icemodel.namelists.cumulativeoutputs())
      method = 4;
   end
end

function tf = matchesField(name, candidates)
   %MATCHESFIELD Match a channel without codegen-unsupported cell set operations.

   % Search the canonical lists directly so fixed-array generated code does not
   % need INTERSECT or ISMEMBER for cell arrays.
   tf = false;
   for n = 1:numel(candidates)
      if strcmp(name, candidates{n})
         tf = true;
         return
      end
   end
end

function [hourly, bin_start, bin_end] = aggregateNativeBins(TT)
   %AGGREGATENATIVEBINS Preserve MATLAB behavior on non-fixed time grids.

   % Ordinary channels use the mean; the canonical budget and cumulative
   % classes then replace it with their required reductions.
   vars = TT.Properties.VariableNames;
   hourly = retime(TT, 'hourly', 'mean');
   hourly = replaceAggregation( ...
      hourly, TT, vars, icemodel.namelists.budgetoutputs('sum'), 'sum');
   hourly = replaceAggregation(hourly, TT, vars, ...
      vars(icemodel.isIncrementChannel(vars)), 'sum');
   hourly = replaceAggregation(hourly, TT, vars, ...
      icemodel.namelists.budgetoutputs('first'), 'firstvalue');
   hourly = replaceAggregation(hourly, TT, vars, ...
      icemodel.namelists.budgetoutputs('last'), 'lastvalue');
   hourly = replaceAggregation(hourly, TT, vars, ...
      icemodel.namelists.cumulativeoutputs(), 'lastvalue');
   hourly = hourly(:, vars);

   % Record contiguous raw bounds for every native hourly label. RETIME may
   % insert an empty hour across a gap; ACCUMARRAY leaves those bounds at zero.
   bin_time = dateshift(TT.Properties.RowTimes, 'start', 'hour');
   [in_output, label] = ismember(bin_time, hourly.Properties.RowTimes);
   rows = find(in_output);
   n_bins = height(hourly);
   bin_start = accumarray(label(in_output), rows, [n_bins, 1], @min, 0);
   bin_end = accumarray(label(in_output), rows, [n_bins, 1], @max, 0);
end

function hourly = replaceAggregation(hourly, TT, vars, candidates, method)
   %REPLACEAGGREGATION Replace one present variable class with METHOD.

   % Ignore absent diagnostic channels so every output profile shares this path.
   names = intersect(candidates, vars, 'stable');
   if isempty(names)
      return
   end

   % Native retiming preserves each variable's numeric class and partial bins.
   values = retime(TT(:, names), 'hourly', method);
   hourly(:, names) = values;
end
