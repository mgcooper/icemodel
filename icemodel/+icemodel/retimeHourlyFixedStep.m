function TT = retimeHourlyFixedStep(TT)
   %RETIMEHOURLYFIXEDSTEP Aggregate aligned 15-minute data to hourly means.
   %
   %  TT = icemodel.retimeHourlyFixedStep(TT)
   %
   % This helper is intended for the postprocess path where the model output
   % is known to be 15-minute data (opts.dt == 900). The reshape-based
   % aggregation requires only that the sample count is divisible by 4 and
   % the first sample is hour-aligned.

   % Return early on empty inputs so callers can stay simple.
   if isempty(TT)
      return
   end

   % Enforce the narrow fixed-step contract used by postprocess.
   if ~isFixedStepHourlyCompatible(TT.Properties.RowTimes)
      error(['retimeHourlyFixedStep requires a sample count divisible by 4 ', ...
         'starting on an hourly boundary'])
   end

   % Collapse each block of four quarter-hour samples into one hourly mean.
   vars = TT.Properties.VariableNames;
   data = TT{:, vars};
   n_hours = size(data, 1) / 4;
   n_vars = size(data, 2);
   data = reshape(mean(reshape(data, 4, n_hours, n_vars), 1), n_hours, n_vars);
   time = TT.Properties.RowTimes(1:4:end);
   TT = array2timetable(data, 'RowTimes', time, 'VariableNames', vars);

   % Match the legacy postprocess behavior that drops any synthetic Feb 29
   % row inserted by timetable-based hourly retiming.
   TT = TT(~(month(TT.Properties.RowTimes) == 2 ...
      & day(TT.Properties.RowTimes) == 29), :);
end

function tf = isFixedStepHourlyCompatible(time)
   %ISFIXEDSTEPHOURLYCOMPATIBLE Check the fixed-step retime contract.
   %
   % The caller (postprocess) guarantees the cadence via opts.dt == 900.
   % This check only enforces the structural requirements for the
   % reshape-based aggregation: divisible by 4 and hour-aligned start.

   tf = numel(time) >= 4 ...
      && mod(numel(time), 4) == 0 ...
      && minute(time(1)) == 0 ...
      && second(time(1)) == 0;
end
