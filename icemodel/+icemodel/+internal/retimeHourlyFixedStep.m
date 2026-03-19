function TT = retimeHourlyFixedStep(TT)
   %RETIMEHOURLYFIXEDSTEP Aggregate aligned 15-minute data to hourly means.
   %
   %  TT = icemodel.internal.retimeHourlyFixedStep(TT)
   %
   % This helper is intended for the postprocess path where the model output
   % is already known to be regular 15-minute data. If the row times are not
   % aligned to exact hourly windows, the caller should fall back to MATLAB's
   % generic RETIME implementation instead.

   % Return early on empty inputs so callers can stay simple.
   if isempty(TT)
      return
   end

   % Enforce the narrow fixed-step contract used by postprocess.
   if ~isFixedStepHourlyCompatible(TT.Properties.RowTimes)
      error(['retimeHourlyFixedStep requires data aligned to exact 15-minute ', ...
         'steps starting on an hourly boundary'])
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
   %ISFIXEDSTEPHOURLYCOMPATIBLE Check the narrow fixed-step retime contract.

   tf = true;
   if numel(time) < 4 || mod(numel(time), 4) ~= 0
      tf = false;
      return
   end

   % Require exact quarter-hour spacing so the reshape-based aggregation is
   % identical to the hourly averaging windows used by RETIME.
   dt = diff(time);
   tf = tf && all(dt == minutes(15));
   tf = tf && minute(time(1)) == 0;
   tf = tf && second(time(1)) == 0;
end
