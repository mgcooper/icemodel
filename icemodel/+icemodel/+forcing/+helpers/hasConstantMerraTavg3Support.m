function tf = hasConstantMerraTavg3Support(T)
   %HASCONSTANTMERRATAVG3SUPPORT True when glacier channels hold each UTC block.
   %
   %  tf = icemodel.forcing.helpers.hasConstantMerraTavg3Support(T)
   %
   % Checks runoff, albedo, snow depth, and SWE over every available row in each
   % UTC 00/03/... three-hour support block. Hourly Data therefore contributes up
   % to three rows per block and 15-minute met contributes up to twelve. Partial
   % leading/trailing blocks are checked over their available rows; paired NaNs
   % compare equal so a mixed finite/missing block cannot masquerade as a hold.

   tf = istimetable(T);
   if ~tf
      return
   end
   channels = intersect(["runoff", "albedo", "snowd", "swe"], ...
      string(T.Properties.VariableNames), 'stable');
   row_times = T.Properties.RowTimes;
   if isempty(channels) || isempty(row_times)
      return
   end
   if ~isdatetime(row_times) || string(row_times.TimeZone) ~= "UTC" ...
         || ~issorted(row_times)
      tf = false;
      return
   end

   % Compare adjacent rows within each UTC block. Requiring a sorted time axis
   % makes this one linear pass rather than rescanning the full table per block.
   block_start = dateshift(row_times, 'start', 'hour') ...
      - hours(mod(hour(row_times), 3));
   same_block = block_start(2:end) == block_start(1:end-1);
   for name = reshape(channels, 1, [])
      values = T.(char(name));
      if ~isnumeric(values)
         tf = false;
         return
      end

      % Paired NaNs compare equal, but a finite/missing transition fails.
      current = values(2:end, :);
      previous = values(1:end-1, :);
      if ~isequaln(current(same_block, :), previous(same_block, :))
         tf = false;
         return
      end
   end
end
