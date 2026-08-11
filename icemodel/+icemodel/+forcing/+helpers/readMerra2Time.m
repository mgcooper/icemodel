function Time = readMerra2Time(filename)
   %READMERRA2TIME Decode only a MERRA-2 file's native UTC time coordinate.
   %
   %  Time = icemodel.forcing.helpers.readMerra2Time(filename)
   %
   % This function reads the `time` coordinate and its units attribute. It does
   % not open a gridded science variable. It returns the native snapshot or
   % interval-center stamps unchanged. The builders do any relabeling.

   arguments
      filename (1, 1) string
   end

   % MERRA-2 stores time as minutes from a timestamp embedded in the units.
   raw = double(ncread(filename, 'time'));
   units = string(ncreadatt(filename, 'time', 'units'));
   token = regexp(units, ...
      'minutes since (\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2})', ...
      'tokens', 'once');
   if isempty(token)
      error('icemodel:forcing:readMerra2Time:badUnits', ...
         'unexpected MERRA time units in %s: %s', filename, units)
   end
   origin = datetime(token{1}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', ...
      'TimeZone', 'UTC');
   Time = origin + minutes(raw(:));
end
