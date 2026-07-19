function time = gcnetTime(value, units)
   %GCNETTIME Convert Vandecrux/GC-Net numeric time to UTC datetimes.
   units = string(units);
   if startsWith(units, "days since ")
      raw = extractAfter(units, "days since ");
      base = datetime(raw, 'InputFormat', 'yyyy-M-d H:m:s', ...
         'TimeZone', 'UTC');
      time = base + days(value);
   elseif startsWith(units, "hours since ")
      raw = extractAfter(units, "hours since ");
      base = datetime(raw, 'InputFormat', 'yyyy-M-d H:m:s', ...
         'TimeZone', 'UTC');
      time = base + hours(value);
   else
      error('icemodel:forcing:helpers:gcnetTime:unsupportedTimeUnits', ...
         'unsupported GC-Net/Vandecrux time units: %s', units)
   end
end
