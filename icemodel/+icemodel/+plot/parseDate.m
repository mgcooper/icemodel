function t = parseDate(value, fallback)
   %PARSEDATE Normalize optional date-like input to a datetime.
   t = icemodel.forcing.helpers.optionalDate(value, fallback);
   if isempty(fallback.TimeZone) && ~isempty(t.TimeZone)
      t.TimeZone = '';
   elseif ~isempty(fallback.TimeZone) && isempty(t.TimeZone)
      t.TimeZone = fallback.TimeZone;
   end
end
