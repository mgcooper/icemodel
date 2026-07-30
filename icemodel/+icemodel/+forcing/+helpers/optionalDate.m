function t = optionalDate(value, fallback)
   %OPTIONALDATE Normalize optional date-like inputs to UTC datetimes.
   if isempty(value) || isequal(value, "")
      t = fallback;
   elseif isa(value, 'datetime')
      if all(isnat(value))
         t = fallback;
      else
         t = value;
      end
   else
      t = datetime(value, 'TimeZone', 'UTC');
   end
   if isempty(t.TimeZone)
      t.TimeZone = 'UTC';
   end
end
