function t = ensureUtc(t)
   %ENSUREUTC Coerce a date/datetime input to a UTC-tagged datetime.
   %
   %  t = icemodel.verification.setup.ensureUtc(t)
   %
   %  Accepts a string / char date specifier or a datetime, and returns a
   %  datetime tagged with TimeZone='UTC'. Existing time-zone information
   %  on a datetime input is preserved (the TimeZone is only set when
   %  empty). Used by the ESM-SnowMIP forcing / observation builders to
   %  normalize startdate / enddate kwargs.

   if isstring(t) || ischar(t)
      t = datetime(t, 'TimeZone', 'UTC');
   elseif isdatetime(t) && isempty(t.TimeZone)
      t.TimeZone = 'UTC';
   end
end
