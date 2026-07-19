function [t1, t2] = periodBounds(period)
   %PERIODBOUNDS Parse a manifest period to UTC datetimes.
   if strlength(string(period.start)) == 0 || strlength(string(period.end)) == 0
      t1 = NaT;
      t2 = NaT;
      return
   end
   t1 = icemodel.verification.setup.ensureUtc(period.start);
   t2 = icemodel.verification.setup.ensureUtc(period.end);
end
