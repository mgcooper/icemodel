function iter = time2iter(time, met)
   
   if ~isdatetime(time)
      [time, tf, dtype] = todatetime(time);
   end
   if isempty(time.TimeZone)
      time.TimeZone = met.Time.TimeZone;
   elseif time.TimeZone ~= met.Time.TimeZone
      warning('time.TimeZone does not equal met.Time.TimeZone')
      time.TimeZone = met.Time.TimeZone;
   end
   iter = find(met.Time == time);
end
