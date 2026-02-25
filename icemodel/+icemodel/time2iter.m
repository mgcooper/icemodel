function step = time2iter(time, met)

   if ~isdatetime(time)
      [time, tf, dtype] = todatetime(time);
   end
   if istimetable(met)
      metTime = met.Time;
   elseif isdatetime(met)
      metTime = met;
   end

   if isempty(time.TimeZone)
      time.TimeZone = metTime.TimeZone;
   elseif time.TimeZone ~= metTime.TimeZone
      warning('time.TimeZone does not equal met.Time.TimeZone')
      time.TimeZone = metTime.TimeZone;
   end
   step = find(metTime == time);
end
