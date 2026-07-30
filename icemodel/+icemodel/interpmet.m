function met = interpmet(met, newdt)
   %INTERPMET Interpolate met data

   if istimetable(met)

      yyyy = year(met.Time(1));

      switch newdt

         case '15m'

            % Match the source timezone so retime never compares zoned and
            % unzoned datetime arrays.
            timezone = met.Time.TimeZone;
            stime = datetime(yyyy, 1, 1, 0, 0, 0, ...
               'TimeZone', timezone);
            etime = datetime(yyyy, 12, 31, 23, 45, 0, ...
               'TimeZone', timezone);
            newtime = stime:minutes(15):etime;
            met = retime(met, newtime, 'linear', 'EndValues', 'extrap');
      end
   else
      warning('interpmet currently only supports timetables')
   end
end
