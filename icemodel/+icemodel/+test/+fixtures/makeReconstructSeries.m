function series = makeReconstructSeries()
   %MAKERECONSTRUCTSERIES One year of hourly synthetic met, smooth channels.
   %
   %  series = icemodel.test.fixtures.makeReconstructSeries()
   %
   % Shared fixture for the reconstruction suites: one leap year of
   % hourly samples with smooth seasonal + diurnal structure, so
   % boundary-jump scales stay sane and every suite exercises the same
   % series shape.

   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + ...
      hours(0:24 * 366 - 1)).';
   doy = day(times, 'dayofyear');
   % Smooth seasonal + diurnal structure keeps boundary-jump scales sane.
   tair = 255 + 15 * sin(2 * pi * (doy - 30) / 366) + ...
      2 * sin(2 * pi * hour(times) / 24);
   rh = 80 + 10 * sin(2 * pi * doy / 366);
   psfc = 90000 + 500 * sin(2 * pi * doy / 366);
   lwd = 230 + 60 * sin(2 * pi * (doy - 30) / 366);
   series = timetable(times, tair, rh, psfc, lwd, ...
      'VariableNames', {'tair', 'rh', 'psfc', 'lwd'});
end
