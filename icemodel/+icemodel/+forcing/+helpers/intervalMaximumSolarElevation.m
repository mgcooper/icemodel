function maximum_elevation = intervalMaximumSolarElevation( ...
      times, latitude, longitude, interval)
   %INTERVALMAXIMUMSOLARELEVATION Maximum sun angle over each time support.
   %
   %  maximum_elevation = ...
   %     icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
   %     times, latitude, longitude, interval)
   %
   % TIMES label interval starts. The returned column is the maximum NOAA
   % solar elevation at the start, quarter points, and end of each support.

   sample_times = times(:) + (0:0.25:1) .* interval;
   elevation = icemodel.forcing.helpers.solarElevation( ...
      sample_times, latitude, longitude);
   maximum_elevation = max(elevation, [], 2);
end
