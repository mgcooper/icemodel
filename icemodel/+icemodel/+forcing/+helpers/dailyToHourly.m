function hourly = dailyToHourly(daily, t_daily, t_hourly, kwargs)
   %DAILYTOHOURLY Interpolate daily data onto an hourly (or finer) time axis.
   %
   %  hourly = icemodel.forcing.helpers.dailyToHourly(daily, t_daily, t_hourly)
   %  hourly = ... dailyToHourly(_, method="linear", extrapolate=false)
   %
   % Gridded sources mix hourly channels with daily ones (e.g. MAR daily
   % snow depth, cloud cover, surface temperature, pressure). This helper
   % interpolates the daily series onto the met-file time axis. With
   % extrapolate=true (default) the final partial day is extended by
   % linear extrapolation, matching the legacy builders' handling of the
   % hours after the last daily sample.
   %
   % Inputs
   %  daily    - daily data, one series per column (a row vector is
   %             treated as a single series)
   %  t_daily  - datetimes of the daily samples
   %  t_hourly - target datetimes
   %
   % Outputs
   %  hourly - interpolated data, numel(t_hourly) rows
   %
   % See also: interp1, icemodel.interpmet

   arguments
      daily double
      t_daily (:, 1) datetime
      t_hourly (:, 1) datetime
      kwargs.method (1, 1) string = "linear"
      kwargs.extrapolate (1, 1) logical = true
   end

   if isrow(daily)
      daily = daily(:);
   end

   if kwargs.extrapolate
      hourly = interp1(t_daily, daily, t_hourly, kwargs.method, 'extrap');
   else
      hourly = interp1(t_daily, daily, t_hourly, kwargs.method);
   end
end
