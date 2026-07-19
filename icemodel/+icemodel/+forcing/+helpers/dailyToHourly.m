function hourly = dailyToHourly(daily, t_daily, t_hourly, kwargs)
   %DAILYTOHOURLY Interpolate daily data onto an hourly (or finer) time axis.
   %
   %  hourly = icemodel.forcing.helpers.dailyToHourly(daily, t_daily, t_hourly)
   %  hourly = ... dailyToHourly(_, method="linear", extrapolate=false)
   %  hourly = ... dailyToHourly(_, bounds=[0 1])
   %
   % Gridded sources mix hourly channels with daily ones (e.g. MAR daily
   % snow depth, cloud cover, surface temperature, pressure). This helper
   % interpolates the daily series onto the met-file time axis. With
   % extrapolate=true (default), targets outside the native support use the
   % nearest endpoint. This preserves the first/last daily sample without
   % inventing an unbounded linear trend at separately processed year edges.
   % Optional bounds validate both source and interpolated finite values.
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
      kwargs.bounds (1, 2) double = [-Inf Inf]
   end

   if isrow(daily)
      daily = daily(:);
   end

   if kwargs.bounds(1) > kwargs.bounds(2)
      error('icemodel:forcing:dailyToHourly:invalidBounds', ...
         'daily-to-hourly bounds must be ordered [minimum maximum]')
   end
   % Conservative spatial averaging can place an otherwise bounded source a
   % few ulps outside its physical interval. Snap only machine-scale roundoff
   % to the boundary; materially invalid source values must still fail closed.
   bound_scale = max([1, abs(kwargs.bounds(isfinite(kwargs.bounds)))]);
   bound_tolerance = 256 * eps(bound_scale);
   finite = isfinite(daily);
   near_lower = finite & daily < kwargs.bounds(1) ...
      & daily >= kwargs.bounds(1) - bound_tolerance;
   near_upper = finite & daily > kwargs.bounds(2) ...
      & daily <= kwargs.bounds(2) + bound_tolerance;
   daily(near_lower) = kwargs.bounds(1);
   daily(near_upper) = kwargs.bounds(2);
   if any(daily(finite) < kwargs.bounds(1) ...
         | daily(finite) > kwargs.bounds(2))
      % Report the extrema so a source-product defect can be distinguished
      % from numerical noise without instrumenting the caller.
      source_min = min(daily(finite));
      source_max = max(daily(finite));
      error('icemodel:forcing:dailyToHourly:sourceOutOfBounds', ...
         ['daily source values [%.17g, %.17g] fall outside the requested ' ...
         'physical bounds [%.17g, %.17g]'], source_min, source_max, ...
         kwargs.bounds(1), kwargs.bounds(2))
   end

   if kwargs.extrapolate
      % Clamp only the query coordinate. Interpolation inside the native
      % support is unchanged, while both unsupported tails hold an observed
      % endpoint instead of extending the last linear slope.
      query_time = t_hourly;
      query_time(query_time < t_daily(1)) = t_daily(1);
      query_time(query_time > t_daily(end)) = t_daily(end);
      hourly = interp1(t_daily, daily, query_time, kwargs.method);
   else
      hourly = interp1(t_daily, daily, t_hourly, kwargs.method);
   end

   % Apply the same roundoff-only snap after interpolation so the advertised
   % postcondition is exact without hiding a substantive overshoot.
   finite = isfinite(hourly);
   near_lower = finite & hourly < kwargs.bounds(1) ...
      & hourly >= kwargs.bounds(1) - bound_tolerance;
   near_upper = finite & hourly > kwargs.bounds(2) ...
      & hourly <= kwargs.bounds(2) + bound_tolerance;
   hourly(near_lower) = kwargs.bounds(1);
   hourly(near_upper) = kwargs.bounds(2);
   if any(hourly(finite) < kwargs.bounds(1) ...
         | hourly(finite) > kwargs.bounds(2))
      error('icemodel:forcing:dailyToHourly:outputOutOfBounds', ...
         'interpolated daily values fall outside the requested physical bounds')
   end
end
