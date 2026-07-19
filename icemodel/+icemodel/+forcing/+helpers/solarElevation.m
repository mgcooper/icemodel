function elevation = solarElevation(Time, latitude, longitude)
   %SOLARELEVATION Approximate NOAA geometric solar elevation in degrees.
   %
   % `Time` is interpreted as UTC when unzoned. Zoned inputs preserve their
   % instants while converting to UTC. Longitude is degrees east.
   arguments
      Time datetime
      latitude (1, 1) double {mustBeFinite}
      longitude (1, 1) double {mustBeFinite}
   end

   % NOAA's fractional-year approximation supplies declination and equation of
   % time without external ephemeris dependencies.
   Time.TimeZone = 'UTC';
   decimal_hour = hour(Time) + minute(Time) ./ 60 + second(Time) ./ 3600;
   year_number = year(Time);
   is_leap_year = mod(year_number, 4) == 0 ...
      & (mod(year_number, 100) ~= 0 | mod(year_number, 400) == 0);
   days_in_year = 365 + is_leap_year;
   gamma = 2 * pi ./ days_in_year .* ...
      (day(Time, 'dayofyear') - 1 + (decimal_hour - 12) ./ 24);
   equation_of_time = 229.18 .* (0.000075 + 0.001868 .* cos(gamma) ...
      - 0.032077 .* sin(gamma) - 0.014615 .* cos(2 .* gamma) ...
      - 0.040849 .* sin(2 .* gamma));
   declination = 0.006918 - 0.399912 .* cos(gamma) ...
      + 0.070257 .* sin(gamma) - 0.006758 .* cos(2 .* gamma) ...
      + 0.000907 .* sin(2 .* gamma) - 0.002697 .* cos(3 .* gamma) ...
      + 0.00148 .* sin(3 .* gamma);

   % UTC has zero timezone offset. Wrap the solar hour angle to [-180, 180]
   % before evaluating the standard solar-zenith relation.
   true_solar_minutes = decimal_hour .* 60 + equation_of_time + 4 .* longitude;
   hour_angle = mod(true_solar_minutes ./ 4, 360) - 180;
   cos_zenith = sind(latitude) .* sin(declination) ...
      + cosd(latitude) .* cos(declination) .* cosd(hour_angle);
   cos_zenith = min(max(cos_zenith, -1), 1);
   elevation = asind(cos_zenith);
end
