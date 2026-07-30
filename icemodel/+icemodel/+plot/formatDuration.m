function label = formatDuration(duration_hours)
   %FORMATDURATION Render a duration in hours as a human-readable label.
   %
   %  label = icemodel.plot.formatDuration(3)      % "3 h"
   %  label = icemodel.plot.formatDuration(32300)  % "3.7 y"
   %
   % Role
   %  Single source of the duration labels report figures and tables
   %  use: hours below two days, days below two years, years above —
   %  never scientific notation.
   %
   % See also: icemodel.plot.compareTimeseries

   arguments
      duration_hours (1, 1) double {mustBeNonnegative}
   end

   if duration_hours < 48
      label = sprintf('%.3g h', duration_hours);
   elseif duration_hours < 2 * 24 * 365.25
      label = sprintf('%.3g d', duration_hours / 24);
   else
      label = sprintf('%.3g y', duration_hours / (24 * 365.25));
   end
end
