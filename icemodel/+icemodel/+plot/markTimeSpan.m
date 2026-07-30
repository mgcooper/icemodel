function h = markTimeSpan(ax, t_start, t_end, kwargs)
   %MARKTIMESPAN Mark a time span on an axes without touching the legend.
   %
   %  h = icemodel.plot.markTimeSpan(ax, t1, t2)
   %
   % Role
   %  Single source of the span-annotation style report figures use to
   %  highlight an interval (a filled gap, an event window): one
   %  boundary line at each end, excluded from the legend so overlay
   %  labels stay clean.
   %
   % Returns
   %  h : the two constant-line handles.
   %
   % See also: icemodel.plot.compareTimeseries, xline

   arguments
      ax (1, 1) matlab.graphics.axis.Axes
      t_start (1, 1) datetime
      t_end (1, 1) datetime
      kwargs.line_style (1, :) char = ':'
      kwargs.color (1, 3) double = [0.4 0.4 0.4]
   end

   h = [ ...
      xline(ax, t_start, kwargs.line_style, 'Color', kwargs.color, ...
      'HandleVisibility', 'off')
      xline(ax, t_end, kwargs.line_style, 'Color', kwargs.color, ...
      'HandleVisibility', 'off')];
end
