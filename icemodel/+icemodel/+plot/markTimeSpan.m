function h = markTimeSpan(ax, t_start, t_end, kwargs)
   %MARKTIMESPAN Mark a time span on an axes without touching the legend.
   %
   %  h = icemodel.plot.markTimeSpan(ax, t1, t2)
   %
   % Role
   %  Draws the span annotation report figures use to highlight an interval
   %  (a filled gap, an event window). The annotation is excluded from the
   %  legend so overlay labels stay clean.
   %
   %  style="lines" draws one boundary line at each end. style="fill" shades
   %  the interval instead, for panels that highlight many spans at once and
   %  would be unreadable with boundary lines.
   %
   % Returns
   %  h : the two constant-line handles for style="lines", or the single
   %      region handle for style="fill".
   %
   % See also: icemodel.plot.compareTimeseries, xline

   arguments
      ax (1, 1) matlab.graphics.axis.Axes
      t_start (1, 1) datetime
      t_end (1, 1) datetime
      kwargs.line_style (1, :) char = ':'
      kwargs.color (1, 3) double = [0.4 0.4 0.4]
      kwargs.style (1, 1) string {mustBeMember(kwargs.style, ...
         ["lines", "fill"])} = "lines"
      kwargs.face_alpha (1, 1) double {mustBeInRange( ...
         kwargs.face_alpha, 0, 1)} = 0.12
   end

   if kwargs.style == "fill"
      h = xregion(ax, t_start, t_end, 'FaceColor', kwargs.color, ...
         'FaceAlpha', kwargs.face_alpha, 'EdgeColor', 'none');
      h.HandleVisibility = 'off';
      return
   end

   h = [ ...
      xline(ax, t_start, kwargs.line_style, 'Color', kwargs.color, ...
      'HandleVisibility', 'off')
      xline(ax, t_end, kwargs.line_style, 'Color', kwargs.color, ...
      'HandleVisibility', 'off')];
end
