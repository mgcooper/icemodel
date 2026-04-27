function h = timeseries(time, values, kwargs)
   %TIMESERIES Plot one or more datetime-indexed series on one axes.
   %
   %  h = icemodel.plot.timeseries(time, values)
   %  h = icemodel.plot.timeseries(time, values, axes=ax)
   %  h = icemodel.plot.timeseries(time, values, display_name="obs")
   %
   % TIME is an N-by-1 datetime vector and VALUES is an N-by-M numeric array.
   % The helper exists so callers do not need to open-code the common
   % "datetime on x, one variable on y" plotting path. Sparse observational
   % series often contain isolated valid samples on a denser model grid; in
   % that case the helper overlays point markers automatically so the data do
   % not disappear when no line segments can be drawn.

   arguments
      time (:, 1) datetime
      values {mustBeNumeric}
      kwargs.axes matlab.graphics.axis.Axes = gca
      kwargs.display_name string = strings(0, 1)
      kwargs.line_style (1, :) char = '-'
      kwargs.color = []
      kwargs.line_width (1, 1) double {mustBePositive} = 1.2
      kwargs.marker_style (1, :) char = 'auto'
      kwargs.marker_size (1, 1) double {mustBePositive} = 12
   end

   if isvector(values)
      values = values(:);
   end
   if size(values, 1) ~= numel(time)
      error('values must have the same number of rows as time')
   end

   n_series = size(values, 2);
   display_name = normalizeNames(kwargs.display_name, n_series);
   h = gobjects(n_series, 1);
   hold_was_on = ishold(kwargs.axes);
   hold(kwargs.axes, 'on')

   for i = 1:n_series
      series = values(:, i);
      plot_args = cell(1, 8);
      plot_args(1:4) = {'LineStyle', kwargs.line_style, ...
         'LineWidth', kwargs.line_width};
      n_args = 4;
      if ~isempty(kwargs.color)
         plot_args(n_args + (1:2)) = {'Color', kwargs.color};
         n_args = n_args + 2;
      end
      if display_name(i) ~= ""
         plot_args(n_args + (1:2)) = {'DisplayName', char(display_name(i))};
         n_args = n_args + 2;
      end
      h(i) = plot(kwargs.axes, time, series, plot_args{1:n_args});
      addSparseMarkers(kwargs.axes, time, series, kwargs, display_name(i));
   end

   if ~hold_was_on
      hold(kwargs.axes, 'off')
   end
end

function addSparseMarkers(ax, time, series, kwargs, display_name)
   %ADDSPARSEMARKERS Overlay markers when a sparse series would otherwise vanish.

   marker_style = resolveMarkerStyle(series, kwargs.marker_style);
   if marker_style == "none"
      return
   end

   marker_args = cell(1, 8);
   marker_args(1:6) = {'LineStyle', 'none', 'Marker', char(marker_style), ...
      'MarkerSize', kwargs.marker_size};
   n_args = 6;
   if ~isempty(kwargs.color)
      marker_args(n_args + (1:2)) = {'Color', kwargs.color};
      n_args = n_args + 2;
   end

   % Keep the main line as the legend representative so one series still
   % contributes only one legend entry.
   marker_args(n_args + (1:2)) = {'HandleVisibility', 'off'};
   n_args = n_args + 2;
   plot(ax, time, series, marker_args{1:n_args});
end

function marker_style = resolveMarkerStyle(series, requested_style)
   %RESOLVEMARKERSTYLE Choose whether sparse samples need explicit markers.

   if ~strcmp(requested_style, 'auto')
      marker_style = string(requested_style);
      return
   end

   finite_mask = isfinite(series);
   n_finite = nnz(finite_mask);
   if n_finite == 0
      marker_style = "none";
      return
   end

    prev_finite = [false; finite_mask(1:end-1)];
    next_finite = [finite_mask(2:end); false];
    isolated_mask = finite_mask & ~prev_finite & ~next_finite;
   finite_fraction = n_finite / numel(series);
   if any(isolated_mask) || finite_fraction <= 0.25
      marker_style = ".";
   else
      marker_style = "none";
   end
end

function display_name = normalizeNames(display_name, n_series)
   %NORMALIZENAMES Expand or validate per-series display names.

   if isempty(display_name)
      display_name = repmat("", n_series, 1);
      return
   end
   display_name = reshape(display_name, [], 1);
   if isscalar(display_name)
      display_name = repmat(display_name, n_series, 1);
      return
   end
   if numel(display_name) ~= n_series
      error('display_name must be scalar or match the number of series')
   end
end
