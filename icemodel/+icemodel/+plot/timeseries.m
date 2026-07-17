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
      kwargs.marker_size (1, 1) double {mustBePositive} = 6
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

   washeld = ishold(kwargs.axes);
   hold(kwargs.axes, 'on')

   for n = 1:n_series
      series = values(:, n);
      [line_time, line_series] = breakOmittedTimeGaps(time, series);
      plot_args = cell(1, 8);
      plot_args(1:4) = {'LineStyle', kwargs.line_style, ...
         'LineWidth', kwargs.line_width};
      n_args = 4;
      if ~isempty(kwargs.color)
         plot_args(n_args + (1:2)) = {'Color', kwargs.color};
         n_args = n_args + 2;
      end
      if display_name(n) ~= ""
         plot_args(n_args + (1:2)) = {'DisplayName', char(display_name(n))};
         n_args = n_args + 2;
      end
      h(n) = plot(kwargs.axes, line_time, line_series, plot_args{1:n_args});
      addSparseMarkers(kwargs.axes, time, series, kwargs, display_name(n));
   end

   if ~washeld
      hold(kwargs.axes, 'off')
   end
end

function [line_time, line_series] = breakOmittedTimeGaps(time, series)
   %BREAKOMITTEDTIMEGAPS Insert NaNs only when cadence is safely inferable.
   line_time = time;
   line_series = series;
   steps_s = seconds(diff(time));
   if numel(steps_s) < 3 || any(steps_s <= 0)
      return
   end

   cadence_s = median(steps_s);
   tolerance_s = max(1e-6, cadence_s * 1e-6);
   nominal = abs(steps_s - cadence_s) <= tolerance_s;
   if nnz(nominal) < 2
      return
   end

   gaps = find(steps_s > 1.5 * cadence_s);
   if isempty(gaps)
      return
   end

   % Insert one unavailable midpoint per omitted interval. Original finite
   % samples remain untouched, so sparse-marker and interval observations keep
   % their values while the line renderer cannot bridge the outage.
   n_rows = numel(time) + numel(gaps);
   line_time = repmat(time(1), n_rows, 1);
   line_series = nan(n_rows, 1);
   source_row = 1;
   output_row = 1;
   for gap = reshape(gaps, 1, [])
      count = gap - source_row + 1;
      rows = output_row:(output_row + count - 1);
      line_time(rows) = time(source_row:gap);
      line_series(rows) = series(source_row:gap);
      output_row = output_row + count;
      line_time(output_row) = time(gap) + (time(gap + 1) - time(gap)) / 2;
      output_row = output_row + 1;
      source_row = gap + 1;
   end
   line_time(output_row:end) = time(source_row:end);
   line_series(output_row:end) = series(source_row:end);
end

function addSparseMarkers(ax, time, series, kwargs, ~)
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
