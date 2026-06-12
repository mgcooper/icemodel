function out = scatterplot(x, y, kwargs)
   %SCATTERPLOT Plot paired data with 1:1 and fitted-line overlays.
   %
   %  out = icemodel.plot.scatterplot(x, y)
   %  out = icemodel.plot.scatterplot({x1, x2}, {y1, y2}, display_name=names)
   %
   % This helper owns the common "candidate versus target" scatter contract:
   % finite-pair filtering, point clouds, one 1:1 reference line, one linear
   % least-squares fit per series, and a legend built from the plotted fit
   % handles. Callers can reuse it without open-coding the same plot mechanics.

   arguments
      x
      y
      kwargs.axes matlab.graphics.axis.Axes = gca
      kwargs.display_name string = strings(0, 1)
      kwargs.colors double = []
      kwargs.x_label (1, 1) string = "target"
      kwargs.y_label (1, 1) string = "candidate"
      kwargs.marker_size (1, 1) double {mustBePositive} = 18
      kwargs.marker_face_alpha (1, 1) double = 0.30
      kwargs.line_width (1, 1) double {mustBePositive} = 1.6
      kwargs.legend_location (1, 1) string = "best"
      kwargs.legend_font_size (1, 1) double {mustBePositive} = 10
      kwargs.show_bias (1, 1) logical = false
   end

   [x_series, y_series] = normalizeSeries(x, y);
   n_series = numel(x_series);
   names = normalizeNames(kwargs.display_name, n_series);
   colors = normalizeColors(kwargs.colors, n_series);

   washeld = ishold(kwargs.axes);
   hold(kwargs.axes, 'on')

   scatter_handles = gobjects(n_series, 1);
   fit_handles = gobjects(n_series, 1);
   stats = repmat(statTemplate(), n_series, 1);
   all_values = cell(n_series, 1);

   for i = 1:n_series
      [x_finite, y_finite] = finitePairs(x_series{i}, y_series{i});
      stats(i).name = names(i);
      stats(i).n = numel(x_finite);
      all_values{i} = [x_finite; y_finite];

      if isempty(x_finite)
         continue
      end

      % Points are hidden from the legend because the fit lines carry the
      % series identity plus R-squared. This matches the THF validation plots.
      scatter_handles(i) = scatter(kwargs.axes, x_finite, y_finite, ...
         kwargs.marker_size, colors(i, :), 'filled', ...
         'MarkerFaceAlpha', kwargs.marker_face_alpha, ...
         'HandleVisibility', 'off');

      if numel(x_finite) < 2
         continue
      end

      [x_fit, y_fit, stats(i).slope, stats(i).intercept, stats(i).r_squared] = ...
         fitLine(x_finite, y_finite);
      fit_handles(i) = plot(kwargs.axes, x_fit, y_fit, '-', ...
         'Color', colors(i, :), ...
         'LineWidth', kwargs.line_width, ...
         'DisplayName', sprintf('%s R^2=%.3f', names(i), ...
         stats(i).r_squared));
      stats(i).bias = mean(y_finite - x_finite);
   end

   one_to_one = gobjects(1);
   limits = dataLimits(all_values);
   if all(isfinite(limits))
      one_to_one = plot(kwargs.axes, limits, limits, 'k--', ...
         'LineWidth', 1.0, 'DisplayName', '1:1');
      xlim(kwargs.axes, limits)
      ylim(kwargs.axes, limits)
   end

   xlabel(kwargs.axes, kwargs.x_label)
   ylabel(kwargs.axes, kwargs.y_label)
   grid(kwargs.axes, 'on')

   legend_handles = fit_handles(isgraphics(fit_handles));
   if isgraphics(one_to_one)
      legend_handles = [legend_handles; one_to_one];
   end
   if ~isempty(legend_handles)
      legend(kwargs.axes, legend_handles, 'Location', kwargs.legend_location, ...
         'FontSize', kwargs.legend_font_size, 'Box', 'off')
   end

   if kwargs.show_bias && n_series == 1 && isfinite(stats(1).bias)
      text(kwargs.axes, 0.05, 0.93, sprintf('bias = %+.3g', stats(1).bias), ...
         'Units', 'normalized', 'FontSize', kwargs.legend_font_size);
   end

   if ~washeld
      hold(kwargs.axes, 'off')
   end

   out = struct( ...
      'scatter', scatter_handles, ...
      'fit', fit_handles, ...
      'one_to_one', one_to_one, ...
      'stats', struct2table(stats));
end

function [x_series, y_series] = normalizeSeries(x, y)
   %NORMALIZESERIES Return paired column-vector cells.

   x_series = toSeriesCells(x);
   y_series = toSeriesCells(y);

   if isscalar(x_series) && numel(y_series) > 1
      x_series = repmat(x_series, numel(y_series), 1);
   elseif isscalar(y_series) && numel(x_series) > 1
      y_series = repmat(y_series, numel(x_series), 1);
   end

   if numel(x_series) ~= numel(y_series)
      error('x and y must contain the same number of series')
   end
end

function series = toSeriesCells(values)
   %TOSERIESCELLS Convert vectors, matrices, or cell arrays into cells.

   if iscell(values)
      series = reshape(values, [], 1);
      for i = 1:numel(series)
         series{i} = series{i}(:);
      end
      return
   end

   if isvector(values)
      series = {values(:)};
      return
   end

   series = num2cell(values, 1)';
   for i = 1:numel(series)
      series{i} = series{i}(:);
   end
end

function names = normalizeNames(names, n_series)
   %NORMALIZENAMES Expand or validate legend names.

   if isempty(names)
      names = "series " + string(1:n_series)';
      return
   end

   names = reshape(names, [], 1);
   if isscalar(names)
      names = repmat(names, n_series, 1);
   elseif numel(names) ~= n_series
      error('display_name must be scalar or match the number of series')
   end
end

function colors = normalizeColors(colors, n_series)
   %NORMALIZECOLORS Return one RGB color per series.

   if isempty(colors)
      colors = lines(n_series);
      return
   end

   if size(colors, 1) == 1 && n_series > 1
      colors = repmat(colors, n_series, 1);
   end
   if size(colors, 1) ~= n_series || size(colors, 2) ~= 3
      error('colors must be an N-by-3 RGB array')
   end
end

function [x_finite, y_finite] = finitePairs(x, y)
   %FINITEPAIRS Drop missing or non-finite values from paired data.

   if numel(x) ~= numel(y)
      error('paired x and y series must have equal lengths')
   end

   keep = isfinite(x) & isfinite(y);
   x_finite = x(keep);
   y_finite = y(keep);
end

function [x_fit, y_fit, slope, intercept, r_squared] = fitLine(x, y)
   %FITLINE Compute a least-squares fit and R-squared.

   coeffs = polyfit(x, y, 1);
   slope = coeffs(1);
   intercept = coeffs(2);
   x_fit = linspace(min(x), max(x), 200)';
   y_fit = polyval(coeffs, x_fit);

   r_squared = NaN;
   if std(x) > 0 && std(y) > 0
      C = corrcoef(x, y);
      r_squared = C(1, 2) ^ 2;
   end
end

function limits = dataLimits(values)
   %DATALIMITS Return padded equal x/y limits from all finite plotted values.

   values = vertcat(values{:});
   if isempty(values)
      limits = [NaN NaN];
      return
   end

   lo = min(values);
   hi = max(values);
   if lo == hi
      pad = max(abs(lo), 1) * 0.05;
      limits = [lo - pad, hi + pad];
   else
      pad = 0.03 * (hi - lo);
      limits = [lo - pad, hi + pad];
   end
end

function row = statTemplate()
   %STATTEMPLATE Return one row of scatter diagnostics.

   row = struct( ...
      'name', "", ...
      'n', 0, ...
      'slope', NaN, ...
      'intercept', NaN, ...
      'r_squared', NaN, ...
      'bias', NaN);
end
