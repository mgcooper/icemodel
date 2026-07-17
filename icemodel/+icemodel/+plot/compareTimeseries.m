function out = compareTimeseries(data, varname, kwargs)
   %COMPARETIMESERIES Overlay one variable from multiple timetables.
   %
   %  out = icemodel.plot.compareTimeseries({tt1, tt2}, "tair", names=...)
   %
   % DATA is a timetable or a cell array of timetables. VARNAME is plotted from
   % each timetable that contains a finite numeric column. The helper owns the
   % common visual-QA overlay path so verification plots do not duplicate
   % retiming, temperature conversion, unit labels, and line-only rendering.
   % SHORTWAVE_WEIGHTED_MEAN uses ALBEDO_SOURCE to preserve the source's
   % physical meaning. RADIOMETER ratios use the ratio of daily reflected and
   % incoming shortwave energy with six hours of valid support. MODEL_STATE
   % values use an exact-grid finite-sample arithmetic mean without a solar
   % support gate. MODEL_RATIO values use the daily energy ratio without the
   % observation-specific six-hour gate. Already-daily products without SWD
   % retain a finite-sample fallback. Daily WDIR means are circular and use
   % WSPD as the vector weight when available.

   arguments
      data
      varname (1, 1) string
      kwargs.axes matlab.graphics.axis.Axes = gca
      kwargs.names string = strings(0, 1)
      kwargs.frequency (1, 1) string ...
         {mustBeMember(kwargs.frequency, ["native", "hourly", "daily"])} ...
         = "native"
      kwargs.aggregation (1, :) string ...
         {mustBeMember(kwargs.aggregation, ...
         ["mean", "mean_omitmissing", "shortwave_weighted_mean", ...
          "daily_total"])} = "mean"
      kwargs.albedo_source (1, :) string ...
         {mustBeMember(kwargs.albedo_source, ...
         ["radiometer", "model_state", "model_ratio"])} = "radiometer"
      kwargs.line_width (1, 1) double {mustBePositive} = 1.2
      kwargs.line_style (1, :) char = '-'
      kwargs.display_unit (1, 1) string ...
         {mustBeMember(kwargs.display_unit, ["", "mm", "mmWE"])} = ""
      kwargs.show_legend (1, 1) logical = true
      kwargs.startdate = ""
      kwargs.enddate = ""
   end

   tables = normalizeTables(data);
   names = normalizeNames(kwargs.names, numel(tables));
   aggregations = reshape(kwargs.aggregation, 1, []);
   if isscalar(aggregations)
      aggregations = repmat(aggregations, 1, numel(tables));
   elseif numel(aggregations) ~= numel(tables)
      error('icemodel:plot:compareTimeseries:badAggregation', ...
         'aggregation must be scalar or match the number of timetables')
   end
   albedo_sources = reshape(kwargs.albedo_source, 1, []);
   if isscalar(albedo_sources)
      albedo_sources = repmat(albedo_sources, 1, numel(tables));
   elseif numel(albedo_sources) ~= numel(tables)
      error('icemodel:plot:compareTimeseries:badAlbedoSource', ...
         'albedo_source must be scalar or match the number of timetables')
   end
   colors = icemodel.plot.sourceColor(names);
   lines_out = gobjects(numel(tables), 1);
   plotted = false(numel(tables), 1);
   y_unit = "";

   % Hold state is preserved so callers can compose this helper inside larger
   % figure layouts without surprising existing axes content.
   washeld = ishold(kwargs.axes);
   hold(kwargs.axes, 'on')

   for k = 1:numel(tables)
      series_kwargs = kwargs;
      series_kwargs.aggregation = aggregations(k);
      series_kwargs.albedo_source = albedo_sources(k);
      [tt, unit] = extractSeries(tables{k}, varname, series_kwargs);
      if isempty(tt)
         continue
      end
      if y_unit == ""
         y_unit = unit;
      end
      lines_out(k) = icemodel.plot.timeseries(tt.Time, tt.value, ...
         axes=kwargs.axes, display_name=names(k), color=colors(k, :), ...
         line_width=kwargs.line_width, line_style=kwargs.line_style, ...
         marker_style='none');
      plotted(k) = true;
   end

   if ~washeld
      hold(kwargs.axes, 'off')
   end

   ylabel(kwargs.axes, ylabelText(varname, y_unit), 'Interpreter', 'none')
   grid(kwargs.axes, 'on')
   if any(plotted) && kwargs.show_legend
      styleLegend(kwargs.axes)
   end

   out = struct('lines', lines_out(plotted), 'plotted', plotted, ...
      'unit', y_unit);
end

function tables = normalizeTables(data)
   %NORMALIZETABLES Convert a timetable or cell array to a row cell array.

   if istimetable(data)
      tables = {icemodel.plot.canonicalTimeDimension(data)};
      return
   end
   if iscell(data)
      tables = reshape(data, 1, []);
      for k = 1:numel(tables)
         if istimetable(tables{k})
            tables{k} = icemodel.plot.canonicalTimeDimension(tables{k});
         end
      end
      return
   end
   error('icemodel:plot:compareTimeseries:badInput', ...
      'data must be a timetable or cell array of timetables')
end

function names = normalizeNames(names, n_tables)
   %NORMALIZENAMES Return one display name per input timetable.

   if isempty(names)
      names = "series " + string(1:n_tables);
      return
   end
   names = reshape(names, 1, []);
   if isscalar(names)
      names = repmat(names, 1, n_tables);
   elseif numel(names) ~= n_tables
      error('icemodel:plot:compareTimeseries:badNames', ...
         'names must be scalar or match the number of timetables')
   end
end

function [tt, unit] = extractSeries(T, varname, kwargs)
   %EXTRACTSERIES Pull one finite numeric variable into a plotting timetable.

   tt = timetable.empty;
   unit = "";
   if ~istimetable(T) || ~ismember(varname, string(T.Properties.VariableNames))
      return
   end

   values = T.(char(varname));
   if ~isnumeric(values) || ~isvector(values) || ~any(isfinite(values))
      return
   end
   unit = icemodel.plot.variableUnit(T, varname);
   values = values(:);
   if shouldConvertTemperature(varname, unit, values)
      values = values - 273.15;
      unit = "degC";
   end

   % Resolve physical daily semantics centrally so every caller of the shared
   % renderer inherits ratio and circular-variable behavior automatically.
   aggregation = dailyAggregation(varname, kwargs.frequency, ...
      kwargs.aggregation, kwargs.albedo_source);
   tt = timetable(values, 'RowTimes', T.Time, 'VariableNames', {'value'});
   if any(aggregation == ...
         ["shortwave_weighted_mean", "shortwave_weighted_model_mean"]) ...
         && varname == "albedo" ...
         && ismember("swd", string(T.Properties.VariableNames))
      % Weight source-screened albedo by measured incoming shortwave. This is
      % the daily reflected/incoming energy ratio when albedo = swu/swd, and it
      % prevents low-sun ratios from dominating accumulation-site diagnostics.
      weight = T.swd;
      if isnumeric(weight) && isvector(weight) ...
            && numel(weight) == height(T)
         tt.weight = weight(:);
      end
   elseif varname == "wdir" ...
         && ismember("wspd", string(T.Properties.VariableNames))
      % A mean wind direction is the direction of the mean wind vector, so
      % stronger winds carry proportionally more directional information.
      weight = T.wspd;
      if isnumeric(weight) && isvector(weight) ...
            && numel(weight) == height(T)
         tt.weight = weight(:);
      end
   end
   tt = icemodel.plot.canonicalTimeDimension(tt);
   % Infer dense support from the unfiltered source axis. A short plot window
   % must not turn one posting from a known dense series into a complete day.
   daily_support = inferDailySupport(utcTimes(tt.Time));
   tt = icemodel.plot.filterDateWindow(tt, kwargs.startdate, kwargs.enddate);
   if isempty(tt)
      return
   end
   switch kwargs.aggregation
      case "daily_total"
         if ~isIntervalObservation(T, varname)
            [tt, unit] = dailyTotal(tt, unit, daily_support);
         end
      case {"mean", "mean_omitmissing", "shortwave_weighted_mean"}
         if kwargs.frequency == "daily"
            tt = dailyMean(tt, daily_support, aggregation);
         elseif kwargs.frequency ~= "native"
            tt = retime(tt, char(kwargs.frequency), 'mean');
         end
   end
   % Presentation scaling happens after temporal aggregation so rates are
   % integrated exactly once and interval observations retain their support.
   [tt, unit] = convertDisplayUnit(tt, unit, kwargs.display_unit);
   % Aggregation can prove that every requested interval is incomplete. Do not
   % create an all-NaN graphics object whose legend promises a nonexistent line.
   if isempty(tt) || ~any(isfinite(tt.value))
      tt = timetable.empty;
   end
end

function aggregation = dailyAggregation(varname, frequency, requested, ...
      albedo_source)
   %DAILYAGGREGATION Resolve variable-aware physical reduction semantics.

   aggregation = requested;
   if frequency ~= "daily"
      return
   end

   % Source radiometers are ratios, native model albedo is a surface state, and
   % a ratio derived from model fluxes remains a ratio but does not inherit the
   % sparse-observation support gate. Keeping these cases distinct prevents
   % valid modeled polar-night state from being censored or a flux ratio from
   % being averaged as though it were a prognostic state.
   if varname == "albedo" ...
         && any(requested == ...
         ["mean", "mean_omitmissing", "shortwave_weighted_mean"])
      if albedo_source == "model_state"
         aggregation = "mean_omitmissing";
      elseif albedo_source == "model_ratio"
         aggregation = "shortwave_weighted_model_mean";
      else
         aggregation = "shortwave_weighted_mean";
      end
   elseif varname == "wdir" && requested == "mean"
      % Direction is circular: 350 and 10 degrees average to north, not south.
      aggregation = "circular_mean_degrees";
   elseif varname == "wdir" && requested == "mean_omitmissing"
      % Preserve omit-missing support while applying the same circular geometry.
      aggregation = "circular_mean_degrees_omitmissing";
   end
end

function [tt, unit] = convertDisplayUnit(tt, unit, requested)
   %CONVERTDISPLAYUNIT Apply explicit water-equivalent display scaling.

   % Source metadata uses several spellings for metre water equivalent. Keep
   % the artifact untouched while presenting one canonical unit label.
   key = lower(regexprep(string(unit), '[\s._-]', ''));
   if any(key == ["mwe", "mwaterequivalent"])
      unit = "mWE";
   end

   if requested == ""
      return
   end
   if requested == "mm" && unit == "m"
      tt.value = tt.value .* 1000;
      unit = "mm";
   elseif requested == "mmWE" && any(unit == ["m", "mWE"])
      % Canonical mass-balance depths are water equivalent even when a source
      % abbreviates the unit as metres.
      tt.value = tt.value .* 1000;
      unit = "mmWE";
   end
end

function tf = shouldConvertTemperature(varname, unit, values)
   %SHOULDCONVERTTEMPERATURE True for Kelvin temperature channels.

   is_temp_name = any(varname == ["tair", "tsfc", "tice10m"]) ...
      || startsWith(varname, "tice") ...
      || contains(lower(varname), "temperature");
   tf = is_temp_name && (unit == "K" || median(values, 'omitnan') > 100);
end

function tf = isIntervalObservation(T, varname)
   %ISINTERVALOBSERVATION True when values already represent interval totals.
   tf = false;
   metadata = T.Properties.UserData;
   if ~isstruct(metadata) || ~isfield(metadata, 'interval_observation') ...
         || ~logical(metadata.interval_observation)
      return
   end
   if isfield(metadata, 'interval_variables')
      tf = any(string(metadata.interval_variables) == varname);
   else
      tf = true;
   end
end

function [tt, unit] = dailyTotal(tt, unit, support)
   %DAILYTOTAL Integrate proven dense rates without relabelling sparse data.

   if height(tt) == 0
      return
   end

   % A sparse amount already represents its native observation support. A
   % sparse rate has no defensible represented duration, so keep its timestamp
   % but mask the requested total instead of inventing an integration interval.
   is_rate = any(unit == ["m s-1", "mWE/h"]);
   if ~support.dense
      if is_rate
         tt.value(:) = NaN;
      end
      return
   end

   values = tt.value;
   if unit == "m s-1"
      values = values .* support.cadence_s;
      unit = "m";
   elseif unit == "mWE/h"
      values = values .* support.cadence_s ./ 3600;
      unit = "mWE";
   end
   tt.value = values;
   tt = aggregateDenseDays(tt, support, "sum");
end

function tt = dailyMean(tt, support, method)
   %DAILYMEAN Aggregate dense series without presenting partial days as means.

   % Native daily/sparse observations must retain their actual timestamps. Only
   % a safely inferred subdaily series enters exact-grid daily aggregation.
   if ~support.dense
      return
   end
   if method == "mean_omitmissing"
      % Some variables (notably albedo) are structurally undefined at otherwise
      % present native timestamps. Their daily plot mean uses available values
      % but still requires the complete timestamp grid for that UTC day.
      method = "mean_omitmissing";
   elseif any(method == ...
         ["shortwave_weighted_mean", "shortwave_weighted_model_mean"])
      if ismember("weight", string(tt.Properties.VariableNames))
         if method == "shortwave_weighted_model_mean"
            method = "weighted_mean_available";
         else
            method = "weighted_mean_omitmissing";
         end
      else
         % Sparse albedo products such as MODIS have no collocated shortwave
         % weights. Preserve their finite-sample reduction as the safe fallback.
         method = "mean_omitmissing";
      end
   end
   tt = aggregateDenseDays(tt, support, method);
end

function support = inferDailySupport(time)
   %INFERDAILYSUPPORT Identify a safely evidenced dense subdaily source axis.

   support = struct('dense', false, 'cadence_s', NaN, 'phase_s', NaN, ...
      'expected_per_day', 0, 'has_reversed', false, 'tolerance_s', NaN);
   if numel(time) < 2 || any(isnat(time))
      return
   end

   % Absolute steps retain cadence evidence for a reversed dense axis, which
   % must fail closed rather than being misclassified as a sparse observation.
   steps_s = seconds(diff(time));
   nonzero_steps = abs(steps_s(isfinite(steps_s) & steps_s ~= 0));
   if isempty(nonzero_steps)
      return
   end
   cadence_s = median(nonzero_steps, 'omitnan');
   tolerance_s = max(1e-6, cadence_s * 1e-9);
   expected = 86400 / cadence_s;
   if ~isfinite(cadence_s) || cadence_s <= 0 || cadence_s >= 86400 ...
         || abs(expected - round(expected)) > tolerance_s
      return
   end

   % At least one day's sample count and repeated nominal steps distinguish a
   % dense source from an isolated set of subdaily observations.
   nominal = abs(nonzero_steps - cadence_s) <= tolerance_s;
   expected = round(expected);
   if numel(time) < expected || nnz(nominal) < 2
      return
   end

   support.dense = true;
   support.cadence_s = cadence_s;
   support.phase_s = dominantDailyPhase(time, cadence_s, tolerance_s);
   support.expected_per_day = expected;
   support.has_reversed = any(steps_s < -tolerance_s);
   support.tolerance_s = tolerance_s;
end

function phase_s = dominantDailyPhase(time, cadence_s, tolerance_s)
   %DOMINANTDAILYPHASE Infer the repeated within-day grid offset.

   day_start = dateshift(time, 'start', 'day');
   phase = mod(seconds(time - day_start), cadence_s);
   phase(abs(phase - cadence_s) <= tolerance_s) = 0;

   % Quantize only at the comparison tolerance so one off-grid posting cannot
   % move an otherwise stable native phase.
   phase = round(phase ./ tolerance_s) .* tolerance_s;
   phase_s = mode(phase);
end

function tt = aggregateDenseDays(tt, support, method)
   %AGGREGATEDENSEDAYS Summarize only exact complete native timestamp grids.

   time = utcTimes(tt.Time);
   first_day = dateshift(min(time), 'start', 'day');
   last_day = dateshift(max(time), 'start', 'day');
   day = (first_day:caldays(1):last_day)';
   value = nan(numel(day), 1);

   % Any reversed step invalidates the declared dense axis. The vectorized path
   % below keeps duplicate, missing, off-grid, and nonfinite failures local to
   % their UTC day without rebuilding a full-length time mask for every day.
   if support.has_reversed
      tt = timetable(value, 'RowTimes', day, 'VariableNames', {'value'});
      return
   end

   % Assign each row to one UTC day and one expected native-grid slot. Every
   % reducer requires each slot exactly once; strict reducers also require every
   % value, while omit-missing reducers require finite measurement support.
   % A daily radiometer albedo based on less than one quarter of the native
   % grid is too sensitive to the few low-sun samples at polar dawn/dusk, so
   % the shortwave-weighted reducer requires at least six represented hours.
   day_start = dateshift(time, 'start', 'day');
   day_index = round(days(day_start - first_day)) + 1;
   offset_s = seconds(time - day_start);
   slot = round((offset_s - support.phase_s) ./ support.cadence_s);
   expected_offset_s = support.phase_s + slot .* support.cadence_s;
   on_grid = slot >= 0 & slot < support.expected_per_day ...
      & abs(offset_s - expected_offset_s) <= support.tolerance_s;

   n_days = numel(day);
   n_slots = support.expected_per_day;
   finite = isfinite(tt.value);
   weighted = any(method == ...
      ["weighted_mean_omitmissing", "weighted_mean_available"]);
   circular = any(method == ...
      ["circular_mean_degrees", "circular_mean_degrees_omitmissing"]);
   circular_weighted = circular ...
      && ismember("weight", string(tt.Properties.VariableNames));
   if weighted
      weight = tt.weight;
      finite = finite & isfinite(weight) & weight > 0;
   elseif circular_weighted
      weight = tt.weight;
      finite = finite & isfinite(weight) & weight >= 0;
   end
   row_count = accumarray(day_index, ones(size(day_index)), [n_days, 1]);
   finite_count = zeros(n_days, 1);
   if any(finite)
      % Empty numeric accumarray indices are not portable inputs, so every
      % optional accumulator is guarded explicitly rather than relying on it.
      finite_count = accumarray(day_index(finite), ones(nnz(finite), 1), ...
         [n_days, 1]);
   end
   slot_index = (day_index(on_grid) - 1) .* n_slots + slot(on_grid) + 1;
   slot_count = zeros(n_days * n_slots, 1);
   if ~isempty(slot_index)
      slot_count = accumarray(slot_index, ones(size(slot_index)), ...
         [n_days * n_slots, 1]);
   end
   unique_grid = all(reshape(slot_count, n_slots, n_days) == 1, 1)';
   bad_step = find(diff(time) <= seconds(0)) + 1;
   bad_order = false(n_days, 1);
   if ~isempty(bad_step)
      bad_order = accumarray(day_index(bad_step), ...
         ones(numel(bad_step), 1), [n_days, 1]) > 0;
   end
   complete_grid = row_count == n_slots & unique_grid & ~bad_order;
   if method == "weighted_mean_available"
      % A modeled flux ratio is defined whenever the exact native day includes
      % at least one positive-SWD sample; do not apply a radiometer-duration
      % gate to deterministic model output.
      complete = complete_grid & finite_count > 0;
   elseif weighted
      minimum_finite = max(1, ceil(n_slots / 4));
      complete = complete_grid & finite_count >= minimum_finite;
   elseif any(method == ["mean_omitmissing", ...
         "circular_mean_degrees_omitmissing"])
      complete = complete_grid & finite_count > 0;
   else
      complete = complete_grid & finite_count == n_slots;
   end

   % Sum finite rows once per day; completeness masks all partial-day results.
   aggregate = zeros(n_days, 1);
   if any(finite)
      samples = tt.value(finite);
      if weighted
         samples = samples .* weight(finite);
      end
      aggregate = accumarray(day_index(finite), samples, ...
         [n_days, 1], @sum, 0);
   end
   if method == "mean"
      aggregate = aggregate ./ n_slots;
   elseif method == "mean_omitmissing"
      aggregate = aggregate ./ finite_count;
   elseif circular
      vector_weight = ones(nnz(finite), 1);
      if circular_weighted
         vector_weight = weight(finite);
      end
      sine = accumarray(day_index(finite), ...
         vector_weight .* sind(tt.value(finite)), ...
         [n_days, 1], @sum, 0);
      cosine = accumarray(day_index(finite), ...
         vector_weight .* cosd(tt.value(finite)), ...
         [n_days, 1], @sum, 0);
      resultant = hypot(sine, cosine);
      aggregate = mod(atan2d(sine, cosine), 360);
      aggregate(abs(aggregate) <= 1e-12) = 360;
      scale = n_slots .* ones(n_days, 1);
      if circular_weighted
         scale = accumarray(day_index(finite), weight(finite), ...
            [n_days, 1], @sum, 0);
      end
      aggregate(resultant <= 1e-12 .* scale) = NaN;
   elseif weighted
      denominator = accumarray(day_index(finite), weight(finite), ...
         [n_days, 1], @sum, 0);
      aggregate = aggregate ./ denominator;
   end
   value(complete) = aggregate(complete);
   tt = timetable(value, 'RowTimes', day, 'VariableNames', {'value'});
end

function time = utcTimes(time)
   %UTCTIMES Attach or convert row times to the daily-summary UTC convention.

   time.TimeZone = 'UTC';
end

function styleLegend(ax)
   %STYLELEGEND Keep legend colors readable under light and dark MATLAB themes.

   lgd = legend(ax, 'Location', 'best', 'Interpreter', 'none', ...
      'FontSize', 7, 'Color', 'w', 'TextColor', 'k', ...
      'EdgeColor', [0.75 0.75 0.75]);
   lgd.NumColumns = min(3, numel(lgd.String));
end

function label = ylabelText(varname, unit)
   %YLABELTEXT Compose a compact variable + unit y label.

   label = varname;
   if unit ~= ""
      label = varname + " [" + unit + "]";
   end
end
