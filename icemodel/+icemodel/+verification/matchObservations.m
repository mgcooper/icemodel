function matches = matchObservations(observations, candidate, kwargs)
   %MATCHOBSERVATIONS Match interval SMB and dated subsurface profiles.
   %
   %  matches = icemodel.verification.matchObservations( ...
   %     observations, candidate)
   %  matches = icemodel.verification.matchObservations( ...
   %     observations, candidate, variables=["density", "smb"], ...
   %     make_plot=true, plot_kind="profile")
   %
   % OBSERVATIONS and CANDIDATE may be canonical verification bundles with a
   % `.data` payload or the payload structs themselves. The returned
   % `intervals` and `profiles` tables retain physical identity, support,
   % observed/modelled values, variables, and units for metrics or focused
   % model-validation plots.

   arguments
      observations (1, 1) struct
      candidate
      kwargs.variables string = strings(0, 1)
      kwargs.make_plot (1, 1) logical = false
      kwargs.plot_kind (1, 1) string ...
         {mustBeMember(kwargs.plot_kind, ["auto", "profile", "interval"])} ...
         = "auto"
      kwargs.plot_id (1, 1) string = ""
      kwargs.plot_variable (1, 1) string = ""
      kwargs.visible (1, 1) string ...
         {mustBeMember(kwargs.visible, ["on", "off"])} = "off"
   end

   % Normalize full verification bundles to their heterogeneous data payloads.
   observation_data = bundleData(observations);
   candidate_data = bundleData(candidate);
   variables = resolveVariables(observation_data, kwargs.variables);

   % Collect variable results without repeatedly growing output tables.
   metric_cells = cell(numel(variables), 1);
   aligned_cells = cell(numel(variables), 1);
   interval_cells = cell(numel(variables), 1);
   profile_cells = cell(numel(variables), 1);
   for n = 1:numel(variables)
      [metric_cells{n}, aligned_cells{n}, interval_cells{n}, ...
         profile_cells{n}] = matchVariable( ...
         observation_data, candidate_data, variables(n));
   end

   % Empty-table prototypes keep the public schema stable when no match exists.
   metrics = concatenateMetrics(metric_cells);
   intervals = concatenateTables(interval_cells, intervalPrototype());
   profiles = concatenateTables(profile_cells, profilePrototype());
   aligned = cell2struct(aligned_cells, cellstr(variables), 1);

   % Focused plotting is opt-in and consumes the returned rows rather than
   % re-reading or re-matching the source bundles.
   figure_handle = gobjects(0);
   if kwargs.make_plot
      figure_handle = plotFocusedMatch(intervals, profiles, kwargs);
   end

   matches = struct( ...
      'intervals', intervals, ...
      'profiles', profiles, ...
      'metrics', metrics, ...
      'aligned', aligned, ...
      'figure', figure_handle);
end

function data = bundleData(bundle)
   %BUNDLEDATA Return the canonical data payload or an already-raw payload.
   data = bundle;
   if isstruct(bundle) && isfield(bundle, 'data')
      data = bundle.data;
   end
end

function variables = resolveVariables(observation_data, requested)
   %RESOLVEVARIABLES Return requested variables or every observation field.
   variables = reshape(string(requested), 1, []);
   variables = variables(strlength(variables) > 0);
   if ~isempty(variables)
      variables = unique(variables, 'stable');
      return
   end
   if ~isstruct(observation_data)
      error('icemodel:verification:matchObservations:badObservations', ...
         'observations must contain a struct data payload')
   end
   variables = reshape(string(fieldnames(observation_data)), 1, []);
   variables = variables(~ismember(variables, ["format", "metadata"]));
end

function [metrics, aligned, intervals, profiles] = matchVariable( ...
      observations, candidate, varname)
   %MATCHVARIABLE Match one observation variable by interval or profile support.
   [intervals, profiles] = deal(intervalPrototype(), profilePrototype());
   aligned = table();
   row = metricRow(varname, "");

   % Missing fields remain soft diagnostics so one unavailable model axis does
   % not abort comparison of the remaining observations.
   if ~isstruct(observations) || ~isfield(observations, varname)
      row.status = "missing_target_variable";
      metrics = struct2table(row);
      return
   end
   [has_candidate, candidate_value] = candidateValue(candidate, varname);
   if ~has_candidate
      row.status = "missing_candidate_variable";
      metrics = struct2table(row);
      return
   end

   observation_value = attachBundleDepth( ...
      observations, varname, observations.(char(varname)));
   candidate_value = attachBundleDepth(candidate, varname, candidate_value);
   if isIntervalTable(observation_value)
      [row, intervals, aligned] = matchIntervals( ...
         observation_value, candidate_value, varname);
      metrics = struct2table(row);
      return
   end

   % Profile variables can contain several identities and dates, so one metric
   % row is returned per physical observation profile.
   [rows, profiles, aligned] = matchProfiles( ...
      observation_value, candidate_value, varname);
   metrics = struct2table(rows);
end

function [tf, value] = candidateValue(candidate, varname)
   %CANDIDATEVALUE Resolve a variable from a struct or table candidate payload.
   tf = false;
   value = [];
   if isstruct(candidate) && isfield(candidate, varname)
      value = candidate.(char(varname));
      tf = true;
   elseif (istable(candidate) || istimetable(candidate)) ...
         && ismember(varname, string(candidate.Properties.VariableNames))
      value = candidate;
      tf = true;
   end
end

function tf = isIntervalTable(value)
   %ISINTERVALTABLE True for explicit start/end observation-support tables.
   tf = istable(value) ...
      && all(ismember(["start_date", "end_date"], ...
      string(value.Properties.VariableNames))) ...
      && isdatetime(value.start_date) && isdatetime(value.end_date);
end

function [row, intervals, aligned] = matchIntervals( ...
      observation_table, candidate_value, varname)
   %MATCHINTERVALS Integrate a model rate over complete observation intervals.
   row = metricRow(varname, "");
   intervals = intervalPrototype();
   aligned = table();
   observation_series = tableSeries(observation_table, [varname, "value"]);
   candidate_series = tableSeries(candidate_value, [varname, "value"]);
   if isempty(candidate_series.values)
      row.status = "missing_candidate_variable";
      return
   end

   [observed, modeled, raw_aligned] = ...
      icemodel.verification.helpers.alignObservationSeries( ...
      observation_series, candidate_series);
   row = populateMetric(row, observed, modeled);
   if isempty(observed)
      row.status = "no_overlap";
      return
   end

   % Recover the source row identity for each retained interval without
   % changing the public half-open interval semantics.
   source_rows = intervalSourceRows(observation_table, raw_aligned);
   observation_id = intervalIds(observation_table, source_rows, varname);
   variable = repmat(varname, numel(observed), 1);
   unit = repmat(observation_series.units, numel(observed), 1);
   intervals = table(observation_id, variable, raw_aligned.start_date, ...
      raw_aligned.end_date, observed, modeled, unit, ...
      'VariableNames', {'observation_id', 'variable', 'start_date', ...
      'end_date', 'observed', 'modeled', 'unit'});
   aligned = intervals;
end

function rows = intervalSourceRows(observation_table, aligned)
   %INTERVALSOURCEROWS Map aligned interval support back to source table rows.
   if ismember("source_row", string(aligned.Properties.VariableNames))
      rows = aligned.source_row;
      return
   end

   % Older aligned tables did not retain their source row explicitly.
   rows = zeros(height(aligned), 1);
   for n = 1:height(aligned)
      rows(n) = find(observation_table.start_date == aligned.start_date(n) ...
         & observation_table.end_date == aligned.end_date(n), 1);
   end
end

function ids = intervalIds(observation_table, rows, varname)
   %INTERVALIDS Preserve source ids and synthesize only when none is supplied.
   names = string(observation_table.Properties.VariableNames);
   preferred = ["observation_id", "measurement_id", "id"];
   hit = preferred(find(ismember(preferred, names), 1));
   if ~isempty(hit)
      ids = string(observation_table.(char(hit))(rows));
      return
   end

   % The variable plus exact interval is deterministic across repeated runs and
   % remains readable in metric/plot selection output.
   ids = varname + "|" ...
      + string(observation_table.start_date(rows), 'yyyy-MM-dd''T''HH:mm:ss') ...
      + "|" ...
      + string(observation_table.end_date(rows), 'yyyy-MM-dd''T''HH:mm:ss');
end

function [rows, profiles, aligned] = matchProfiles( ...
      observation_value, candidate_value, varname)
   %MATCHPROFILES Pair physical profile groups before depth interpolation.
   observation_groups = profileGroupList(observation_value, varname, ...
      "observation");
   candidate_groups = profileGroupList(candidate_value, varname, "candidate");
   if isempty(observation_groups)
      row = metricRow(varname, "");
      row.status = "missing_target_variable";
      rows = row;
      profiles = profilePrototype();
      aligned = table();
      return
   end

   % One metric row records each requested physical observation profile.
   rows = repmat(metricRow(varname, ""), numel(observation_groups), 1);
   profile_cells = cell(numel(observation_groups), 1);
   for n = 1:numel(observation_groups)
      group = observation_groups(n);
      rows(n).experiment = profileMatchId(group.id, group.datetime);
      [candidate_group, status] = selectCandidateGroup( ...
         group, candidate_groups);
      if status ~= "ok"
         rows(n).status = status;
         profile_cells{n} = profilePrototype();
         continue
      end

      observation_series = tableSeries(group.data, [varname, "value"]);
      candidate_series = tableSeries(candidate_group.data, [varname, "value"]);
      [observed, modeled, one_aligned] = ...
         icemodel.verification.helpers.alignObservationSeries( ...
         observation_series, candidate_series);
      rows(n) = populateMetric(rows(n), observed, modeled);
      if isempty(observed)
         rows(n).status = "no_overlap";
         profile_cells{n} = profilePrototype();
         continue
      end

      % Public rows retain observation identity/date and canonical target units.
      n_depths = numel(observed);
      match_id = repmat(profileMatchId(group.id, group.datetime), ...
         n_depths, 1);
      profile_id = repmat(group.id, n_depths, 1);
      variable = repmat(varname, n_depths, 1);
      profile_time = repmat(group.datetime, n_depths, 1);
      unit = repmat(observation_series.units, n_depths, 1);
      if ismember("depth", string(one_aligned.Properties.VariableNames))
         matched_depth = one_aligned.depth;
      else
         % Legacy profile bundles can carry a temperature timetable without a
         % depth column. Preserve its time-aligned values with unknown depth.
         matched_depth = NaN(n_depths, 1);
      end
      profile_cells{n} = table(match_id, profile_id, variable, profile_time, ...
         matched_depth, observed, modeled, unit, ...
         'VariableNames', {'match_id', 'profile_id', 'variable', ...
         'datetime', 'depth', 'observed', 'modeled', 'unit'});
   end

   profiles = concatenateTables(profile_cells, profilePrototype());
   aligned = profiles;
end

function groups = profileGroupList(value, varname, role)
   %PROFILEGROUPLIST Normalize dated tables and legacy numeric depth profiles.
   if istable(value) || istimetable(value)
      names = string(value.Properties.VariableNames);
      has_depth = depthVariable(names) ~= "";
      if ~has_depth || isLegacyDepthTimetable(value, names)
         % A no-depth table is a time series stored inside a legacy profile
         % bundle. A depth timetable with one distinct row time per depth is a
         % legacy single profile, not a dated tall profile history.
         groups = singleProfileGroup(value, role);
         return
      end
      if role == "candidate"
         % Candidate histories retain each model timestamp as one profile so
         % repeated depths on the same UTC date can never be pooled.
         groups = icemodel.verification.helpers.profileGroups( ...
            value, time_resolution="timestamp");
      else
         groups = icemodel.verification.helpers.profileGroups(value);
      end
      return
   end
   if isnumeric(value)
      % A legacy numeric profile has no identity/date; it remains one explicit
      % group and may be matched only as one profile.
      values = double(value);
      if ~isvector(values)
         error('icemodel:verification:matchObservations:undatedMatrix', ...
            '%s %s matrix requires dated table metadata', role, varname)
      end
      data = table((0:numel(values) - 1)', values(:), ...
         'VariableNames', {'depth', char(varname)});
      groups = struct('id', role, ...
         'datetime', NaT(1, 1, 'TimeZone', 'UTC'), 'data', data);
      return
   end
   groups = repmat(struct('id', "", ...
      'datetime', NaT(1, 1, 'TimeZone', 'UTC'), 'data', table()), 0, 1);
end

function tf = isLegacyDepthTimetable(value, names)
   %ISLEGACYDEPTHTIMETABLE True for one row-time-tagged undated depth profile.
   identity_names = ["profile_id", "name_key", "name"];
   has_identity = any(ismember(identity_names, names)) ...
      || all(ismember(["latitude", "longitude"], names));
   tf = istimetable(value) && ~has_identity ...
      && numel(unique(value.Properties.RowTimes)) == height(value);
end

function groups = singleProfileGroup(value, role)
   %SINGLEPROFILEGROUP Wrap one legacy time series or undated depth profile.
   groups = struct('id', role, ...
      'datetime', NaT(1, 1, 'TimeZone', 'UTC'), 'data', value);
end

function [selected, status] = selectCandidateGroup(observation, candidates)
   %SELECTCANDIDATEGROUP Resolve one unambiguous candidate identity and date.
   selected = struct('id', "", ...
      'datetime', NaT(1, 1, 'TimeZone', 'UTC'), 'data', table());
   status = "missing_candidate_date";
   if isempty(candidates)
      status = "missing_candidate_variable";
      return
   end

   % Dated observations require the same UTC calendar date. A single undated
   % legacy profile remains a valid reusable candidate for backward compatibility.
   candidate_dates = NaT(numel(candidates), 1, 'TimeZone', 'UTC');
   for n = 1:numel(candidates)
      candidate_dates(n) = dateshift(candidates(n).datetime, 'start', 'day');
   end
   observation_date = dateshift(observation.datetime, 'start', 'day');
   if isnat(observation_date)
      hits = true(numel(candidates), 1);
   else
      hits = candidate_dates == observation_date;
      if ~any(hits) && isscalar(candidates) ...
            && isnat(candidates(1).datetime)
         hits = true;
      end
   end
   candidates = candidates(hits);
   if isempty(candidates)
      return
   end
   if isscalar(candidates)
      selected = candidates;
      status = "ok";
      return
   end

   % Multiple physical identities on one date require an exact identity match.
   % A model history normally has one reusable identity, so keep every profile
   % with that identity for timestamp selection below.
   candidate_ids = reshape(string({candidates.id}), [], 1);
   if numel(unique(candidate_ids)) > 1
      exact_identity = candidate_ids == observation.id;
      if ~any(exact_identity)
         status = "ambiguous_candidate_profile";
         return
      end
      candidates = candidates(exact_identity);
   end
   if isscalar(candidates)
      selected = candidates;
      status = "ok";
      return
   end

   % Full model histories can contain several profiles on the accepted UTC
   % date. Prefer the exact observation timestamp, otherwise the unique nearest
   % timestamp; never concatenate repeated depths from separate model states.
   if isnat(observation.datetime)
      status = "ambiguous_candidate_profile";
      return
   end
   candidate_times = reshape([candidates.datetime], [], 1);
   exact_time = candidate_times == observation.datetime;
   if nnz(exact_time) == 1
      selected = candidates(exact_time);
      status = "ok";
      return
   end
   distance = abs(seconds(candidate_times - observation.datetime));
   nearest = distance == min(distance);
   if nnz(nearest) == 1
      selected = candidates(nearest);
      status = "ok";
   else
      status = "ambiguous_candidate_profile";
   end
end

function id = profileMatchId(profile_id, profile_time)
   %PROFILEMATCHID Compose a stable selector from identity and UTC date.
   id = string(profile_id);
   if ~isnat(profile_time)
      id = id + "|" + string(profile_time, 'yyyy-MM-dd');
   end
end

function series = tableSeries(value, aliases)
   %TABLESERIES Extract one numeric variable and its declared comparison axis.
   series = emptySeries();
   if ~(istable(value) || istimetable(value))
      return
   end
   names = string(value.Properties.VariableNames);
   hit = aliases(find(ismember(aliases, names), 1));
   if isempty(hit) || ~isnumeric(value.(char(hit)))
      return
   end

   series.variable = string(hit);
   series.units = seriesUnit(value, hit);
   series.values = double(value.(char(hit))(:));
   depth_hit = depthVariable(names);
   if depth_hit ~= "" && isnumeric(value.(char(depth_hit))) ...
         && numel(value.(char(depth_hit))) == numel(series.values)
      series.axis = double(value.(char(depth_hit))(:));
      series.axis_kind = "depth";
      return
   end
   if istimetable(value)
      series.axis = value.Properties.RowTimes;
      series.axis_kind = "time";
      return
   end
   if isIntervalTable(value)
      series.axis = [value.start_date(:), value.end_date(:)];
      series.axis_kind = "interval";
      return
   end
   datetime_hit = datetimeVariable(value, names);
   if datetime_hit ~= ""
      series.axis = value.(char(datetime_hit));
      series.axis_kind = "time";
   end
end

function value = attachBundleDepth(bundle, varname, value)
   %ATTACHBUNDLEDEPTH Preserve legacy outer depth metadata for numeric vectors.
   if ~isnumeric(value) || ~isvector(value) || ~isstruct(bundle) ...
         || ~isfield(bundle, 'depth') || numel(bundle.depth) ~= numel(value)
      return
   end

   % Convert only the legacy vector representation; dated tables retain their
   % own explicit depth and identity metadata.
   value = table(double(bundle.depth(:)), double(value(:)), ...
      'VariableNames', {'depth', char(varname)});
end

function series = emptySeries()
   %EMPTYSERIES Return one normalized observation/model series prototype.
   series = struct('values', [], 'axis', [], 'axis_kind', "", ...
      'variable', "", 'units', "");
end

function hit = datetimeVariable(value, names)
   %DATETIMEVARIABLE Return the first explicit datetime table variable.
   hit = "";
   preferred = ["datetime", "time", "Time", "date"];
   for name = preferred
      if ismember(name, names) && isdatetime(value.(char(name)))
         hit = name;
         return
      end
   end
end

function hit = depthVariable(names)
   %DEPTHVARIABLE Return the first recognized profile-depth column.
   preferred = ["depth", "midpoint", "z", "depth_m", "midpoint_m"];
   hit = preferred(find(ismember(preferred, names), 1));
   if isempty(hit)
      hit = "";
   end
end

function unit = seriesUnit(value, varname)
   %SERIESUNIT Return the table unit without inventing missing metadata.
   unit = "";
   if ~(istable(value) || istimetable(value))
      return
   end
   names = string(value.Properties.VariableNames);
   idx = find(names == varname, 1);
   if ~isempty(idx) && numel(value.Properties.VariableUnits) >= idx
      unit = string(value.Properties.VariableUnits{idx});
   end
end

function row = metricRow(varname, experiment)
   %METRICROW Return one canonical comparecase metric row.
   [names, defaults] = icemodel.verification.helpers.metricRowSchema();
   row = cell2struct(defaults, cellstr(names), 1);
   row.variable = string(varname);
   row.experiment = string(experiment);
end

function row = populateMetric(row, observed, modeled)
   %POPULATEMETRIC Compute canonical diagnostics from finite matched values.
   row.n = uint64(numel(observed));
   if isempty(observed)
      row.status = "no_overlap";
      return
   end
   delta = modeled - observed;
   row.bias = mean(delta);
   row.rmse = sqrt(mean(delta .^ 2));
   if numel(observed) > 1 && std(observed) > 0 && std(modeled) > 0
      C = corrcoef(observed, modeled);
      row.correlation = C(1, 2);
   end
   row.peak_target = max(observed);
   row.peak_candidate = max(modeled);
   row.peak_error = row.peak_candidate - row.peak_target;
   row.status = "ok";
end

function metrics = concatenateMetrics(cells)
   %CONCATENATEMETRICS Concatenate variable metric tables once.
   nonempty = ~cellfun(@isempty, cells);
   if ~any(nonempty)
      [names, defaults] = icemodel.verification.helpers.metricRowSchema();
      metrics = struct2table(cell2struct(defaults, cellstr(names), 1));
      metrics(1, :) = [];
      return
   end
   metrics = vertcat(cells{nonempty});
end

function result = concatenateTables(cells, prototype)
   %CONCATENATETABLES Concatenate nonempty tables while preserving one schema.
   nonempty = ~cellfun(@isempty, cells);
   if any(nonempty)
      result = vertcat(cells{nonempty});
   else
      result = prototype;
   end
end

function T = intervalPrototype()
   %INTERVALPROTOTYPE Return the stable public interval-match table schema.
   T = table(strings(0, 1), strings(0, 1), ...
      NaT(0, 1, 'TimeZone', 'UTC'), NaT(0, 1, 'TimeZone', 'UTC'), ...
      zeros(0, 1), zeros(0, 1), strings(0, 1), ...
      'VariableNames', {'observation_id', 'variable', 'start_date', ...
      'end_date', 'observed', 'modeled', 'unit'});
end

function T = profilePrototype()
   %PROFILEPROTOTYPE Return the stable public profile-match table schema.
   T = table(strings(0, 1), strings(0, 1), strings(0, 1), ...
      NaT(0, 1, 'TimeZone', 'UTC'), zeros(0, 1), zeros(0, 1), ...
      zeros(0, 1), strings(0, 1), ...
      'VariableNames', {'match_id', 'profile_id', 'variable', 'datetime', ...
      'depth', 'observed', 'modeled', 'unit'});
end

function f = plotFocusedMatch(intervals, profiles, kwargs)
   %PLOTFOCUSEDMATCH Plot one returned profile or interval through shared tools.
   kind = kwargs.plot_kind;
   if kind == "auto"
      if ~isempty(profiles)
         kind = "profile";
      elseif ~isempty(intervals)
         kind = "interval";
      else
         error('icemodel:verification:matchObservations:noPlottableMatch', ...
            'no matched profile or interval is available to plot')
      end
   end

   % Select exactly one returned match so focused validation never overlays
   % unrelated physical profiles or accumulation intervals.
   if kind == "profile"
      rows = selectProfileRows(profiles, kwargs.plot_id, ...
         kwargs.plot_variable);
      f = plotProfileRows(rows, kwargs.visible);
   else
      row = selectIntervalRow(intervals, kwargs.plot_id, ...
         kwargs.plot_variable);
      f = plotIntervalRow(row, kwargs.visible);
   end
end

function rows = selectProfileRows(profiles, requested_id, requested_variable)
   %SELECTPROFILEROWS Select one profile match id and variable.
   rows = profiles;
   if requested_variable ~= ""
      rows = rows(rows.variable == requested_variable, :);
   end
   if isempty(rows)
      error('icemodel:verification:matchObservations:noProfileMatch', ...
         'no matched profile satisfies the requested variable')
   end
   if requested_id == ""
      requested_id = rows.match_id(1);
   end
   rows = rows(rows.match_id == requested_id, :);
   if isempty(rows)
      error('icemodel:verification:matchObservations:noProfileMatch', ...
         'profile match id %s was not found', requested_id)
   end
   if numel(unique(rows.variable)) ~= 1
      error('icemodel:verification:matchObservations:ambiguousProfilePlot', ...
         'plot_variable is required when a profile id has multiple variables')
   end
end

function row = selectIntervalRow(intervals, requested_id, requested_variable)
   %SELECTINTERVALROW Select one interval observation id and variable.
   rows = intervals;
   if requested_variable ~= ""
      rows = rows(rows.variable == requested_variable, :);
   end
   if isempty(rows)
      error('icemodel:verification:matchObservations:noIntervalMatch', ...
         'no matched interval satisfies the requested variable')
   end
   if requested_id == ""
      row = rows(1, :);
      return
   end
   rows = rows(rows.observation_id == requested_id, :);
   if height(rows) ~= 1
      error('icemodel:verification:matchObservations:noIntervalMatch', ...
         'interval observation id %s did not select exactly one row', ...
         requested_id)
   end
   row = rows;
end

function f = plotProfileRows(rows, visible)
   %PLOTPROFILEROWS Render observed and modeled values on their matched depths.
   varname = rows.variable(1);
   observed = table(rows.depth, rows.observed, ...
      'VariableNames', {'depth', char(varname)});
   modeled = table(rows.depth, rows.modeled, ...
      'VariableNames', {'depth', char(varname)});
   observed.Properties.VariableUnits = {'m', char(rows.unit(1))};
   modeled.Properties.VariableUnits = {'m', char(rows.unit(1))};
   f = figure('Visible', char(visible), 'Color', 'white');
   ax = axes(f);
   icemodel.plot.profile({observed, modeled}, varname, ...
      names=["observed", "modeled"], axes=ax);
   title(ax, rows.match_id(1), 'Interpreter', 'none')
end

function f = plotIntervalRow(row, visible)
   %PLOTINTERVALROW Render observed and modeled accumulated interval values.
   Time = [row.start_date; row.end_date];
   observed = timetable(Time, repmat(row.observed, 2, 1), ...
      'VariableNames', {char(row.variable)});
   modeled = timetable(Time, repmat(row.modeled, 2, 1), ...
      'VariableNames', {char(row.variable)});
   observed.Properties.VariableUnits = {char(row.unit)};
   modeled.Properties.VariableUnits = {char(row.unit)};
   f = figure('Visible', char(visible), 'Color', 'white');
   ax = axes(f);
   icemodel.plot.compareTimeseries({observed, modeled}, row.variable, ...
      names=["observed", "modeled"], axes=ax);
   title(ax, row.observation_id, 'Interpreter', 'none')
end
