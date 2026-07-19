function result = comparecase(case_id, kwargs)
   %COMPARECASE Compare one staged verification target against a candidate.
   %
   %  result = icemodel.verification.comparecase("cdp")
   %  result = icemodel.verification.comparecase("wfj", make_plot=false)
   %  result = icemodel.verification.comparecase("colbeck1976", ...
   %     artifact_dir=fullfile(icemodel.getpath('test'), 'artifacts', 'tmp'))
   %
   % Inputs
   %  case_id                    Staged verification case id.
%  evaluation_data_root       Base evaluation-data root. When blank, the
%                             repo-local data/eval tree is used.
%  input_data_root            Optional paired input-data root for staged
%                             met/userdata artifacts.
%  icemodel_config_casename   Config casename used to resolve the default
%                             evaluation-data root without mutating config.
%  dataset_family             Optional family filter for shared case ids.
%  artifact_dir               Optional folder where metrics, aligned series,
   %                             and figures are written.
   %  candidate                  Optional in-memory candidate bundle.
   %  candidate_file             Optional MAT file containing `candidate` or
   %                             `reference`.
   %  make_plot                  Whether to create a comparison figure.
   %  save_plot                  Whether to export the comparison figure when
   %                             artifact_dir is provided.
   %  plot_visible               Figure visibility used for generated plots.
   %
   % Outputs
   %  result   Struct containing the resolved manifest, metric table, and any
   %           written artifact paths.
   %
   % Role
   %  This is a normal verification workflow entry point. It compares staged
   %  targets with either a supplied model candidate or the staged smoke
   %  reference.

   arguments
      case_id (1, :) string
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
      kwargs.dataset_family (1, 1) string = ""
      kwargs.artifact_dir (1, 1) string = ""
      kwargs.candidate = []
      kwargs.candidate_file (1, 1) string = ""
      kwargs.make_plot (1, 1) logical = true
      kwargs.save_plot (1, 1) logical = true
      kwargs.plot_visible (1, 1) string = "off"
   end

   % Resolve the case manifest first; all downstream artifact paths and
   % comparison variables come from this single contract.
   manifest = icemodel.verification.loadmanifest(case_id, ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "input_data_root", kwargs.input_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename, ...
      "dataset_family", kwargs.dataset_family);

   % Load the target and candidate bundles. With no candidate supplied, the
   % staged smoke reference is used so the suite runs before a snow model exists.
   % The eval target is a forcing-agnostic, data-only observations.mat bundle
   % referenced via evaluation_file/evaluation_path (ESM-SnowMIP, SUMup, and
   % freshly staged PROMICE cases); Snow/Colbeck cases load their evaluation.mat
   % the same way. Only PROMICE fixtures staged before the observations.mat
   % contract lack that file; they fall back to reconstituting the PROMICE-obs
   % target on demand from the staged per-year userdata files the manifest
   % declares.
   if isfield(manifest, 'evaluation_path') ...
         && strlength(string(manifest.evaluation_path)) > 0 ...
         && isfile(manifest.evaluation_path)
      targets = icemodel.verification.helpers.loadArtifact( ...
         manifest.evaluation_path, "targets");
   else
      targets = icemodel.verification.helpers.loadColocatedData( ...
         manifest, "promice");
   end

   % Some cases stage multiple target sources keyed under one evaluation.mat
   % (e.g. colbeck1976 carries numerical_summa and analytical_clark2017).
   % comparecase consumes one bundle at a time; pick the default
   % numerical_summa source. Use compareSolutions for the multi-axis driver.
   if ~isfield(targets, 'format') && isfield(targets, 'numerical_summa')
      targets = targets.numerical_summa;
   end

   candidate = icemodel.verification.helpers.resolveCandidateBundle(manifest, ...
      "candidate", kwargs.candidate, "candidate_file", kwargs.candidate_file);

   % Dispatch by artifact format rather than dataset family so future
   % verification families can reuse the same comparison entry point.
   switch targets.format
      case "timeseries"
         [metrics, aligned] = compareTimeseriesBundle( ...
            targets.data, candidate.data, manifest.comparison_variables);
      case "experiment_bundle"
         [metrics, aligned] = compareExperimentBundle( ...
            targets.experiments, candidate.experiments, ...
            manifest.comparison_variables);
      case "subsurface_profile_bundle"
         matches = icemodel.verification.matchObservations( ...
            targets.data, candidate.data, ...
            variables=string(manifest.comparison_variables));
         metrics = matches.metrics;
         aligned = matches.aligned;
      case "retmip_protocol_bundle"
         [metrics, aligned] = compareRetmipProtocolBundle( ...
            targets.data, candidate.data, manifest.comparison_variables);
      otherwise
         error('unsupported verification target format: %s', targets.format)
   end

   % Acceptance-gate mode is keyed by case_type, not by family. comparecase
   % itself only emits diagnostic metrics (per-variable status + bias/RMSE);
   % the hard PASS/FAIL tolerance gate lives in colbeck.compareSolutions for
   % synthetic_process verification cases. firn_observational and esm_site are
   % SOFT (diagnostic only) - they report metrics and never hard-fail here, so
   % a missing or noisy firn observation lane cannot break the suite.
   gate_mode = acceptanceGateMode(manifest);

   % Artifact writing is optional so the same function can serve interactive
   % checks and batch smoke-suite runs.
   [artifact_dir, metrics_path, figure_path, scatter_figure_path] = ...
      deal("", "", "", "");
   if ~isblanktext(kwargs.artifact_dir)
      artifact_dir = fullfile(kwargs.artifact_dir, case_id);
      icemodel.helpers.ensureDirExists(artifact_dir);

      metrics_path = fullfile(artifact_dir, "metrics.csv");
      writetable(metrics, metrics_path);

      save(fullfile(artifact_dir, "result.mat"), 'metrics', 'aligned', 'manifest');
   end

   % Plotting has two independent policies. make_plot controls whether figures
   % are created, save_plot controls PNG export, and plot_visible controls
   % whether figures remain open for interactive review.
   if kwargs.make_plot
      if kwargs.save_plot && ~isblanktext(artifact_dir)
         figure_path = fullfile(artifact_dir, "comparison.png");
         if targets.format == "timeseries"
            scatter_figure_path = fullfile(artifact_dir, "scatter.png");
         end
      end

      should_make_visible_plot = kwargs.plot_visible ~= "off";
      should_make_saved_plot = kwargs.save_plot && ~isblanktext(figure_path);

      if should_make_visible_plot || should_make_saved_plot

         f = icemodel.verification.plotcase( ...
            case_id, ...
            "evaluation_data_root", kwargs.evaluation_data_root, ...
            "input_data_root", kwargs.input_data_root, ...
            "icemodel_config_casename", kwargs.icemodel_config_casename, ...
            "dataset_family", kwargs.dataset_family, ...
            "source", "compare", ...
            "candidate", candidate, ...
            "visible", kwargs.plot_visible, ...
            "output_file", figure_path);

         if kwargs.plot_visible == "off"
            close(f)
         end
      end

      % Scatter plots are useful for site time series, but the current Colbeck
      % process benchmark is better reviewed as time-series/process panels until
      % analytical profile or flux targets are staged explicitly.
      should_make_saved_scatter = kwargs.save_plot ...
         && ~isblanktext(scatter_figure_path);
      if targets.format == "timeseries" ...
            && (should_make_visible_plot || should_make_saved_scatter)

         f = icemodel.verification.plotscatter( ...
            case_id, ...
               "evaluation_data_root", kwargs.evaluation_data_root, ...
               "input_data_root", kwargs.input_data_root, ...
               "icemodel_config_casename", kwargs.icemodel_config_casename, ...
               "dataset_family", kwargs.dataset_family, ...
               "candidate", candidate, ...
            "visible", kwargs.plot_visible, ...
            "output_file", scatter_figure_path);

         if kwargs.plot_visible == "off"
            close(f)
         end
      end
   end

   % Return paths even when blank so callers do not need field-existence checks.
   result = struct( ...
      'case_id', case_id, ...
      'manifest', manifest, ...
      'metrics', metrics, ...
      'aligned', aligned, ...
      'gate_mode', gate_mode, ...
      'artifact_dir', artifact_dir, ...
      'metrics_path', metrics_path, ...
      'figure_path', figure_path, ...
      'scatter_figure_path', scatter_figure_path);
end

function gate_mode = acceptanceGateMode(manifest)
   %ACCEPTANCEGATEMODE Classify a case's acceptance gate by case_type.
   %
   % "hard"        synthetic_process - formal RMSE tolerances apply (Colbeck;
   %               gated in colbeck.compareSolutions, not here).
   % "soft"        esm_site / firn_observational - diagnostic metrics only,
   %               no hard PASS/FAIL. This is the firn observational contract:
   %               report bias/RMSE/correlation, never gate the suite.
   case_type = "";
   if isfield(manifest, "case_type")
      case_type = string(manifest.case_type);
   end
   switch case_type
      case "synthetic_process"
         gate_mode = "hard";
      otherwise
         gate_mode = "soft";
   end
end

function [metrics, aligned] = compareTimeseriesBundle(target_tt, candidate_tt, ...
      variable_names)
   %COMPARETIMESERIESBUNDLE Compare one target/candidate timetable pair.

   % Preallocate one metric row and one aligned-series slot per variable.
   n_vars = numel(variable_names);
   rows = repmat(metricRowTemplate(), n_vars, 1);
   aligned_cells = cell(n_vars, 1);

   % Compare variables independently so missing-variable status is reported per
   % variable instead of aborting the whole case.
   for n = 1:n_vars
      [rows(n), aligned_cells{n}] = compareOneVariable(target_tt, candidate_tt, ...
         variable_names(n), "");
   end

   metrics = struct2table(rows);
   aligned = cell2struct(aligned_cells, variable_names, 1);
end

function [metrics, aligned] = compareExperimentBundle(targets, candidates, ...
      variable_names)
   %COMPAREEXPERIMENTBUNDLE Compare a struct of experiment timetables.

   % Experiment bundles compare matching named experiments, such as Colbeck
   % exp1/exp2/exp3. The field order is treated as part of the contract.
   [target_names, target_values] = deal(fieldnames(targets), ...
      struct2cell(targets));
   [candidate_names, candidate_values] = deal(fieldnames(candidates), ...
      struct2cell(candidates));

   if ~isequal(target_names, candidate_names)
      error('candidate experiment names do not match target experiment names')
   end

   % Preallocate the rectangular experiment-by-variable metric table.
   n_rows = numel(target_names) * numel(variable_names);
   rows = repmat(metricRowTemplate(), n_rows, 1);
   aligned_groups = cell(numel(target_names), 1);
   row_idx = 0;

   % Keep the nested loop explicit so the metric table order matches the
   % plotted experiment grid and remains easy to inspect.
   for n = 1:numel(target_names)
      sync_rows = cell(numel(variable_names), 1);
      for m = 1:numel(variable_names)
         row_idx = row_idx + 1;
         [rows(row_idx), sync_rows{m}] = compareOneVariable( ...
            target_values{n}, candidate_values{n}, variable_names(m), ...
            target_names{n});
      end
      aligned_groups{n} = cell2struct(sync_rows, variable_names, 1);
   end

   metrics = struct2table(rows);
   aligned = cell2struct(aligned_groups, target_names, 1);
end

function [metrics, aligned] = compareRetmipProtocolBundle(targets, candidates, ...
      variable_names)
   %COMPARERETMIPPROTOCOLBUNDLE Compare RetMIP protocol surface/profile axes.

   n_vars = numel(variable_names);
   rows = repmat(metricRowTemplate(), n_vars, 1);
   aligned_cells = cell(n_vars, 1);
   for n = 1:n_vars
      [rows(n), aligned_cells{n}] = compareOneRetmipProtocolVariable( ...
         targets, candidates, variable_names(n));
   end

   metrics = struct2table(rows);
   aligned = cell2struct(aligned_cells, variable_names, 1);
end

function [row, aligned] = compareOneRetmipProtocolVariable(targets, ...
      candidates, varname)
   %COMPAREONERETMIPPROTOCOLVARIABLE Compare one RetMIP protocol axis.
   row = metricRowTemplate();
   row.variable = varname;
   row.experiment = "";
   aligned = table();

   target_series = retmipProtocolSeries(targets, varname);
   if isempty(target_series.values)
      row.status = "missing_target_variable";
      return
   end
   candidate_series = retmipProtocolSeries(candidates, varname);
   if isempty(candidate_series.values)
      row.status = "missing_candidate_variable";
      return
   end

   [target, candidate, aligned] = alignComparableSeries( ...
      target_series, candidate_series);

   row.n = uint64(numel(target));
   if isempty(target)
      row.status = "no_overlap";
      return
   end

   delta = candidate - target;
   row.bias = mean(delta);
   row.rmse = sqrt(mean(delta .^ 2));
   if numel(target) > 1 && std(target) > 0 && std(candidate) > 0
      C = corrcoef(target, candidate);
      row.correlation = C(1, 2);
   end
   row.peak_target = max(target);
   row.peak_candidate = max(candidate);
   row.peak_error = row.peak_candidate - row.peak_target;
   row.status = "ok";
end

function series = retmipProtocolSeries(bundle, varname)
   %RETMIPPROTOCOLSERIES Extract one RetMIP axis with its comparison coordinate.
   series = emptyComparableSeries();
   if ~isstruct(bundle)
      return
   end
   if isfield(bundle, 'surface')
      series = tableVariableSeries(bundle.surface, retmipAliases(varname));
   end
   if ~isempty(series.values) || ~isfield(bundle, 'profiles') ...
         || ~isstruct(bundle.profiles)
      return
   end
   profiles = bundle.profiles;
   names = fieldnames(profiles);
   for k = 1:numel(names)
      series = tableVariableSeries(profiles.(names{k}), retmipAliases(varname));
      if ~isempty(series.values)
         return
      end
   end
end

function aliases = retmipAliases(varname)
   %RETMIPALIASES Return RetMIP protocol/profile column aliases.
   switch string(varname)
      case "density"
         aliases = ["density", "density_kgm3", "rho"];
      case "subsurface_temperature"
         aliases = ["subsurface_temperature", "temp", "T"];
      case "lwc"
         aliases = ["lwc", "LWC", "liquid_water_content"];
      otherwise
         aliases = string(varname);
   end
end

function series = tableVariableSeries(tbl, aliases)
   %TABLEVARIABLESERIES Extract a numeric table column and its coordinate axis.
   series = emptyComparableSeries();
   if ~(istable(tbl) || istimetable(tbl))
      return
   end
   names = string(tbl.Properties.VariableNames);
   hit = aliases(find(ismember(aliases, names), 1));
   if isempty(hit) || ~isnumeric(tbl.(char(hit)))
      return
   end
   idx = find(names == hit, 1, 'first');
   series.variable = string(hit);
   if ~isempty(idx) && numel(tbl.Properties.VariableUnits) >= idx
      series.units = string(tbl.Properties.VariableUnits{idx});
   end
   series.values = double(tbl.(char(hit))(:));
   depth_hit = depthTableVariable(names);
   if strlength(depth_hit) > 0 && isnumeric(tbl.(char(depth_hit))) ...
         && numel(tbl.(char(depth_hit))) == numel(series.values)
      % Profile bundles can be timetables with a depth column; compare those
      % by profile depth, not by coincidental observation timestamps.
      series.axis = double(tbl.(char(depth_hit))(:));
      series.axis_kind = "depth";
      return
   end
   if istimetable(tbl)
      series.axis = tbl.Properties.RowTimes;
      series.axis_kind = "time";
      return
   end
   interval_hit = intervalTableVariables(tbl, names);
   if all(strlength(interval_hit) > 0)
      series.axis = [tbl.(char(interval_hit(1)))(:), ...
         tbl.(char(interval_hit(2)))(:)];
      series.axis_kind = "interval";
      return
   end
   time_hit = datetimeTableVariable(tbl, names);
   if strlength(time_hit) > 0
      series.axis = tbl.(char(time_hit));
      series.axis_kind = "time";
      return
   end
end

function series = emptyComparableSeries()
   %EMPTYCOMPARABLESERIES Prototype for aligned soft-comparison values.
   series = struct('values', [], 'axis', [], 'axis_kind', "", ...
      'variable', "", 'units', "");
end

function [target, candidate, aligned] = alignComparableSeries( ...
      target_series, candidate_series)
   %ALIGNCOMPARABLESERIES Align values by their declared comparison axis.
   [target, candidate, aligned] = ...
      icemodel.verification.helpers.alignObservationSeries( ...
      target_series, candidate_series);
end

function hit = intervalTableVariables(tbl, names)
   %INTERVALTABLEVARIABLES Return start/end datetime columns for intervals.
   hit = strings(1, 2);
   if all(ismember(["start_date", "end_date"], names)) ...
         && isdatetime(tbl.start_date) && isdatetime(tbl.end_date)
      hit = ["start_date", "end_date"];
   end
end

function hit = datetimeTableVariable(tbl, names)
   %DATETIMETABLEVARIABLE Return the first datetime table variable.
   hit = "";
   preferred = ["time", "Time"];
   preferred_hit = preferred(find(ismember(preferred, names), 1));
   if ~isempty(preferred_hit) && isdatetime(tbl.(char(preferred_hit)))
      hit = preferred_hit;
      return
   end
   for n = 1:numel(names)
      if isdatetime(tbl.(char(names(n))))
         hit = names(n);
         return
      end
   end
end

function hit = depthTableVariable(names)
   %DEPTHTABLEVARIABLE Return the first recognized profile-depth column.
   depth_names = ["depth", "midpoint", "z", "depth_m", "midpoint_m"];
   hit = depth_names(find(ismember(depth_names, names), 1));
   if isempty(hit)
      hit = "";
   end
end

function [row, sync_tt] = compareOneVariable(target_tt, candidate_tt, ...
      varname, experiment)
   %COMPAREONEVARIABLE Compute one metric row and aligned series pair.

   % Start with the canonical row schema so every exit path returns the same
   % table columns.
   row = metricRowTemplate();
   row.experiment = string(experiment);
   row.variable = varname;

   % Missing variables are a per-variable diagnostic, not a fatal case error.
   if ~ismember(varname, target_tt.Properties.VariableNames)
      row.status = "missing_target_variable";
      sync_tt = timetable();
      return
   end
   if ~ismember(varname, candidate_tt.Properties.VariableNames)
      row.status = "missing_candidate_variable";
      sync_tt = timetable();
      return
   end

   % Align on common timestamps and drop non-finite values before computing
   % scalar diagnostics.
   sync_tt = synchronize(target_tt(:, varname), candidate_tt(:, varname), ...
      'intersection', 'first');
   sync_tt.Properties.VariableNames = {'target', 'candidate'};
   valid = isfinite(sync_tt.target) & isfinite(sync_tt.candidate);
   sync_tt = sync_tt(valid, :);

   row.n = uint64(height(sync_tt));
   if height(sync_tt) == 0
      row.status = "no_overlap";
      return
   end

   % Core all-case metrics: bias, RMSE, optional correlation, peak amplitude,
   % and peak timing.
   delta = sync_tt.candidate - sync_tt.target;
   row.bias = mean(delta);
   row.rmse = sqrt(mean(delta .^ 2));
   if height(sync_tt) > 1 && std(sync_tt.target) > 0 && std(sync_tt.candidate) > 0
      C = corrcoef(sync_tt.target, sync_tt.candidate);
      row.correlation = C(1, 2);
   end

   [row.peak_target, idx_target] = max(sync_tt.target);
   [row.peak_candidate, idx_candidate] = max(sync_tt.candidate);
   row.peak_error = row.peak_candidate - row.peak_target;
   row.peak_time_error_hours = hours(sync_tt.Time(idx_candidate) - ...
      sync_tt.Time(idx_target));

   % Snow-depth and SWE series get two additional timing diagnostics: when
   % the series first rises above an onset threshold (snow accumulation
   % start), and when the post-peak series first returns below the same
   % threshold (melt-out). Both errors are signed (candidate - target) in
   % hours; an unobserved onset / melt-out leaves the metric as NaN.
   if any(varname == ["snow_depth_m", "swe_kg_m2"])
      threshold = meltThreshold(varname);
      onset_target = firstAboveThreshold(sync_tt.target, sync_tt.Time, threshold);
      onset_candidate = firstAboveThreshold(sync_tt.candidate, sync_tt.Time, threshold);
      if ~isnat(onset_target) && ~isnat(onset_candidate)
         row.snow_onset_time_error_hours = hours(onset_candidate - onset_target);
      end

      melt_target = firstBelowThresholdAfterPeak(sync_tt.target, ...
         sync_tt.Time, threshold);
      melt_candidate = firstBelowThresholdAfterPeak(sync_tt.candidate, ...
         sync_tt.Time, threshold);
      if ~isnat(melt_target) && ~isnat(melt_candidate)
         row.melt_out_time_error_hours = hours(melt_candidate - melt_target);
      end
   end
end

function row = metricRowTemplate()
   %METRICROWTEMPLATE Return one empty metric row with canonical fields.

   [names, defaults] = icemodel.verification.helpers.metricRowSchema();
   row = cell2struct(defaults, cellstr(names), 1);
end

function threshold = meltThreshold(varname)
   %MELTTHRESHOLD Return a variable-specific melt-out threshold.

   switch varname
      case 'snow_depth_m'
         threshold = 0.01;
      case 'swe_kg_m2'
         threshold = 1.0;
      otherwise
         threshold = 0.0;
   end
end

function t = firstBelowThresholdAfterPeak(values, time, threshold)
   %FIRSTBELOWTHRESHOLDAFTERPEAK Return the first post-peak time below threshold.
   % This is the canonical melt-out diagnostic.

   t = NaT;
   if isempty(values)
      return
   end

   [~, peak_idx] = max(values);
   idx = find(values(peak_idx:end) <= threshold, 1, 'first');
   if isempty(idx)
      return
   end
   t = time(peak_idx + idx - 1);
end

function t = firstAboveThreshold(values, time, threshold)
   %FIRSTABOVETHRESHOLD Return the first time the series rises above threshold.
   % This is the canonical snow-onset diagnostic.

   t = NaT;
   if isempty(values)
      return
   end
   idx = find(values > threshold, 1, 'first');
   if isempty(idx)
      return
   end
   t = time(idx);
end
