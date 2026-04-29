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
   %                             path is resolved from icemodel.config.
   %  icemodel_config_casename   Config casename used to resolve the default
   %                             evaluation-data root without mutating config.
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
      kwargs.icemodel_config_casename (1, 1) string = "test"
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
      "icemodel_config_casename", kwargs.icemodel_config_casename);

   % Load the target and candidate bundles. With no candidate supplied, the
   % staged smoke reference is used so the suite runs before a snow model exists.
   targets = icemodel.verification.helpers.loadArtifact( ...
      manifest.evaluation_path, "targets");

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
      otherwise
         error('unsupported verification target format: %s', targets.format)
   end

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
            "icemodel_config_casename", kwargs.icemodel_config_casename, ...
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
            "icemodel_config_casename", kwargs.icemodel_config_casename, ...
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
      'artifact_dir', artifact_dir, ...
      'metrics_path', metrics_path, ...
      'figure_path', figure_path, ...
      'scatter_figure_path', scatter_figure_path);
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
   for i = 1:n_vars
      [rows(i), aligned_cells{i}] = compareOneVariable(target_tt, candidate_tt, ...
         variable_names(i), "");
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
   for i = 1:numel(target_names)
      sync_rows = cell(numel(variable_names), 1);
      for j = 1:numel(variable_names)
         row_idx = row_idx + 1;
         [rows(row_idx), sync_rows{j}] = compareOneVariable( ...
            target_values{i}, candidate_values{i}, variable_names(j), ...
            target_names{i});
      end
      aligned_groups{i} = cell2struct(sync_rows, variable_names, 1);
   end

   metrics = struct2table(rows);
   aligned = cell2struct(aligned_groups, target_names, 1);
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

   % Snow-depth and SWE series get one additional timing metric: when the
   % post-peak series first returns to a near-zero melt-out threshold.
   if any(varname == ["snow_depth_m", "swe_kg_m2"])
      threshold = meltThreshold(varname);
      target_time = firstBelowThreshold(sync_tt.target, sync_tt.Time, threshold);
      candidate_time = firstBelowThreshold(sync_tt.candidate, sync_tt.Time, threshold);
      if ~isnat(target_time) && ~isnat(candidate_time)
         row.melt_out_time_error_hours = hours(candidate_time - target_time);
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

function t = firstBelowThreshold(values, time, threshold)
   %FIRSTBELOWTHRESHOLD Return the first post-peak time below threshold.

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
