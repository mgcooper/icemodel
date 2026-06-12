function result = compareSolutions(kwargs)
   %COMPARESOLUTIONS Compare cached + computed Colbeck solutions side-by-side.
   %
   %  result = icemodel.verification.colbeck.compareSolutions()
   %  result = icemodel.verification.colbeck.compareSolutions(make_plot=true)
   %
   %  Single comparison driver for the Colbeck 1976 case. Loads the two cached
   %  target sources (numerical_summa, analytical_clark2017) from the unified
   %  colbeck1976/evaluation.mat, runs the IceModel candidate in both numerical
   %  and analytical kinds, and computes the four pairwise comparison rows per
   %  experiment. Two pairs are formal verification axes; the other two are
   %  diagnostic.
   %
   %  Pairwise axes (per experiment, per variable):
   %    formal:
   %      numerical_summa       vs numerical_icemodel
   %      analytical_clark2017  vs numerical_icemodel
   %    diagnostic:
   %      numerical_summa       vs analytical_icemodel
   %      analytical_clark2017  vs analytical_icemodel
   %
   %  Outputs
   %    result   Struct with fields:
   %      targets        Struct keyed numerical_summa / analytical_clark2017
   %      candidates     Struct keyed numerical_icemodel / analytical_icemodel
   %      metrics_table  Long-format table (axis_role, target_source,
   %                     candidate_source, experiment, variable, n, bias,
   %                     rmse, correlation, status)
   %      summary        Pass/fail printout from formal axes
   %
   %  See also: icemodel.verification.colbeck.runCase,
   %    icemodel.verification.colbeck.analyticalSolution

   arguments
      kwargs.artifact_dir (1, 1) string = ""
      kwargs.experiment_names (1, :) string = ["exp1", "exp2", "exp3"]
      kwargs.make_plot (1, 1) logical = true
      kwargs.save_plot (1, 1) logical = false
      kwargs.plot_visible (1, 1) string = "off"
      kwargs.tolerance_storage_m (1, 1) double = 5e-3    % m
      kwargs.tolerance_outflow_mps (1, 1) double = 5e-7  % m s-1 (5% of q_top=1e-5)
   end

   % Per-variable RMSE tolerances drive the formal pass/fail decision.
   tolerance = struct( ...
      'snow_liquid_water_storage_m', kwargs.tolerance_storage_m, ...
      'bottom_outflow_mps',          kwargs.tolerance_outflow_mps);

   % Load the case manifest and the cached target bundle (multi-source).
   manifest = icemodel.verification.loadmanifest("colbeck1976");
   targets = loadTargets(manifest);

   % Run the IceModel candidate in both numerical and analytical kinds.
   candidates = runCandidates(manifest, kwargs.experiment_names);

   % Four pairwise comparison axes: (target_source, candidate_source, role).
   pairs = { ...
      "numerical_summa",      "numerical_icemodel",  "formal"
      "analytical_clark2017", "numerical_icemodel",  "formal"
      "numerical_summa",      "analytical_icemodel", "diagnostic"
      "analytical_clark2017", "analytical_icemodel", "diagnostic"};

   % Long-format metrics table over (pair, experiment, variable).
   variables = string(manifest.comparison_variables);
   metrics_table = computeMetricTable(targets, candidates, pairs, ...
      variables, kwargs.experiment_names, manifest.case_id);

   % Pass/fail printout from the two formal axes.
   summary = printSummary(metrics_table, tolerance);

   % Optional CSV artifact alongside the figure.
   csv_path = "";
   if ~isblanktext(kwargs.artifact_dir)
      icemodel.helpers.ensureDirExists(kwargs.artifact_dir);
      csv_path = fullfile(kwargs.artifact_dir, "compareSolutions_metrics.csv");
      writetable(metrics_table, csv_path);
   end

   % Optional 4-way overlay figure.
   figure_path = "";
   if kwargs.make_plot
      [f, figure_path] = plotSolutions(targets, candidates, ...
         kwargs.experiment_names, kwargs.artifact_dir, ...
         kwargs.save_plot, kwargs.plot_visible);
      if kwargs.plot_visible == "off"
         close(f)
      end
   end

   result = struct( ...
      'targets',       targets, ...
      'candidates',    candidates, ...
      'metrics_table', metrics_table, ...
      'summary',       summary, ...
      'figure_path',   figure_path, ...
      'csv_path',      csv_path);
end

%% Targets: load the cached multi-source bundle from evaluation.mat.
function targets = loadTargets(manifest)

   % The staged artifact stores a struct keyed by target_source name; load and
   % return it as-is so downstream code can index by source.
   targets = icemodel.verification.helpers.loadArtifact( ...
      manifest.evaluation_path, "targets");
end

%% Candidates: run IceModel for both numerical and analytical kinds.
function candidates = runCandidates(manifest, experiment_names)

   % Each kind produces an experiment_bundle with the same schema, so the
   % comparison driver treats them uniformly downstream.
   candidates = struct( ...
      'numerical_icemodel',  icemodel.verification.colbeck.runCase( ...
      manifest, kind="numerical", experiment_names=experiment_names), ...
      'analytical_icemodel', icemodel.verification.colbeck.runCase( ...
      manifest, kind="analytical", experiment_names=experiment_names));
end

%% Metric table: preallocated long-format rows over (pair, experiment, var).
function metrics_table = computeMetricTable(targets, candidates, pairs, ...
      variables, experiment_names, case_id)

   % Preallocate the row block so the inner loops can assign by index.
   n_pairs = size(pairs, 1);
   n_exp   = numel(experiment_names);
   n_var   = numel(variables);
   n_rows  = n_pairs * n_exp * n_var;

   status           = strings(n_rows, 1);
   n                = zeros(n_rows, 1);
   bias             = nan(n_rows, 1);
   rmse             = nan(n_rows, 1);
   correlation      = nan(n_rows, 1);
   target_source    = strings(n_rows, 1);
   candidate_source = strings(n_rows, 1);
   axis_role        = strings(n_rows, 1);
   case_id_col      = strings(n_rows, 1);
   experiment       = strings(n_rows, 1);
   variable         = strings(n_rows, 1);

   % Fill the preallocated columns row-by-row.
   r = 0;
   for p = 1:n_pairs
      ts = string(pairs{p, 1});
      cs = string(pairs{p, 2});
      ar = string(pairs{p, 3});
      target_bundle    = targets.(char(ts));
      candidate_bundle = candidates.(char(cs));
      for i = 1:n_exp
         exp_name = experiment_names(i);
         target_tt    = target_bundle.experiments.(char(exp_name));
         candidate_tt = candidate_bundle.experiments.(char(exp_name));
         for v = 1:n_var
            varname = variables(v);
            row = computeOnePair(target_tt, candidate_tt, varname);

            r = r + 1;
            status(r)           = row.status;
            n(r)                = row.n;
            bias(r)             = row.bias;
            rmse(r)             = row.rmse;
            correlation(r)      = row.correlation;
            target_source(r)    = ts;
            candidate_source(r) = cs;
            axis_role(r)        = ar;
            case_id_col(r)      = case_id;
            experiment(r)       = exp_name;
            variable(r)         = varname;
         end
      end
   end

   % Assemble the long-format table with the column order callers expect.
   metrics_table = table(status, n, bias, rmse, correlation, ...
      target_source, candidate_source, axis_role, case_id_col, ...
      experiment, variable, ...
      'VariableNames', {'status', 'n', 'bias', 'rmse', 'correlation', ...
      'target_source', 'candidate_source', 'axis_role', 'case_id', ...
      'experiment', 'variable'});
end

%% Per-pair statistics: synchronize then compute bias / rmse / correlation.
function row = computeOnePair(target_tt, candidate_tt, varname)

   % Bail out if the variable is missing from either timetable.
   if ~ismember(varname, string(target_tt.Properties.VariableNames)) || ...
         ~ismember(varname, string(candidate_tt.Properties.VariableNames))
      row = struct('status', "not_applicable", 'n', 0, ...
         'bias', NaN, 'rmse', NaN, 'correlation', NaN);
      return
   end

   % Align target and candidate on the target's time grid; nearest-neighbor
   % matching handles candidates sampled on slightly different cadences.
   sync = synchronize(target_tt(:, varname), candidate_tt(:, varname), ...
      'first', 'nearest');
   target_vals = sync{:, 1};
   cand_vals   = sync{:, 2};

   % Drop sample pairs where either side is non-finite.
   ok_pair = isfinite(target_vals) & isfinite(cand_vals);
   target_vals = target_vals(ok_pair);
   cand_vals   = cand_vals(ok_pair);
   n = numel(target_vals);

   if n == 0
      row = struct('status', "not_applicable", 'n', 0, ...
         'bias', NaN, 'rmse', NaN, 'correlation', NaN);
      return
   end

   % Standard scalar comparison metrics.
   bias = mean(cand_vals - target_vals);
   rmse = sqrt(mean((cand_vals - target_vals) .^ 2));
   if std(target_vals) > 0 && std(cand_vals) > 0 && n > 1
      cc = corrcoef(target_vals, cand_vals);
      correlation = cc(1, 2);
   else
      correlation = NaN;
   end
   row = struct('status', "ok", 'n', n, ...
      'bias', bias, 'rmse', rmse, 'correlation', correlation);
end

%% Summary: pass/fail printout from the formal axes only.
function summary = printSummary(metrics_table, tolerance)

   % Restrict to formal axes; diagnostic rows do not gate pass/fail.
   formal = metrics_table(metrics_table.axis_role == "formal", :);
   tol = arrayfun(@(v) tolerance.(char(v)), formal.variable);
   pass = formal.status == "ok" & formal.rmse < tol;

   % Header with the active per-variable RMSE tolerances.
   fprintf("\n=== Colbeck 1976 verification summary ===\n");
   fprintf("Formal axes (per-variable RMSE tolerance):\n");
   fprintf("  storage  < %.3g m\n", tolerance.snow_liquid_water_storage_m);
   fprintf("  outflow  < %.3g m s-1\n", tolerance.bottom_outflow_mps);

   % One line per formal row with PASS/FAIL marker and metric values.
   for n = 1:height(formal)
      mark = "FAIL";
      if pass(n)
         mark = "PASS";
      end
      fprintf("  [%s] %s vs %s @ %s %s: rmse=%.3e bias=%.3e n=%d\n", ...
         mark, formal.target_source(n), formal.candidate_source(n), ...
         formal.experiment(n), formal.variable(n), ...
         formal.rmse(n), formal.bias(n), formal.n(n));
   end
   summary = struct('formal', formal, 'all_pass', all(pass));
   fprintf("Overall: %s\n\n", string(ternary(summary.all_pass, "PASS", "FAIL")));
end

function out = ternary(cond, a, b)
   if cond
      out = a;
   else
      out = b;
   end
end

%% Plot: 4-way overlay of target + candidate series per (experiment, variable).
function [f, figure_path] = plotSolutions(targets, candidates, ...
      experiment_names, artifact_dir, save_plot, plot_visible)

   variables = ["snow_liquid_water_storage_m"; "bottom_outflow_mps"];
   n_exp = numel(experiment_names);
   n_var = numel(variables);

   % One tile per (experiment, variable) pair.
   f = figure('Visible', plot_visible, 'Position', [100, 100, 1200, 800]);
   tl = tiledlayout(f, n_exp, n_var, 'TileSpacing', 'compact', ...
      'Padding', 'compact');

   for n = 1:n_exp
      name = char(experiment_names(n));
      for m = 1:n_var
         v = char(variables(m));
         ax = nexttile(tl);
         hold(ax, 'on');

         % Overlay the two cached targets and the two IceModel candidates.
         plotOneSeries(ax, targets.numerical_summa.experiments.(name), ...
            v, 'numerical\_summa', '-.');
         plotOneSeries(ax, targets.analytical_clark2017.experiments.(name), ...
            v, 'analytical\_clark2017', '--');
         plotOneSeries(ax, candidates.numerical_icemodel.experiments.(name), ...
            v, 'numerical\_icemodel', '-');
         plotOneSeries(ax, candidates.analytical_icemodel.experiments.(name), ...
            v, 'analytical\_icemodel', ':');
         hold(ax, 'off');

         title(ax, sprintf('%s - %s', name, strrep(v, '_', '\_')));
         xlabel(ax, 'Time');
         ylabel(ax, '');
         if n == 1 && m == 1
            legend(ax, 'show', 'Location', 'best');
         end
         grid(ax, 'on');
      end
   end
   sgtitle(tl, 'Colbeck 1976: cached vs computed solutions (4-way overlay)');

   % Save artifact if a destination was requested.
   figure_path = "";
   if save_plot && ~isblanktext(artifact_dir)
      icemodel.helpers.ensureDirExists(artifact_dir);
      figure_path = fullfile(artifact_dir, "compareSolutions.png");
      exportgraphics(f, figure_path, 'Resolution', 150);
   end
end

function plotOneSeries(ax, tt, varname, display_name, line_style)

   % Skip variables that are missing or fully non-finite for this bundle.
   if ~ismember(varname, tt.Properties.VariableNames)
      return
   end
   values = tt.(varname);
   if all(~isfinite(values))
      return
   end
   icemodel.plot.timeseries(tt.Time, values, axes=ax, ...
      display_name=string(display_name), line_style=line_style);
end
