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

   % Per-variable RMSE tolerances. The verification suite is the real-time
   % correctness signal: failures here mean the snow model is not reproducing
   % the reference accurately.
   tolerance = struct( ...
      'snow_liquid_water_storage_m', kwargs.tolerance_storage_m, ...
      'bottom_outflow_mps',          kwargs.tolerance_outflow_mps);

   manifest = icemodel.verification.loadmanifest("colbeck1976");
   targets = icemodel.verification.helpers.loadArtifact( ...
      manifest.evaluation_path, "targets");

   candidates = struct( ...
      'numerical_icemodel',  icemodel.verification.colbeck.runCase( ...
      manifest, kind="numerical", ...
      experiment_names=kwargs.experiment_names), ...
      'analytical_icemodel', icemodel.verification.colbeck.runCase( ...
      manifest, kind="analytical", ...
      experiment_names=kwargs.experiment_names));

   % Four pairwise axes: (target_source, candidate_source, axis_role).
   pairs = { ...
      "numerical_summa",      "numerical_icemodel",  "formal"
      "analytical_clark2017", "numerical_icemodel",  "formal"
      "numerical_summa",      "analytical_icemodel", "diagnostic"
      "analytical_clark2017", "analytical_icemodel", "diagnostic"};

   metric_rows = table();
   variables = string(manifest.comparison_variables);
   for p = 1:size(pairs, 1)
      target_source    = string(pairs{p, 1});
      candidate_source = string(pairs{p, 2});
      axis_role        = string(pairs{p, 3});
      target_bundle    = targets.(char(target_source));
      candidate_bundle = candidates.(char(candidate_source));
      for i = 1:numel(kwargs.experiment_names)
         experiment_name = kwargs.experiment_names(i);
         for v = 1:numel(variables)
            varname = variables(v);
            row = computeOnePair( ...
               target_bundle.experiments.(char(experiment_name)), ...
               candidate_bundle.experiments.(char(experiment_name)), ...
               varname);
            row.target_source    = target_source;
            row.candidate_source = candidate_source;
            row.axis_role        = axis_role;
            row.case_id          = manifest.case_id;
            row.experiment       = experiment_name;
            row.variable         = varname;
            metric_rows = [metric_rows; row]; %#ok<AGROW>
         end
      end
   end
   metrics_table = metric_rows;

   summary = printSummary(metrics_table, tolerance);

   figure_path = "";
   csv_path = "";
   if ~isblanktext(kwargs.artifact_dir)
      icemodel.helpers.ensureDirExists(kwargs.artifact_dir);
      csv_path = fullfile(kwargs.artifact_dir, "compareSolutions_metrics.csv");
      writetable(metrics_table, csv_path);
   end

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

function row = computeOnePair(target_tt, candidate_tt, varname)
   if ~ismember(varname, string(target_tt.Properties.VariableNames)) || ...
         ~ismember(varname, string(candidate_tt.Properties.VariableNames))
      row = table("not_applicable", 0, NaN, NaN, NaN, ...
         'VariableNames', ["status","n","bias","rmse","correlation"]);
      return
   end

   sync = synchronize(target_tt(:, varname), candidate_tt(:, varname), ...
      'first', 'nearest');
   target_vals = sync{:, 1};
   cand_vals = sync{:, 2};
   ok_pair = isfinite(target_vals) & isfinite(cand_vals);
   target_vals = target_vals(ok_pair);
   cand_vals = cand_vals(ok_pair);
   n = numel(target_vals);

   if n == 0
      row = table("not_applicable", 0, NaN, NaN, NaN, ...
         'VariableNames', ["status","n","bias","rmse","correlation"]);
      return
   end

   bias = mean(cand_vals - target_vals);
   rmse = sqrt(mean((cand_vals - target_vals) .^ 2));
   if std(target_vals) > 0 && std(cand_vals) > 0 && n > 1
      cc = corrcoef(target_vals, cand_vals);
      correlation = cc(1, 2);
   else
      correlation = NaN;
   end
   row = table("ok", n, bias, rmse, correlation, ...
      'VariableNames', ["status","n","bias","rmse","correlation"]);
end

function summary = printSummary(metrics_table, tolerance)
   formal = metrics_table(metrics_table.axis_role == "formal", :);
   tol = arrayfun(@(v) tolerance.(char(v)), formal.variable);
   pass = formal.status == "ok" & formal.rmse < tol;
   fprintf("\n=== Colbeck 1976 verification summary ===\n");
   fprintf("Formal axes (per-variable RMSE tolerance):\n");
   fprintf("  storage  < %.3g m\n", tolerance.snow_liquid_water_storage_m);
   fprintf("  outflow  < %.3g m s-1\n", tolerance.bottom_outflow_mps);
   for i = 1:height(formal)
      mark = "FAIL";
      if pass(i)
         mark = "PASS";
      end
      fprintf("  [%s] %s vs %s @ %s %s: rmse=%.3e bias=%.3e n=%d\n", ...
         mark, formal.target_source(i), formal.candidate_source(i), ...
         formal.experiment(i), formal.variable(i), ...
         formal.rmse(i), formal.bias(i), formal.n(i));
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

function [f, figure_path] = plotSolutions(targets, candidates, ...
      experiment_names, artifact_dir, save_plot, plot_visible)
   variables = ["snow_liquid_water_storage_m"; "bottom_outflow_mps"];
   n_exp = numel(experiment_names);
   n_var = numel(variables);

   f = figure('Visible', plot_visible, 'Position', [100, 100, 1200, 800]);
   tl = tiledlayout(f, n_exp, n_var, 'TileSpacing', 'compact', ...
      'Padding', 'compact');

   for i = 1:n_exp
      name = char(experiment_names(i));
      for j = 1:n_var
         v = char(variables(j));
         ax = nexttile(tl);
         hold(ax, 'on');
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
         if i == 1 && j == 1
            legend(ax, 'show', 'Location', 'best');
         end
         grid(ax, 'on');
      end
   end
   sgtitle(tl, 'Colbeck 1976: cached vs computed solutions (4-way overlay)');

   figure_path = "";
   if save_plot && ~isblanktext(artifact_dir)
      icemodel.helpers.ensureDirExists(artifact_dir);
      figure_path = fullfile(artifact_dir, "compareSolutions.png");
      exportgraphics(f, figure_path, 'Resolution', 150);
   end
end

function plotOneSeries(ax, tt, varname, display_name, line_style)
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
