function report_file = buildTestSuiteReport(suite_kind, results, kwargs)
   %BUILDTESTSUITEREPORT Render numerical or performance suite results.
   %
   %  report_file = icemodel.verification.report.buildTestSuiteReport( ...
   %     "regression", results)
   %  report_file = icemodel.verification.report.buildTestSuiteReport( ...
   %     "performance", results, render=false, output_dir=tempdir)
   %
   % The saved suite result struct is the only report input. The renderer writes
   % a compact CSV, MATLAB-exported figures, generated QMD, and a self-contained
   % Quarto HTML report beside the suite MAT artifacts.

   arguments
      suite_kind (1, 1) string {mustBeMember(suite_kind, ...
         ["regression", "performance"])}
      results (1, 1) struct
      kwargs.render (1, 1) logical = true
      kwargs.output_dir (1, 1) string = ""
      kwargs.baseline_root (1, 1) string = ""
   end

   % Resolve one shared artifact directory so multi-model runs produce one
   % combined report rather than one report per saved MAT file.
   output_dir = resolveOutputDir(results, kwargs.output_dir);
   icemodel.helpers.ensureDirExists(output_dir)
   asset_dir = fullfile(output_dir, "report-assets");
   icemodel.helpers.ensureDirExists(asset_dir)

   % Normalize suite-specific result fields into one compact report table.
   if suite_kind == "regression"
      summary = regressionSummary(results);
      assets = regressionFigures(summary, asset_dir);
      title_text = "Numerical regression suite report";
   else
      summary = performanceSummary(results);
      baseline_root = kwargs.baseline_root;
      if baseline_root == ""
         baseline_root = fullfile(icemodel.getpath("test"), "baselines");
      end
      assets = performanceFigures(summary, asset_dir, baseline_root);
      title_text = "Performance suite report";
   end

   % Preserve the compact machine-readable summary beside the human report.
   summary_file = fullfile(output_dir, suite_kind + "-suite-summary.csv");
   writetable(summary, summary_file)

   % Generate plain QMD so the report remains inspectable without executing code.
   qmd_file = fullfile(output_dir, suite_kind + "-suite-report.qmd");
   report_file = fullfile(output_dir, suite_kind + "-suite-report.html");
   lines = reportMarkdown(title_text, suite_kind, results, summary, assets, ...
      summary_file, report_file);
   writelines(lines, qmd_file)

   % Rendering is optional only for focused tests and report-source inspection.
   if kwargs.render
      command = "quarto render " ...
         + icemodel.shellQuote(qmd_file);
      [status, output] = system(command);
      if status ~= 0
         error("icemodel:verification:report:quartoFailed", ...
            "Quarto failed to render %s:\n%s", qmd_file, output)
      end
      assert(isfile(report_file), ...
         "icemodel:verification:report:missingHtml", ...
         "Quarto did not create the expected report: %s", report_file)
   end
end

function output_dir = resolveOutputDir(results, requested)
   %RESOLVEOUTPUTDIR Resolve the common directory for combined suite artifacts.

   if requested ~= ""
      output_dir = requested;
      return
   end

   % Default report placement follows the artifact paths emitted by the runners.
   assert(isfield(results, "artifact_file") && ~isempty(results.artifact_file), ...
      "icemodel:verification:report:missingArtifactPath", ...
      "Provide output_dir when results do not contain artifact_file.")
   files = string(results.artifact_file(:));
   dirs = arrayfun(@(file) string(fileparts(char(file))), files);
   dirs = unique(dirs);
   assert(isscalar(dirs), ...
      "icemodel:verification:report:multipleArtifactDirs", ...
      "Combined suite artifacts must share one output directory.")
   output_dir = dirs;
end

function summary = regressionSummary(results)
   %REGRESSIONSUMMARY Build the compact numerical comparison table.

   assert(isfield(results, "report") && istable(results.report), ...
      "icemodel:verification:report:missingRegressionTable", ...
      "Regression results must contain a report table.")
   report = results.report;
   vars = string(report.Properties.VariableNames);

   % Add the two iteration deltas that make solver-cost changes visible.
   if all(ismember(["mean_Tice_numiter", ...
         "baseline_mean_Tice_numiter"], vars))
      report.mean_iteration_delta = report.mean_Tice_numiter ...
         - report.baseline_mean_Tice_numiter;
   end
   if all(ismember(["max_Tice_numiter", ...
         "baseline_max_Tice_numiter"], vars))
      report.max_iteration_delta = report.max_Tice_numiter ...
         - report.baseline_max_Tice_numiter;
   end

   % Retain identity, physical deltas, convergence, and pass/fail fields only.
   keep = intersect(["case_id", "tier", "smbmodel", "sitename", "solver", ...
      "runoff_pct_delta", "melt_pct_delta", "runoff_eval_pct_delta", ...
      "melt_eval_pct_delta", "mean_iteration_delta", ...
      "max_iteration_delta", "n_not_converged", "closure_seb_rmse", ...
      "baseline_closure_seb_rmse", "closure_seb_max_abs", ...
      "baseline_closure_seb_max_abs", "passed"], ...
      string(report.Properties.VariableNames), 'stable');
   summary = report(:, keep);
end

function summary = performanceSummary(results)
   %PERFORMANCESUMMARY Build the compact timing comparison table.

   assert(isfield(results, "case_summary") && istable(results.case_summary), ...
      "icemodel:verification:report:missingPerformanceTable", ...
      "Performance results must contain a case_summary table.")
   report = results.case_summary;
   vars = string(report.Properties.VariableNames);

   % A normalized ratio makes cases with different absolute runtimes comparable.
   if all(ismember(["median_wall_s", "ref_wall_s"], vars))
      report.runtime_ratio = report.median_wall_s ./ report.ref_wall_s;
   end
   keep = intersect(["case_id", "tier", "smbmodel", "sitename", "solver", ...
      "median_wall_s", "ref_wall_s", "floor_wall_s", "gate_wall_s", ...
      "runtime_ratio", "baseline_compatible", "passed_perf", ...
      "compare_reason"], ...
      string(report.Properties.VariableNames), 'stable');
   summary = report(:, keep);
end

function assets = regressionFigures(summary, asset_dir)
   %REGRESSIONFIGURES Plot physical percent changes and iteration deltas.

   assets = strings(0, 1);
   vars = string(summary.Properties.VariableNames);
   delta_vars = intersect(["runoff_pct_delta", "melt_pct_delta", ...
      "runoff_eval_pct_delta", "melt_eval_pct_delta"], vars, 'stable');
   delta_vars = delta_vars(arrayfun(@(name) ...
      any(isfinite(summary.(name)), "all"), delta_vars));

   % Group all physical metrics in one figure so case-scale changes are obvious.
   if ~isempty(delta_vars) && height(summary) > 0
      file = fullfile(asset_dir, "regression-percent-deltas.png");
      [fig, ax] = newReportFigure(height(summary));
      barh(ax, summary{:, delta_vars})
      grid(ax, "on")
      xlabel(ax, "Difference from accepted baseline (%; + means more)")
      title(ax, "Change in runoff and melt")
      icemodel.verification.report.configureCategoryAxis(ax, summary.case_id)
      labels = replace(delta_vars, ...
         ["runoff_pct_delta", "melt_pct_delta", ...
          "runoff_eval_pct_delta", "melt_eval_pct_delta"], ...
         ["Full-run runoff", "Full-run melt", ...
          "Evaluation-window runoff", "Evaluation-window melt"]);
      lgd = legend(ax, labels, Location="eastoutside");
      formatReportLegend(lgd)
      icemodel.verification.report.exportAndClose(fig, file)
      assets(end + 1) = "report-assets/regression-percent-deltas.png";
   end

   % Plot iteration changes separately because they are counts, not percentages.
   iteration_vars = intersect(["mean_iteration_delta", ...
      "max_iteration_delta"], vars, 'stable');
   iteration_vars = iteration_vars(arrayfun(@(name) ...
      any(isfinite(summary.(name)), "all"), iteration_vars));
   if ~isempty(iteration_vars) && height(summary) > 0
      file = fullfile(asset_dir, "regression-iteration-deltas.png");
      [fig, ax] = newReportFigure(height(summary));
      barh(ax, summary{:, iteration_vars})
      xline(ax, 0, "k-")
      grid(ax, "on")
      xlabel(ax, "Difference from accepted baseline (+ means more work)")
      title(ax, "Change in solver iterations")
      icemodel.verification.report.configureCategoryAxis(ax, summary.case_id)
      labels = replace(iteration_vars, ...
         ["mean_iteration_delta", "max_iteration_delta"], ...
         ["Mean iterations", "Maximum iterations"]);
      lgd = legend(ax, labels, Location="eastoutside");
      formatReportLegend(lgd)
      icemodel.verification.report.exportAndClose(fig, file)
      assets(end + 1) = "report-assets/regression-iteration-deltas.png";
   end

   % Surface-residual statistics expose coupling and convergence problems
   % that cumulative runoff and melt totals can hide.
   rmse_vars = intersect(["closure_seb_rmse", ...
      "baseline_closure_seb_rmse"], vars, 'stable');
   rmse_vars = rmse_vars(arrayfun(@(name) ...
      any(isfinite(summary.(name)), "all"), rmse_vars));
   max_vars = intersect(["closure_seb_max_abs", ...
      "baseline_closure_seb_max_abs"], vars, 'stable');
   max_vars = max_vars(arrayfun(@(name) ...
      any(isfinite(summary.(name)), "all"), max_vars));
   if (~isempty(rmse_vars) || ~isempty(max_vars)) && height(summary) > 0
      file = fullfile(asset_dir, "regression-seb-closure.png");
      [fig, placeholder] = newReportFigure(height(summary));
      delete(placeholder)
      n_panels = double(~isempty(rmse_vars)) + double(~isempty(max_vars));
      fig.Position(3) = 900 * n_panels;
      layout = tiledlayout(fig, 1, n_panels, ...
         TileSpacing="compact", Padding="compact");

      % Give typical and worst residuals separate scales so neither is hidden.
      if ~isempty(rmse_vars)
         ax = nexttile(layout);
         barh(ax, summary{:, rmse_vars})
         grid(ax, "on")
         xlabel(ax, "Residual RMSE (W m^{-2})")
         title(ax, "Typical residual")
         icemodel.verification.report.configureCategoryAxis(ax, summary.case_id)
         labels = replace(rmse_vars, ...
            ["closure_seb_rmse", "baseline_closure_seb_rmse"], ...
            ["Current", "Accepted baseline"]);
         lgd = legend(ax, labels, ...
            Location="southoutside", Orientation="horizontal");
         formatReportLegend(lgd)
      end

      % The maximum absolute residual exposes short-lived closure outliers.
      if ~isempty(max_vars)
         ax = nexttile(layout);
         barh(ax, summary{:, max_vars})
         grid(ax, "on")
         xlabel(ax, "Maximum absolute residual (W m^{-2})")
         title(ax, "Worst residual")
         icemodel.verification.report.configureCategoryAxis(ax, summary.case_id)
         labels = replace(max_vars, ...
            ["closure_seb_max_abs", "baseline_closure_seb_max_abs"], ...
            ["Current", "Accepted baseline"]);
         lgd = legend(ax, labels, ...
            Location="southoutside", Orientation="horizontal");
         formatReportLegend(lgd)
      end
      title(layout, "Surface-energy residuals (smaller is better)", ...
         Color="k")
      icemodel.verification.report.exportAndClose(fig, file)
      assets(end + 1) = "report-assets/regression-seb-closure.png";
   end
end

function assets = performanceFigures(summary, asset_dir, baseline_root)
   %PERFORMANCEFIGURES Plot current timing, ratios, and accepted history.

   assets = strings(0, 1);
   vars = string(summary.Properties.VariableNames);

   % Show absolute current/reference timings for operational runtime context.
   if all(ismember(["median_wall_s", "ref_wall_s"], vars)) ...
         && height(summary) > 0
      file = fullfile(asset_dir, "performance-current-reference.png");
      [fig, ax] = newReportFigure(height(summary));
      barh(ax, [summary.median_wall_s, summary.ref_wall_s])
      grid(ax, "on")
      xlabel(ax, "Wall time (s)")
      title(ax, "Current and accepted median runtime")
      icemodel.verification.report.configureCategoryAxis(ax, summary.case_id)
      lgd = legend(ax, ["current", "accepted"], Location="eastoutside");
      formatReportLegend(lgd)
      icemodel.verification.report.exportAndClose(fig, file)
      assets(end + 1) = "report-assets/performance-current-reference.png";
   end

   % Normalize each case to its own reference so the performance gate is legible.
   if ismember("runtime_ratio", vars) && height(summary) > 0
      file = fullfile(asset_dir, "performance-runtime-ratio.png");
      [fig, ax] = newReportFigure(height(summary));
      barh(ax, summary.runtime_ratio)
      hold(ax, "on")
      xline(ax, 1, "k-")
      if all(ismember(["floor_wall_s", "gate_wall_s", ...
             "ref_wall_s"], vars))
          scatter(ax, summary.floor_wall_s ./ summary.ref_wall_s, ...
             1:height(summary), 45, "filled", Marker="<")
          scatter(ax, summary.gate_wall_s ./ summary.ref_wall_s, ...
             1:height(summary), 45, "filled", Marker=">")
          lgd = legend(ax, ["current/reference", "accepted", ...
             "lower gate", "upper gate"], Location="eastoutside");
       else
          lgd = legend(ax, ["current/reference", "accepted"], ...
             Location="eastoutside");
       end
      formatReportLegend(lgd)
      hold(ax, "off")
      grid(ax, "on")
      xlabel(ax, "Runtime ratio")
      title(ax, "Runtime relative to accepted baseline")
      icemodel.verification.report.configureCategoryAxis(ax, summary.case_id)
      icemodel.verification.report.exportAndClose(fig, file)
      assets(end + 1) = "report-assets/performance-runtime-ratio.png";
   end

   % Accepted rolling archives provide the only durable longitudinal history.
   history = loadPerformanceHistory(baseline_root, string(summary.case_id));
   if ~isempty(history)
      file = fullfile(asset_dir, "performance-history.png");
      cases = unique(history.case_id, 'stable');
      n_cols = 2;
      n_rows = ceil(numel(cases) / n_cols);
      fig = icemodel.plot.newFigure(width=1300, ...
         height=max(720, 280 * n_rows));
      layout = tiledlayout(fig, n_rows, n_cols, ...
         TileSpacing="compact", Padding="compact");
      for k = 1:numel(cases)
         rows = history(history.case_id == cases(k), :);
         ax = nexttile(layout);
         plot(ax, rows.timestamp_utc, rows.median_wall_s, "-o", ...
            LineWidth=1.3, MarkerSize=5)
         grid(ax, "on")
         ylabel(ax, "Median wall time (s)")
         title(ax, cases(k), Interpreter="none")
         icemodel.verification.report.formatReportAxes(ax)
      end
      title(layout, "Rolling performance baseline history", Color="k")
      icemodel.verification.report.exportAndClose(fig, file)
      assets(end + 1) = "report-assets/performance-history.png";
   end
end

function history = loadPerformanceHistory(baseline_root, requested_cases)
   %LOADPERFORMANCEHISTORY Load accepted rolling baseline timing rows.

   history = table();
   if ~isfolder(baseline_root)
      return
   end

   % Load only managed rolling baseline MAT files, excluding profiler sidecars.
   files = dir(fullfile(baseline_root, "**", ...
      "perf_baseline_*_rolling_*.mat"));
   collected = cell(numel(files), 1);
   for k = 1:numel(files)
      filename = fullfile(files(k).folder, files(k).name);
      saved = load(filename, "PerfBaseline");
      if ~isfield(saved, "PerfBaseline") || ~istable(saved.PerfBaseline)
         continue
      end
      baseline = saved.PerfBaseline;
      vars = string(baseline.Properties.VariableNames);
      if ~all(ismember(["case_id", "median_wall_s"], vars))
         continue
      end

      % Prefer the accepted timestamp stored with each row; legacy files fall
      % back to their filesystem modification time.
      if ismember("last_updated_utc", vars)
         timestamp = baseline.last_updated_utc;
         if ~isdatetime(timestamp)
            timestamp = datetime(string(timestamp), TimeZone="UTC");
         elseif isempty(timestamp.TimeZone)
            timestamp.TimeZone = "UTC";
         end
      else
         timestamp = repmat(datetime(files(k).datenum, ...
            ConvertFrom="datenum", TimeZone="UTC"), height(baseline), 1);
      end
      rows = table(timestamp, string(baseline.case_id), ...
         double(baseline.median_wall_s), ...
         VariableNames=["timestamp_utc", "case_id", "median_wall_s"]);
      collected{k} = rows;
   end
   history = [history; vertcat(collected{:})];

   % Restrict the plot to the cases in the current comparison and valid samples.
   if ~isempty(history)
      history = history(ismember(history.case_id, requested_cases) ...
         & isfinite(history.median_wall_s) & ~isnat(history.timestamp_utc), :);
      history = sortrows(history, ["case_id", "timestamp_utc"]);
   end
end

function lines = reportMarkdown(title_text, suite_kind, results, summary, ...
      assets, summary_file, report_file)
   %REPORTMARKDOWN Build generated QMD using only saved result values.

   passed = isfield(results, "passed") && logical(results.passed);
   outcome = "FAILED";
   if passed
      outcome = "PASSED";
   end
   generated = string(datetime("now", TimeZone="UTC", ...
      Format="yyyy-MM-dd HH:mm:ss 'UTC'"));

   % Keep report configuration self-contained and free of executable code.
   [~, output_name, output_ext] = fileparts(report_file);
   lines = [ ...
      "---"
      "title: """ + title_text + """"
      "date: """ + generated + """"
      "format:"
      "  html:"
      "    embed-resources: true"
      "    toc: true"
      "output-file: """ + output_name + output_ext + """"
      "---"
      ""
      "## Outcome"
      ""
      "**" + outcome + "** — " + string(height(summary)) + ...
         " formal case(s)."
      ""
      reportMetadata(summary, suite_kind, results)];

   % Regression reports lead with the decision-relevant facts in plain language.
   if suite_kind == "regression"
      lines = [lines
         ""
         regressionExplanation(summary)];
   end
   lines = [lines
      ""
      "## Visual summary"
      ""];

   % Put the visual evidence before the compact table, so reviewers see it
   % first.
   if isempty(assets)
      lines(end + 1) = "No plottable suite metrics were present.";
   else
      asset_lines = strings(2 * numel(assets), 1);
      for k = 1:numel(assets)
         asset_lines(2 * k - 1) = "";
         asset_lines(2 * k) = "![" + assetCaption(assets(k)) + "](" ...
            + assets(k) + ")";
      end
      lines = [lines; asset_lines];
   end

   % Finish with a compact exact-value table and links to machine-readable data.
   [~, csv_name, csv_ext] = fileparts(summary_file);
   lines = [lines
      ""
      "## Summary table"
      ""
      "The field names match the downloadable CSV, which retains full stored " ...
         + "numeric precision."
      ""
      icemodel.verification.report.markdownTable(summary)
      ""
      "[Download the compact CSV](" + csv_name + csv_ext + ")"];

   % Explain only the longitudinal evidence present for the selected suite.
   if suite_kind == "performance"
      lines = [lines
         ""
         "Performance history, when shown, contains accepted rolling baselines. " ...
            + "Changes can reflect code, MATLAB, host, or suite-contract revisions."];
   end
   lines(end + 1) = "";
end

function lines = reportMetadata(summary, suite_kind, results)
   %REPORTMETADATA Format the small run identity block.

   vars = string(summary.Properties.VariableNames);
   lines = strings(0, 1);
   lines(end + 1, 1) = "- Suite: `" + suite_kind + "`";
   if ismember("tier", vars)
      lines(end + 1, 1) = "- Tier: " ...
         + icemodel.verification.report.markdownCode( ...
         join(unique(string(summary.tier), 'stable'), ", "));
   end
   if ismember("smbmodel", vars)
      lines(end + 1, 1) = "- Models: " ...
         + icemodel.verification.report.markdownCode( ...
         join(unique(string(summary.smbmodel), 'stable'), ", "));
   end
   % Combined model reports carry one metadata struct per model; run identity
   % fields are shared, so the first struct is the report-level source.
   meta = struct();
   if isfield(results, "meta") && isstruct(results.meta) ...
         && ~isempty(results.meta)
      meta = results.meta(1);
   end
   if isfield(meta, "baseline_tag") ...
         && strlength(string(meta.baseline_tag)) > 0
      lines(end + 1, 1) = "- Accepted baseline: " ...
         + icemodel.verification.report.markdownCode(string(meta.baseline_tag));
   elseif isfield(meta, "baseline_type") ...
         && strlength(string(meta.baseline_type)) > 0
      lines(end + 1, 1) = "- Accepted baseline: " ...
         + icemodel.verification.report.markdownCode(string(meta.baseline_type));
   end
   if isfield(meta, "input_path") ...
         && strlength(string(meta.input_path)) > 0
      lines(end + 1, 1) = "- Current input root: " ...
         + icemodel.verification.report.markdownCode(string(meta.input_path));
   end
end

function lines = regressionExplanation(summary)
   %REGRESSIONEXPLANATION State the comparison result without test jargon.

   lines = [ ...
      "## What this means"
      ""
      "**The outcome above is a strict software regression check. A FAILED " ...
         + "outcome means values exceeded the accepted baseline's tolerances " ...
         + "or that a formal case could not complete. Scientific validity " ...
         + "must be assessed separately against observations and physical " ...
         + "constraints.**"
      ""
      "Positive runoff or melt percentages mean the current run produced " ...
         + "more than the accepted baseline; negative percentages mean less. " ...
         + "The figure shows whichever full-run or evaluation-window " ...
         + "percentage changes are available and omits unavailable series."
      ""
      "This report is read-only. It did not accept or replace any baseline."
      ""];
   vars = string(summary.Properties.VariableNames);

   % Summarize each model and site group, so that repeated solvers do not hide
   % the pattern.
   if all(ismember(["smbmodel", "sitename"], vars))
      metric_names = ["runoff_pct_delta", "melt_pct_delta", ...
         "runoff_eval_pct_delta", "melt_eval_pct_delta"];
      metric_labels = ["full-run runoff", "full-run melt", ...
         "evaluation-window runoff", "evaluation-window melt"];
      pairs = unique([string(summary.smbmodel), string(summary.sitename)], ...
         'rows', 'stable');
      group_lines = strings(size(pairs, 1), 1);
      for k = 1:size(pairs, 1)
         use = string(summary.smbmodel) == pairs(k, 1) ...
            & string(summary.sitename) == pairs(k, 2);
         model = replace(pairs(k, 1), ...
            ["icemodel", "skinmodel"], ["IceModel", "SkinModel"]);
         site = replace(upper(pairs(k, 2)), ...
            ["KANM", "KANL"], ["KAN_M", "KAN_L"]);
         n_cases = sum(use);
         if n_cases == 1
            case_text = "1 formal case";
            count_label = "case";
         else
            case_text = string(n_cases) + " formal cases";
            count_label = "cases";
         end

         % Treat every saved physical metric independently and name absences.
         metric_text = strings(size(metric_names));
         for m = 1:numel(metric_names)
            values = [];
            if ismember(metric_names(m), vars)
               values = summary.(metric_names(m));
               values = values(use);
               values = values(isfinite(values));
            end
            if isempty(values)
               value_range = "unavailable";
            elseif isscalar(values)
               value_range = string(sprintf("%+.2f%%", values));
            else
               value_range = string(sprintf("%+.2f%% to %+.2f%%", ...
                  min(values), max(values)));
            end
            metric_text(m) = metric_labels(m) + " " + value_range ...
               + " (" + string(numel(values)) + "/" + string(n_cases) ...
               + " " + count_label + ")";
         end
         group_lines(k) = "- " ...
            + icemodel.verification.report.markdownCode(model) + " at " ...
            + icemodel.verification.report.markdownCode(site) ...
            + ": " + join(metric_text, "; ") ...
            + " across " + case_text + ".";
      end
      lines = [lines; group_lines];
   end

   % Report convergence and solver effort separately from physical totals.
   if ismember("n_not_converged", vars)
      counts = summary.n_not_converged;
      saved = isfinite(counts);
      if ~any(saved)
         lines(end + 1, 1) = "- **Non-converged timesteps:** unavailable.";
      else
         lines(end + 1, 1) = "- **Non-converged timesteps:** " ...
            + string(sum(counts(saved))) + " across cases with saved values " ...
            + "(coverage " + string(sum(saved)) + "/" ...
            + string(numel(counts)) + "; missing " ...
            + string(sum(~saved)) + ").";
      end
   end
   if all(ismember(["mean_iteration_delta", ...
         "max_iteration_delta"], vars))
      mean_delta = summary.mean_iteration_delta;
      max_delta = summary.max_iteration_delta;
      mean_delta = mean_delta(isfinite(mean_delta));
      max_delta = max_delta(isfinite(max_delta));
      if ~isempty(mean_delta) && ~isempty(max_delta)
         lines(end + 1, 1) = "- **Solver effort:** mean iterations changed " ...
            + string(sprintf("%+.2f to %+.2f", min(mean_delta), ...
            max(mean_delta))) + "; maximum iterations changed " ...
            + string(sprintf("%+.2f to %+.2f", min(max_delta), ...
            max(max_delta))) + ".";
      end
   end

   % Separate the current closure quality from comparison evidence, which can
   % be absent.
   current_closure = ["closure_seb_rmse", "closure_seb_max_abs"];
   if all(ismember(current_closure, vars))
      rmse = summary.closure_seb_rmse;
      max_abs = summary.closure_seb_max_abs;
      rmse = rmse(isfinite(rmse));
      max_abs = max_abs(isfinite(max_abs));
      if ~isempty(rmse) && ~isempty(max_abs)
         lines(end + 1, 1) = "- **Current energy-balance residuals:** " ...
            + "RMSE reached " + string(sprintf("%.2f", max(rmse))) ...
            + " W m^-2; the worst absolute residual reached " ...
            + string(sprintf("%.2f", max(max_abs))) + " W m^-2. " ...
            + "Smaller is better.";
      end
   end
   baseline_closure = ["baseline_closure_seb_rmse", ...
      "baseline_closure_seb_max_abs"];
   if all(ismember(baseline_closure, vars))
      baseline_values = [summary.baseline_closure_seb_rmse; ...
         summary.baseline_closure_seb_max_abs];
      if ~any(isfinite(baseline_values))
         lines(end + 1, 1) = "- **Closure comparison unavailable:** the " ...
            + "accepted baseline has no saved surface-energy residuals, so " ...
            + "the closure figure describes current runs only.";
      end
   end
end





function caption = assetCaption(asset)
   %ASSETCAPTION Return a concise caption from one stable asset name.

   [~, name] = fileparts(asset);
   caption = replace(string(name), ["-", "_"], " ");
end


function [fig, ax] = newReportFigure(n_rows)
   %NEWREPORTFIGURE Create a light report figure sized for its case rows.

   height_px = max(520, 58 * n_rows + 180);
   fig = icemodel.plot.newFigure(width=1300, height=height_px);
   ax = axes(fig);
   icemodel.verification.report.formatReportAxes(ax)
end



function formatReportLegend(lgd)
   %FORMATREPORTLEGEND Isolate exported legends from interactive themes.

   lgd.Color = "w";
   lgd.TextColor = "k";
   lgd.EdgeColor = [0.75 0.75 0.75];
end
