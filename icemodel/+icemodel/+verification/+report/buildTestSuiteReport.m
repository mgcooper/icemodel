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
      command = "quarto render " + shellQuote(qmd_file);
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

   % Group all physical metrics in one figure so case-scale changes are obvious.
   if ~isempty(delta_vars) && height(summary) > 0
      file = fullfile(asset_dir, "regression-percent-deltas.png");
      [fig, ax] = newReportFigure(height(summary));
      barh(ax, summary{:, delta_vars})
      grid(ax, "on")
      xlabel(ax, "Change from accepted baseline (%)")
      title(ax, "Numerical changes by formal case")
      configureCaseAxis(ax, summary.case_id)
      lgd = legend(ax, replace(delta_vars, "_pct_delta", ""), ...
         Interpreter="none", Location="eastoutside");
      formatReportLegend(lgd)
      exportAndClose(fig, file)
      assets(end + 1) = "report-assets/regression-percent-deltas.png";
   end

   % Plot iteration changes separately because they are counts, not percentages.
   iteration_vars = intersect(["mean_iteration_delta", ...
      "max_iteration_delta"], vars, 'stable');
   if ~isempty(iteration_vars) && height(summary) > 0
      file = fullfile(asset_dir, "regression-iteration-deltas.png");
      [fig, ax] = newReportFigure(height(summary));
      barh(ax, summary{:, iteration_vars})
      xline(ax, 0, "k-")
      grid(ax, "on")
      xlabel(ax, "Iteration-count change")
      title(ax, "Solver iteration changes by formal case")
      configureCaseAxis(ax, summary.case_id)
      lgd = legend(ax, replace(iteration_vars, "_", " "), ...
         Interpreter="none", Location="eastoutside");
      formatReportLegend(lgd)
      exportAndClose(fig, file)
      assets(end + 1) = "report-assets/regression-iteration-deltas.png";
   end

   % Surface-residual statistics expose coupling and convergence problems
   % that cumulative runoff and melt totals can hide.
   closure_vars = intersect(["closure_seb_rmse", ...
      "baseline_closure_seb_rmse", "closure_seb_max_abs", ...
      "baseline_closure_seb_max_abs"], vars, 'stable');
   if ~isempty(closure_vars) && height(summary) > 0
      file = fullfile(asset_dir, "regression-seb-closure.png");
      [fig, ax] = newReportFigure(height(summary));
      barh(ax, summary{:, closure_vars})
      grid(ax, "on")
      xlabel(ax, "Surface-energy residual (W m^{-2})")
      title(ax, "Surface-energy closure by formal case")
      configureCaseAxis(ax, summary.case_id)
      lgd = legend(ax, replace(closure_vars, "_", " "), ...
         Interpreter="none", Location="eastoutside");
      formatReportLegend(lgd)
      exportAndClose(fig, file)
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
      configureCaseAxis(ax, summary.case_id)
      lgd = legend(ax, ["current", "accepted"], Location="eastoutside");
      formatReportLegend(lgd)
      exportAndClose(fig, file)
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
      configureCaseAxis(ax, summary.case_id)
      exportAndClose(fig, file)
      assets(end + 1) = "report-assets/performance-runtime-ratio.png";
   end

   % Accepted rolling archives provide the only durable longitudinal history.
   history = loadPerformanceHistory(baseline_root, string(summary.case_id));
   if ~isempty(history)
      file = fullfile(asset_dir, "performance-history.png");
      cases = unique(history.case_id, 'stable');
      n_cols = 2;
      n_rows = ceil(numel(cases) / n_cols);
      fig = figure(Visible="off", Color="w", ...
         Position=[100 100 1300 max(720, 280 * n_rows)]);
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
         formatReportAxes(ax)
      end
      title(layout, "Rolling performance baseline history", Color="k")
      exportAndClose(fig, file)
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
      history = [history; rows]; %#ok<AGROW>
   end

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
      reportMetadata(summary, suite_kind)
      ""
      "## Visual summary"
      ""];

   % Put visual evidence before the compact table, as the primary review surface.
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
      "## Compact results"
      ""
      markdownTable(summary)
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

function lines = reportMetadata(summary, suite_kind)
   %REPORTMETADATA Format the small run identity block.

   vars = string(summary.Properties.VariableNames);
   lines = strings(0, 1);
   lines(end + 1, 1) = "- Suite: `" + suite_kind + "`";
   if ismember("tier", vars)
      lines(end + 1, 1) = "- Tier: `" ...
         + join(unique(string(summary.tier), 'stable'), ", ") + "`";
   end
   if ismember("smbmodel", vars)
      lines(end + 1, 1) = "- Models: `" ...
         + join(unique(string(summary.smbmodel), 'stable'), ", ") + "`";
   end
end

function lines = markdownTable(values)
   %MARKDOWNTABLE Convert a compact scalar table to Markdown.

   vars = string(values.Properties.VariableNames);
   header = "| " + join(replace(vars, "_", " "), " | ") + " |";
   divider = "| " + join(repmat("---", size(vars)), " | ") + " |";
   lines = [header; divider];
   for row = 1:height(values)
      cells = strings(size(vars));
      for col = 1:numel(vars)
         column = values.(vars(col));
         cells(col) = formatValue(column(row, :));
      end
      lines(end + 1) = "| " + join(cells, " | ") + " |"; %#ok<AGROW>
   end
end

function text = formatValue(value)
   %FORMATVALUE Format one scalar table value for Markdown.

   if iscell(value)
      value = value{1};
   end
   if isdatetime(value)
      if isnat(value)
         text = "NA";
      else
         text = string(value, "yyyy-MM-dd HH:mm:ss z");
      end
   elseif isnumeric(value)
      if isempty(value) || ~isfinite(value)
         text = "NA";
      else
         text = string(sprintf("%.5g", value));
      end
   elseif islogical(value)
      text = string(value);
   else
      text = join(string(value), ", ");
   end
   text = replace(replace(text, "|", "\|"), newline, " ");
end

function caption = assetCaption(asset)
   %ASSETCAPTION Return a concise caption from one stable asset name.

   [~, name] = fileparts(asset);
   caption = replace(string(name), ["-", "_"], " ");
end

function exportAndClose(fig, filename)
   %EXPORTANDCLOSE Export one report figure and release its graphics state.

   cleanup = onCleanup(@() close(fig));
   exportgraphics(fig, filename, Resolution=160)
end

function [fig, ax] = newReportFigure(n_rows)
   %NEWREPORTFIGURE Create a light report figure sized for its case rows.

   height_px = max(520, 58 * n_rows + 180);
   fig = figure(Visible="off", Color="w", ...
      Position=[100 100 1300 height_px]);
   ax = axes(fig);
   formatReportAxes(ax)
end

function configureCaseAxis(ax, case_ids)
   %CONFIGURECASEAXIS Label horizontal case bars without clipping identifiers.

   yticks(ax, 1:numel(case_ids))
   yticklabels(ax, string(case_ids))
   ax.YDir = "reverse";
   ax.TickLabelInterpreter = "none";
   formatReportAxes(ax)
end

function formatReportAxes(ax)
   %FORMATREPORTAXES Isolate exported figures from interactive theme defaults.

   ax.Color = "w";
   ax.XColor = "k";
   ax.YColor = "k";
   ax.GridColor = [0.65 0.65 0.65];
   ax.GridAlpha = 0.25;
   ax.FontSize = 11;
   ax.Box = "off";
   ax.Title.Color = "k";
   ax.XLabel.Color = "k";
   ax.YLabel.Color = "k";
end

function formatReportLegend(lgd)
   %FORMATREPORTLEGEND Isolate exported legends from interactive themes.

   lgd.Color = "w";
   lgd.TextColor = "k";
   lgd.EdgeColor = [0.75 0.75 0.75];
end

function quoted = shellQuote(value)
   %SHELLQUOTE Quote one local path for the Quarto shell command.

   assert(~contains(value, char(34)), ...
      "icemodel:verification:report:unsupportedPath", ...
      "Report paths containing double quotes are unsupported.")
   quoted = string(sprintf('"%s"', char(value)));
end
