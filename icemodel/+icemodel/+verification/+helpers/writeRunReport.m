function report_path = writeRunReport(run_dir, case_results, cases, kwargs)
   %WRITERUNREPORT Write a concise markdown report for a verification run.
   %
   %  report_path = icemodel.verification.helpers.writeRunReport(...
   %     run_dir, case_results, cases)
   %  report_path = icemodel.verification.helpers.writeRunReport(...
   %     run_dir, case_results, cases, run_icemodel=true)
   %
   %  Produces a single human-readable report at <run_dir>/report.md
   %  that ties the run's case-level metrics, comparison figures, and
   %  scatter figures together. Designed for the snow-model
   %  development feedback loop: an agent or developer can open this
   %  one file after a run and immediately see (1) whether all cases
   %  ran, (2) per-case metric headlines, and (3) where to find the
   %  underlying figures and per-case CSVs.
   %
   %  Inputs
   %    run_dir       Run artifact directory.
   %    case_results  Cell array of per-case result structs from
   %                  icemodel.verification.comparecase. Each carries
   %                  a metric table; the report iterates over them
   %                  directly so the long-format summary table held
   %                  by the runner does not need to be passed in.
   %    cases         Resolved case manifest array from listcases.
   %
   %  Name-value
   %    run_name : string (default = basename of run_dir)
   %    run_icemodel : logical (default false)
   %        Mark the report as a synthetic-candidate run.
   %    plotted : logical (default false)
   %        True when comparison / scatter figures were produced and
   %        saved alongside the report. Controls whether per-case
   %        figure links are included.
   %
   %  Returns
   %    report_path : string
   %        Absolute path to the written report.

   arguments
      run_dir       (1, 1) string
      case_results  cell
      cases         struct
      kwargs.run_name     (1, 1) string = ""
      kwargs.run_icemodel (1, 1) logical = false
      kwargs.plotted      (1, 1) logical = false
   end

   if kwargs.run_name == ""
      [~, run_name] = fileparts(char(run_dir));
      kwargs.run_name = string(run_name);
   end

   report_path = fullfile(run_dir, "report.md");
   fid = fopen(report_path, 'w');
   if fid < 0
      error('icemodel:verification:writeRunReport:openFailed', ...
         'unable to open %s for writing', report_path);
   end
   cleanup = onCleanup(@() fclose(fid));

   % --- Header ---------------------------------------------------------
   fprintf(fid, "# Snow verification run %s\n\n", kwargs.run_name);
   fprintf(fid, "Run directory: `%s`\n\n", run_dir);
   if kwargs.run_icemodel
      fprintf(fid, "**Candidate source:** synthetic-snow hook ");
      fprintf(fid, "(`icemodel.verification.runIcemodelSnowCandidate`).\n");
      fprintf(fid, "Metrics include the deliberate perturbations:\n");
      fprintf(fid, "snow_depth +0.02 m, swe x 1.05, surface_temp +0.25 K, ");
      fprintf(fid, "liquid_water x 1.05. ");
      fprintf(fid, "The Colbeck case uses real IceModel candidates and is unaffected.\n\n");
   else
      fprintf(fid, "**Candidate source:** staged smoke reference (no model run).\n\n");
   end

   % --- Per-case headline summary -------------------------------------
   fprintf(fid, "## Per-case headline\n\n");
   fprintf(fid, "| Case | Family | Window | Variables | OK | NA |\n");
   fprintf(fid, "|------|--------|--------|-----------|----|----|\n");
   for i = 1:numel(case_results)
      cr = case_results{i};
      m = cr.metrics;
      ok = nnz(m.status == "ok");
      na = nnz(m.status == "not_applicable");
      family = caseFamily(cases, cr.case_id);
      window = caseWindow(cases, cr.case_id);
      fprintf(fid, "| %s | %s | %s | %d | %d | %d |\n", ...
         cr.case_id, family, window, height(m), ok, na);
   end
   fprintf(fid, "\n");

   % --- Per-variable metrics by case ----------------------------------
   fprintf(fid, "## Metrics\n\n");
   fprintf(fid, "Bias = mean(candidate - target). RMSE / correlation are ");
   fprintf(fid, "evaluated on aligned finite pairs. Peak / melt-out timing ");
   fprintf(fid, "errors are signed offsets in hours.\n\n");
   for i = 1:numel(case_results)
      cr = case_results{i};
      fprintf(fid, "### `%s`\n\n", cr.case_id);
      m = cr.metrics;
      if isempty(m)
         fprintf(fid, "_no metrics_\n\n");
         continue
      end
      fprintf(fid, "| Variable | Status | n | Bias | RMSE | Corr | Peak Err | Peak ΔT (h) | Onset ΔT (h) | Melt-out ΔT (h) |\n");
      fprintf(fid, "|----------|--------|---|------|------|------|----------|-------------|--------------|------------------|\n");
      for k = 1:height(m)
         row = m(k, :);
         fprintf(fid, "| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s |\n", ...
            string(getOr(row, 'variable', '')), ...
            string(row.status), ...
            row.n, ...
            fmtNum(getOr(row, 'bias', NaN)), ...
            fmtNum(getOr(row, 'rmse', NaN)), ...
            fmtNum(getOr(row, 'correlation', NaN)), ...
            fmtNum(getOr(row, 'peak_error', NaN)), ...
            fmtNum(getOr(row, 'peak_time_error_hours', NaN)), ...
            fmtNum(getOr(row, 'snow_onset_time_error_hours', NaN)), ...
            fmtNum(getOr(row, 'melt_out_time_error_hours', NaN)));
      end
      fprintf(fid, "\n");

      if kwargs.plotted
         fig_paths = collectFigurePaths(run_dir, cr.case_id);
         if ~isempty(fig_paths)
            fprintf(fid, "**Figures:**\n\n");
            for p = fig_paths(:)'
               fprintf(fid, "- `%s`\n", p);
            end
            fprintf(fid, "\n");
         end
      end
   end

   % --- Footer ---------------------------------------------------------
   fprintf(fid, "## Artifacts\n\n");
   fprintf(fid, "- `summary.csv` (long-format metric table)\n");
   fprintf(fid, "- `summary.mat` (`summary` + per-case `case_results` cell array)\n");
   for i = 1:numel(case_results)
      cr = case_results{i};
      if ~isblanktext(cr.metrics_path)
         [~, base, ext] = fileparts(char(cr.metrics_path));
         fprintf(fid, "- `%s/%s%s`\n", cr.case_id, base, ext);
      end
   end
   fprintf(fid, "\n");
end

% =====================================================================
% Local helpers
% =====================================================================

function family = caseFamily(cases, case_id)
   match = string({cases.case_id}) == string(case_id);
   if any(match)
      c = cases(find(match, 1));
      family = string(c.dataset_family);
   else
      family = "?";
   end
end

function window = caseWindow(cases, case_id)
   match = string({cases.case_id}) == string(case_id);
   if any(match)
      c = cases(find(match, 1));
      if isfield(c, 'comparison_window') && isstruct(c.comparison_window)
         window = sprintf("%s..%s", string(c.comparison_window.start), ...
            string(c.comparison_window.end));
      else
         window = "?";
      end
   else
      window = "?";
   end
end

function paths = collectFigurePaths(run_dir, case_id)
   case_dir = fullfile(run_dir, char(case_id));
   if exist(case_dir, 'dir') ~= 7
      paths = strings(0, 1);
      return
   end
   d = dir(fullfile(case_dir, '*.png'));
   if isempty(d)
      paths = strings(0, 1);
      return
   end
   paths = string({d.name})';
   paths = arrayfun(@(p) fullfile(string(case_id), p), paths);
end

function v = getOr(row, name, default)
   if ismember(name, row.Properties.VariableNames)
      v = row.(name);
   else
      v = default;
   end
end

function s = fmtNum(v)
   if isnan(v)
      s = "—";
      return
   end
   if abs(v) < 1e-3 && v ~= 0
      s = string(sprintf('%.3e', v));
   else
      s = string(sprintf('%.3g', v));
   end
end
