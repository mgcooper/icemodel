classdef IcemodelRegressionTest < matlab.unittest.TestCase
   %ICEMODELREGRESSIONTEST Numerical regression checks for formal point runs.
   %
   % This class runs the compact formal regression matrix, compares the core
   % scalar outputs to the requested baseline, and writes one artifact table
   % summarizing the modeled runoff/melt context for each case.

   properties
      rel_tol_scalar (1, 1) double = 1e-6
      abs_tol_scalar (1, 1) double = 1e-9
      rel_tol_runoff_m3 (1, 1) double = 1e-4
      abs_tol_runoff_m3 (1, 1) double = 1.0
   end

   methods (TestMethodSetup)
      function configurePaths(~)
         % Configure the public demo/data paths used by the formal suite.
         IcemodelRegressionTest.ensureModelConfigPaths();
      end
   end

   methods (Test)
      function testCoreRegression(testCase)
         % Read the requested case matrix/baseline from the runner env vars.
         tier = IcemodelRegressionTest.getenvRequired('ICEMODEL_TEST_TIER');
         smbmodel = IcemodelRegressionTest.getenvRequired( ...
            'ICEMODEL_TEST_SMBMODEL_FILTER');
         solver = IcemodelRegressionTest.getenvAsNumbers( ...
            'ICEMODEL_TEST_SOLVER_FILTER');
         simyear = str2double(IcemodelRegressionTest.getenvRequired( ...
            'ICEMODEL_TEST_SIMYEAR_FILTER'));
         smoke_sites = IcemodelRegressionTest.getenvAsStrings( ...
            'ICEMODEL_TEST_SMOKE_SITES');
         full_sites = IcemodelRegressionTest.getenvAsStrings( ...
            'ICEMODEL_TEST_FULL_SITES');
         baseline_tag = getenv('ICEMODEL_REGRESSION_BASELINE');
         run_name = getenv('ICEMODEL_TEST_RUN_NAME');

         runoff_ref = icemodel.test.helpers.loadRunoffReference();
         cases = icemodel.test.helpers.getRegressionCaseMatrix( ...
            tier=tier, smbmodel=smbmodel, solver=solver, simyear=simyear, ...
            smoke_sites=smoke_sites, full_sites=full_sites);
         testCase.assertNotEmpty(cases, ...
            sprintf('no regression cases matched tier=%s smbmodel=%s', ...
            tier, smbmodel));
         baseline = icemodel.test.helpers.loadRegressionBaseline(baseline_tag, ...
            smbmodel);

         report_rows = struct([]);
         case_opts = struct([]);
         r = 0;

         % Run each formal case, compare to baseline, and collect a report row.
         for icase = 1:height(cases)
            c = cases(icase, :);
            opts = icemodel.test.helpers.setModelOptsForCase(c);

            [ice1, ice2] = icemodel.test.helpers.runSmbModel(opts);
            [ice1, ice2] = icemodel.postprocess( ...
               ice1, ice2, opts, opts.output_years); %#ok<ASGLU>

            ridx = icemodel.test.helpers.findRunoffReferenceRow(runoff_ref, c);
            met = icemodel.test.helpers.loadProcessedMetForOutputYears(opts);
            if isempty(ridx)
               refrow = [];
            else
               refrow = runoff_ref(ridx, :);
            end
            S = icemodel.test.helpers.summarizeIce1Metrics(ice1, met, refrow);

            % Init values
            [base_runoff_final, base_melt_final, base_runoff_eval, ...
               base_melt_eval, runoff_delta, melt_delta, runoff_eval_delta, ...
               melt_eval_delta, runoff_pct_delta, melt_pct_delta, ...
               runoff_eval_pct_delta, melt_eval_pct_delta, ...
               base_mean_Tice_numiter, base_max_Tice_numiter, ...
               base_n_not_converged] = deal(nan);

            bid = icemodel.test.helpers.findCaseRow(baseline, string(c.case_id));
            if ~isempty(bid)
               base_runoff_final = ...
                  testCase.getBaselineValue(baseline, bid, 'runoff_final');
               base_melt_final = ...
                  testCase.getBaselineValue(baseline, bid, 'melt_final');
               base_runoff_eval = ...
                  testCase.getBaselineValue(baseline, bid, 'runoff_eval');
               base_melt_eval = ...
                  testCase.getBaselineValue(baseline, bid, 'melt_eval');
               base_mean_Tice_numiter = ...
                  testCase.getBaselineValue(baseline, bid, 'mean_Tice_numiter');
               base_max_Tice_numiter = ...
                  testCase.getBaselineValue(baseline, bid, 'max_Tice_numiter');
               base_n_not_converged = ...
                  testCase.getBaselineValue(baseline, bid, 'n_not_converged');
               metric_names = string(fieldnames(S));
               for imetric = 1:numel(metric_names)
                  this_metric = metric_names(imetric);
                  actual = S.(char(this_metric));
                  if isfinite(actual)
                     testCase.checkAgainstBaseline(actual, ...
                        baseline, bid, char(this_metric));
                  end
               end
            end

            if isfinite(base_runoff_final)
               runoff_delta = S.runoff_final - base_runoff_final;
               runoff_pct_delta = IcemodelRegressionTest.computePctDelta( ...
                  S.runoff_final, base_runoff_final);
            end
            if isfinite(base_melt_final)
               melt_delta = S.melt_final - base_melt_final;
               melt_pct_delta = IcemodelRegressionTest.computePctDelta( ...
                  S.melt_final, base_melt_final);
            end
            if isfinite(base_runoff_eval)
               runoff_eval_delta = ...
                  S.runoff_eval - base_runoff_eval;
               runoff_eval_pct_delta = ...
                  IcemodelRegressionTest.computePctDelta( ...
                  S.runoff_eval, base_runoff_eval);
            end
            if isfinite(base_melt_eval)
               melt_eval_delta = ...
                  S.melt_eval - base_melt_eval;
               melt_eval_pct_delta = ...
                  IcemodelRegressionTest.computePctDelta( ...
                  S.melt_eval, base_melt_eval);
            end

            row = struct(); %#ok<*AGROW>
            row.case_id = string(c.case_id);
            row.tier = string(c.tier);
            row.baseline_tag = string(baseline_tag);
            row.smbmodel = string(c.smbmodel);
            row.sitename = string(c.sitename);
            row.forcings = string(c.forcings);
            row.simyear = c.simyear;
            row.solver = c.solver;
            row = testCase.copyMetricFields(row, S);
            row.baseline_runoff_final = base_runoff_final;
            row.baseline_melt_final = base_melt_final;
            row.baseline_runoff_eval = base_runoff_eval;
            row.baseline_melt_eval = base_melt_eval;
            row.baseline_mean_Tice_numiter = base_mean_Tice_numiter;
            row.baseline_max_Tice_numiter = base_max_Tice_numiter;
            row.baseline_n_not_converged = base_n_not_converged;
            row.runoff_delta = runoff_delta;
            row.melt_delta = melt_delta;
            row.runoff_eval_delta = runoff_eval_delta;
            row.melt_eval_delta = melt_eval_delta;
            row.runoff_pct_delta = runoff_pct_delta;
            row.melt_pct_delta = melt_pct_delta;
            row.runoff_eval_pct_delta = runoff_eval_pct_delta;
            row.melt_eval_pct_delta = melt_eval_pct_delta;
            row.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
            report_rows = vertcat(report_rows, row);

            r = r + 1;

            case_opts(r).case_id = string(c.case_id);
            case_opts(r).case = table2struct(c);
            case_opts(r).opts = opts;
         end

         report = struct2table(report_rows);
         meta = IcemodelRegressionTest.reportMeta( ...
            tier, smbmodel, solver, simyear, smoke_sites, full_sites, ...
            baseline_tag, run_name);
         % Save one artifact for this regression compare run.
         testCase.logReport(report, case_opts, meta);
      end

   end

   methods (Access = private)
      function row = copyMetricFields(~, row, S)
         names = string(fieldnames(S));
         for i = 1:numel(names)
            name = char(names(i));
            row.(name) = S.(name);
         end
      end

      function checkAgainstBaseline(testCase, actual, baseline, row, varname)
         % Compare one scalar regression metric against the selected baseline.
         if ~ismember(varname, baseline.Properties.VariableNames)
            return
         end
         expected = baseline.(varname)(row);
         if ~isfinite(expected)
            return
         end
         if endsWith(string(varname), "_m3")
            tol = max(testCase.abs_tol_runoff_m3, ...
               testCase.rel_tol_runoff_m3 * abs(expected));
         else
            tol = max(testCase.abs_tol_scalar, ...
               testCase.rel_tol_scalar * abs(expected));
         end
         testCase.verifyLessThanOrEqual(abs(actual - expected), tol, ...
            sprintf('baseline mismatch var=%s', varname));
      end

      function x = getBaselineValue(~, baseline, row, varname)
         if ismember(varname, baseline.Properties.VariableNames)
            x = baseline.(varname)(row);
         else
            x = nan;
         end
      end

      function logReport(~, report, case_opts, meta)
         % Save the full regression artifact, then print only the baseline-vs-
         % current summary that is useful during development.
         rootdir = IcemodelRegressionTest.repoRoot();
         outdir = fullfile(rootdir, 'test', 'artifacts', char(meta.run_name));
         if exist(outdir, 'dir') ~= 7
            mkdir(outdir);
         end
         if meta.baseline_tag == ""
            baseline_tag = 'nobaseline';
         else
            baseline_tag = char(icemodel.test.helpers.sanitizeTag(meta.baseline_tag));
         end
         model_tag = IcemodelRegressionTest.smbmodelLabel(meta.smbmodel_filter);
         solver_tag = IcemodelRegressionTest.solverLabel(meta.solver_filter);
         outfile = fullfile(outdir, ...
            sprintf('regression_report_%s%s_%s.mat', ...
            char(meta.tier), char(model_tag + solver_tag), baseline_tag));
         save(outfile, 'report', 'case_opts', 'meta');
         disp(outfile)
         disp(makeDisplayTable(report, ...
            'runoff_final', 'baseline_runoff_final', 'runoff_pct_delta', ...
            'runoff', 'baseline', 'pct_delta'))
         fprintf('\n')
         disp(makeDisplayTable(report, ...
            'melt_final', 'baseline_melt_final', 'melt_pct_delta', ...
            'melt', 'baseline', 'pct_delta'))
         fprintf('\n')
         disp(makeDisplayTable(report, ...
            'runoff_eval', 'baseline_runoff_eval', ...
            'runoff_eval_pct_delta', ...
            'runoff_eval', 'baseline_eval', 'pct_delta'))
         fprintf('\n')
         disp(makeDisplayTable(report, ...
            'melt_eval', 'baseline_melt_eval', ...
            'melt_eval_pct_delta', ...
            'melt_eval', 'baseline_eval', 'pct_delta'))
      end
   end

   methods (Static, Access = private)
      function ensureModelConfigPaths()
         rootdir = IcemodelRegressionTest.repoRoot();
         addpath(fullfile(rootdir, 'test'));
         icemodel.test.helpers.configureModelPaths(rootdir);
      end

      function rootdir = repoRoot()
         rootdir = icemodel.internal.fullpath();
      end

      function s = getenvRequired(name)
         s = getenv(name);
         if isempty(s)
            error('missing required regression env var: %s', name)
         end
      end

      function meta = reportMeta(tier, smbmodel, solver, simyear, ...
            smoke_sites, full_sites, baseline_tag, run_name)

         [run_date, run_id, run_name] = ...
            icemodel.test.helpers.resolveRunStamp(run_name);

         meta = struct();
         meta.tier = string(tier);
         meta.smbmodel_filter = string(smbmodel);
         meta.solver_filter = solver;
         meta.simyear = simyear;
         meta.smoke_sites = smoke_sites;
         meta.full_sites = full_sites;
         meta.baseline_tag = string(baseline_tag);
         meta.run_date = run_date;
         meta.run_id = run_id;
         meta.run_name = run_name;
         meta.rel_tol_scalar = IcemodelRegressionTest.scalarRelTol();
         meta.abs_tol_scalar = IcemodelRegressionTest.scalarAbsTol();
         meta.rel_tol_runoff_m3 = IcemodelRegressionTest.runoffRelTol();
         meta.abs_tol_runoff_m3 = IcemodelRegressionTest.runoffAbsTol();
         meta.case_builder = "icemodel.test.helpers.setModelOptsForCase";
         meta.opts_source = "icemodel.setopts defaults";
         meta.reset_fields = "solver";
         meta.input_path = string(getenv('ICEMODEL_INPUT_PATH'));
         meta.output_path = string(getenv('ICEMODEL_OUTPUT_PATH'));
         meta.suite_file = string(mfilename('fullpath'));
         meta.matlab_version = string(version);
         meta.host = string(computer);
         meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
      end

      function label = smbmodelLabel(smbmodel)
         smbmodel = string(smbmodel);
         if any(strcmpi(smbmodel, "all"))
            label = "";
         else
            label = "_" + icemodel.test.helpers.smbmodelTag(smbmodel);
         end
      end

      function x = getenvAsNumbers(name)
         s = getenv(name);
         if isempty(s)
            x = [];
         else
            x = str2double(split(string(s), ','));
            x = x(isfinite(x));
         end
      end

      function x = getenvAsStrings(name)
         s = getenv(name);
         if isempty(s)
            x = strings(0, 1);
         else
            x = split(string(s), ',');
            x = x(x ~= "");
         end
      end

      function label = solverLabel(solver)
         if isempty(solver)
            label = "";
         else
            label = "_s" + join(string(solver), '-');
         end
      end

      function x = scalarRelTol()
         x = 1e-6;
      end

      function x = scalarAbsTol()
         x = 1e-9;
      end

      function x = runoffRelTol()
         x = 1e-4;
      end

      function x = runoffAbsTol()
         x = 1.0;
      end
   end

   methods (Static, Access = private)
      function pct = computePctDelta(actual, baseline)
         pct = nan;
         if isfinite(actual) && isfinite(baseline) && baseline ~= 0
            pct = 100 * (actual - baseline) / baseline;
         end
      end
   end
end

function T = makeDisplayTable(report, current_var, baseline_var, pct_var, ...
      current_label, baseline_label, pct_label)

   T = table(report.case_id, report.(current_var), report.(baseline_var), ...
      report.(pct_var), 'VariableNames', ...
      {'case_id', current_label, baseline_label, pct_label});
end
