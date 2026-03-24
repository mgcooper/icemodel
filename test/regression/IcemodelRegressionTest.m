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
      env_cleanup
   end

   methods (TestMethodSetup)
      function configurePaths(testCase)
         % Install the canonical test config so direct class runs and
         % runner-based runs see the same environment.
         [~, ~, ~, ~, testCase.env_cleanup] = ...
            icemodel.test.helpers.bootstrapTestEnvironment();
      end
   end

   methods (TestMethodTeardown)
      function restorePaths(testCase)
         % Release the setup cleanup handle after each regression case.
         testCase.env_cleanup = [];
      end
   end

   methods (Test)
      function testCoreRegression(testCase)
         % Read the requested case matrix/baseline selectors from the
         % runner-owned environment variables.
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

         % Resolve the retained formal cases and the matched accepted baseline
         % before entering the per-case compare loop.
         runoff_ref = icemodel.test.helpers.loadRunoffReference();
         cases = icemodel.test.helpers.getRegressionCaseMatrix( ...
            tier=tier, smbmodel=smbmodel, solver=solver, simyear=simyear, ...
            smoke_sites=smoke_sites, full_sites=full_sites);

         % Confirm the case exists.
         testCase.assertNotEmpty(cases, ...
            sprintf('no regression cases matched tier=%s smbmodel=%s', ...
            tier, smbmodel));

         % Load the baseline.
         baseline = icemodel.test.helpers.loadRegressionBaseline( ...
            baseline_tag, smbmodel);

         % Accumulate the compare report and resolved opts for one saved
         % artifact per regression run.
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

            % Summarize the retained output years against the matched runoff
            % reference row, if one exists for this formal case.
            idx = icemodel.test.helpers.findRunoffReferenceRow(runoff_ref, c);
            met = icemodel.test.helpers.loadProcessedMetForOutputYears(opts);
            if isempty(idx)
               refrow = [];
            else
               refrow = runoff_ref(idx, :);
            end
            S = icemodel.test.helpers.summarizeIce1Metrics(ice1, met, refrow);

            % Metrics with baseline lookup and delta/pct-delta pairs.
            delta_specs = {
               'runoff_final',  'runoff_delta',      'runoff_pct_delta'
               'melt_final',    'melt_delta',        'melt_pct_delta'
               'runoff_eval',   'runoff_eval_delta', 'runoff_eval_pct_delta'
               'melt_eval',     'melt_eval_delta',   'melt_eval_pct_delta'
               };

            % Metrics with baseline lookup only (no delta pair).
            baseline_only = ["mean_Tice_numiter", "max_Tice_numiter", ...
               "n_not_converged"];
            all_baseline_fields = [ ...
               string(delta_specs(:, 1)); baseline_only(:)];

            % Initialize fields nan
            base = struct();
            for f = all_baseline_fields'
               base.(f) = nan;
            end

            % Load the accepted baseline values for this case.
            bid = icemodel.test.helpers.findCaseRow(baseline, string(c.case_id));

            % Populate the fields
            if ~isempty(bid)
               for f = all_baseline_fields'
                  base.(f) = testCase.getBaselineValue(baseline, bid, f);
               end
               metric_names = string(fieldnames(S));
               for imetric = 1:numel(metric_names)
                  actual = S.(metric_names(imetric));
                  if isfinite(actual)
                     testCase.checkAgainstBaseline(actual, ...
                        baseline, bid, metric_names(imetric));
                  end
               end
            end

            % Build the report row with case identity, current metrics,
            % baseline values, and computed deltas.
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

            for f = all_baseline_fields'
               row.("baseline_" + f) = base.(f);
            end

            % Compute absolute and percent deltas for each tracked metric pair
            % using the delta_specs table. computeDelta returns NaN when the
            % baseline value is NaN (case not found) or zero (undefined pct).
            for i = 1:size(delta_specs, 1)
               [row.(delta_specs{i, 2}), row.(delta_specs{i, 3})] = ...
                  IcemodelRegressionTest.computeDelta( ...
                  S.(delta_specs{i, 1}), base.(delta_specs{i, 1}));
            end
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

         % Save the artifact and display the comparison summary.
         artifact_file = testCase.saveArtifacts(report, case_opts, meta);
         disp(artifact_file)
         icemodel.test.helpers.displayRegressionReport(report)
      end

   end

   methods (Access = private)
      function row = copyMetricFields(~, row, S)
         % Copy the scalar metric struct into one artifact row.
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
         tol = testCase.metricTolerance(char(varname), expected);
         testCase.verifyLessThanOrEqual(abs(actual - expected), tol, ...
            sprintf('baseline mismatch var=%s', varname));
      end

      function x = getBaselineValue(~, baseline, row, varname)
         % Read one baseline value when the saved table carries that column.
         if ismember(varname, baseline.Properties.VariableNames)
            x = baseline.(varname)(row);
         else
            x = nan;
         end
      end

      function tol = metricTolerance(testCase, varname, expected)
         % Return one metric-specific scalar tolerance.

         if endsWith(string(varname), "_m3")
            tol = max(testCase.abs_tol_runoff_m3, ...
               testCase.rel_tol_runoff_m3 * abs(expected));
            return
         end

         if any(startsWith(string(varname), ["runoff_", "melt_"]))
            tol = max(1e-4, 1e-4 * abs(expected));
            return
         end

         if contains(varname, "numiter")
            tol = 0.5;
            return
         end

         if contains(varname, "not_converged")
            tol = 1.0;
            return
         end

         if startsWith(varname, "closure_") || startsWith(varname, "gof_")
            tol = 5e-3;
            return
         end

         tol = max(testCase.abs_tol_scalar, ...
            testCase.rel_tol_scalar * abs(expected));
      end

      function artifact_file = saveArtifacts(~, report, case_opts, meta)
         %saveArtifacts Save the regression comparison artifact for one run.

         % Create the run-specific artifact folder before saving the report.
         testdir = icemodel.getpath('test');
         outdir = fullfile(testdir, 'artifacts', char(meta.run_name));
         if exist(outdir, 'dir') ~= 7
            mkdir(outdir);
         end

         % Format the baseline/model/solver tags used by the saved filename.
         if meta.baseline_tag == ""
            baseline_tag = 'nobaseline';
         else
            baseline_tag = char( ...
               icemodel.test.helpers.sanitizeTag(meta.baseline_tag));
         end
         model_tag = IcemodelRegressionTest.smbmodelLabel( ...
            meta.smbmodel_filter);
         solver_tag = IcemodelRegressionTest.solverLabel( ...
            meta.solver_filter);

         % Build the artifact filename.
         artifact_file = fullfile(outdir, ...
            sprintf('regression_report_%s%s_%s.mat', ...
            char(meta.tier), char(model_tag + solver_tag), baseline_tag));

         % Save the artifact file.
         save(artifact_file, 'report', 'case_opts', 'meta');
      end
   end

   methods (Static, Access = private)
      function s = getenvRequired(name)
         %GETENVREQUIRED Read one required regression env var or error cleanly.
         %
         % The test bootstrap installs config paths, but the runner still
         % owns the concrete ICEMODEL_TEST_* selectors for one compare run.
         s = getenv(name);
         if isempty(s)
            error('missing required regression env var: %s', name)
         end
      end

      function meta = reportMeta(tier, smbmodel, solver, simyear, ...
            smoke_sites, full_sites, baseline_tag, run_name)
         %REPORTMETA Build the saved metadata struct for one regression artifact.

         % Resolve the shared run stamp before filling the saved metadata.
         [run_date, run_id, run_name] = ...
            icemodel.test.helpers.resolveRunStamp(run_name);

         % Record the compare selector, tolerances, and resolved path context.
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
         meta.metric_tolerance_policy = ...
            "metric-specific tolerances for runoff/melt, iteration, closure, and GOF metrics";
         meta.case_builder = "icemodel.test.helpers.setModelOptsForCase";
         meta.opts_source = "icemodel.setopts defaults";
         meta.spinup_policy = ...
            "formal regression runs include the canonical leading spinup year";
         meta.reset_fields = "solver";
         meta.input_path = string(icemodel.getpath('input'));
         meta.output_path = string(icemodel.getpath('output'));
         meta.suite_file = string(mfilename('fullpath'));
         meta.matlab_version = string(version);
         meta.host = string(computer);
         meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
      end

      function label = smbmodelLabel(smbmodel)
         %SMBMODELLABEL Format the smbmodel selector for artifact filenames.
         smbmodel = string(smbmodel);
         if any(strcmpi(smbmodel, "all"))
            label = "";
         else
            label = "_" + icemodel.test.helpers.smbmodelTag(smbmodel);
         end
      end

      function x = getenvAsNumbers(name)
         %GETENVASNUMBERS Parse a numeric vector env var.
         s = getenv(name);
         if isempty(s)
            x = [];
         else
            x = str2double(split(string(s), ','));
            x = x(isfinite(x));
         end
      end

      function x = getenvAsStrings(name)
         %GETENVASSTRINGS Parse a comma-delimited string env var.
         s = getenv(name);
         if isempty(s)
            x = strings(0, 1);
         else
            x = split(string(s), ',');
            x = x(x ~= "");
         end
      end

      function label = solverLabel(solver)
         %SOLVERLABEL Format the solver filter for artifact filenames.
         if isempty(solver)
            label = "";
         else
            label = "_s" + join(string(solver), '-');
         end
      end

      function x = scalarRelTol()
         %SCALARRELTOL Return the scalar relative regression tolerance.
         x = 1e-6;
      end

      function x = scalarAbsTol()
         %SCALARABSTOL Return the scalar absolute regression tolerance.
         x = 1e-9;
      end

      function x = runoffRelTol()
         %RUNOFFRELTOL Return the runoff relative regression tolerance.
         x = 1e-4;
      end

      function x = runoffAbsTol()
         %RUNOFFABSTOL Return the runoff absolute regression tolerance.
         x = 1.0;
      end
   end

   methods (Static, Access = private)
      function [abs_delta, pct_delta] = computeDelta(actual, baseline)
         %COMPUTEDELTA Compute absolute and percent change against a baseline.
         abs_delta = nan;
         pct_delta = nan;
         if isfinite(actual) && isfinite(baseline)
            abs_delta = actual - baseline;
            if baseline ~= 0
               pct_delta = 100 * abs_delta / baseline;
            end
         end
      end
   end
end
