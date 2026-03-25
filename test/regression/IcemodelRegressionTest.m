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
      caseinfo
   end

   methods (TestMethodSetup)
      function configureCases(testCase)
         % Install the canonical test config so direct class runs and
         % runner-based runs see the same environment.
         [~, ~, ~, ~, testCase.env_cleanup] = ...
            icemodel.test.helpers.bootstrapTestEnvironment();

         % Resolve the case matrix, baseline, and runoff reference once
         % before the per-case compare loop.
         testCase.caseinfo = buildRegressionCaseInfo();
      end
   end

   methods (TestMethodTeardown)
      function restoreConfig(testCase)
         % Release the setup cleanup handle after each regression case.
         testCase.env_cleanup = [];
      end
   end

   methods (Test)
      function testCoreRegression(testCase)
         % Unpack the pre-resolved case matrix, baseline, and runoff
         % reference from the caseinfo struct built in configureCases.
         cases = testCase.caseinfo.cases;
         baseline = testCase.caseinfo.baseline;
         runoff_ref = testCase.caseinfo.runoff_ref;
         baseline_tag = testCase.caseinfo.baseline_tag;

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
         meta = IcemodelRegressionTest.reportMeta(testCase.caseinfo);

         % Save the artifact for the runner-level display and results struct.
         testCase.saveArtifacts(report, case_opts, meta);
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

      function saveArtifacts(~, report, case_opts, meta)
         %saveArtifacts Save the regression comparison artifact for one run.

         % Build the canonical artifact path.
         artifact_file = icemodel.test.helpers.artifactFilePath( ...
            "regression", tier=meta.tier, smbmodel=meta.smbmodel_filter, ...
            solver=meta.solver_filter, ...
            baseline_tag=meta.baseline_tag, run_name=meta.run_name);

         % Create the run-specific artifact folder before saving.
         outdir = fileparts(artifact_file);
         if exist(outdir, 'dir') ~= 7
            mkdir(outdir);
         end

         % Save the artifact file.
         save(artifact_file, 'report', 'case_opts', 'meta');

         % Export the artifact path so the runner can load it after the
         % unittest framework returns control.
         setenv('ICEMODEL_REGRESSION_ARTIFACT_FILE', artifact_file);
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

      function meta = reportMeta(info)
         %REPORTMETA Build the metadata struct for one regression artifact.
         %
         % Accepts the caseinfo struct built by buildRegressionCaseInfo, which
         % carries all the selector and site fields needed for the artifact
         % metadata.

         % Resolve the shared run stamp before filling the saved metadata.
         [run_date, run_id, run_name] = ...
            icemodel.test.helpers.resolveRunStamp(info.run_name);

         % Record the compare selector, tolerances, and resolved path context.
         meta = struct();
         meta.tier = string(info.tier);
         meta.smbmodel_filter = string(info.smbmodel);
         meta.solver_filter = info.solver;
         meta.simyear = info.simyear;
         meta.smoke_sites = info.smoke_sites;
         meta.full_sites = info.full_sites;
         meta.baseline_tag = string(info.baseline_tag);
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

function info = buildRegressionCaseInfo()
   %BUILDREGRESSIONCASEINFO Resolve the full regression case matrix and baseline.
   %
   % Reads the suite-selection env vars installed by the runner (or by
   % bootstrapTestEnvironment for direct class runs), builds the case matrix,
   % loads the baseline and runoff reference, and returns a struct carrying
   % everything testCoreRegression needs.

   % Read the suite selector from environment variables.
   tier = getenvOrDefault('ICEMODEL_TEST_TIER', 'smoke');
   smbmodel = getenvOrDefault('ICEMODEL_TEST_SMBMODEL_FILTER', 'all');
   solver = parseNumericEnv('ICEMODEL_TEST_SOLVER_FILTER');
   simyear = str2double(getenvOrDefault('ICEMODEL_TEST_SIMYEAR_FILTER', '2016'));
   smoke_sites = parseStringEnv('ICEMODEL_TEST_SMOKE_SITES', "kanm");
   full_sites = parseStringEnv('ICEMODEL_TEST_FULL_SITES', ["kanm"; "kanl"]);
   baseline_tag = getenvOrDefault('ICEMODEL_REGRESSION_BASELINE', 'rolling');
   run_name = string(getenv('ICEMODEL_TEST_RUN_NAME'));

   % Build the formal regression case matrix from the resolved selectors.
   cases = icemodel.test.helpers.getRegressionCaseMatrix( ...
      tier=tier, smbmodel=smbmodel, solver=solver, simyear=simyear, ...
      smoke_sites=smoke_sites, full_sites=full_sites);

   assert(~isempty(cases), ...
      'regression case matrix is empty for tier=%s smbmodel=%s solver=[%s]', ...
      tier, smbmodel, join(string(solver), ','));

   % Load the accepted baseline and runoff reference for comparison.
   baseline = icemodel.test.helpers.loadBaseline("regression", ...
      smbmodel=smbmodel, baseline_tag=baseline_tag);
   runoff_ref = icemodel.test.helpers.loadReference("runoff");

   % Pack everything into a single struct for testCoreRegression.
   info = struct();
   info.tier = string(tier);
   info.smbmodel = string(smbmodel);
   info.solver = solver;
   info.simyear = simyear;
   info.smoke_sites = smoke_sites;
   info.full_sites = full_sites;
   info.baseline_tag = string(baseline_tag);
   info.run_name = run_name;
   info.cases = cases;
   info.baseline = baseline;
   info.runoff_ref = runoff_ref;
end

function s = getenvOrDefault(name, default)
   %GETENVORDEFAULT Read an env var with a fallback default.
   s = getenv(name);
   if isempty(s)
      s = default;
   end
end

function x = parseNumericEnv(name)
   %PARSENUMERICENV Parse a comma-delimited numeric env var.
   s = getenv(name);
   if isempty(s)
      x = [];
   else
      x = str2double(split(string(s), ','));
      x = x(isfinite(x));
   end
end

function x = parseStringEnv(name, default)
   %PARSESTRINGENV Parse a comma-delimited string env var with a fallback.
   s = getenv(name);
   if isempty(s)
      x = default;
   else
      x = split(string(s), ',');
      x = x(x ~= "");
   end
end
