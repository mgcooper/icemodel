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
         tier = IcemodelRegressionTest.getenvOrDefault( ...
            'ICEMODEL_TEST_TIER', 'smoke');
         smbmodel = IcemodelRegressionTest.getenvOrDefault( ...
            'ICEMODEL_TEST_SMBMODEL_FILTER', 'all');
         baseline_tag = IcemodelRegressionTest.getenvOrDefault( ...
            'ICEMODEL_REGRESSION_BASELINE', '');
         run_name = IcemodelRegressionTest.getenvOrDefault( ...
            'ICEMODEL_TEST_RUN_NAME', '');
         runoff_ref = test.helpers.loadRunoffReference();
         cases = test.helpers.getRegressionCaseMatrix(tier, smbmodel);
         testCase.assertNotEmpty(cases, ...
            sprintf('no regression cases matched tier=%s smbmodel=%s', ...
            tier, smbmodel));
         baseline = test.helpers.loadRegressionBaseline(baseline_tag, ...
            smbmodel);

         report_rows = struct([]);
         case_opts = struct([]);
         r = 0;

         % Run each formal case, compare to baseline, and collect a report row.
         for icase = 1:height(cases)
            c = cases(icase, :);
            opts = test.helpers.buildFormalCaseOpts(c);

            [ice1, ice2] = runModel(opts);
            [ice1, ice2] = POSTPROC(ice1, ice2, opts, c.simyear); %#ok<ASGLU>
            S = test.helpers.summarizeIce1Metrics(ice1);

            bid = test.helpers.findCaseRow(baseline, string(c.case_id));
            if ~isempty(bid)
               testCase.checkAgainstBaseline(S.runoff_final, ...
                  baseline, bid, 'runoff_final');
               testCase.checkAgainstBaseline(S.melt_final, ...
                  baseline, bid, 'melt_final');
               testCase.checkAgainstBaseline(S.mean_Tice_numiter, ...
                  baseline, bid, 'mean_Tice_numiter');
               testCase.checkAgainstBaseline(S.max_Tice_numiter, ...
                  baseline, bid, 'max_Tice_numiter');
               testCase.checkAgainstBaseline(S.n_not_converged, ...
                  baseline, bid, 'n_not_converged');
            end

            ridx = test.helpers.findRunoffReferenceRow(runoff_ref, c);

            obs_final = nan;
            mar_final = nan;
            merra_final = nan;
            racmo_final = nan;
            icemodel_final_m3 = nan;
            obs_diff_m3 = nan;
            mar_diff_m3 = nan;
            merra_diff_m3 = nan;
            racmo_diff_m3 = nan;

            if ~isempty(ridx)
               obs_final = runoff_ref.obs_final_m3(ridx);
               mar_final = runoff_ref.mar_final_m3(ridx);
               merra_final = runoff_ref.merra_final_m3(ridx);
               racmo_final = runoff_ref.racmo_final_m3(ridx);
               icemodel_final_m3 = test.helpers.computeCatchmentRunoffFinal( ...
                  ice1, runoff_ref.area_med_m2(ridx), runoff_ref.t1(ridx), ...
                  runoff_ref.t2(ridx));

               if ~isempty(bid) && ismember('icemodel_final_m3', ...
                     baseline.Properties.VariableNames)

                  ref_icemodel = baseline.icemodel_final_m3(bid);

                  if isfinite(ref_icemodel)
                     tol = max(testCase.abs_tol_runoff_m3, ...
                        testCase.rel_tol_runoff_m3 * abs(ref_icemodel));

                     testCase.verifyLessThanOrEqual( ...
                        abs(icemodel_final_m3 - ref_icemodel), tol, ...
                        sprintf('runoff m3 mismatch case=%s', c.case_id));
                  end
               end

               obs_diff_m3 = icemodel_final_m3 - obs_final;
               mar_diff_m3 = icemodel_final_m3 - mar_final;
               merra_diff_m3 = icemodel_final_m3 - merra_final;
               racmo_diff_m3 = icemodel_final_m3 - racmo_final;
            end

            r = r + 1;
            report_rows(r).case_id = string(c.case_id); %#ok<*AGROW>
            report_rows(r).tier = string(c.tier);
            report_rows(r).baseline_tag = string(baseline_tag);
            report_rows(r).smbmodel = string(c.smbmodel);
            report_rows(r).sitename = string(c.sitename);
            report_rows(r).forcings = string(c.forcings);
            report_rows(r).simyear = c.simyear;
            report_rows(r).solver = c.solver;
            report_rows(r).runoff_final = S.runoff_final;
            report_rows(r).melt_final = S.melt_final;
            report_rows(r).mean_Tice_numiter = S.mean_Tice_numiter;
            report_rows(r).max_Tice_numiter = S.max_Tice_numiter;
            report_rows(r).n_not_converged = S.n_not_converged;
            report_rows(r).obs_final_m3 = obs_final;
            report_rows(r).mar_final_m3 = mar_final;
            report_rows(r).merra_final_m3 = merra_final;
            report_rows(r).racmo_final_m3 = racmo_final;
            report_rows(r).icemodel_final_m3 = icemodel_final_m3;
            report_rows(r).obs_diff_m3 = obs_diff_m3;
            report_rows(r).mar_diff_m3 = mar_diff_m3;
            report_rows(r).merra_diff_m3 = merra_diff_m3;
            report_rows(r).racmo_diff_m3 = racmo_diff_m3;
            report_rows(r).timestamp_utc = datetime('now', 'TimeZone', 'UTC');

            case_opts(r).case_id = string(c.case_id);
            case_opts(r).case = table2struct(c);
            case_opts(r).opts = opts;
         end

         report = struct2table(report_rows);
         meta = IcemodelRegressionTest.reportMeta( ...
            tier, smbmodel, baseline_tag, run_name);
         % Save one artifact for this regression compare run.
         testCase.logReport(report, case_opts, meta);
      end
   end

   methods (Access = private)
      function checkAgainstBaseline(testCase, actual, baseline, row, varname)
         % Compare one scalar regression metric against the selected baseline.
         if ~ismember(varname, baseline.Properties.VariableNames)
            return
         end
         expected = baseline.(varname)(row);
         if ~isfinite(expected)
            return
         end
         tol = max(testCase.abs_tol_scalar, ...
            testCase.rel_tol_scalar * abs(expected));
         testCase.verifyLessThanOrEqual(abs(actual - expected), tol, ...
            sprintf('baseline mismatch var=%s', varname));
      end

      function logReport(~, report, case_opts, meta)
         % Save the regression artifact and print a compact summary table.
         rootdir = IcemodelRegressionTest.repoRoot();
         outdir = fullfile(rootdir, 'test', 'artifacts', char(meta.run_name));
         if exist(outdir, 'dir') ~= 7
            mkdir(outdir);
         end
         if meta.baseline_tag == ""
            baseline_tag = 'nobaseline';
         else
            baseline_tag = char(test.helpers.sanitizeTag(meta.baseline_tag));
         end
         model_tag = IcemodelRegressionTest.smbmodelLabel(meta.smbmodel_filter);
         outfile = fullfile(outdir, ...
            sprintf('regression_report_%s%s_%s.mat', ...
            char(meta.tier), char(model_tag), baseline_tag));
         save(outfile, 'report', 'case_opts', 'meta');
         disp(outfile)
         disp(report(:, {'case_id', 'runoff_final', 'melt_final', ...
            'icemodel_final_m3'}))
      end
   end

   methods (Static, Access = private)
      function ensureModelConfigPaths()
         rootdir = IcemodelRegressionTest.repoRoot();
         addpath(fullfile(rootdir, 'test'));
         test.helpers.configureModelPaths(rootdir);
      end

      function rootdir = repoRoot()
         testdir = fileparts(mfilename('fullpath'));
         rootdir = fileparts(testdir);
      end

      function s = getenvOrDefault(name, default)
         s = getenv(name);
         if isempty(s)
            s = default;
         end
      end

      function meta = reportMeta(tier, smbmodel, baseline_tag, run_name)
         [run_date, run_id, run_name] = test.helpers.resolveRunStamp(run_name);
         meta = struct();
         meta.tier = string(tier);
         meta.smbmodel_filter = string(smbmodel);
         meta.baseline_tag = string(baseline_tag);
         meta.run_date = run_date;
         meta.run_id = run_id;
         meta.run_name = run_name;
         meta.rel_tol_scalar = IcemodelRegressionTest.scalarRelTol();
         meta.abs_tol_scalar = IcemodelRegressionTest.scalarAbsTol();
         meta.rel_tol_runoff_m3 = IcemodelRegressionTest.runoffRelTol();
         meta.abs_tol_runoff_m3 = IcemodelRegressionTest.runoffAbsTol();
         meta.case_builder = "test.helpers.buildFormalCaseOpts";
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
            label = "_" + test.helpers.smbmodelTag(smbmodel);
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
end

function varargout = runModel(opts)
   % Dispatch to the requested model kernel.
   switch opts.smbmodel
      case 'icemodel'
         [varargout{1:nargout}] = icemodel(opts);
      case 'skinmodel'
         [varargout{1:nargout}] = skinmodel(opts);
      otherwise
         error('unsupported smbmodel: %s', opts.smbmodel)
   end
end
