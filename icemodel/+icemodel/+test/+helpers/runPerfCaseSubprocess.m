function runPerfCaseSubprocess(spec_file)
   %RUNPERFCASESUBPROCESS Run one formal perf case in this fresh session.
   %
   %  icemodel.test.helpers.runPerfCaseSubprocess('/path/to/case_spec.mat')
   %
   % This is the child entry point for run_perf_suite(isolation="process").
   % The parent writes one spec MAT per case and launches
   % `matlab -batch` on this function, so every measured case starts from
   % a fresh process: no inherited JIT state, persistents, or heap layout.
   % The spec carries the case row, sampling controls, and the test
   % environment the parent resolved. The result MAT carries the same
   % struct icemodel.test.helpers.runPerfCase returns.
   %
   % Spec fields:
   %  case_struct - one case row from getPerfCaseMatrix, as a struct
   %  n_runs      - fixed sample count for the experiment
   %  config_case - icemodel_config_casename for bootstrapTestEnvironment
   %  data_root   - resolved formal data root ("" selects the default)
   %  result_file - absolute path this function writes
   %
   % See also: icemodel.test.helpers.runPerfCase, run_perf_suite

   spec = load(spec_file);

   % Configure the same formal test environment the parent resolved. The
   % cleanup handle restores the environment on any error exit; the happy
   % path deletes it explicitly after the result is saved.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment( ...
      icemodel_config_casename=spec.config_case, ...
      data_root=spec.data_root);

   % Formal wall-clock timings must never inherit an interactive profiler.
   profile off

   % Build the same experiment shape the in-session runner uses.
   testdir = icemodel.getpath('test');
   suite = testsuite(fullfile(testdir, 'regression', 'IcemodelPerfTest.m'));
   experiment = matlab.perftest.TimeExperiment.withFixedSampleSize( ...
      spec.n_runs, 'NumWarmups', 1);

   % Run the one case and persist the normalized result for the parent.
   c = struct2table(spec.case_struct, 'AsArray', true);
   perf_data = icemodel.test.helpers.runPerfCase(experiment, suite, c);
   save(spec.result_file, '-struct', 'perf_data');

   % Restore the environment now that the result is on disk.
   delete(suite_cleanup)
end
