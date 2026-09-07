function tests = test_run_unit_suite_progress
   %TEST_RUN_UNIT_SUITE_PROGRESS Verify the runner's progress-log contract.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Bootstrap once so run_unit_suite (at the test root) resolves when
   % this file runs standalone; the runner also bootstraps internally.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
end

function teardownOnce(testCase)
   % Release the bootstrap cleanup handle after the last test.
   testCase.TestData.cleanup = [];
end

function setup(testCase)
   % Each test gets its own scratch directory for log files.
   testCase.TestData.tmp = tempname;
   mkdir(testCase.TestData.tmp);
end

function teardown(testCase)
   % Remove the scratch directory.
   if isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
end

function test_unwritable_progress_log_fails_fast(testCase)
   % An unopenable log path must error before any test runs, because a
   % silently missing log defeats the hang-diagnosis purpose.
   bad_log = fullfile(testCase.TestData.tmp, 'missing-parent', 'p.log');
   testCase.verifyError(@() run_unit_suite( ...
      selector="test_verification_setup_helpers", progress_log=bad_log), ...
      'icemodel:test:runUnitSuite:progressLogOpenFailed');
end

function test_progress_log_records_file_boundaries(testCase)
   % A selector run must leave one flushed line before and one after the
   % file, each naming the file, so a hung run names its last file.
   logfile = fullfile(testCase.TestData.tmp, 'progress.log');
   % Capture the runner's own console progress so this check stays quiet.
   evalc(['results = run_unit_suite(' ...
      'selector="test_verification_setup_helpers", ' ...
      sprintf('progress_log="%s");', logfile)]);
   testCase.verifyTrue(all([results.Passed]));
   lines = readlines(logfile);
   lines = lines(strlength(lines) > 0);
   testCase.verifyEqual(numel(lines), 2);
   testCase.verifyTrue(all(contains(lines, ...
      "test_verification_setup_helpers")));
   testCase.verifySubstring(lines(1), "[1/1]");
end
