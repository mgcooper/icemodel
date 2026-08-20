function tests = test_perf_measurement_controls
   %TEST_PERF_MEASUREMENT_CONTROLS Verify the formal-timing measurement gates.
   %
   % The perf suite's formal verdicts depend on these controls: the sample
   % dispersion validity gate, the session-activity record, the
   % contaminated-session refusal, the invalid-measurement retry policy,
   % the ambient-anchor verdict, and the baseline compatibility check.
   % These tests pin each control at the function level. The two
   % measurement-branch tests below drive measurePerfCase directly: the
   % process branch launches one real subprocess (spec MAT, command
   % construction, child bootstrap, measurement, result MAT), and the
   % session branch measures in this session. They are the slowest tests
   % in this file by design because each runs a real model case.
   %
   % See also: icemodel.test.helpers.perfSampleValidity,
   %  icemodel.test.helpers.assertCleanPerfSession,
   %  icemodel.test.helpers.runPerfCaseSubprocess, run_perf_suite
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Every test may mutate the session-activity record, so save the
   % caller's value and start each test from a clean session.
   testCase.TestData.prior_activity = ...
      getenv('ICEMODEL_TEST_SESSION_ACTIVITY');
   setenv('ICEMODEL_TEST_SESSION_ACTIVITY', '');
end

function teardown(testCase)
   % Restore the caller's session-activity record after each test.
   setenv('ICEMODEL_TEST_SESSION_ACTIVITY', ...
      testCase.TestData.prior_activity);
end

function test_clean_samples_are_valid(testCase)
   % A tight sample set with a valid framework flag supports a verdict.
   [returned, reason, dispersion] = ...
      icemodel.test.helpers.perfSampleValidity([50.1; 50.4; 49.9], true);
   testCase.verifyTrue(returned);
   testCase.verifyEqual(reason, "");
   testCase.verifyLessThan(dispersion, 1.5);
end

function test_dispersion_above_the_gate_is_invalid(testCase)
   % One slow outlier marks interference: max/median = 2 exceeds 1.5.
   [returned, reason, dispersion] = ...
      icemodel.test.helpers.perfSampleValidity([50; 50; 100], true);
   testCase.verifyFalse(returned);
   testCase.verifySubstring(reason, "dispersion");
   testCase.verifyEqual(dispersion, 2.0, 'RelTol', 1e-12);
end

function test_nonfinite_or_nonpositive_samples_are_invalid(testCase)
   % Timings must be finite and positive to mean anything.
   returned = icemodel.test.helpers.perfSampleValidity([50; NaN; 50], true);
   testCase.verifyFalse(returned);
   returned = icemodel.test.helpers.perfSampleValidity([50; -1; 50], true);
   testCase.verifyFalse(returned);
   returned = icemodel.test.helpers.perfSampleValidity( ...
      zeros(0, 1), true);
   testCase.verifyFalse(returned);
end

function test_framework_invalid_flag_overrides_clean_samples(testCase)
   % The perf framework's own Valid flag overrides the sample checks for
   % protocol errors, so clean-looking samples cannot rescue an invalid
   % result.
   [returned, reason] = ...
      icemodel.test.helpers.perfSampleValidity([50; 50; 50], false);
   testCase.verifyFalse(returned);
   testCase.verifySubstring(reason, "framework");
end

function test_custom_dispersion_gate_is_honored(testCase)
   % A caller-supplied gate replaces the 1.5 default in both directions.
   returned = icemodel.test.helpers.perfSampleValidity( ...
      [50; 50; 100], true, 2.5);
   testCase.verifyTrue(returned);
   returned = icemodel.test.helpers.perfSampleValidity( ...
      [50; 50; 55], true, 1.05);
   testCase.verifyFalse(returned);
end

function test_session_activity_records_runners_in_order(testCase)
   % The record starts empty, then appends one label per runner call.
   returned = icemodel.test.helpers.testSessionActivity();
   testCase.verifyEmpty(returned);

   icemodel.test.helpers.markTestSessionDirty("run_unit_suite");
   icemodel.test.helpers.markTestSessionDirty("run_perf_suite");
   returned = icemodel.test.helpers.testSessionActivity();
   expected = ["run_unit_suite", "run_perf_suite"];
   testCase.verifyEqual(returned, expected);
end

function test_session_mode_refuses_a_dirty_session(testCase)
   % An in-session formal run in a dirty session is a refused verdict,
   % not a warning, because a contaminated verdict is worse than none.
   icemodel.test.helpers.markTestSessionDirty("run_unit_suite");
   testCase.verifyError( ...
      @() icemodel.test.helpers.assertCleanPerfSession("session"), ...
      'icemodel:test:perf:contaminatedSession');
end

function test_session_mode_accepts_a_clean_session(testCase)
   % A clean session passes without error.
   icemodel.test.helpers.assertCleanPerfSession("session");
   testCase.verifyTrue(true);
end

function test_process_mode_ignores_session_history(testCase)
   % Process isolation is immune: every case starts a fresh subprocess.
   icemodel.test.helpers.markTestSessionDirty("run_unit_suite");
   icemodel.test.helpers.assertCleanPerfSession("process");
   testCase.verifyTrue(true);
end

function test_process_branch_measures_through_a_fresh_subprocess(testCase)
   % The process branch of measurePerfCase owns spec creation, command
   % construction, the matlab -batch launch, and result loading; the
   % child runs runPerfCaseSubprocess (spec deserialization, environment
   % bootstrap, suite discovery, one real measurement, result
   % serialization). One call with one sample covers that whole path.
   [experiment, suite, c, scratch] = measurementFixture(testCase);

   returned = icemodel.test.helpers.measurePerfCase(experiment, suite, ...
      c, "process", "verification", "", 1, scratch, 1);
   verifyPerfDataSchema(testCase, returned, 1);

   % The spec and result MAT files must stay auditable in the artifact
   % folder under the attempt-numbered names.
   testCase.verifyEqual(exist(fullfile(scratch, sprintf( ...
      'isolated_%s_attempt1_spec.mat', c.case_id)), 'file'), 2);
   testCase.verifyEqual(exist(fullfile(scratch, sprintf( ...
      'isolated_%s_attempt1_result.mat', c.case_id)), 'file'), 2);
end

function test_session_branch_measures_in_this_session(testCase)
   % The session branch measures in this MATLAB session after
   % clear-functions hygiene and returns the same result shape as the
   % process branch.
   [experiment, suite, c, scratch] = measurementFixture(testCase);

   returned = icemodel.test.helpers.measurePerfCase(experiment, suite, ...
      c, "session", "verification", "", 1, scratch, 1);
   verifyPerfDataSchema(testCase, returned, 1);
end

function [experiment, suite, c, scratch] = measurementFixture(testCase)
   % Shared setup for the two measurement-branch tests: the cheapest
   % formal case, a one-sample experiment, and one fixture-owned scratch
   % directory in the system temporary folder (teardown registers before
   % the first write, and the folder never touches the worktree).
   cases = icemodel.test.helpers.getPerfCaseMatrix( ...
      tier="smoke", smbmodel="icemodel", solver=2);
   testCase.assertGreaterThanOrEqual(height(cases), 1);
   c = cases(1, :);
   testdir = icemodel.getpath('test');
   suite = testsuite(fullfile(testdir, 'regression', ...
      'IcemodelPerfTest.m'));
   experiment = matlab.perftest.TimeExperiment.withFixedSampleSize( ...
      1, 'NumWarmups', 1);
   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   scratch = fixture.Folder;
end

function verifyPerfDataSchema(testCase, returned, n_runs)
   % Both measurement branches must return the runPerfCase fields with
   % finite positive samples and a usable validity flag.
   expected = ["samples"; "activity"; "sample_times"; ...
      "activity_times"; "valid"; "n_warmups"];
   testCase.verifyTrue(all(ismember(expected, fieldnames(returned))));
   testCase.verifyNumElements(returned.sample_times, n_runs);
   testCase.verifyTrue(all(isfinite(returned.sample_times) ...
      & returned.sample_times > 0));
   testCase.verifyClass(returned.valid, 'logical');
end

function test_retry_remeasures_an_invalid_sample_set(testCase)
   % One dispersed sample set earns one automatic re-measure; a clean
   % second set supports a verdict.
   sets = {[50; 50; 100], [50.1; 50.4; 49.9]};
   measure = @(attempt) struct( ...
      'sample_times', sets{attempt}, 'valid', true);
   [returned, valid, reason, ~, attempts] = ...
      icemodel.test.helpers.retryInvalidMeasurement(measure);
   testCase.verifyTrue(valid);
   testCase.verifyEqual(reason, "");
   testCase.verifyEqual(attempts, 2);
   testCase.verifyEqual(returned.sample_times, sets{2});
end

function test_retry_returns_the_second_invalid_set_as_invalid(testCase)
   % Two invalid sets fail the case as "measurement invalid" rather than
   % pass or fail on contaminated numbers.
   measure = @(~) struct('sample_times', [50; 50; 100], 'valid', true);
   [~, valid, reason, dispersion, attempts] = ...
      icemodel.test.helpers.retryInvalidMeasurement(measure);
   testCase.verifyFalse(valid);
   testCase.verifySubstring(reason, "dispersion");
   testCase.verifyEqual(dispersion, 2.0, 'RelTol', 1e-12);
   testCase.verifyEqual(attempts, 2);
end

function test_retry_absorbs_one_transient_launch_failure(testCase)
   % The first subprocess launch failure earns one retry; the second
   % attempt's clean measurement supports a verdict.
   [~, valid, ~, ~, attempts] = ...
      icemodel.test.helpers.retryInvalidMeasurement( ...
      @(attempt) throwingMeasure(attempt, 1, ...
      'icemodel:test:perf:subprocessFailed'));
   testCase.verifyTrue(valid);
   testCase.verifyEqual(attempts, 2);
end

function test_retry_rethrows_a_second_launch_failure(testCase)
   % A second launch failure is real and must stop the run.
   testCase.verifyError( ...
      @() icemodel.test.helpers.retryInvalidMeasurement( ...
      @(attempt) throwingMeasure(attempt, inf, ...
      'icemodel:test:perf:subprocessFailed')), ...
      'icemodel:test:perf:subprocessFailed');
end

function test_retry_rethrows_other_errors_immediately(testCase)
   % Only the transient launch identifier earns a retry; any other error
   % propagates on the first attempt.
   testCase.verifyError( ...
      @() icemodel.test.helpers.retryInvalidMeasurement( ...
      @(attempt) throwingMeasure(attempt, inf, ...
      'icemodel:test:someOtherError')), 'icemodel:test:someOtherError');
end

function perf_data = throwingMeasure(attempt, n_failures, identifier)
   % Stub measurement that throws IDENTIFIER on the first N_FAILURES
   % attempts and returns a clean sample set after them. The output is
   % assigned before the throw so every path sets the declared return.
   perf_data = struct('sample_times', [50; 50; 50], 'valid', true);
   if attempt <= n_failures
      error(identifier, 'stub failure on attempt %d', attempt)
   end
end

function test_anchor_certifies_a_stable_run(testCase)
   % An anchor within tolerance of its own first measurement certifies
   % the run's timings.
   [returned, ratio] = icemodel.test.helpers.ambientAnchorVerdict( ...
      50.0, [51; 52; 51], true, 0.15);
   testCase.verifyTrue(returned);
   testCase.verifyEqual(ratio, 51 / 50, 'RelTol', 1e-12);
end

function test_anchor_drift_invalidates_the_run(testCase)
   % A drifted anchor marks every verdict ambient-invalid.
   returned = icemodel.test.helpers.ambientAnchorVerdict( ...
      50.0, [70; 70; 70], true, 0.15);
   testCase.verifyFalse(returned);
end

function test_an_invalid_anchor_cannot_certify_stability(testCase)
   % An anchor whose own sample set failed the validity gate cannot
   % certify ambient stability, even when its median lands in tolerance.
   returned = icemodel.test.helpers.ambientAnchorVerdict( ...
      50.0, [50; 50; 50], false, 0.15);
   testCase.verifyFalse(returned);
end

function test_a_nonfinite_anchor_ratio_invalidates_the_run(testCase)
   % An empty or all-NaN anchor sample set produces a NaN ratio, which
   % must read as unstable rather than pass through the tolerance test.
   [returned, ratio] = icemodel.test.helpers.ambientAnchorVerdict( ...
      50.0, nan(3, 1), true, 0.15);
   testCase.verifyFalse(returned);
   testCase.verifyFalse(isfinite(ratio));
end

function test_matching_protocol_and_environment_are_compatible(testCase)
   % A baseline from this environment and protocol supports a timing gate.
   meta = struct('matlab_version', string(version), ...
      'host', string(computer), 'isolation', "process");
   [returned, reason] = ...
      icemodel.test.helpers.perfBaselineCompatibility(meta, "process");
   testCase.verifyTrue(returned);
   testCase.verifyEqual(reason, "");
end

function test_isolation_mismatch_is_incompatible(testCase)
   % Cross-protocol timings are not comparable, whatever the environment.
   meta = struct('matlab_version', string(version), ...
      'host', string(computer), 'isolation', "session");
   [returned, reason] = ...
      icemodel.test.helpers.perfBaselineCompatibility(meta, "process");
   testCase.verifyFalse(returned);
   testCase.verifySubstring(reason, "isolation");
end

function test_a_baseline_without_the_isolation_field_is_session(testCase)
   % A baseline saved before the isolation field exists counts as
   % session-protocol: incompatible with process, compatible with session.
   meta = struct('matlab_version', string(version), ...
      'host', string(computer));
   returned = ...
      icemodel.test.helpers.perfBaselineCompatibility(meta, "process");
   testCase.verifyFalse(returned);
   returned = ...
      icemodel.test.helpers.perfBaselineCompatibility(meta, "session");
   testCase.verifyTrue(returned);
end

function test_missing_environment_metadata_is_incompatible(testCase)
   % A baseline without environment metadata cannot support a timing gate.
   [returned, reason] = ...
      icemodel.test.helpers.perfBaselineCompatibility(struct(), "process");
   testCase.verifyFalse(returned);
   testCase.verifySubstring(reason, "metadata");
end
