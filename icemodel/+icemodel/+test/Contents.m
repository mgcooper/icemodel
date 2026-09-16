% +TEST
%
%   Contents file for +TEST and its subfolders.
%
%   +TEST
%   README.md
%
%   +TEST/+FIXTURES
%   icemodel.test.fixtures.cleanupSyntheticWorkspace         - Restore env vars and remove a temp workspace
%   icemodel.test.fixtures.makeReconstructSeries             - One year of hourly synthetic met, smooth channels
%   icemodel.test.fixtures.makeSyntheticColumnState          - Build a resolved synthetic column kernel state
%   icemodel.test.fixtures.makeSyntheticMetFile              - Build a simple synthetic forcing timetable for tests
%   icemodel.test.fixtures.makeSyntheticWorkspace            - Create an isolated icemodel test workspace
%   icemodel.test.fixtures.writeSyntheticMetFile             - Write a synthetic met file for tests
%   icemodel.test.fixtures.writeSyntheticUserdataFile        - Write a synthetic yearly userdata timetable
%
%   +TEST/+HELPERS
%   icemodel.test.helpers.ambientAnchorVerdict               - Decide whether ambient conditions held for a run
%   icemodel.test.helpers.archiveManagedBaseline             - Archive a rolling baseline before overwrite
%   icemodel.test.helpers.artifactFilePath                   - Return the canonical artifact file path
%   icemodel.test.helpers.assertAmbientBaselineAcceptance    - Validate the final baseline anchor
%   icemodel.test.helpers.assertCleanPerfSession             - Refuse an in-session formal run in a dirty session
%   icemodel.test.helpers.assertCommonBaselineRevision       - Require one source revision for a model set
%   icemodel.test.helpers.assertFormalBaselineCandidate      - Reject incomplete state before publication
%   icemodel.test.helpers.assertFormalBaselineForcing        - Verify a baseline's registered forcing identity
%   icemodel.test.helpers.assertFormalBenchmarkCandidate     - Reject invalid component timing evidence
%   icemodel.test.helpers.assertNewReleaseBaselineTarget     - Reject an existing immutable release file
%   icemodel.test.helpers.baselineFilePath                   - Return the canonical baseline file path
%   icemodel.test.helpers.baselineProfilerDir                - Return the profiler-artifact folder for a baseline file
%   icemodel.test.helpers.benchmarkSuiteSignature            - Hash the managed core benchmark suite
%   icemodel.test.helpers.bootstrapTestEnvironment           - Add test paths and install one scoped data config
%   icemodel.test.helpers.buildSyntheticOpts                 - Build resolved OPTS for synthetic unit-test runs
%   icemodel.test.helpers.buildThfValidationCases            - Build focused real-case THF validation cases
%   icemodel.test.helpers.captureBaselineProfile             - Save a build-time profiler report alongside a baseline
%   icemodel.test.helpers.captureExpectedWarning             - Run fcn once, verify its warning, and capture output
%   icemodel.test.helpers.commitBaselineProfilePublication   - Remove a retained sidecar backup
%   icemodel.test.helpers.displayPerfResults                 - Display compact performance results from run_perf_suite
%   icemodel.test.helpers.displayPerfSummary                 - Display a compact perf case summary and benchmark table
%   icemodel.test.helpers.displayRegressionResults           - Display compact results from run_regression_suite
%   icemodel.test.helpers.displayRegressionSummary           - Display compact regression compare summaries
%   icemodel.test.helpers.findCaseRow                        - Find the first baseline/report row matching CASE_ID
%   icemodel.test.helpers.findRunoffReferenceRow             - Resolve runoff reference row for one formal case
%   icemodel.test.helpers.formalBaselinePolicy               - Return forcing and default-root policy for a baseline
%   icemodel.test.helpers.formalPerformanceVerdict           - Evaluate one formal timing comparison row
%   icemodel.test.helpers.formalRegressionMetricEvidence     - Check one saved metric is comparable
%   icemodel.test.helpers.getFormalForcing                   - Return the forcing identity for one formal baseline
%   icemodel.test.helpers.getFormalTestSuiteCases            - Return the canonical formal test-suite cases
%   icemodel.test.helpers.getPerfCaseMatrix                  - Return the canonical formal performance case matrix
%   icemodel.test.helpers.getRegressionCaseMatrix            - Return the canonical formal regression case matrix
%   icemodel.test.helpers.getRunoffSite                      - Map formal station cases to runoff-validation catchments
%   icemodel.test.helpers.loadArtifact                       - Load a saved test artifact
%   icemodel.test.helpers.loadBaseline                       - Load a rolling or release baseline table
%   icemodel.test.helpers.loadProcessedMetForOutputYears     - Load processed met limited to output years
%   icemodel.test.helpers.loadReference                      - Load a test reference table
%   icemodel.test.helpers.loadSavedTable                     - Load a saved table-like object from a MAT file
%   icemodel.test.helpers.machineHostname                    - Return the current machine's trimmed hostname
%   icemodel.test.helpers.makeFormalCaseId                   - Return canonical formal-suite identifier for one model run
%   icemodel.test.helpers.managedBaselineSiblings            - Return the baseline files one build writes
%   icemodel.test.helpers.markTestSessionDirty               - Record that this MATLAB session ran a test suite
%   icemodel.test.helpers.measurePerfCase                    - Measure one case under the selected isolation protocol
%   icemodel.test.helpers.normalizeFormalCaseId              - Normalize legacy formal-suite identifiers
%   icemodel.test.helpers.perfBaselineCompatibility          - Decide whether wall-time comparison is fair
%   icemodel.test.helpers.perfMeasurementPolicy              - Formal timing gate thresholds
%   icemodel.test.helpers.performanceGate                    - Compare one runtime to a two-sided accepted band
%   icemodel.test.helpers.perfSampleValidity                 - Decide whether one case's timing samples are usable
%   icemodel.test.helpers.prepareBaselineBuild               - Resolve shared setup for perf/regression baseline builds
%   icemodel.test.helpers.printFilePath                      - Print a file path truncated to the test/ directory
%   icemodel.test.helpers.publishBaselineBundleSet           - Publish a complete model baseline set
%   icemodel.test.helpers.publishBaselineProfile             - Replace one managed profiler sidecar
%   icemodel.test.helpers.referenceFilePath                  - Return the canonical reference file path
%   icemodel.test.helpers.regressionCaseGates                - Evaluate every regression gate for one formal case
%   icemodel.test.helpers.regressionFailures                 - List the failed cases and gates of one regression run
%   icemodel.test.helpers.removeBaselineProfileStage         - Remove one owned profiler staging directory
%   icemodel.test.helpers.removeReleaseSnapshotArtifacts     - Remove a snapshot MAT file and sidecar
%   icemodel.test.helpers.resolveBaselineBuild               - Resolve baseline type/tag and default output file
%   icemodel.test.helpers.resolveBaselineSelector            - Parse rolling vs release baseline selectors
%   icemodel.test.helpers.resolveBootstrapRelease            - Load or create registered release baselines
%   icemodel.test.helpers.resolveReleaseDataRoots            - Resolve model and fixture roots for a baseline
%   icemodel.test.helpers.resolveRequestedSmbmodels          - Expand one requested formal smbmodel selector
%   icemodel.test.helpers.resolveRunStamp                    - Resolve shared batch run identifiers for test artifacts
%   icemodel.test.helpers.retryInvalidMeasurement            - Measure once; re-measure once if invalid
%   icemodel.test.helpers.rollbackBaselineProfilePublication - Restore the prior profiler sidecar
%   icemodel.test.helpers.runBenchmarkDiagnostics            - Run and compare managed component benchmarks
%   icemodel.test.helpers.runModelCase                       - Resolve, execute, and postprocess one supported model case
%   icemodel.test.helpers.runPerfCase                        - Run one formal performance case and normalize the result
%   icemodel.test.helpers.runPerfCaseSubprocess              - Run one formal perf case in this fresh session
%   icemodel.test.helpers.runSmbModel                        - Dispatch to the requested core SMB model kernel
%   icemodel.test.helpers.runThfValidationCase               - Run one real-case THF validation scenario
%   icemodel.test.helpers.sanitizeTag                        - Replace punctuation and whitespace for filename-safe tags
%   icemodel.test.helpers.setModelOptsForCase                - Build resolved model OPTS for one case
%   icemodel.test.helpers.smbmodelTag                        - Return canonical smbmodel tag for filenames and identifiers
%   icemodel.test.helpers.snapshotBaseline                   - Save a release snapshot from the rolling test baseline
%   icemodel.test.helpers.sourceRevisionGuard                - Capture or verify the current source revision
%   icemodel.test.helpers.summarizeIce1Metrics               - Extract formal regression metrics from output and refs
%   icemodel.test.helpers.testSessionActivity                - Return the suites this MATLAB session has run
%   icemodel.test.helpers.transactionalSnapshotSet           - Create and validate an aggregate release set
%   icemodel.test.helpers.worktreeRevision                   - Return a content-sensitive Git source identity
%
%   +TEST/+VERIFY
%   icemodel.test.verify.verifyEqualNested                   - Recursively compare nested structs, timetables, and arrays
%   icemodel.test.verify.verifyProcessedOutputBounds         - Verify physical bounds of the processed output
%
%   updatecontents.m generated this file on 14 Sep 2026 at 20:14:00.
