function results = run_perf_suite(kwargs)
   %RUN_PERF_SUITE Run formal icemodel performance benchmark suite.
   %
   %  results = run_perf_suite()
   %  results = run_perf_suite(tier="smoke")
   %  results = run_perf_suite(isolation="process")
   %  results = run_perf_suite(tier="smoke", smbmodel="skinmodel")
   %  results = run_perf_suite(tier="smoke", smbmodel="icemodel", solver=2)
   %  results = run_perf_suite(tier="smoke", smbmodel="icemodel", solver=[2 3])
   %  results = run_perf_suite(tier="smoke", smbmodel="icemodel", solver=[1 3])
   %  results = run_perf_suite(simyear=2017, smoke_sites="kanm", ...
   %     full_sites=["kanm"; "kanl"])
   %  results = run_perf_suite(tier="full", baseline="v1.1")
   %  results = run_perf_suite(data_root="/path/to/test/data")
   %  results = run_perf_suite(fixture_root="/path/to/provisioned/data")
   %  results = run_perf_suite(artifact_root="/path/to/artifacts")
   %  results = run_perf_suite(tier="smoke", measure_anchor=false)
   %
   % Use this for normal performance comparisons against an existing rolling or
   % release baseline.
   %
   % This function does not update baselines; it only runs the formal cases,
   % compares runtime to the requested baseline, and writes one artifact under
   % test/artifacts/<run_name>/ by default. ARTIFACT_ROOT selects a different
   % root. The function also writes a Quarto HTML summary for the combined run.
   %
   % By default it also runs the managed core benchmark suite and saves those
   % diagnostic timings alongside the formal perf artifact.
   %
   % Formal perf cases run the same way regression does: one leading
   % spinup year plus one retained output year when the case matrix
   % carries only SIMYEAR.
   %
   % The optional solver filter accepts any subset of
   % icemodel.namelists.solver().
   % DATA_ROOT overrides the default test case for isolated fixture comparisons.
   %
   % ISOLATION selects the measurement protocol. "process" (default) runs
   % every case in a fresh `matlab -batch` subprocess and is the required
   % mode for formal accept/reject verdicts on refactors. The opt-in
   % "session" mode times every case in this MATLAB session for quick
   % diagnostics; it refuses to start when this session already ran
   % another suite, because inherited JIT state and persistents
   % contaminate formal timings. Both modes randomize the
   % case order (the seed is recorded in the artifact) and gate each
   % case's samples on a dispersion check: an invalid sample set is
   % re-measured once, then fails as "measurement invalid" instead of
   % producing an incorrect verdict.
   %
   % SMOKE_SITES and FULL_SITES are advanced overrides for the site lists
   % used by each formal tier. Most callers should leave them at the
   % defaults and select only TIER.
   %
   % ARTIFACT_ROOT writes run folders outside test/artifacts. Use it for
   % temporary or externally managed comparison artifacts.
   %
   % MEASURE_ANCHOR controls the ambient anchor, one more subprocess that
   % re-measures the first executed case at the end of the run. It is the
   % only check that sees a load or thermal shift across cases, so it stays
   % on by default. Pass false for a quick read-only diagnostic comparison;
   % the artifact then records anchor_measured=false and fails the
   % measurement quality conditions, which a comparison never needs.
   %
   % results.passed is the comparison verdict: every case inside its
   % tolerance band. results.quality is the measurement quality verdict of
   % perfMeasurementQuality. The two are separate: a run whose timing moved
   % outside the band still writes its artifact and still reports quality,
   % and a run measured on a loaded machine can pass its comparison. Every
   % artifact records a machine-state attestation (foreign MATLAB processes,
   % one-minute load average, AC power) sampled at run start, after each
   % case, and at run end.
   %
   % CLI entrypoint:
   %  matlab -batch "run('/ABS/PATH/icemodel/test/run_perf_suite.m')"

   arguments (Input)

      kwargs.tier (1, :) string ...
         {icemodel.validators.mustBeTestTierName(kwargs.tier)} ...
         = "smoke"

      kwargs.smbmodel (1, :) string ...
         {icemodel.validators.mustBeTestSmbmodelSelector(kwargs.smbmodel)} ...
         = "all"

      kwargs.solver {icemodel.validators.mustBeSolverFilter(kwargs.solver)} ...
         = []

      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} ...
         = 2016

      kwargs.smoke_sites string ...
         = "kanm"

      kwargs.full_sites string ...
         = ["kanm"; "kanl"]

      kwargs.n_runs (1, 1) double {mustBeInteger, mustBePositive} ...
         = 3

      kwargs.tol_perf (1, 1) double {mustBePositive} ...
         = icemodel.test.helpers.perfMeasurementPolicy().tol_perf

      kwargs.include_benchmarks (1, 1) logical ...
         = true

      kwargs.benchmark_sampling_profile (1, :) string ...
         {icemodel.validators.mustBeBenchmarkSamplingProfileName( ...
         kwargs.benchmark_sampling_profile)} ...
         = "default"

      kwargs.baseline (1, :) string ...
         = "rolling"

      kwargs.run_name string ...
         = string.empty()

      kwargs.artifact_root (1, 1) string ...
         = ""

      kwargs.build_report (1, 1) logical ...
         = true

      kwargs.data_root (1, 1) string ...
         = ""

      % FIXTURE_ROOT selects where a release's required fixture capabilities
      % are verified. A release that runs the model from its own provisioned
      % data (v1.2 onward) rejects a DATA_ROOT that names a different tree
      % with icemodel:test:releaseDataRootMismatch.
      kwargs.fixture_root (1, 1) string ...
         = ""

      % Process isolation is the default because it is the only protocol
      % that supports a formal accept/reject verdict. Session mode is
      % opt-in for quick diagnostics.
      kwargs.isolation (1, 1) string ...
         {mustBeMember(kwargs.isolation, ["session", "process"])} ...
         = "process"

      % The anchor is on by default because it is the only check that sees
      % a shift across cases. build_perf_baseline has no opt-out.
      kwargs.measure_anchor (1, 1) logical ...
         = true
   end

   % Deal out arguments.
   [tier, smbmodel, solver, simyear, smoke_sites, full_sites, n_runs, ...
      tol_perf, include_benchmarks, benchmark_sampling_profile, ...
      baseline_selector, run_name, artifact_root, build_report, isolation, ...
      measure_anchor] = deal( ...
      kwargs.tier, kwargs.smbmodel, kwargs.solver, kwargs.simyear, ...
      reshape(kwargs.smoke_sites, [], 1), reshape(kwargs.full_sites, [], 1), ...
      kwargs.n_runs, kwargs.tol_perf, kwargs.include_benchmarks, ...
      kwargs.benchmark_sampling_profile, ...
      kwargs.baseline, kwargs.run_name, kwargs.artifact_root, ...
      kwargs.build_report, ...
      kwargs.isolation, kwargs.measure_anchor);

   % Resolve full path to the test/ dir.
   testdir = icemodel.getpath('test');

   % Resolve the requested baseline and shared batch run identifier.
   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector(baseline_selector);
   baseline_policy = ...
      icemodel.test.helpers.formalBaselinePolicy(baseline_selector);

   [data_root, fixture_root] = ...
      icemodel.test.helpers.resolveReleaseDataRoots( ...
      baseline_policy, kwargs.data_root, kwargs.fixture_root);

   [run_date, run_id, run_name] = ...
      icemodel.test.helpers.resolveRunStamp(run_name);

   % Bootstrap the source/test trees once, then configure the formal paths.
   % Keep the cleanup handle in scope so the caller's config is restored
   % when this entrypoint returns.
   [~, input_path, output_path, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment( ...
      icemodel_config_casename=baseline_policy.config_case, ...
      data_root=data_root);

   % Carry the configured root through the unittest class's nested setup. A
   % caller-supplied root keeps precedence over the resolved verification root.
   if isblanktext(data_root)
      data_root = string(fileparts(input_path));
   end

   % Verify release assets separately from the formal model input tree.
   if ~isempty(baseline_policy.required_fixture_capabilities)
      icemodel.verification.setup.fetchFixtures( ...
         baseline_policy.baseline_tag, ...
         capabilities=baseline_policy.required_fixture_capabilities, ...
         root=fixture_root, download=false);
   end
   % Hold the env restore for the rest of this function; cleared at the end.
   data_root_cleanup = configurePerfDataRootEnv(data_root);

   % Refuse a contaminated in-session formal run now, immediately before
   % measurement and after the bootstrap and capability seams, so the
   % no-run probes in test_baseline_contracts stop at their intended
   % earlier checks. Record this run in the session activity for later
   % runners; the activity at start rides the artifact metadata.
   icemodel.test.helpers.assertCleanPerfSession(isolation);
   session_activity = icemodel.test.helpers.testSessionActivity();
   icemodel.test.helpers.markTestSessionDirty("run_perf_suite");

   % Formal wall-clock timings must never inherit an interactive profiler.
   profile off

   % Expand the requested formal model selector once at the entrypoint.
   models = icemodel.test.helpers.resolveRequestedSmbmodels(smbmodel);

   % Build the MATLAB perf experiment once, then reuse it for each
   % single-model perf workflow below.
   suite = testsuite(fullfile(testdir, 'regression', ...
      'IcemodelPerfTest.m'));
   experiment = matlab.perftest.TimeExperiment.withFixedSampleSize( ...
      n_runs, 'NumWarmups', 1);

   % Capture the source identity before the first measurement. The value is
   % git describe --always --dirty, so a tracked edit made before the run
   % marks every artifact of the run with the -dirty suffix.
   revision_at_start = icemodel.test.helpers.worktreeRevision();

   % Run the single-model workflow for each requested model and merge the
   % saved compare summaries into one returned struct.
   per_model = arrayfun(@(mdl) runSingleModelPerfSuite( ...
      input_path, output_path, testdir, experiment, suite, tier, ...
      mdl, solver, simyear, smoke_sites, full_sites, baseline_type, ...
      baseline_tag, run_date, run_id, run_name, n_runs, tol_perf, ...
      include_benchmarks, benchmark_sampling_profile, isolation, ...
      baseline_policy.config_case, data_root, session_activity, ...
      revision_at_start, artifact_root, measure_anchor), ...
      models, 'UniformOutput', false);

   % Combine results into a common struct.
   results = combinePerfResults(per_model);

   % Display the results.
   icemodel.test.helpers.displayPerfResults(results)

   % Render the saved comparison with
   % icemodel.verification.report.buildTestSuiteReport unless the caller
   % requested an artifact-only run.
   results.report_file = "";
   if build_report
      results.report_file = ...
         icemodel.verification.report.buildTestSuiteReport("performance", results);
   end

   % Restore the data-root environment now that every timing has been taken.
   % An early error still restores it, because the object dies with the scope.
   delete(data_root_cleanup)

   % Restore the caller's config last, after every path-dependent cleanup has
   % already run against the configured test environment.
   delete(suite_cleanup)
end

function results = runSingleModelPerfSuite(input_path, output_path, ...
      thisdir, experiment, suite, tier, smbmodel, solver, simyear, ...
      smoke_sites, full_sites, baseline_type, baseline_tag, run_date, ...
      run_id, run_name, n_runs, tol_perf, include_benchmarks, ...
      benchmark_sampling_profile, isolation, config_case, data_root, ...
      session_activity, revision_at_start, artifact_root, measure_anchor)
   %RUNSINGLEMODELPERFSUITE Run the formal perf workflow for one smbmodel.

   % Build the deterministic case list and load the matching managed baseline.
   cases = icemodel.test.helpers.getPerfCaseMatrix( ...
      tier=tier, smbmodel=smbmodel, solver=solver, simyear=simyear, ...
      baseline=baseline_tag, smoke_sites=smoke_sites, full_sites=full_sites);

   if isempty(cases)
      error('no performance cases matched tier=%s smbmodel=%s', tier, smbmodel)
   end

   % Perf baselines and benchmark baselines are keyed by one comparison
   % year, so check that before loading any baseline.
   benchmark_year = unique(cases.simyear);
   assert(isscalar(benchmark_year), ...
      'formal perf suite expects exactly one benchmark year')

   % Load the accepted baseline that matches this concrete formal model.
   [baseline, baseline_meta] = icemodel.test.helpers.loadBaseline("perf", ...
      smbmodel=smbmodel, baseline_tag=baseline_tag, simyear=benchmark_year);
   icemodel.test.helpers.assertFormalBaselineForcing( ...
      baseline, baseline_tag);
   [baseline_compatible, compare_reason] = ...
      icemodel.test.helpers.perfBaselineCompatibility( ...
      baseline_meta, isolation);

   % Accumulate the measured sample/activity rows for the saved artifact.
   [sample_rows, activity_rows, case_rows, case_opts] = deal(struct([]));
   [r_sample, r_activity, r_case] = deal(0);

   % Randomize the case order so no case always inherits the same
   % predecessor's session state. The seed rides the artifact so the
   % exact order is reproducible.
   case_order_seed = randi(2^31 - 2);
   rng_prior = rng(case_order_seed, 'twister');
   case_order = randperm(height(cases));
   rng(rng_prior);

   % Resolve the artifact folder now: subprocess spec/result files live
   % beside the saved comparison so the raw measurements stay auditable.
   artifact_dir = fileparts(icemodel.test.helpers.artifactFilePath( ...
      "perf", tier=tier, smbmodel=smbmodel, solver=solver, ...
      baseline_type=baseline_type, baseline_tag=baseline_tag, ...
      run_name=run_name, artifact_root=artifact_root));
   if isolation == "process" && exist(artifact_dir, 'dir') ~= 7
      mkdir(artifact_dir);
   end

   % Sample the machine state at run start, after each case, and after the
   % anchor, so the artifact can say whether the machine was quiet.
   n_state_samples = height(cases) + 1 + measure_anchor;
   state_samples = repmat( ...
      icemodel.test.helpers.sampleMachineState(), n_state_samples, 1);
   n_sampled = 1;

   % Run the per-case performance experiment and compare to baseline.
   for iorder = 1:height(cases)
      icase = case_order(iorder);

      % Measure the case under the selected isolation protocol, with one
      % automatic re-measure when the sample set fails the validity gate.
      c = cases(icase, :);
      [perf_data, valid, measurement_reason, dispersion, ...
         n_measure_attempts] = ...
         icemodel.test.helpers.retryInvalidMeasurement(@(attempt) ...
         icemodel.test.helpers.measurePerfCase(experiment, suite, c, ...
         isolation, config_case, data_root, n_runs, artifact_dir, ...
         attempt));

      % Extract results for this case.
      samples = perf_data.samples;
      activity = perf_data.activity;
      sample_times = perf_data.sample_times;
      activity_times = perf_data.activity_times;

      % Compare the measured runtime to the accepted baseline row.
      case_compare_reason = compare_reason;
      if isempty(baseline) ...
            || ~ismember('case_id', baseline.Properties.VariableNames)
         bid = [];
      else
         bid = find(string(baseline.case_id) == string(c.case_id));
      end
      [passed_perf, ref_wall, floor_wall, gate_wall, ...
         case_compare_reason] = ...
         icemodel.test.helpers.formalPerformanceVerdict( ...
         valid, median(sample_times, 'omitnan'), baseline, bid, ...
         baseline_compatible, tol_perf, case_compare_reason);

      % Save the per-sample timings for later inspection.
      for i = 1:height(samples)
         r_sample = r_sample + 1;
         sample_rows(r_sample).case_id = string(c.case_id);
         sample_rows(r_sample).tier = string(c.tier);
         sample_rows(r_sample).smbmodel = string(c.smbmodel);
         sample_rows(r_sample).sitename = string(c.sitename);
         sample_rows(r_sample).forcings = string(c.forcings);
         sample_rows(r_sample).simyear = c.simyear;
         sample_rows(r_sample).solver = c.solver;
         sample_rows(r_sample).baseline_type = baseline_type;
         sample_rows(r_sample).baseline_tag = baseline_tag;
         sample_rows(r_sample).sample = i;
         sample_rows(r_sample).wall_s = sample_times(i);
      end

      % Save the per-activity timings that MATLAB reports for this case.
      for i = 1:height(activity)
         r_activity = r_activity + 1;
         activity_rows(r_activity).case_id = string(c.case_id);
         activity_rows(r_activity).tier = string(c.tier);
         activity_rows(r_activity).smbmodel = string(c.smbmodel);
         activity_rows(r_activity).sitename = string(c.sitename);
         activity_rows(r_activity).forcings = string(c.forcings);
         activity_rows(r_activity).simyear = c.simyear;
         activity_rows(r_activity).solver = c.solver;
         if ismember('Objective', activity.Properties.VariableNames)
            activity_rows(r_activity).objective = string(activity.Objective(i));
         else
            activity_rows(r_activity).objective = "";
         end
         activity_rows(r_activity).wall_s = activity_times(i);
      end

      % Save the compact per-case summary and the resolved opts struct.
      r_case = r_case + 1;
      case_rows(r_case).case_id = string(c.case_id);
      case_rows(r_case).tier = string(c.tier);
      case_rows(r_case).smbmodel = string(c.smbmodel);
      case_rows(r_case).sitename = string(c.sitename);
      case_rows(r_case).forcings = string(c.forcings);
      case_rows(r_case).simyear = c.simyear;
      case_rows(r_case).solver = c.solver;
      case_rows(r_case).baseline_type = baseline_type;
      case_rows(r_case).baseline_tag = baseline_tag;
      case_rows(r_case).n_runs = n_runs;
      case_rows(r_case).n_warmups = perf_data.n_warmups;
      case_rows(r_case).tol_perf = tol_perf;
      case_rows(r_case).median_wall_s = median(sample_times, 'omitnan');
      case_rows(r_case).mean_wall_s = mean(sample_times, 'omitnan');
      case_rows(r_case).min_wall_s = min(sample_times, [], 'omitnan');
      case_rows(r_case).max_wall_s = max(sample_times, [], 'omitnan');
      case_rows(r_case).ref_wall_s = ref_wall;
      case_rows(r_case).floor_wall_s = floor_wall;
      case_rows(r_case).gate_wall_s = gate_wall;
      case_rows(r_case).baseline_compatible = baseline_compatible;
      case_rows(r_case).compare_reason = case_compare_reason;
      case_rows(r_case).valid = valid;
      case_rows(r_case).isolation = isolation;
      case_rows(r_case).measurement_reason = measurement_reason;
      case_rows(r_case).dispersion = dispersion;
      case_rows(r_case).n_measure_attempts = n_measure_attempts;
      case_rows(r_case).passed_perf = passed_perf;
      % The machine state right after this case, so a reader can tell
      % which case ran under load.
      n_sampled = n_sampled + 1;
      state_samples(n_sampled) = icemodel.test.helpers.sampleMachineState();
      case_rows(r_case).load_average_1min = ...
         state_samples(n_sampled).load_average_1min;
      case_rows(r_case).foreign_matlab_processes = ...
         state_samples(n_sampled).foreign_matlab_processes;
      case_rows(r_case).last_updated_utc = datetime('now', 'TimeZone', 'UTC');

      case_opts(r_case).case_id = string(c.case_id);
      case_opts(r_case).case = table2struct(c);
      case_opts(r_case).opts = icemodel.test.helpers.setModelOptsForCase(c);
   end

   % Ambient anchor: re-measure the first executed case once at the end
   % of the run and compare against its own first measurement. The
   % per-case dispersion gate cannot see load or scheduling shifts that
   % are steady WITHIN each case but different ACROSS cases (measured on
   % this host as a ~35 percent case-level swing under a constant
   % background load). A drifted anchor records ambient_stable=false once
   % for the run; every case keeps its own comparison verdict and numbers,
   % because a reader needs them to diagnose the run, and the run fails
   % the measurement quality conditions instead.
   % The anchor tolerance comes from icemodel.test.helpers.perfMeasurementPolicy.
   anchor_tol = icemodel.test.helpers.perfMeasurementPolicy().anchor_tol;
   anchor_ratio = nan;
   ambient_stable = false;
   if height(cases) > 0 && measure_anchor
      c_anchor = cases(case_order(1), :);
      first_median = case_rows([case_rows.case_id] == ...
         string(c_anchor.case_id)).median_wall_s;

      % The anchor sample set passes the same validity gate and single
      % re-measure as a formal case: an invalid anchor cannot certify
      % ambient stability. The attempt offset keeps its subprocess spec
      % and result files clear of the case's own attempt files.
      [anchor_data, anchor_valid] = ...
         icemodel.test.helpers.retryInvalidMeasurement(@(attempt) ...
         icemodel.test.helpers.measurePerfCase(experiment, suite, ...
         c_anchor, isolation, config_case, data_root, n_runs, ...
         artifact_dir, 98 + attempt));
      [ambient_stable, anchor_ratio] = ...
         icemodel.test.helpers.ambientAnchorVerdict(first_median, ...
         anchor_data.sample_times, anchor_valid, anchor_tol);
      n_sampled = n_sampled + 1;
      state_samples(n_sampled) = icemodel.test.helpers.sampleMachineState();
   end
   attestation = icemodel.test.helpers.summarizeMachineState( ...
      state_samples(1:n_sampled));

   % The failed list is the comparison verdict only. The string conversion
   % keeps an all-pass run's empty list a string array, because the empty
   % struct-field concatenation is numeric.
   failed_cases = reshape(string( ...
      [case_rows(~[case_rows.passed_perf]).case_id]), [], 1);

   % Build the saved artifact tables for this concrete formal model. The
   % rows accumulated in randomized execution order; sort the summary by
   % case id so displays and diffs stay stable across runs.
   sample_detail = struct2table(sample_rows);
   activity_detail = struct2table(activity_rows);
   case_summary = struct2table(case_rows);
   case_summary = sortrows(case_summary, 'case_id');

   % Record the compare metadata for this model-specific artifact.
   meta = struct();
   meta.tier = tier;
   meta.smbmodel_filter = smbmodel;
   meta.solver_filter = solver;
   meta.baseline_type = baseline_type;
   meta.baseline_tag = baseline_tag;
   meta.run_date = run_date;
   meta.run_id = run_id;
   meta.run_name = run_name;
   meta.simyear = benchmark_year;
   meta.smoke_sites = smoke_sites;
   meta.full_sites = full_sites;
   meta.baseline_file = icemodel.test.helpers.baselineFilePath("perf", ...
      smbmodel=smbmodel, baseline_type=baseline_type, ...
      baseline_tag=baseline_tag, simyear=benchmark_year);
   meta.case_builder = "icemodel.test.helpers.setModelOptsForCase";
   meta.opts_source = "icemodel.setopts defaults";
   meta.spinup_policy = ...
      "formal perf runs include the canonical leading spinup year";
   meta.reset_fields = "solver";
   meta.n_runs = n_runs;
   meta.n_warmups = 1;
   meta.tol_perf = tol_perf;
   meta.include_benchmarks = include_benchmarks;
   meta.benchmark_sampling_profile = benchmark_sampling_profile;
   meta.isolation = isolation;
   meta.case_order_seed = case_order_seed;
   meta.case_order = case_order;
   meta.session_activity_at_start = session_activity;
   meta.anchor_measured = measure_anchor && height(cases) > 0;
   meta.anchor_ratio = anchor_ratio;
   meta.ambient_stable = ambient_stable;
   % A comparison run never accepts a drift override; the field is present
   % so the quality conditions read the same fields as a baseline file.
   meta.ambient_drift_accepted = false;
   meta.attestation = attestation;
   meta.experiment = "matlab.perftest.TimeExperiment.withFixedSampleSize";
   meta.timing_scope = "IcemodelPerfTest.testCoreRuntime (runSmbModel only)";
   meta.timing_notes = sprintf([ ...
      'median_wall_s is the median of %d timed samples (wall-clock seconds). ' ...
      '%d warmup run(s) are executed first and excluded from all summary ' ...
      'statistics. Setup, teardown, and reporting overhead are outside the ' ...
      'measured region. The timed model call includes spinup and output years.'], ...
      n_runs, 1);
   meta.input_path = string(input_path);
   meta.output_path = string(output_path);
   meta.suite_file = string(fullfile(thisdir, 'regression', ...
      'IcemodelPerfTest.m'));
   meta.matlab_version = string(version);
   meta.host = string(computer);
   meta.hostname = icemodel.test.helpers.machineHostname();
   meta.artifact_root = artifact_root;

   % The A/A gate certifies two runs of the same code, so record the
   % identity captured before the first measurement. A tracked edit made
   % during the run is not detected; a single developer knows about such an
   % edit and answers it with a rerun.
   meta.git_revision = revision_at_start;

   % Runs that measured different input trees are not comparable; the
   % A/A gate compares this resolved root.
   meta.data_root = string(data_root);

   meta.baseline_meta = baseline_meta;
   meta.baseline_compatible = baseline_compatible;
   meta.compare_reason = compare_reason;
   meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');

   % Run the supporting component benchmarks for this saved compare artifact.
   benchmark = icemodel.test.helpers.runBenchmarkDiagnostics( ...
      benchmark_year, baseline_tag, smbmodel, ...
      include_benchmarks, benchmark_sampling_profile);

   % Save the artifacts file.
   artifact_file = saveArtifacts(sample_detail, activity_detail, ...
      case_summary, case_opts, benchmark, meta);

   % Assemble the results.
   results = struct();
   results.case_summary = case_summary;
   results.sample_detail = sample_detail;
   results.activity_detail = activity_detail;
   results.case_opts = case_opts;
   results.benchmark = benchmark;
   results.meta = meta;
   results.artifact_file = string(artifact_file);
   results.failed_cases = failed_cases;
   results.passed = isempty(failed_cases);
   % The quality verdict is separate from the comparison verdict; a release
   % gate blocks on quality and reads the comparison as information.
   results.quality = ...
      icemodel.test.helpers.perfMeasurementQuality(case_summary, meta);
   results.quality.smbmodels = string(smbmodel);
end

function results = combinePerfResults(per_model)
   %COMBINEPERFRESULTS Merge one-or-more single-model perf run results.

   if isscalar(per_model)
      results = per_model{1};
      return
   end

   % Extract each returned field once, then concatenate the per-model pieces.
   % The aggregate struct carries both table-like fields and scalar pass/fail
   % metadata, which a single vertcat cannot collapse.
   case_summary = cellfun(@(s) s.case_summary, per_model, ...
      'UniformOutput', false);
   sample_detail = cellfun(@(s) s.sample_detail, per_model, ...
      'UniformOutput', false);
   activity_detail = cellfun(@(s) s.activity_detail, per_model, ...
      'UniformOutput', false);
   case_opts = cellfun(@(s) s.case_opts(:), per_model, ...
      'UniformOutput', false);
   benchmark = cellfun(@(s) s.benchmark, per_model, ...
      'UniformOutput', false);
   meta = cellfun(@(s) s.meta, per_model, ...
      'UniformOutput', false);
   artifact_file = cellfun(@(s) string(s.artifact_file(:)), per_model, ...
      'UniformOutput', false);
   failed_cases = cellfun(@(s) string(s.failed_cases(:)), per_model, ...
      'UniformOutput', false);
   pass_flags = cellfun(@(s) s.passed, per_model);
   quality_flags = cellfun(@(s) s.quality.passed, per_model);
   quality_conditions = cellfun(@(s) s.quality.conditions, per_model);
   quality_reasons = cellfun(@(s) reshape(s.quality.reasons, [], 1), ...
      per_model, 'UniformOutput', false);
   quality_models = cellfun(@(s) s.quality.smbmodels, per_model);

   results = struct();
   results.case_summary = vertcat(case_summary{:});
   results.sample_detail = vertcat(sample_detail{:});
   results.activity_detail = vertcat(activity_detail{:});
   results.case_opts = vertcat(case_opts{:});
   results.benchmark = vertcat(benchmark{:});
   results.meta = vertcat(meta{:});
   results.artifact_file = vertcat(artifact_file{:});
   results.failed_cases = vertcat(failed_cases{:});
   results.passed = all(pass_flags);
   % Keep the single-model shape: one conditions struct and one model name
   % per element, so a consumer can see which model failed which condition.
   results.quality = struct('passed', all(quality_flags), ...
      'conditions', reshape(quality_conditions, [], 1), ...
      'reasons', vertcat(quality_reasons{:}), ...
      'smbmodels', reshape(quality_models, [], 1));
end

function artifact_file = saveArtifacts(sample_detail, ...
      activity_detail, case_summary, case_opts, benchmark, meta)
   %saveArtifacts Save the perf comparison artifact bundle for one run.

   % Build the artifact path.
   artifact_file = icemodel.test.helpers.artifactFilePath("perf", ...
      tier=meta.tier, smbmodel=meta.smbmodel_filter, ...
      solver=meta.solver_filter, baseline_type=meta.baseline_type, ...
      baseline_tag=meta.baseline_tag, run_name=meta.run_name, ...
      artifact_root=meta.artifact_root);

   % Create the run-specific artifact folder before saving.
   outdir = fileparts(artifact_file);
   if exist(outdir, 'dir') ~= 7
      mkdir(outdir);
   end

   % Save the benchmark struct as-is so the artifact content matches the
   % results struct returned by run_perf_suite.
   save(artifact_file, 'sample_detail', 'activity_detail', ...
      'case_summary', 'case_opts', 'benchmark', 'meta');

   % Print the loaded filename to the console.
   icemodel.test.helpers.printFilePath(artifact_file, "save");
end

function cleanup = configurePerfDataRootEnv(data_root)
   %CONFIGUREPERFDATAROOTENV Scope the runner-selected formal data root.

   % The performance TestCase bootstraps once per measured case, so expose the
   % outer runner selection explicitly and restore any interactive prior value.
   name = 'ICEMODEL_TEST_DATA_ROOT';
   old_value = getenv(name);
   cleanup = onCleanup(@() setenv(name, old_value));
   setenv(name, char(data_root));
end
