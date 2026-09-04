function PerfBaseline = build_perf_baseline(kwargs)
   %BUILD_PERF_BASELINE Build rolling or versioned model performance baselines.
   %
   %  PerfBaseline = build_perf_baseline(baseline="rolling")
   %  PerfBaseline = build_perf_baseline(baseline_tag="v1.1")
   %  PerfBaseline = build_perf_baseline(baseline_tag="v1.1", tier="full", ...
   %     smbmodel="skinmodel")
   %  PerfBaseline = build_perf_baseline(baseline="rolling", ...
   %     smbmodel="icemodel", solver=2)
   %  PerfBaseline = build_perf_baseline(baseline="rolling", ...
   %     smbmodel="icemodel", solver=[1 3])
   %  PerfBaseline = build_perf_baseline(simyear=2017, smoke_sites="kanm", ...
   %     full_sites=["kanm"; "kanl"])
   %  PerfBaseline = build_perf_baseline(data_root="/path/to/test/data")
   %
   % Use this when you want to accept new runtime measurements as a rolling
   % or versioned perf baseline. This writes baseline files only; it does not
   % produce compare artifacts or evaluate pass/fail against an older baseline.
   %
   % Formal perf cases always run one leading spinup year plus one
   % retained output year when the case matrix carries only SIMYEAR.
   %
   % ISOLATION selects the measurement protocol and defaults to "process",
   % matching run_perf_suite: every case measures in a fresh
   % `matlab -batch` subprocess, and the saved metadata records the
   % protocol so the compatibility check can pair baseline and run. The
   % opt-in "session" mode refuses a session that already ran another
   % suite.
   %
   % The saved MAT file also carries the managed core benchmark timings.
   % Profiler artifacts are an opt-in, single-model diagnostic.
   %
   % A custom OUTPUT_FILE is supported only when SMBMODEL resolves to one
   % concrete formal model. Multi-model requests write the managed per-model
   % baseline files under test/baselines/.
   %
   % The optional solver filter accepts any subset of [1 2 3].
   % The formal benchmark year and smoke/full site selections are explicit
   % here rather than hidden in the case-matrix helper.

   arguments (Input)

      kwargs.baseline (1, :) string ...
         {icemodel.validators.mustBeRollingBaselineName(kwargs.baseline)} ...
         = "rolling"

      kwargs.baseline_tag string ...
         = string.empty()

      kwargs.tier (1, :) string ...
         {icemodel.validators.mustBeTestTierName(kwargs.tier)} ...
         = "full"

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

      kwargs.include_profile_artifacts (1, 1) logical ...
         = false

      kwargs.profile_history_size (1, 1) double {mustBeInteger, ...
         mustBePositive} ...
         = 25000000

      kwargs.output_file string ...
         = string.empty()

      kwargs.data_root (1, 1) string ...
         = ""

      % Process isolation is the default because it is the only protocol
      % that supports a formal accept/reject verdict against this
      % baseline. Session mode is opt-in for quick diagnostics.
      kwargs.isolation (1, 1) string ...
         {mustBeMember(kwargs.isolation, ["session", "process"])} ...
         = "process"
   end

   % Resolve the baseline-owned default tree before installing scoped config.
   baseline_selector = kwargs.baseline_tag;
   if isblanktext(baseline_selector)
      baseline_selector = kwargs.baseline;
   end
   baseline_policy = ...
      icemodel.test.helpers.formalBaselinePolicy(baseline_selector);

   % Keep the cleanup handle in scope so the caller's config is restored
   % when this entrypoint returns.
   [~, input_path, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment( ...
      icemodel_config_casename=baseline_policy.config_case, ...
      data_root=kwargs.data_root);

   % The perf TestCase bootstraps each measured case, so retain the configured
   % root through the environment the runner sets for the class. A caller-
   % supplied root keeps precedence over the resolved verification root.
   data_root = kwargs.data_root;
   if isblanktext(data_root)
      data_root = string(fileparts(input_path));
   end
   data_root_name = 'ICEMODEL_TEST_DATA_ROOT';
   prior_data_root = getenv(data_root_name);
   data_root_cleanup = onCleanup( ...
      @() setenv(data_root_name, prior_data_root));
   setenv(data_root_name, char(data_root));

   % Baseline timings must never inherit an interactive profiler.
   profile off

   % The component benchmark suite runs in this session before the model
   % cases, so under session isolation it would warm the session the
   % model medians are then measured in. Refuse the combination before
   % any session-state check: it is invalid no matter the session state.
   if kwargs.isolation == "session" && kwargs.include_benchmarks
      error('icemodel:test:perf:benchmarksWarmSession', ...
         ['include_benchmarks=true warms the session before the model ', ...
         'cases. Use isolation="process" (the default) or pass ', ...
         'include_benchmarks=false for an in-session build.'])
   end

   % An in-session baseline build in a dirty session records contaminated
   % timings; process isolation is immune.
   icemodel.test.helpers.assertCleanPerfSession(kwargs.isolation);

   % Record this build in the session activity so a later in-session
   % formal run refuses the session this matrix already warmed.
   icemodel.test.helpers.markTestSessionDirty("build_perf_baseline");

   % Unpack the parsed inputs.
   [baseline, baseline_tag, tier, smbmodel, solver, simyear, smoke_sites, ...
      full_sites, n_runs, tol_perf, include_benchmarks, ...
      benchmark_sampling_profile, ...
      include_profile_artifacts, ...
      profile_history_size, output_file] ...
      = deal( ...
      kwargs.baseline, kwargs.baseline_tag, kwargs.tier, kwargs.smbmodel, ...
      kwargs.solver, kwargs.simyear, reshape(kwargs.smoke_sites, [], 1), ...
      reshape(kwargs.full_sites, [], 1), kwargs.n_runs, kwargs.tol_perf, ...
      kwargs.include_benchmarks, kwargs.benchmark_sampling_profile, ...
      kwargs.include_profile_artifacts, ...
      kwargs.profile_history_size, ...
      kwargs.output_file);

   % Expand the requested formal model selector once at the entrypoint.
   models = icemodel.test.helpers.resolveRequestedSmbmodels(smbmodel);

   % Profiling more than one model in one process would warm code and data
   % before the later model's timing pass. Keep that diagnostic isolated.
   if include_profile_artifacts && numel(models) > 1
      error("icemodel:test:profileRequiresSingleModel", ...
         ['include_profile_artifacts=true requires one concrete smbmodel. ', ...
         'Run each profiled model in a fresh MATLAB process.'])
   end

   % A custom output file is only coherent for one concrete model build.
   if numel(models) > 1 && ~isblanktext(output_file)
      error(['output_file overrides only one managed baseline file. Omit ', ...
         'it when smbmodel expands to more than one formal model.'])
   end

   % Measure the managed component benchmarks once for the whole build.
   % They are model-independent, so measuring them inside the per-model
   % builder would repeat the same suite for every model.
   BenchmarkBaseline = table();
   benchmark_meta = struct();
   if include_benchmarks
      [BenchmarkBaseline, benchmark_meta] = buildBenchmarkBaseline( ...
         sampling_profile=benchmark_sampling_profile);
      icemodel.test.helpers.assertFormalBenchmarkCandidate(BenchmarkBaseline);
      benchmark_meta.source = "run_benchmark_suite";
   end

   % Build the baselines.
   baselines = arrayfun(@(mdl) buildSingleModelPerfBaseline( ...
      baseline, baseline_tag, tier, mdl, solver, simyear, ...
      smoke_sites, full_sites, n_runs, tol_perf, ...
      include_benchmarks, benchmark_sampling_profile, ...
      BenchmarkBaseline, benchmark_meta, ...
      include_profile_artifacts, ...
      profile_history_size, output_file, kwargs.isolation, ...
      baseline_policy.config_case, data_root), ...
      models, 'UniformOutput', false);

   % Collapse to a single table.
   PerfBaseline = vertcat(baselines{:});

   % Restore the caller config now that this entrypoint is done.
   delete(suite_cleanup)
end

function PerfBaseline = buildSingleModelPerfBaseline(baseline, ...
      baseline_tag, tier, smbmodel, solver, simyear, smoke_sites, ...
      full_sites, n_runs, tol_perf, include_benchmarks, ...
      benchmark_sampling_profile, BenchmarkBaseline, benchmark_meta, ...
      include_profile_artifacts, ...
      profile_history_size, output_file, isolation, config_case, ...
      data_root)
   %BUILDSINGLEMODELPERFBASELINE Build one perf baseline file.

   % Resolve the baseline target, configure paths, and load formal cases.
   [baseline_type, baseline_tag, output_file, input_path, output_path, ...
      cases] = icemodel.test.helpers.prepareBaselineBuild( ...
      "perf", baseline, baseline_tag, tier, smbmodel, output_file, simyear, ...
      solver, smoke_sites, full_sites);
   testdir = icemodel.getpath('test');

   % Accepted wall-clock measurements must not inherit profiler state from an
   % interactive session or a preceding diagnostic call.
   profile off

   % Set up the formal perf class and one fixed-sample experiment that will
   % be reused across the accepted case matrix.
   suite = testsuite(fullfile(testdir, 'regression', ...
      'IcemodelPerfTest.m'));
   experiment = matlab.perftest.TimeExperiment.withFixedSampleSize( ...
      n_runs, 'NumWarmups', 1);

   % Subprocess spec/result files land beside the run artifacts so the raw
   % accepted measurements stay auditable.
   artifact_dir = fileparts(icemodel.test.helpers.artifactFilePath( ...
      "perf", tier=tier, smbmodel=smbmodel, solver=solver, ...
      baseline_type=baseline_type, baseline_tag=baseline_tag, ...
      run_name="baseline_build"));
   if isolation == "process" && exist(artifact_dir, 'dir') ~= 7
      mkdir(artifact_dir);
   end

   % Preallocate row containers for the accepted baseline summary and opts.
   rows = struct([]);
   case_opts = struct([]);
   k = 0;

   % Randomize the case order so no case always inherits the same
   % predecessor's state. The metadata records the seed so the order is
   % reproducible. The saved rows are keyed by case_id, so the
   % comparison never depends on row order.
   case_order_seed = randi(2^31 - 2);
   rng_prior = rng(case_order_seed, 'twister');
   case_order = randperm(height(cases));
   rng(rng_prior);

   % Measure each formal case and save the accepted timing summary.
   for iorder = 1:height(cases)
      c = cases(case_order(iorder), :);
      fprintf('Perf baseline case %d/%d: %s\n', ...
         iorder, height(cases), c.case_id)
      [perf_data, valid_gate] = ...
         icemodel.test.helpers.retryInvalidMeasurement(@(attempt) ...
         icemodel.test.helpers.measurePerfCase(experiment, suite, c, ...
         isolation, config_case, data_root, n_runs, artifact_dir, ...
         attempt));
      sample_times = perf_data.sample_times;

      k = k + 1;
      rows(k).case_id = string(c.case_id);
      rows(k).tier = string(c.tier);
      rows(k).smbmodel = string(c.smbmodel);
      rows(k).sitename = string(c.sitename);
      rows(k).forcings = string(c.forcings);
      rows(k).simyear = c.simyear;
      rows(k).solver = c.solver;
      rows(k).baseline_type = baseline_type;
      rows(k).baseline_tag = baseline_tag;
      rows(k).smbmodel_filter = smbmodel;
      rows(k).n_runs = n_runs;
      rows(k).n_warmups = perf_data.n_warmups;
      rows(k).tol_perf = tol_perf;
      rows(k).median_wall_s = median(sample_times, 'omitnan');
      rows(k).mean_wall_s = mean(sample_times, 'omitnan');
      rows(k).min_wall_s = min(sample_times, [], 'omitnan');
      rows(k).max_wall_s = max(sample_times, [], 'omitnan');
      rows(k).ref_wall_s = nan;
      rows(k).gate_wall_s = nan;
      rows(k).valid = valid_gate;
      rows(k).passed_perf = valid_gate;
      rows(k).last_updated_utc = datetime('now', 'TimeZone', 'UTC');

      case_opts(k).case_id = string(c.case_id);
      case_opts(k).case = table2struct(c);
      case_opts(k).opts = icemodel.test.helpers.setModelOptsForCase(c);
   end

   % The per-case dispersion gate cannot see load shifts that are steady
   % WITHIN each case but different ACROSS cases. Re-measure the first
   % executed case (rows(1), by construction) and refuse to accept a
   % drifted or uncertifiable baseline; comparisons against it would gate
   % on biased medians.
   anchor_tol = icemodel.test.helpers.perfMeasurementPolicy().anchor_tol;
   c_anchor = cases(case_order(1), :);
   [anchor_data, anchor_valid] = ...
      icemodel.test.helpers.retryInvalidMeasurement(@(attempt) ...
      icemodel.test.helpers.measurePerfCase(experiment, suite, ...
      c_anchor, isolation, config_case, data_root, n_runs, ...
      artifact_dir, 98 + attempt));
   [ambient_stable, anchor_ratio] = ...
      icemodel.test.helpers.ambientAnchorVerdict( ...
      rows(1).median_wall_s, anchor_data.sample_times, anchor_valid, ...
      anchor_tol);
   if ~ambient_stable
      error('icemodel:test:perf:ambientDrift', ...
         ['ambient conditions shifted during the baseline build, or ' ...
         'the anchor re-measurement was invalid (anchor ratio %.3f); ' ...
         'a drifted baseline must not be accepted'], anchor_ratio)
   end

   % Convert the accepted case rows into the saved baseline table.
   PerfBaseline = struct2table(rows);
   selector = baseline_tag;
   if baseline_type == "rolling"
      selector = "rolling";
   end
   icemodel.test.helpers.assertFormalBaselineCandidate( ...
      "perf", PerfBaseline, cases, selector);

   % Record the build metadata alongside the accepted baseline values.
   meta = struct();
   meta.tier = tier;
   meta.smbmodel_filter = smbmodel;
   meta.simyear = simyear;
   meta.smoke_sites = smoke_sites;
   meta.full_sites = full_sites;
   meta.baseline_type = baseline_type;
   meta.baseline_tag = baseline_tag;
   meta.case_builder = "icemodel.test.helpers.setModelOptsForCase";
   meta.opts_source = "icemodel.setopts defaults";
   meta.spinup_policy = ...
      "formal perf runs include the canonical leading spinup year";
   meta.reset_fields = "solver";
   meta.n_runs = n_runs;
   meta.n_warmups = 1;
   meta.tol_perf = tol_perf;
   meta.isolation = isolation;
   meta.case_order_seed = case_order_seed;
   meta.anchor_ratio = anchor_ratio;
   meta.timing_scope = "IcemodelPerfTest.testCoreRuntime (runSmbModel only)";
   meta.timing_notes = sprintf([ ...
      'median_wall_s is the median of %d timed samples (wall-clock seconds). ' ...
      '%d warmup run(s) are executed first and excluded from all summary ' ...
      'statistics. Setup, teardown, and reporting overhead are outside the ' ...
      'measured region. The timed model call includes spinup and output years.'], ...
      n_runs, 1);
   meta.include_benchmarks = include_benchmarks;
   meta.benchmark_sampling_profile = benchmark_sampling_profile;
   meta.include_profile_artifacts = include_profile_artifacts;
   meta.profile_history_size = profile_history_size;
   meta.input_path = string(input_path);
   meta.output_path = string(output_path);
   meta.suite_file = string(fullfile(testdir, 'regression', ...
      'IcemodelPerfTest.m'));
   meta.matlab_version = string(version);
   meta.host = string(computer);
   meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');

   % Attach the managed component benchmark baseline (measured once at
   % the entrypoint) to the same file so the accepted end-to-end timings
   % and their supporting kernel diagnostics stay linked.
   if ~isempty(BenchmarkBaseline)
      n_rows = height(BenchmarkBaseline);
      BenchmarkBaseline.baseline_type = repmat(baseline_type, n_rows, 1);
      BenchmarkBaseline.baseline_tag = repmat(baseline_tag, n_rows, 1);
      BenchmarkBaseline.last_updated_utc = repmat( ...
         datetime('now', 'TimeZone', 'UTC'), n_rows, 1);
      benchmark_meta.baseline_type = baseline_type;
      benchmark_meta.baseline_tag = baseline_tag;
   end

   % Save profiler diagnostics in a separate rerun so the accepted timing
   % pass above stays focused on the managed perf measurements.
   profile_summary = table();
   profile_meta = struct();
   profile_artifacts = struct();
   if include_profile_artifacts
      [profile_summary, profile_meta, profile_artifacts] = ...
         icemodel.test.helpers.captureBaselineProfile( ...
         "perf", cases, output_file, history_size=profile_history_size);
   end

   % Archive only after the end-to-end, benchmark, and optional profile
   % candidates have all passed their own validation and completed.
   if baseline_type == "rolling"
      icemodel.test.helpers.archiveManagedBaseline(output_file, "perf");
   end

   % Save the rolling or release perf baseline file.
   outdir = fileparts(char(output_file));
   if exist(outdir, 'dir') ~= 7
      mkdir(outdir);
   end
   save(char(output_file), 'PerfBaseline', 'case_opts', 'meta', ...
      'BenchmarkBaseline', 'benchmark_meta', 'profile_summary', ...
      'profile_meta', 'profile_artifacts');
end

function [BenchmarkBaseline, meta] = buildBenchmarkBaseline(kwargs)
   %BUILDBENCHMARKBASELINE Measure the managed core benchmark suite once.

   arguments
      kwargs.sampling_profile (1, :) string ...
         {icemodel.validators.mustBeBenchmarkSamplingProfileName( ...
         kwargs.sampling_profile)} = "default"
   end

   % Reuse the public benchmark runner so the managed baseline reflects the
   % exact component suite developers can run directly.
   results = run_benchmark_suite( ...
      sampling_profile=kwargs.sampling_profile, show_summary=false);

   [suite_signature, suite_files] = ...
      icemodel.test.helpers.benchmarkSuiteSignature();

   BenchmarkBaseline = sampleSummary(results);
   if ~isempty(BenchmarkBaseline) ...
         && ismember('Name', BenchmarkBaseline.Properties.VariableNames)
      BenchmarkBaseline.Name = string(BenchmarkBaseline.Name);
   end
   BenchmarkBaseline.Valid = reshape(logical([results.Valid]), [], 1);

   % Record the benchmark experiment settings with the saved table.
   meta = struct();
   meta.sampling_profile = kwargs.sampling_profile;
   meta.include_subfolders = false;
   meta.suite_signature = suite_signature;
   meta.suite_files = suite_files;
   meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   meta.matlab_version = string(version);
   meta.host = string(computer);
end
