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
   %
   % Use this when you want to accept new runtime measurements as a rolling
   % or versioned perf baseline. This writes baseline files only; it does not
   % produce compare artifacts or evaluate pass/fail against an older baseline.
   %
   % Formal perf cases use the canonical suite runtime contract: one leading
   % spinup year plus one retained output year when the case matrix carries
   % only SIMYEAR.
   %
   % The saved MAT file also carries the managed core benchmark timings and,
   % by default, a saved profiler report from a separate diagnostic rerun.
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
         = 0.20

      kwargs.include_benchmarks (1, 1) logical ...
         = true

      kwargs.benchmark_sampling_profile (1, :) string ...
         {icemodel.validators.mustBeBenchmarkSamplingProfileName( ...
         kwargs.benchmark_sampling_profile)} ...
         = "default"

      kwargs.include_profile_artifacts (1, 1) logical ...
         = true

      kwargs.profile_history_size (1, 1) double {mustBeInteger, ...
         mustBePositive} ...
         = 25000000

      kwargs.output_file string ...
         = string.empty()
   end

   % Keep the cleanup handle in scope so the caller's config is restored
   % when this entrypoint returns.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Deal out inputs.
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

   % A custom output file is only coherent for one concrete model build.
   if numel(models) > 1 && ~isblanktext(output_file)
      error(['output_file overrides only one managed baseline file. Omit ', ...
         'it when smbmodel expands to more than one formal model.'])
   end

   % Build the baselines.
   baselines = arrayfun(@(mdl) buildSingleModelPerfBaseline( ...
      baseline, baseline_tag, tier, mdl, solver, simyear, ...
      smoke_sites, full_sites, n_runs, tol_perf, include_benchmarks, ...
      benchmark_sampling_profile, ...
      include_profile_artifacts, ...
      profile_history_size, output_file), ...
      models, 'UniformOutput', false);

   % Collapse to a single table.
   PerfBaseline = vertcat(baselines{:});
end

function PerfBaseline = buildSingleModelPerfBaseline(baseline, ...
      baseline_tag, tier, smbmodel, solver, simyear, smoke_sites, ...
      full_sites, n_runs, tol_perf, include_benchmarks, ...
      benchmark_sampling_profile, ...
      include_profile_artifacts, ...
      profile_history_size, output_file)
   %BUILDSINGLEMODELPERFBASELINE Build one canonical perf baseline file.

   % Resolve the baseline target, configure paths, and load formal cases.
   [baseline_type, baseline_tag, output_file, input_path, output_path, ...
      cases] = icemodel.test.helpers.prepareBaselineBuild( ...
      "perf", baseline, baseline_tag, tier, smbmodel, output_file, simyear, ...
      solver, smoke_sites, full_sites);
   testdir = icemodel.getpath('test');

   % Set up the formal perf class and one fixed-sample experiment that will
   % be reused across the accepted case matrix.
   suite = testsuite(fullfile(testdir, 'regression', ...
      'IcemodelPerfTest.m'));
   experiment = matlab.perftest.TimeExperiment.withFixedSampleSize( ...
      n_runs, 'NumWarmups', 1);

   % Preallocate row containers for the accepted baseline summary and opts.
   rows = struct([]);
   case_opts = struct([]);
   k = 0;

   % Measure each formal case and save the accepted timing summary.
   for icase = 1:height(cases)
      c = cases(icase, :);
      fprintf('Perf baseline case %d/%d: %s\n', ...
         icase, height(cases), c.case_id)
      perf_data = icemodel.test.helpers.runPerfCase(experiment, suite, c);
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
      rows(k).valid = perf_data.valid;
      rows(k).passed_perf = perf_data.valid;
      rows(k).last_updated_utc = datetime('now', 'TimeZone', 'UTC');

      case_opts(k).case_id = string(c.case_id);
      case_opts(k).case = table2struct(c);
      case_opts(k).opts = icemodel.test.helpers.setModelOptsForCase(c);
   end

   % Convert the accepted case rows into the saved baseline table.
   PerfBaseline = struct2table(rows);

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

   % Rolling baselines are acceptance targets. Archive the prior managed
   % state before overwriting it so older accepted timings remain available.
   if baseline_type == "rolling"
      icemodel.test.helpers.archiveManagedBaseline(output_file, "perf");
   end

   % Attach the managed component benchmark baseline to the same file so the
   % accepted end-to-end timings and their supporting kernel diagnostics stay
   % linked.
   BenchmarkBaseline = table();
   benchmark_meta = struct();
   if include_benchmarks
      [BenchmarkBaseline, benchmark_meta] = buildBenchmarkBaseline( ...
         sampling_profile=benchmark_sampling_profile);
      if ~isempty(BenchmarkBaseline)
         n_rows = height(BenchmarkBaseline);
         BenchmarkBaseline.baseline_type = repmat(baseline_type, n_rows, 1);
         BenchmarkBaseline.baseline_tag = repmat(baseline_tag, n_rows, 1);
         BenchmarkBaseline.last_updated_utc = repmat( ...
            datetime('now', 'TimeZone', 'UTC'), n_rows, 1);
      end
      benchmark_meta.baseline_type = baseline_type;
      benchmark_meta.baseline_tag = baseline_tag;
      benchmark_meta.source = "run_benchmark_suite";
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
