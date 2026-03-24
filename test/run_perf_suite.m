function results = run_perf_suite(kwargs)
   %RUN_PERF_SUITE Run formal icemodel performance benchmark suite.
   %
   %  results = run_perf_suite()
   %  results = run_perf_suite(tier="smoke")
   %  results = run_perf_suite(tier="smoke", smbmodel="skinmodel")
   %  results = run_perf_suite(tier="smoke", smbmodel="icemodel", solver=2)
   %  results = run_perf_suite(tier="smoke", smbmodel="icemodel", solver=[2 3])
   %  results = run_perf_suite(tier="smoke", smbmodel="icemodel", solver=[1 3])
   %  results = run_perf_suite(simyear=2017, smoke_sites="kanm", ...
   %     full_sites=["kanm"; "kanl"])
   %  results = run_perf_suite(tier="full", baseline="v1.1")
   %
   % Use this for normal performance comparisons against an existing rolling or
   % release baseline.
   %
   % This function does not update baselines; it only runs the formal cases,
   % compares runtime to the requested baseline, and writes one artifact under
   % test/artifacts/<run_name>/.
   %
   % By default it also runs the managed core benchmark suite and saves those
   % diagnostic timings alongside the formal perf artifact.
   %
   % Formal perf cases use the same canonical runtime contract as regression:
   % one leading spinup year plus one retained output year when the case
   % matrix carries only SIMYEAR.
   %
   % The optional solver filter accepts any subset of [1 2 3].
   %
   % SMOKE_SITES and FULL_SITES are advanced overrides for the site lists
   % used by each formal tier. Most callers should leave them at the
   % defaults and select only TIER.
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
         = 0.20

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
   end

   % Deal out arguments.
   [tier, smbmodel, solver, simyear, smoke_sites, full_sites, n_runs, ...
      tol_perf, include_benchmarks, benchmark_sampling_profile, ...
      baseline_selector, run_name] = deal( ...
      kwargs.tier, kwargs.smbmodel, kwargs.solver, kwargs.simyear, ...
      reshape(kwargs.smoke_sites, [], 1), reshape(kwargs.full_sites, [], 1), ...
      kwargs.n_runs, kwargs.tol_perf, kwargs.include_benchmarks, ...
      kwargs.benchmark_sampling_profile, ...
      kwargs.baseline, kwargs.run_name);

   % Resolve full path to the test/ dir.
   testdir = icemodel.getpath('test');

   % Resolve the requested baseline and shared batch run identifier.
   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector(baseline_selector);

   [run_date, run_id, run_name] = ...
      icemodel.test.helpers.resolveRunStamp(run_name);

   % Bootstrap the source/test trees once, then configure the formal paths.
   % Keep the cleanup handle in scope so the caller's config is restored
   % when this entrypoint returns.
   [~, input_path, output_path, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Expand the requested formal model selector once at the entrypoint.
   models = icemodel.test.helpers.resolveRequestedSmbmodels(smbmodel);

   % Build the MATLAB perf experiment once, then reuse it for each
   % single-model perf workflow below.
   suite = testsuite(fullfile(testdir, 'regression', ...
      'IcemodelPerfTest.m'));
   experiment = matlab.perftest.TimeExperiment.withFixedSampleSize( ...
      n_runs, 'NumWarmups', 1);

   % Run the canonical single-model workflow for each requested model and
   % merge the saved compare summaries into one returned struct.
   per_model = arrayfun(@(mdl) runSingleModelPerfSuite( ...
      input_path, output_path, testdir, experiment, suite, tier, ...
      mdl, solver, simyear, smoke_sites, full_sites, baseline_type, ...
      baseline_tag, run_date, run_id, run_name, n_runs, tol_perf, ...
      include_benchmarks, benchmark_sampling_profile), ...
      models, 'UniformOutput', false);

   % Combine results into a common struct.
   results = combinePerfResults(per_model);

   % TODO: Display the results
   % icemodel.test.helpers.displayPerfResults(results)
end

function results = runSingleModelPerfSuite(input_path, output_path, ...
      thisdir, experiment, suite, tier, smbmodel, solver, simyear, ...
      smoke_sites, full_sites, baseline_type, baseline_tag, run_date, ...
      run_id, run_name, n_runs, tol_perf, include_benchmarks, ...
      benchmark_sampling_profile)
   %RUNSINGLEMODELPERFSUITE Run the formal perf workflow for one smbmodel.

   % Build the deterministic case list and load the matching managed baseline.
   cases = icemodel.test.helpers.getPerfCaseMatrix( ...
      tier=tier, smbmodel=smbmodel, solver=solver, simyear=simyear, ...
      smoke_sites=smoke_sites, full_sites=full_sites);

   if isempty(cases)
      error('no performance cases matched tier=%s smbmodel=%s', tier, smbmodel)
   end

   % Perf baselines and benchmark baselines are keyed by one canonical
   % comparison year, so guard that contract before loading any baseline.
   benchmark_year = unique(cases.simyear);
   assert(isscalar(benchmark_year), ...
      'formal perf suite expects exactly one benchmark year')

   % Load the accepted baseline that matches this concrete formal model.
   [baseline, baseline_meta] = icemodel.test.helpers.loadPerfBaseline( ...
      benchmark_year, baseline_tag, smbmodel);
   [baseline_compatible, compare_reason] = perfBaselineCompatibility( ...
      baseline_meta);

   % Accumulate the measured sample/activity rows for the saved artifact.
   [sample_rows, activity_rows, case_rows, case_opts] = deal(struct([]));
   [r_sample, r_activity, r_case] = deal(0);
   failed_cases = strings(0, 1);

   % Run the per-case performance experiment and compare to baseline.
   for icase = 1:height(cases)

      % Run the case.
      c = cases(icase, :);
      perf_data = icemodel.test.helpers.runPerfCase(experiment, suite, c);

      % Extract results for this case.
      samples = perf_data.samples;
      activity = perf_data.activity;
      sample_times = perf_data.sample_times;
      activity_times = perf_data.activity_times;
      valid = perf_data.valid;

      % Compare the measured runtime to the accepted baseline row.
      ref_wall = nan;
      gate_wall = nan;
      passed_perf = valid;
      case_compare_reason = compare_reason;
      bid = icemodel.test.helpers.findCaseRow(baseline, string(c.case_id));
      if ~baseline_compatible
         passed_perf = valid;
      elseif ~isempty(bid)
         ref_wall = baseline.median_wall_s(bid);
         if isfinite(ref_wall) && ref_wall > 0
            if ismember('tol_perf', baseline.Properties.VariableNames) ...
                  && isfinite(baseline.tol_perf(bid)) ...
                  && baseline.tol_perf(bid) > 0
               tol_case = baseline.tol_perf(bid);
            else
               tol_case = tol_perf;
            end
            gate_wall = ref_wall * (1 + tol_case);
            passed_perf = passed_perf ...
               && median(sample_times, 'omitnan') <= gate_wall;
         end
      else
         case_compare_reason = "case not found in compatible perf baseline";
      end

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

      % Save the compact per-case summary and the resolved opts contract.
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
      case_rows(r_case).gate_wall_s = gate_wall;
      case_rows(r_case).baseline_compatible = baseline_compatible;
      case_rows(r_case).compare_reason = case_compare_reason;
      case_rows(r_case).valid = valid;
      case_rows(r_case).passed_perf = passed_perf;
      case_rows(r_case).last_updated_utc = datetime('now', 'TimeZone', 'UTC');

      case_opts(r_case).case_id = string(c.case_id);
      case_opts(r_case).case = table2struct(c);
      case_opts(r_case).opts = icemodel.test.helpers.setModelOptsForCase(c);

      if ~passed_perf
         failed_cases(end+1, 1) = string(c.case_id); %#ok<AGROW>
      end
   end

   % Build the saved artifact tables for this concrete formal model.
   sample_detail = struct2table(sample_rows);
   activity_detail = struct2table(activity_rows);
   case_summary = struct2table(case_rows);

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
   meta.baseline_file = icemodel.test.helpers.defaultBaselinePath( ...
      "perf", baseline_type, baseline_tag, smbmodel, benchmark_year);
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

   % Report if whole-model perf comparison was skipped
   if ~meta.baseline_compatible && ~isblanktext(meta.compare_reason)
      fprintf('Whole-model perf comparison skipped: %s\n', ...
         char(meta.compare_reason));
   end

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

   % Display the perf report for this case.
   % TODO: pass results to renamed displayPerfResults in main function
   icemodel.test.helpers.displayPerfSummary(case_summary, benchmark)
end

function results = combinePerfResults(per_model)
   %COMBINEPERFRESULTS Merge one-or-more single-model perf run results.

   if isscalar(per_model)
      results = per_model{1};
      return
   end

   % Extract each returned field once, then concatenate the per-model pieces.
   % This helper exists because the aggregate struct carries both table-like
   % fields and scalar pass/fail metadata that cannot be collapsed with a
   % single blind vertcat call.
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
end

function artifact_file = saveArtifacts(sample_detail, ...
      activity_detail, case_summary, case_opts, benchmark, meta)
   %saveArtifacts Save the perf comparison artifact bundle for one run.

   % Create the run-specific artifact folder before saving the report.
   testdir = icemodel.getpath('test');
   outdir = fullfile(testdir, 'artifacts', char(meta.run_name));
   if exist(outdir, 'dir') ~= 7
      mkdir(outdir);
   end

   % Format the baseline/model/solver tags used by the saved filename.
   if meta.baseline_type == "rolling"
      baseline_label = 'vs_rolling';
   else
      baseline_label = "vs_" + icemodel.test.helpers.sanitizeTag(meta.baseline_tag);
   end
   model_label = smbmodelLabel(meta.smbmodel_filter);
   solver_label = solverLabel(meta.solver_filter);
   artifact_file = fullfile(outdir, ...
      sprintf('perf_results_%s%s_%s.mat', ...
      char(meta.tier), char(model_label + solver_label), char(baseline_label)));
   benchmark_summary = benchmark.summary;
   benchmark_comparison = benchmark.comparison;
   benchmark_meta = benchmark.meta;

   save(artifact_file, 'sample_detail', 'activity_detail', ...
      'case_summary', 'case_opts', 'benchmark_summary', ...
      'benchmark_comparison', 'benchmark_meta', 'meta');
   disp(artifact_file)
end

function [compatible, reason] = perfBaselineCompatibility(baseline_meta)
   %PERFBASELINECOMPATIBILITY Decide whether wall-time comparison is fair.

   compatible = false;
   reason = "";

   if ~isstruct(baseline_meta) || isempty(fieldnames(baseline_meta))
      reason = "perf baseline metadata not found";
      return
   end

   if ~isfield(baseline_meta, 'matlab_version') || ...
         ~isfield(baseline_meta, 'host') || ...
         isblanktext(baseline_meta.matlab_version) || ...
         isblanktext(baseline_meta.host)
      reason = "perf baseline predates environment metadata";
      return
   end

   current_version = string(version);
   current_host = string(computer);
   baseline_version = string(baseline_meta.matlab_version);
   baseline_host = string(baseline_meta.host);
   compatible = current_version == baseline_version && ...
      current_host == baseline_host;

   if ~compatible
      reason = sprintf([ ...
         'baseline built under MATLAB %s on %s; current environment is ', ...
         'MATLAB %s on %s'], char(baseline_version), char(baseline_host), ...
         char(current_version), char(current_host));
   end
end

function label = smbmodelLabel(smbmodel)
   %SMBMODELLABEL Format the smbmodel selector for perf artifact filenames.
   label = "_" + icemodel.test.helpers.smbmodelTag(string(smbmodel));
end

function label = solverLabel(solver)
   %SOLVERLABEL Format the solver filter for perf artifact filenames.
   if isempty(solver)
      label = "";
   else
      label = "_s" + join(string(solver), '-');
   end
end
