function results = run_perf_suite(kwargs)
   %RUN_PERF_SUITE Run formal icemodel performance benchmark suite.
   %
   %  results = run_perf_suite()
   %  results = run_perf_suite(tier="smoke")
   %  results = run_perf_suite(tier="smoke", smbmodel="skinmodel")
   %  results = run_perf_suite(tier="full", baseline="v1.1")
   %
   % Use this for normal performance comparisons against an existing rolling
   % or release baseline. This function does not update baselines; it only
   % runs the formal cases, compares runtime to the requested baseline, and
   % writes one artifact under test/artifacts/<run_name>/.
   %
   % CLI entrypoint:
   %  matlab -batch "run('/ABS/PATH/icemodel/test/run_perf_suite.m')"

   arguments (Input)
      kwargs.tier (1, :) string {mustBeMember(kwargs.tier, ...
         ["smoke", "full", "all"])} = "smoke"
      kwargs.smbmodel (1, :) string {mustBeMember(kwargs.smbmodel, ...
         ["all", "icemodel", "skinmodel"])} = "all"
      kwargs.n_runs (1, 1) double {mustBeInteger, mustBePositive} = 3
      kwargs.tol_perf (1, 1) double {mustBePositive} = 0.20
      kwargs.baseline (1, :) string = "rolling"
      kwargs.run_name string = string.empty()
   end
   [tier, smbmodel, n_runs, tol_perf, baseline_tag, run_name] = deal( ...
      kwargs.tier, kwargs.smbmodel, kwargs.n_runs, kwargs.tol_perf, ...
      kwargs.baseline, kwargs.run_name);

   % Resolve the requested baseline and shared batch run identifier.
   [baseline_type, baseline_tag] = test.helpers.resolveBaselineSelector( ...
      baseline_tag);
   [run_date, run_id, run_name] = test.helpers.resolveRunStamp(run_name);

   % Ensure test/ is on path and configure the model data paths.
   thisdir = fileparts(mfilename('fullpath'));
   rootdir = fileparts(thisdir);
   addpath(thisdir);
   [input_path, output_path] = test.helpers.configureModelPaths(rootdir);

   % Build the deterministic case list and load the requested baseline.
   cases = test.helpers.getCaseMatrix(tier, smbmodel);
   if isempty(cases)
      error('no performance cases matched tier=%s smbmodel=%s', tier, smbmodel)
   end
   suite = testsuite(fullfile(thisdir, 'IcemodelPerfTest.m'));
   experiment = matlab.perftest.TimeExperiment.withFixedSampleSize( ...
      n_runs, 'NumWarmups', 1);
   baseline = test.helpers.loadPerfBaseline(cases.simyear(1), baseline_tag, ...
      smbmodel);

   % Initial values.
   [sample_rows, activity_rows, case_rows, case_opts] = deal(struct([]));
   [r_sample, r_activity, r_case] = deal(0);
   failed_cases = strings(0, 1);

   % Run the per-case performance experiment and compare to baseline.
   for icase = 1:height(cases)

      % Run the test
      c = cases(icase, :);
      perf_data = test.helpers.runPerfCase(experiment, suite, c);

      % Extract results
      samples = perf_data.samples;
      activity = perf_data.activity;
      sample_times = perf_data.sample_times;
      activity_times = perf_data.activity_times;
      valid = perf_data.valid;

      % Compare to baseline
      ref_wall = nan;
      gate_wall = nan;
      passed_perf = valid;
      bid = test.helpers.findCaseRow(baseline, string(c.case_id));
      if ~isempty(bid)
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
      end

      % Construct rows to build new results tables
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
      case_rows(r_case).valid = valid;
      case_rows(r_case).passed_perf = passed_perf;
      case_rows(r_case).last_updated_utc = datetime('now', 'TimeZone', 'UTC');

      case_opts(r_case).case_id = string(c.case_id);
      case_opts(r_case).case = table2struct(c);
      case_opts(r_case).opts = test.helpers.buildFormalCaseOpts(c);

      if ~passed_perf
         failed_cases(end+1, 1) = string(c.case_id); %#ok<AGROW>
      end
   end

   % Build new results tables from rows.
   sample_detail = struct2table(sample_rows);
   activity_detail = struct2table(activity_rows);
   case_summary = struct2table(case_rows);

   % Record the run metadata and save the compare artifact.
   meta = struct();
   meta.tier = tier;
   meta.smbmodel_filter = smbmodel;
   meta.baseline_type = baseline_type;
   meta.baseline_tag = baseline_tag;
   meta.run_date = run_date;
   meta.run_id = run_id;
   meta.run_name = run_name;
   meta.baseline_file = perfBaselineFile( ...
      rootdir, cases.simyear(1), baseline_type, baseline_tag, smbmodel);
   meta.case_builder = "test.helpers.buildFormalCaseOpts";
   meta.opts_source = "icemodel.setopts defaults";
   meta.reset_fields = "solver";
   meta.n_runs = n_runs;
   meta.n_warmups = 1;
   meta.tol_perf = tol_perf;
   meta.experiment = "matlab.perftest.TimeExperiment.withFixedSampleSize";
   meta.input_path = string(input_path);
   meta.output_path = string(output_path);
   meta.suite_file = string(fullfile(thisdir, 'IcemodelPerfTest.m'));
   meta.matlab_version = string(version);
   meta.host = string(computer);
   meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');

   artifact_file = logArtifacts(rootdir, sample_detail, activity_detail, ...
      case_summary, case_opts, meta);

   results = struct();
   results.case_summary = case_summary;
   results.sample_detail = sample_detail;
   results.activity_detail = activity_detail;
   results.case_opts = case_opts;
   results.meta = meta;
   results.artifact_file = string(artifact_file);
   results.failed_cases = failed_cases;
   results.passed = isempty(failed_cases);
end

function artifact_file = logArtifacts(rootdir, sample_detail, ...
      activity_detail, case_summary, case_opts, meta)
   % Save one performance artifact for this compare run.

   outdir = fullfile(rootdir, 'test', 'artifacts', char(meta.run_name));
   if exist(outdir, 'dir') ~= 7
      mkdir(outdir);
   end
   if meta.baseline_type == "rolling"
      baseline_label = 'vs_rolling';
   else
      baseline_label = "vs_" + test.helpers.sanitizeTag(meta.baseline_tag);
   end
   model_label = smbmodelLabel(meta.smbmodel_filter);
   artifact_file = fullfile(outdir, ...
      sprintf('perf_results_%s%s_%s.mat', ...
      char(meta.tier), char(model_label), char(baseline_label)));
   save(artifact_file, 'sample_detail', 'activity_detail', ...
      'case_summary', 'case_opts', 'meta');
   disp(artifact_file)
   disp(case_summary(:, {'case_id', 'median_wall_s', 'ref_wall_s', ...
      'gate_wall_s', 'passed_perf'}))
end

function pathname = perfBaselineFile(rootdir, simyear, ...
      baseline_type, baseline_tag, smbmodel)

   if smbmodel == "all"
      models = test.helpers.formalSmbmodels();
      pathname = strings(numel(models), 1);
      for i = 1:numel(models)
         pathname(i) = perfBaselineFile(rootdir, simyear, baseline_type, ...
            baseline_tag, models(i));
      end
      return
   end

   model_tag = test.helpers.smbmodelTag(smbmodel);

   if baseline_type == "rolling"
      pathname = fullfile(rootdir, 'test', 'baselines', ...
         sprintf('perf_baseline_%d_rolling_%s.mat', simyear, model_tag));
   else
      pathname = fullfile(rootdir, 'test', 'baselines', ...
         sprintf('perf_baseline_%d_%s_%s.mat', simyear, ...
         test.helpers.sanitizeTag(baseline_tag), model_tag));
   end
end

function label = smbmodelLabel(smbmodel)
   smbmodel = string(smbmodel);
   if any(strcmpi(smbmodel, "all"))
      label = "";
   else
      label = "_" + test.helpers.smbmodelTag(smbmodel);
   end
end
