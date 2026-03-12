function PerfBaseline = build_perf_baseline(kwargs)
   %BUILD_PERF_BASELINE Build rolling or versioned model performance baselines.
   %
   %  PerfBaseline = build_perf_baseline(baseline="rolling")
   %  PerfBaseline = build_perf_baseline(baseline_tag="v1.1")
   %  PerfBaseline = build_perf_baseline(baseline="v1.1", tier="full", ...
   %     smbmodel="skinmodel")
   %  PerfBaseline = build_perf_baseline(baseline="rolling", ...
   %     smbmodel="icemodel", solver=2)
   %  PerfBaseline = build_perf_baseline(baseline="rolling", ...
   %     smbmodel="icemodel", solver=[1 3])
   %
   % Use this when you want to accept new runtime measurements as a rolling
   % or versioned perf baseline. This writes baseline files only; it does not
   % produce compare artifacts or evaluate pass/fail against an older baseline.
   % The optional solver filter accepts any subset of [1 2 3].

   arguments (Input)
      kwargs.baseline string = string.empty()
      kwargs.baseline_tag string = string.empty()
      kwargs.tier (1, :) string {mustBeMember(kwargs.tier, ...
         ["smoke", "full", "all"])} = "full"
      kwargs.smbmodel (1, :) string {mustBeMember(kwargs.smbmodel, ...
         ["all", "icemodel", "skinmodel"])} = "all"
      kwargs.solver {mustBeValidSolverFilter(kwargs.solver)} = []
      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      kwargs.n_runs (1, 1) double {mustBeInteger, mustBePositive} = 3
      kwargs.tol_perf (1, 1) double {mustBePositive} = 0.20
      kwargs.output_file string = string.empty()
   end
   [baseline, baseline_tag, tier, smbmodel, solver, simyear, n_runs, ...
      tol_perf, output_file] = deal(kwargs.baseline, kwargs.baseline_tag, ...
      kwargs.tier, kwargs.smbmodel, kwargs.solver, kwargs.simyear, ...
      kwargs.n_runs, kwargs.tol_perf, kwargs.output_file);

   % Expand smbmodel="all" into one rolling/release baseline file per model.
   if smbmodel == "all"
      if ~isblanktext(output_file)
         error('output_file is only supported for a single smbmodel')
      end
      models = test.helpers.formalSmbmodels();
      baselines = cell(numel(models), 1);
      for i = 1:numel(models)
         baselines{i} = build_perf_baseline( ...
            baseline=baseline, baseline_tag=baseline_tag, tier=tier, ...
            smbmodel=models(i), solver=solver, simyear=simyear, n_runs=n_runs, ...
            tol_perf=tol_perf);
      end
      PerfBaseline = vertcat(baselines{:});
      return
   end

   % Resolve the baseline target, configure paths, and load formal cases.
   [baseline_type, baseline_tag, output_file, rootdir, input_path, ...
      output_path, cases] = test.helpers.prepareBaselineBuild("perf", ...
      baseline, baseline_tag, tier, smbmodel, output_file, simyear, solver);

   % Set up the per-case MATLAB performance experiment once.
   suite = testsuite(fullfile(rootdir, 'test', 'IcemodelPerfTest.m'));
   experiment = matlab.perftest.TimeExperiment.withFixedSampleSize( ...
      n_runs, 'NumWarmups', 1);

   % Preallocate row containers for the accepted baseline summary and opts.
   rows = struct([]);
   case_opts = struct([]);
   k = 0;

   % Measure each formal case and save the accepted timing summary.
   for icase = 1:height(cases)
      c = cases(icase, :);
      perf_data = test.helpers.runPerfCase(experiment, suite, c);
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
      case_opts(k).opts = test.helpers.buildFormalCaseOpts(c);
   end

   % Convert the accepted case rows into the saved baseline table.
   PerfBaseline = struct2table(rows);

   % Record the build metadata alongside the accepted baseline values.
   meta = struct();
   meta.tier = tier;
   meta.smbmodel_filter = smbmodel;
   meta.baseline_type = baseline_type;
   meta.baseline_tag = baseline_tag;
   meta.case_builder = "test.helpers.buildFormalCaseOpts";
   meta.opts_source = "icemodel.setopts defaults";
   meta.reset_fields = "solver";
   meta.n_runs = n_runs;
   meta.n_warmups = 1;
   meta.tol_perf = tol_perf;
   meta.input_path = string(input_path);
   meta.output_path = string(output_path);
   meta.suite_file = string(fullfile(rootdir, 'test', 'IcemodelPerfTest.m'));
   meta.matlab_version = string(version);
   meta.host = string(computer);
   meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');

   % Save the rolling or release perf baseline file.
   outdir = fileparts(char(output_file));
   if exist(outdir, 'dir') ~= 7
      mkdir(outdir);
   end
   save(char(output_file), 'PerfBaseline', 'case_opts', 'meta');
end

function mustBeValidSolverFilter(x)
   if isempty(x)
      return
   end
   mustBeNumeric(x)
   mustBeInteger(x)
   if any(~ismember(x, [1 2 3]))
      error('solver must be empty or a subset of [1 2 3]')
   end
end
