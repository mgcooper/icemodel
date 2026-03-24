function RegressionBaseline = build_regression_baseline(kwargs)
   %BUILD_REGRESSION_BASELINE Build rolling or versioned icemodel regression
   %baselines.
   %
   %  RegressionBaseline = build_regression_baseline(baseline="rolling")
   %  RegressionBaseline = build_regression_baseline(baseline_tag="v1.1")
   %  RegressionBaseline = build_regression_baseline(baseline_tag="v1.01", ...
   %     tier="full", smbmodel="skinmodel")
   %  RegressionBaseline = build_regression_baseline(baseline="rolling", ...
   %     smbmodel="icemodel", solver=2)
   %  RegressionBaseline = build_regression_baseline(baseline="rolling", ...
   %     smbmodel="icemodel", solver=[1 3])
   %  RegressionBaseline = build_regression_baseline(simyear=2017, ...
   %     smoke_sites="kanm", full_sites=["kanm"; "kanl"])
   %
   % Use this when you want to accept new modeled outputs as a rolling or
   % versioned regression baseline. This writes baseline files only; it does
   % not produce compare artifacts or evaluate pass/fail against an older
   % baseline. By default it also saves a profiler report from a separate
   % diagnostic rerun so the accepted baseline and the timing diagnostics are
   % archived together.
   %
   % Formal regression cases use the canonical suite runtime contract: one
   % leading spinup year plus one retained output year when the case matrix
   % carries only SIMYEAR.
   %
   % A custom OUTPUT_FILE is supported only when SMBMODEL resolves to one
   % concrete formal model. Multi-model requests write the managed per-model
   % baseline files under test/baselines/.
   %
   % The optional solver filter accepts any subset of [1 2 3].
   % The formal comparison year and smoke/full site selections are explicit
   % here rather than buried in the regression case matrix helper.

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

   % Deal out arguments.
   [baseline, baseline_tag, tier, smbmodel, solver, simyear, smoke_sites, ...
      full_sites, include_profile_artifacts, ...
      profile_history_size, output_file] = deal( ...
      kwargs.baseline, kwargs.baseline_tag, kwargs.tier, kwargs.smbmodel, ...
      kwargs.solver, kwargs.simyear, reshape(kwargs.smoke_sites, [], 1), ...
      reshape(kwargs.full_sites, [], 1), ...
      kwargs.include_profile_artifacts, kwargs.profile_history_size, ...
      kwargs.output_file);

   % Expand the requested formal model selector once at the entrypoint.
   models = icemodel.test.helpers.resolveRequestedSmbmodels(smbmodel);

   % A custom output file is only coherent for one concrete model build.
   if numel(models) > 1 && ~isblanktext(output_file)
      error(['output_file overrides only one managed baseline file. Omit ', ...
         'it when smbmodel expands to more than one formal model.'])
   end

   % Build the baselines.
   baselines = arrayfun(@(mdl) buildSingleModelRegressionBaseline( ...
      baseline, baseline_tag, tier, mdl, solver, simyear, ...
      smoke_sites, full_sites, include_profile_artifacts, ...
      profile_history_size, output_file), ...
      models, 'UniformOutput', false);

   % Collapse to a single table.
   RegressionBaseline = vertcat(baselines{:});
end

function RegressionBaseline = buildSingleModelRegressionBaseline( ...
      baseline, baseline_tag, tier, smbmodel, solver, simyear, ...
      smoke_sites, full_sites, ...
      include_profile_artifacts, profile_history_size, output_file)
   %BUILDSINGLEMODELREGRESSIONBASELINE Build one canonical regression baseline.

   % Resolve the baseline target, configure paths, and load formal cases.
   [baseline_type, baseline_tag, output_file, input_path, output_path, ...
      cases] = icemodel.test.helpers.prepareBaselineBuild( ...
      "regression", baseline, baseline_tag, tier, smbmodel, output_file, ...
      simyear, solver, smoke_sites, full_sites);

   % Load the static runoff reference once before entering the case loop.
   runoff_ref = icemodel.test.helpers.loadRunoffReference();

   % Preallocate row containers for the accepted baseline summary and opts.
   row_cells = cell(height(cases), 1);
   case_opts = struct([]);

   % Run each formal case and save the accepted scalar regression state.
   for icase = 1:height(cases)
      c = cases(icase, :);
      fprintf('Regression baseline case %d/%d: %s\n', ...
         icase, height(cases), c.case_id)
      opts_run = icemodel.test.helpers.setModelOptsForCase(c);

      % Run the case and postprocess the retained output years into the
      % scalar metrics carried by the regression baseline.
      [ice1, ice2] = icemodel.test.helpers.runSmbModel(opts_run);
      [ice1, ~] = icemodel.postprocess( ...
         ice1, ice2, opts_run, opts_run.output_years);

      % Load the matched runoff reference row, if one exists, before
      % summarizing the retained yearly outputs.
      idx = icemodel.test.helpers.findRunoffReferenceRow(runoff_ref, c);
      met = icemodel.test.helpers.loadProcessedMetForOutputYears(opts_run);
      if isempty(idx)
         refrow = [];
      else
         refrow = runoff_ref(idx, :);
      end
      S = icemodel.test.helpers.summarizeIce1Metrics(ice1, met, refrow);

      row = struct();
      row.case_id = string(c.case_id);
      row.tier = string(c.tier);
      row.baseline_type = baseline_type;
      row.baseline_tag = baseline_tag;
      row.smbmodel = string(c.smbmodel);
      row.sitename = string(c.sitename);
      row.forcings = string(c.forcings);
      row.simyear = c.simyear;
      row.solver = c.solver;
      row = copyMetricFields(row, S);
      row.last_updated_utc = datetime('now', 'TimeZone', 'UTC');
      row_cells{icase} = row;

      case_opts(icase).case_id = string(c.case_id);
      case_opts(icase).case = table2struct(c);
      case_opts(icase).opts = opts_run;
   end

   % Convert the accepted case rows into the saved baseline table.
   RegressionBaseline = struct2table(vertcat(row_cells{:}));

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
      "formal regression runs include the canonical leading spinup year";
   meta.reset_fields = "solver";
   meta.include_profile_artifacts = include_profile_artifacts;
   meta.profile_history_size = profile_history_size;
   meta.input_path = string(input_path);
   meta.output_path = string(output_path);
   meta.matlab_version = string(version);
   meta.host = string(computer);
   meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');

   % Rolling baselines are acceptance targets. Archive the prior managed
   % state before overwriting it so older accepted outputs remain available.
   if baseline_type == "rolling"
      icemodel.test.helpers.archiveManagedBaseline(output_file, "regression");
   end

   % Save profiler diagnostics in a separate rerun so the accepted baseline
   % state above stays independent from the profiling instrumentation.
   profile_summary = table();
   profile_meta = struct();
   profile_artifacts = struct();
   if include_profile_artifacts
      [profile_summary, profile_meta, profile_artifacts] = ...
         icemodel.test.helpers.captureBaselineProfile( ...
         "regression", cases, output_file, history_size=profile_history_size);
   end

   % Save the rolling or release regression baseline file.
   outdir = fileparts(char(output_file));
   if exist(outdir, 'dir') ~= 7
      mkdir(outdir);
   end
   save(char(output_file), 'RegressionBaseline', 'case_opts', 'meta', ...
      'profile_summary', 'profile_meta', 'profile_artifacts');
end

function row = copyMetricFields(row, S)
   %COPYMETRICFIELDS Copy the scalar metric struct into one baseline row.

   names = string(fieldnames(S));
   for i = 1:numel(names)
      name = char(names(i));
      row.(name) = S.(name);
   end
end
