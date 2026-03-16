function RegressionBaseline = build_regression_baseline(kwargs)
   %BUILD_REGRESSION_BASELINE Build rolling or versioned icemodel regression
   %baselines.
   %
   %  RegressionBaseline = build_regression_baseline(baseline="rolling")
   %  RegressionBaseline = build_regression_baseline(baseline_tag="v1.1")
   %  RegressionBaseline = build_regression_baseline(baseline="v1.01", ...
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
   % baseline.
   % The optional solver filter accepts any subset of [1 2 3].
   % The formal comparison year and smoke/full site selections are explicit
   % here rather than buried in the regression case matrix helper.

   arguments (Input)
      kwargs.baseline string = string.empty()
      kwargs.baseline_tag string = string.empty()
      kwargs.tier (1, :) string {mustBeMember(kwargs.tier, ...
         ["smoke", "full", "all"])} = "full"
      kwargs.smbmodel (1, :) string {mustBeMember(kwargs.smbmodel, ...
         ["all", "icemodel", "skinmodel"])} = "all"
      kwargs.solver {icemodel.validators.mustBeSolverFilter(kwargs.solver)} = []
      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      kwargs.smoke_sites string = "kanm"
      kwargs.full_sites string = ["kanm"; "kanl"]
      kwargs.output_file string = string.empty()
   end
   [baseline, baseline_tag, tier, smbmodel, solver, simyear, smoke_sites, ...
      full_sites, output_file] = deal( ...
      kwargs.baseline, kwargs.baseline_tag, kwargs.tier, ...
      kwargs.smbmodel, kwargs.solver, kwargs.simyear, ...
      string(kwargs.smoke_sites(:)), string(kwargs.full_sites(:)), ...
      kwargs.output_file);

   % Expand smbmodel="all" into one rolling/release baseline file per model.
   if smbmodel == "all"
      if ~isblanktext(output_file)
         error('output_file is only supported for a single smbmodel')
      end
      models = icemodel.test.helpers.formalSmbmodels();
      baselines = cell(numel(models), 1);
      for i = 1:numel(models)
         baselines{i} = build_regression_baseline( ...
            baseline=baseline, baseline_tag=baseline_tag, tier=tier, ...
            smbmodel=models(i), solver=solver, simyear=simyear, ...
            smoke_sites=smoke_sites, full_sites=full_sites);
      end
      RegressionBaseline = vertcat(baselines{:});
      return
   end

   % Resolve the baseline target, configure paths, and load formal cases.
   [baseline_type, baseline_tag, output_file, ~, input_path, ...
      output_path, cases] = icemodel.test.helpers.prepareBaselineBuild( ...
      "regression", baseline, baseline_tag, tier, smbmodel, output_file, ...
      simyear, solver, smoke_sites, full_sites);

   % Load the static runoff reference used to derive catchment totals.
   runoff_ref = icemodel.test.helpers.loadRunoffReference();

   % Preallocate row containers for the accepted baseline summary and opts.
   row_cells = cell(height(cases), 1);
   case_opts = struct([]);

   % Run each formal case and save the accepted scalar regression state.
   for icase = 1:height(cases)
      c = cases(icase, :);
      opts_run = icemodel.test.helpers.setModelOptsForCase(c);

      [ice1, ice2] = icemodel.test.helpers.runSmbModel(opts_run);
      [ice1, ~] = icemodel.postprocess( ...
         ice1, ice2, opts_run, opts_run.output_years);
      ridx = icemodel.test.helpers.findRunoffReferenceRow(runoff_ref, c);
      met = icemodel.test.helpers.loadProcessedMetForOutputYears(opts_run);
      if isempty(ridx)
         refrow = [];
      else
         refrow = runoff_ref(ridx, :);
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
   meta.reset_fields = "solver";
   meta.input_path = string(input_path);
   meta.output_path = string(output_path);
   meta.matlab_version = string(version);
   meta.host = string(computer);
   meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');

   % Save the rolling or release regression baseline file.
   outdir = fileparts(char(output_file));
   if exist(outdir, 'dir') ~= 7
      mkdir(outdir);
   end
   save(char(output_file), 'RegressionBaseline', 'case_opts', 'meta');
end

function row = copyMetricFields(row, S)

   names = string(fieldnames(S));
   for i = 1:numel(names)
      name = char(names(i));
      row.(name) = S.(name);
   end
end
