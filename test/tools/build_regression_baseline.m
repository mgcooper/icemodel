function RegressionBaseline = build_regression_baseline(kwargs)
   %BUILD_REGRESSION_BASELINE Build rolling or versioned icemodel regression
   %baselines.
   %
   %  RegressionBaseline = build_regression_baseline(baseline="rolling")
   %  RegressionBaseline = build_regression_baseline(baseline_tag="v1.1")
   %  RegressionBaseline = build_regression_baseline(baseline="v1.01", ...
   %     tier="full", smbmodel="skinmodel")
   %
   % Use this when you want to accept new modeled outputs as a rolling or
   % versioned regression baseline. This writes baseline files only; it does
   % not produce compare artifacts or evaluate pass/fail against an older
   % baseline.

   arguments (Input)
      kwargs.baseline string = string.empty()
      kwargs.baseline_tag string = string.empty()
      kwargs.tier (1, :) string {mustBeMember(kwargs.tier, ...
         ["smoke", "full", "all"])} = "full"
      kwargs.smbmodel (1, :) string {mustBeMember(kwargs.smbmodel, ...
         ["all", "icemodel", "skinmodel"])} = "all"
      kwargs.output_file string = string.empty()
   end
   [baseline, baseline_tag, tier, smbmodel, output_file] = deal( ...
      kwargs.baseline, kwargs.baseline_tag, kwargs.tier, ...
      kwargs.smbmodel, kwargs.output_file);

   % Expand smbmodel="all" into one rolling/release baseline file per model.
   if smbmodel == "all"
      if ~isblanktext(output_file)
         error('output_file is only supported for a single smbmodel')
      end
      models = test.helpers.formalSmbmodels();
      baselines = cell(numel(models), 1);
      for i = 1:numel(models)
         baselines{i} = build_regression_baseline( ...
            baseline=baseline, baseline_tag=baseline_tag, tier=tier, ...
            smbmodel=models(i));
      end
      RegressionBaseline = vertcat(baselines{:});
      return
   end

   % Resolve the baseline target, configure paths, and load formal cases.
   [baseline_type, baseline_tag, output_file, rootdir, input_path, ...
      output_path, cases] = test.helpers.prepareBaselineBuild( ...
      "regression", baseline, baseline_tag, tier, smbmodel, output_file);

   % Load the static runoff reference used to derive catchment totals.
   runoff_ref = test.helpers.loadRunoffReference();

   % Preallocate row containers for the accepted baseline summary and opts.
   rows = struct([]);
   case_opts = struct([]);
   k = 0;

   % Run each formal case and save the accepted scalar regression state.
   for icase = 1:height(cases)
      c = cases(icase, :);
      opts_run = test.helpers.buildFormalCaseOpts(c);

      [ice1, ice2] = runModel(opts_run);
      [ice1, ice2] = POSTPROC(ice1, ice2, opts_run, c.simyear); %#ok<ASGLU>
      S = test.helpers.summarizeIce1Metrics(ice1);
      ridx = test.helpers.findRunoffReferenceRow(runoff_ref, c);

      icemodel_final_m3 = nan;
      if ~isempty(ridx)
         icemodel_final_m3 = test.helpers.computeCatchmentRunoffFinal( ...
            ice1, runoff_ref.area_med_m2(ridx), runoff_ref.t1(ridx), ...
            runoff_ref.t2(ridx));
      end

      k = k + 1;
      rows(k).case_id = string(c.case_id);
      rows(k).tier = string(c.tier);
      rows(k).baseline_type = baseline_type;
      rows(k).baseline_tag = baseline_tag;
      rows(k).smbmodel = string(c.smbmodel);
      rows(k).sitename = string(c.sitename);
      rows(k).forcings = string(c.forcings);
      rows(k).simyear = c.simyear;
      rows(k).solver_mode = c.solver_mode;
      rows(k).runoff_final = S.runoff_final;
      rows(k).melt_final = S.melt_final;
      rows(k).mean_Tice_numiter = S.mean_Tice_numiter;
      rows(k).max_Tice_numiter = S.max_Tice_numiter;
      rows(k).n_not_converged = S.n_not_converged;
      rows(k).icemodel_final_m3 = icemodel_final_m3;
      rows(k).last_updated_utc = datetime('now', 'TimeZone', 'UTC');

      case_opts(k).case_id = string(c.case_id);
      case_opts(k).case = table2struct(c);
      case_opts(k).opts = opts_run;
   end

   % Convert the accepted case rows into the saved baseline table.
   RegressionBaseline = struct2table(rows);

   % Record the build metadata alongside the accepted baseline values.
   meta = struct();
   meta.tier = tier;
   meta.smbmodel_filter = smbmodel;
   meta.baseline_type = baseline_type;
   meta.baseline_tag = baseline_tag;
   meta.case_builder = "test.helpers.buildFormalCaseOpts";
   meta.opts_source = "icemodel.setopts defaults";
   meta.reset_fields = "solver_mode -> bc_type";
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

function varargout = runModel(opts)
   % Dispatch to the requested core model kernel.
   switch opts.smbmodel
      case 'icemodel'
         [varargout{1:nargout}] = icemodel(opts);
      case 'skinmodel'
         [varargout{1:nargout}] = skinmodel(opts);
      otherwise
         error('unsupported smbmodel: %s', opts.smbmodel)
   end
end
