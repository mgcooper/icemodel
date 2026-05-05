function results = run_snow_verification_suite(kwargs)
   %RUN_SNOW_VERIFICATION_SUITE Run the broader snow-verification smoke lane.
   %
   %  results = run_snow_verification_suite()
   %  results = run_snow_verification_suite(tier="full", make_plots=true)
   %  results = run_snow_verification_suite(plot_visible="on")
   %  results = run_snow_verification_suite(run_icemodel=true)
   %  results = run_snow_verification_suite(candidate_provider=@myProvider)
   %  results = run_snow_verification_suite(cases=["cdp", "wfj"])
   %  results = run_snow_verification_suite(write_artifacts=true)
   %
   % This runner is the verification-suite entry point used by agents and
   % interactive development. It reads staged data and writes artifacts; it
   % does not import or refresh setup data.
   %
   % Artifact policy (opt-in)
   % ------------------------
   % Defaults are tuned for interactive development: no figures created, no
   % artifacts written. The runner returns a result struct (summary table +
   % per-case results) so agents and developers can inspect metrics without
   % producing on-disk side effects.
   %
   %   make_plots       = false   create comparison/scatter figures?
   %   save_plots       = false   export those figures as PNG?
   %   write_artifacts  = false   write summary.csv / summary.mat / report.md
   %                              and per-case metrics.csv / result.mat?
   %   plot_visible     = "off"   figure visibility ("on" implies make_plots=true)
   %
   % Pass write_artifacts=true (or any of the plotting flags) to opt in to
   % the persisted-snapshot workflow. The runner only creates the
   % <test>/artifacts/snow-verification/<run_name>/ directory when at least
   % one artifact would be written.

   arguments
      kwargs.tier (1, 1) string ...
         {icemodel.verification.validators.mustBeTierName} ...
         = "smoke"
      kwargs.cases (1, :) string ...
         {icemodel.verification.validators.mustBeCaseIdSubset} = strings(0, 1)
      kwargs.run_name (1, 1) string = ""
      kwargs.make_plots (1, 1) logical = false
      kwargs.save_plots (1, 1) logical = false
      kwargs.plot_visible (1, 1) string = "off"
      kwargs.write_artifacts (1, 1) logical = false
      kwargs.run_icemodel (1, 1) logical = false
      kwargs.candidate_provider = []
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.artifact_root (1, 1) string = ""
      kwargs.startdate = NaT('TimeZone', 'UTC')
      kwargs.enddate   = NaT('TimeZone', 'UTC')
   end

   if kwargs.run_icemodel && ~isempty(kwargs.candidate_provider)
      error('run_icemodel and candidate_provider are mutually exclusive')
   end

   % Visible plots imply we want to create the figure.
   if kwargs.plot_visible ~= "off" && ~kwargs.make_plots
      kwargs.make_plots = true;
   end

   % Persisted figures need both make_plots and an artifact directory; treat
   % save_plots=true as opting in to artifact writing too.
   if kwargs.save_plots
      kwargs.make_plots = true;
      kwargs.write_artifacts = true;
   end

   % Print a clear banner when the synthetic-candidate path is engaged
   % so the metrics that follow cannot be misread as real model output.
   % The Colbeck verification path is excluded from this banner because
   % it produces real numerical / analytical IceModel candidates via
   % icemodel.column.infiltration and does not route through the
   % synthetic-snow hook. Retirement of the synthetic hook is tracked
   % under icemodel-tk6.7.
   if kwargs.run_icemodel
      fprintf('\n');
      fprintf('=== SYNTHETIC-CANDIDATE PATH ENGAGED ===\n');
      fprintf('  run_icemodel=true routes ESM-SnowMIP cases (cdp, wfj)\n');
      fprintf('  through icemodel.verification.syntheticSnowModelRun,\n');
      fprintf('  which perturbs the staged targets by:\n');
      fprintf('    snow_depth_offset_m   = +0.02 (m)\n');
      fprintf('    swe_scale             = x 1.05\n');
      fprintf('    surface_temp_offset_C = +0.25 (K)\n');
      fprintf('    liquid_water_scale    = x 1.05\n');
      fprintf('  The candidate metrics that follow are NOT real model\n');
      fprintf('  output until production snow physics retires this hook.\n');
      fprintf('  Colbeck cases use real IceModel candidates and are\n');
      fprintf('  unaffected by these perturbations.\n');
      fprintf('=========================================\n\n');
   end

   % Install the same test/demo config used by unit tests so fresh-clone runs
   % resolve the committed demo/data verification assets.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Resolve the run-artifact directory only when at least one artifact will
   % be written. Otherwise leave it empty so callers using only interactive
   % figures or returned data do not accumulate empty <run_name> folders.
   write_any_artifacts = kwargs.write_artifacts ...
      || (kwargs.make_plots && kwargs.save_plots);

   if write_any_artifacts
      if isblanktext(kwargs.artifact_root)
         artifact_root = fullfile(string(icemodel.getpath('test')), ...
            "artifacts", "snow-verification");
      else
         artifact_root = kwargs.artifact_root;
      end
      icemodel.helpers.ensureDirExists(artifact_root);

      [~, ~, run_name] = icemodel.test.helpers.resolveRunStamp(kwargs.run_name);
      run_dir = fullfile(artifact_root, run_name);
      icemodel.helpers.ensureDirExists(run_dir);
   else
      [~, ~, run_name] = icemodel.test.helpers.resolveRunStamp(kwargs.run_name);
      run_dir = "";
   end

   % Select cases through the public catalog entry point so the runner stays
   % dataset-family agnostic.
   cases = icemodel.verification.listcases( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "tier", selectTier(kwargs.tier));
   if isempty(cases)
      error('no snow-verification cases available for tier %s', kwargs.tier)
   end

   % Apply the optional explicit case filter after tier selection.
   if ~isempty(kwargs.cases)
      keep = ismember([cases.case_id], kwargs.cases);
      cases = cases(keep);
   end
   if isempty(cases)
      error('requested snow-verification cases were not found')
   end

   % Compare each case, collect per-case artifacts, and concatenate metric rows
   % into one summary table.
   case_results = cell(numel(cases), 1);
   summary_rows = cell(numel(cases), 1);
   for i = 1:numel(cases)

      % Agents can pass a candidate provider that receives the resolved case
      % manifest row and returns a candidate bundle with the same format as the
      % staged reference. With no provider, comparecase uses the smoke reference.
      compare_args = { ...
         "evaluation_data_root", kwargs.evaluation_data_root, ...
         "artifact_dir", run_dir, ...
         "make_plot", kwargs.make_plots, ...
         "save_plot", kwargs.save_plots, ...
         "plot_visible", kwargs.plot_visible};

      if kwargs.run_icemodel
         candidate = icemodel.verification.runIcemodelSnowCandidate(cases(i), ...
            startdate=kwargs.startdate, enddate=kwargs.enddate);
         case_result = icemodel.verification.comparecase( ...
            cases(i).case_id, compare_args{:}, "candidate", candidate);
      elseif isempty(kwargs.candidate_provider)
         case_result = icemodel.verification.comparecase( ...
            cases(i).case_id, compare_args{:});
      else
         if ~isa(kwargs.candidate_provider, 'function_handle')
            error('candidate_provider must be a function handle')
         end
         candidate = kwargs.candidate_provider(cases(i));
         case_result = icemodel.verification.comparecase( ...
            cases(i).case_id, compare_args{:}, "candidate", candidate);
      end

      case_results{i} = case_result;
      metrics = case_result.metrics;
      metrics.case_id = repmat(cases(i).case_id, height(metrics), 1);
      summary_rows{i} = movevars(metrics, 'case_id', 'Before', 1);
   end

   summary = vertcat(summary_rows{:});

   % Persist the human-readable CSV, the richer MATLAB result bundle, and
   % the run report when the caller opted in to artifacts.
   report_path = "";
   if write_any_artifacts
      writetable(summary, fullfile(run_dir, 'summary.csv'));
      save(fullfile(run_dir, 'summary.mat'), 'summary', 'case_results');
      report_path = icemodel.verification.helpers.writeRunReport( ...
         run_dir, case_results, cases, ...
         run_name=string(run_name), ...
         run_icemodel=kwargs.run_icemodel, ...
         plotted=(kwargs.make_plots && kwargs.save_plots));
   end

   results = struct( ...
      'run_name', run_name, ...
      'artifact_dir', run_dir, ...
      'summary', summary, ...
      'case_results', {case_results}, ...
      'cases', {cases}, ...
      'report_path', report_path);
   clear cleanup
end

function tier = selectTier(requested)
   %SELECTTIER Map the runner tier selector onto the current case catalog.

   % The "full" tier currently means no tier filter. Future full-data imports
   % can add explicit full cases without changing the runner call contract.
   switch requested
      case "smoke"
         tier = "smoke";
      case "full"
         tier = "";
      otherwise
         error('unsupported snow-verification tier: %s', requested)
   end
end
