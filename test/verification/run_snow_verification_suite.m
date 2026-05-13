function results = run_snow_verification_suite(kwargs)
   %RUN_SNOW_VERIFICATION_SUITE Run the snow-verification suite.
   %
   %  results = run_snow_verification_suite()
   %  results = run_snow_verification_suite(cases=["cdp", "wfj"])
   %  results = run_snow_verification_suite(make_plots=true)
   %  results = run_snow_verification_suite(plot_visible="on")
   %  results = run_snow_verification_suite(run_icemodel=true)
   %  results = run_snow_verification_suite(candidate_provider=@myProvider)
   %  results = run_snow_verification_suite(write_artifacts=true)
   %  results = run_snow_verification_suite( ...
   %     cases="wfj", startdate="1997-10-01", enddate="1998-09-30")
   %  results = run_snow_verification_suite( ...
   %     cases=icemodel.verification.namelists.snowmipsite, ...
   %     run_icemodel=true, plot_visible="on", save_plots=false)
   %
   % This runner is the verification-suite entry point used by agents and
   % interactive development. It reads staged data and writes artifacts; it
   % does not import or refresh setup data.
   %
   % Default behaviour
   %   With no arguments, runs Col de Porte (cdp) over its
   %   default_smoke_window. CDP is the most canonical / widely-cited
   %   ESM-SnowMIP snow verification site (Menard 2019 ESSD), and a single
   %   site / single year keeps interactive runs fast. Override via cases=
   %   to select other sites and / or startdate / enddate to narrow the
   %   comparison window.
   %
   % Artifact policy (opt-in)
   %   Defaults are tuned for interactive development: no figures created,
   %   no artifacts written. The runner returns a result struct (summary
   %   table + per-case results) so agents and developers can inspect
   %   metrics without producing on-disk side effects.
   %
   %     make_plots       = false   create comparison/scatter figures?
   %     save_plots       = false   export those figures as PNG?
   %     write_artifacts  = false   write summary.csv / summary.mat /
   %                                report.md and per-case
   %                                metrics.csv / result.mat?
   %     plot_visible     = "off"   figure visibility
   %                                ("on" implies make_plots=true)
   %
   %   Pass write_artifacts=true (or any of the plotting flags) to opt
   %   in to the persisted-snapshot workflow. The runner only creates the
   %   <test>/artifacts/snow-verification/<run_name>/ directory when at
   %   least one artifact would be written.

   arguments
      kwargs.cases (1, :) string ...
         {icemodel.verification.validators.mustBeCaseIdSubset} = "cdp"
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
      kwargs.enddate = NaT('TimeZone', 'UTC')
   end

   [run_name, run_dir, write_any_artifacts, kwargs, cleanup] = ...
      resolveRunContext(kwargs);

   % Select cases through the public catalog entry point so the runner stays
   % dataset-family agnostic.
   cases = icemodel.verification.listcases( ...
      "evaluation_data_root", kwargs.evaluation_data_root);
   if isempty(cases)
      error('no snow-verification cases available')
   end

   % Apply the explicit case filter (defaults to "cdp" when caller did not
   % override). Callers wanting "all staged cases" pass case_ids matching
   % every case in the catalog.
   keep = ismember([cases.case_id], kwargs.cases);
   cases = cases(keep);
   if isempty(cases)
      error('requested snow-verification cases were not found')
   end

   % Compare each case, collect per-case artifacts, and concatenate metric rows
   % into one summary table.
   case_results = cell(numel(cases), 1);
   summary_rows = cell(numel(cases), 1);
   for n = 1:numel(cases)

      % Static compare-arg bundle shared by every dispatch branch below.
      compare_args = { ...
         "evaluation_data_root", kwargs.evaluation_data_root, ...
         "artifact_dir", run_dir, ...
         "make_plot", kwargs.make_plots, ...
         "save_plot", kwargs.save_plots, ...
         "plot_visible", kwargs.plot_visible};

      if kwargs.run_icemodel
         % Produce a real (or synthetic-stub) IceModel candidate.
         candidate = icemodel.verification.runIcemodelSnowCandidate(cases(n), ...
            startdate=kwargs.startdate, enddate=kwargs.enddate);

         % Compare the IceModel candidate against the staged reference for this
         % case.
         case_result = icemodel.verification.comparecase( ...
            cases(n).case_id, compare_args{:}, "candidate", candidate);

      elseif isempty(kwargs.candidate_provider)
         % No candidate supplied: comparecase reuses the staged smoke reference
         % as the candidate.
         case_result = icemodel.verification.comparecase( ...
            cases(n).case_id, compare_args{:});
      else
         if ~isa(kwargs.candidate_provider, 'function_handle')
            error('candidate_provider must be a function handle')
         end
         % Caller-provided provider returns a candidate bundle for the resolved
         % case manifest row.
         candidate = kwargs.candidate_provider(cases(n));

         % Compare the provider-built candidate against the staged reference for
         % this case.
         case_result = icemodel.verification.comparecase( ...
            cases(n).case_id, compare_args{:}, "candidate", candidate);
      end

      case_results{n} = case_result;
      metrics = case_result.metrics;
      metrics.case_id = repmat(cases(n).case_id, height(metrics), 1);
      summary_rows{n} = movevars(metrics, 'case_id', 'Before', 1);
   end

   summary = vertcat(summary_rows{:});

   % Persist the human-readable CSV, the richer MATLAB result bundle, and
   % the run report when the caller opted in to artifacts.
   report_path = "";
   if write_any_artifacts

      % Write the flat per-case metrics CSV alongside the MATLAB result bundle.
      writetable(summary, fullfile(run_dir, 'summary.csv'));

      % Bundle summary + per-case structs into a single .mat snapshot.
      save(fullfile(run_dir, 'summary.mat'), 'summary', 'case_results');

      % Render the human-readable Markdown run report into run_dir.
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

%%
function [run_name, run_dir, write_any_artifacts, kwargs, cleanup] = ...
      resolveRunContext(kwargs)
   %RESOLVERUNCONTEXT Resolve comparison window, flags, banner, env, and run_dir.
   %
   % Do the pre-loop setup work: validate flag combinations, default the
   % comparison window for single-site ESM-SnowMIP runs, set implied flags
   % (visible plots, save_plots), print the synthetic-candidate banner when
   % applicable, bootstrap the test environment, and resolve the run-artifact
   % directory.

   if kwargs.run_icemodel && ~isempty(kwargs.candidate_provider)
      error('run_icemodel and candidate_provider are mutually exclusive')
   end

   % Resolve the runtime comparison window. With no explicit dates and a
   % single ESM-SnowMIP site, default to that site's default_smoke_window
   % so interactive runs match the staged smoke fixture without the
   % caller having to know per-site dates. With multiple cases or
   % non-ESM-SnowMIP cases (e.g. colbeck1976), let comparecase use the
   % staged window directly.
   if isnat(kwargs.startdate) && isnat(kwargs.enddate) ...
         && isscalar(kwargs.cases) ...
         && ismember(kwargs.cases, ...
         icemodel.verification.namelists.snowmipsite())
      [kwargs.startdate, kwargs.enddate] = ...
         icemodel.verification.helpers.default_smoke_window(kwargs.cases);
   end

   % Visible plots imply we want to create the figure.
   if kwargs.plot_visible ~= "off" && ~kwargs.make_plots
      kwargs.make_plots = true;
   end

   % Figures need both make_plots and an artifact directory; treat
   % save_plots=true as opting in to artifact writing too.
   if kwargs.save_plots
      kwargs.make_plots = true;
      kwargs.write_artifacts = true;
   end

   % Print a banner when the synthetic-candidate path is chosen so the metrics
   % that follow are not misread as real model output. The Colbeck verification
   % path is excluded from this banner because it produces real numerical /
   % analytical IceModel candidates via icemodel.column.infiltration and does
   % not route through the synthetic-snow hook. Retirement of the synthetic hook
   % is tracked under icemodel-tk6.7.
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
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();

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
end
