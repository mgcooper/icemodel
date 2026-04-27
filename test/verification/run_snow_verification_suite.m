function results = run_snow_verification_suite(kwargs)
   %RUN_SNOW_VERIFICATION_SUITE Run the broader snow-verification smoke lane.
   %
   %  results = run_snow_verification_suite()
   %  results = run_snow_verification_suite(tier="full", make_plots=false)
   %  results = run_snow_verification_suite(plot_visible="on")
   %  results = run_snow_verification_suite(candidate_provider=@runSnowCase)
   %  results = run_snow_verification_suite(cases=["cdp", "wfj"])
   %
   % This runner is the verification-suite entry point used by agents and
   % interactive development. It reads staged data and writes artifacts; it does
   % not import or refresh setup data.

   arguments
      kwargs.tier (1, 1) string ...
         {icemodel.verification.validators.mustBeTierName} ...
         = "smoke"
      kwargs.cases (1, :) string ...
         {icemodel.verification.validators.mustBeCaseIdSubset} = strings(0, 1)
      kwargs.run_name (1, 1) string = ""
      kwargs.make_plots (1, 1) logical = true
      kwargs.save_plots (1, 1) logical = true
      kwargs.plot_visible (1, 1) string = "off"
      kwargs.candidate_provider = []
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.artifact_root (1, 1) string = ""
   end

   % Install the same test/demo config used by unit tests so fresh-clone runs
   % resolve the committed demo/data verification assets.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Resolve and create the run artifact folder before any cases execute.
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
      if isempty(kwargs.candidate_provider)
         case_result = icemodel.verification.comparecase( ...
            cases(i).case_id, ...
            "evaluation_data_root", kwargs.evaluation_data_root, ...
            "artifact_dir", run_dir, ...
            "make_plot", kwargs.make_plots, ...
            "save_plot", kwargs.save_plots, ...
            "plot_visible", kwargs.plot_visible);
      else
         if ~isa(kwargs.candidate_provider, 'function_handle')
            error('candidate_provider must be a function handle')
         end
         candidate = kwargs.candidate_provider(cases(i));
         case_result = icemodel.verification.comparecase( ...
            cases(i).case_id, ...
            "evaluation_data_root", kwargs.evaluation_data_root, ...
            "artifact_dir", run_dir, ...
            "make_plot", kwargs.make_plots, ...
            "save_plot", kwargs.save_plots, ...
            "plot_visible", kwargs.plot_visible, ...
            "candidate", candidate);
      end
 
      case_results{i} = case_result;
      metrics = case_result.metrics;
      metrics.case_id = repmat(cases(i).case_id, height(metrics), 1);
      summary_rows{i} = movevars(metrics, 'case_id', 'Before', 1);
   end

   summary = vertcat(summary_rows{:});

   % Persist both a human-readable CSV and the richer MATLAB result bundle.
   writetable(summary, fullfile(run_dir, 'summary.csv'));
   save(fullfile(run_dir, 'summary.mat'), 'summary', 'case_results');
   results = struct( ...
      'run_name', run_name, ...
      'artifact_dir', run_dir, ...
      'summary', summary, ...
      'case_results', {case_results}, ...
      'cases', {cases});
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
