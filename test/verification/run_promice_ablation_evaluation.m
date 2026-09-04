function results = run_promice_ablation_evaluation(kwargs)
   %RUN_PROMICE_ABLATION_EVALUATION Evaluate PROMICE ablation case-years.
   %
   %  results = run_promice_ablation_evaluation()
   %  results = run_promice_ablation_evaluation( ...
   %     case_ids=["kanm", "kanl"], years=2019)
   %  results = run_promice_ablation_evaluation( ...
   %     case_ids="kanm", write_artifacts=true)
   %
   % The readiness ledger is the only source of forcing identity and window
   % admission. An empty CASE_IDS selection performs readiness only; callers
   % must explicitly select case ids (or "all") before any model is executed.
   % Production execution uses runModelCase with the producer-pinned
   % promice_filled artifact. MODEL_PROVIDER is a focused-test seam whose
   % output must match the postprocessed diagnostic timetable. Each
   % run initializes on January 1 and retains a snow-aware June--October result
   % while comparison metrics use the readiness-selected summer interval.

   arguments
      kwargs.case_ids (1, :) string = string.empty(1, 0)
      kwargs.years (1, :) double {mustBeInteger} = double.empty(1, 0)
      kwargs.data_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "verification"
      kwargs.artifact_root (1, 1) string = ""
      kwargs.run_name (1, 1) string = ""
      kwargs.write_artifacts (1, 1) logical = false
      kwargs.model_provider = []
   end

   % Readiness-only runs are valid. Writing artifacts for one would create a
   % populated run directory and a renderable report that describes zero
   % site-years, which looks like a completed evaluation. Reject the
   % combination before this code creates anything on disk.
   if kwargs.write_artifacts && isempty(kwargs.case_ids)
      error('icemodel:verification:promiceAblationEvaluation:emptySelection', ...
         ['write_artifacts=true requires a non-empty case selection; ' ...
         'pass case_ids="all" for the full cohort or name the cases'])
   end

   % Resolve one paired data tree and one run directory before readiness so
   % every path in the result has the same owner.
   [evaluation_data_root, input_data_root] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      data_root=kwargs.data_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      input_data_root=kwargs.input_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
   [run_name, run_dir, cleanup] = resolveRunDirectory(kwargs);
   readiness = ...
      icemodel.verification.setup.writePromiceAblationReadiness( ...
      evaluation_data_root=evaluation_data_root, ...
      input_data_root=input_data_root, output_dir=run_dir);
   policy = icemodel.verification.namelists.promiceAblationPolicy();

   % Validate the optional test seam once rather than turning every selected
   % row into an identical unavailable result.
   if ~isempty(kwargs.model_provider) ...
         && ~isa(kwargs.model_provider, 'function_handle')
      error('icemodel:verification:promiceAblationEvaluation:modelProvider', ...
         'model_provider must be a function handle')
   end

   % Explicit selection is an execution safety boundary. All ledger rows are
   % retained even when they are excluded or were not selected for execution.
   requested = requestedRows(readiness.rows, kwargs.case_ids, kwargs.years);

   % Resolve the physics stamp before the loop, not at the point it is saved.
   % Resolving it runs icemodel.setopts, which asserts that the workspace
   % exists. A workspace that disappears during a multi-hour cohort must not
   % throw between the last completed site-year and the save that keeps it.
   %
   % A run that executes no model stamps nothing. It writes no artifacts and
   % must not require a workspace it never reads. Gate on a requested row
   % that is also admitted: selecting only readiness-excluded rows requests
   % work the loop below skips before it runs any case.
   %
   % Guard the call. physicsFingerprint resolves its reference through
   % icemodel.setopts, which asserts ICEMODEL_INPUT_PATH exists, and its own
   % docstring says a caller that must not fail on a missing workspace has to
   % guard it. A cohort that supplies input_data_root and
   % evaluation_data_root reaches its forcing through the manifest root that
   % setModelOptsForCase injects. The global path is therefore irrelevant to
   % whether those cases can run. Without this guard such a run dies here before it
   % executes anything. An unstamped run is already a supported state: the
   % report reads a missing stamp as absent and says so without warning.
   % Unavailable is the other state, for a stamped cohort whose current
   % defaults will not resolve.
   physics_fingerprint = struct.empty();
   if any(requested(:) & logical(readiness.rows.admitted(:)))
      try
         physics_fingerprint = ...
            icemodel.verification.helpers.physicsFingerprint();
      catch
         physics_fingerprint = struct.empty();
      end
   end
   n_rows = height(readiness.rows);
   site_year_results = repmat(emptySiteYearResult(), n_rows, 1);
   summary_rows = repmat(emptySummaryRecord(), n_rows, 1);
   nested_by_row = repmat({emptyNestedTable()}, n_rows, 1);
   perturbations_by_row = repmat( ...
      {emptyEndpointPerturbationTable()}, n_rows, 1);
   for k = 1:n_rows
      row = readiness.rows(k, :);
      site_year_results(k) = emptySiteYearResult(row);
      summary_rows(k) = emptySummaryRecord( ...
         row, requested(k), readiness.policy);

      % Admission failures take precedence over selection because they are the
      % scientific reason the case-year cannot be evaluated.
      if ~row.admitted
         site_year_results(k).status = "excluded";
         site_year_results(k).reason = row.exclusion_reason;
         summary_rows(k).status = "excluded";
         summary_rows(k).reason = row.exclusion_reason;
         continue
      end
      if ~requested(k)
         site_year_results(k).status = "not_selected";
         site_year_results(k).reason = "case-year was not explicitly selected";
         summary_rows(k).status = "not_selected";
         summary_rows(k).reason = "case-year was not explicitly selected";
         continue
      end

      % Keep a cohort run recoverable: one unavailable case-year is recorded
      % without discarding completed or excluded rows from the same ledger.
      nested = predeclaredNestedDiagnostics(row, policy);
      perturbations = predeclaredEndpointPerturbations(row, policy);
      site_year_results(k).model_metadata = ...
         plannedRunProvenance(row, policy);
      try
         [site_year_results(k), nested, perturbations] = evaluateRow( ...
            row, evaluation_data_root, input_data_root, ...
            kwargs.model_provider, policy);
         summary_rows(k) = completedSummaryRecord( ...
            row, site_year_results(k), true, readiness.policy);
      catch err
         [failure_reason, error_identifier] = failureDetails(err);
         site_year_results(k).status = "unavailable";
         site_year_results(k).reason = failure_reason;
         site_year_results(k).error_identifier = error_identifier;
         nested.reason(:) = failure_reason;
         nested.error_identifier(:) = error_identifier;
         site_year_results(k).nested_windows = nested;
         perturbations.reason(:) = failure_reason;
         perturbations.error_identifier(:) = error_identifier;
         site_year_results(k).endpoint_perturbations = perturbations;
         summary_rows(k).status = "unavailable";
         summary_rows(k).reason = failure_reason;
         summary_rows(k).error_identifier = error_identifier;
      end
      nested_by_row{k} = nested;
      perturbations_by_row{k} = perturbations;
   end
   nested_windows = emptyNestedTable();
   endpoint_perturbations = emptyEndpointPerturbationTable();
   if n_rows > 0
      nested_windows = vertcat(nested_by_row{:});
      endpoint_perturbations = vertcat(perturbations_by_row{:});
   end
   summary = struct2table(summary_rows);

   % Persist only caller-requested evaluation outputs. The readiness writer is
   % still exercised for transient runs, then its temporary files are removed.
   paths = resultPaths(run_dir, readiness.files, kwargs.write_artifacts);
   results = struct( ...
      'run_name', run_name, ...
      'readiness', readiness, ...
      'site_year_results', site_year_results, ...
      'summary', summary, ...
      'nested_windows', nested_windows, ...
      'endpoint_perturbations', endpoint_perturbations, ...
      'policy', policy, ...
      'paths', paths);

   % Save the digest of the model's default option values and the CITATION.cff
   % version, resolved before the run started. The report recomputes it and
   % warns when it moved.
   %
   % The digest covers a fixed reference resolution, not this cohort's
   % per-site options. It also does not move for a source change that touches
   % no option and no release. A match is not proof of reproducibility.
   %
   % Assign the field separately, so a readiness-only run leaves it out. A
   % present but empty stamp would read as a damaged one at report time.
   if ~isempty(physics_fingerprint)
      results.physics_fingerprint = physics_fingerprint;
   end
   if kwargs.write_artifacts
      writetable(summary, paths.summary_csv)
      writetable(nested_windows, paths.nested_windows_csv)
      writetable(endpoint_perturbations, paths.endpoint_perturbations_csv)
      % Summary and readiness retain the full cohort; only executed cases carry
      % per-site payloads in the MAT artifact to avoid duplicating empty tables.
      artifact = struct('results', results);
      status = string({artifact.results.site_year_results.status});
      keep = ismember(status, ["completed", "unavailable"]);
      artifact.results.site_year_results = ...
         artifact.results.site_year_results(keep);
      save(paths.results_mat, '-struct', 'artifact')
   else
      % Do not return paths to files that the transient cleanup removes.
      results.readiness.files = struct('csv', "", 'json', "");
      results.paths = emptyPaths();
   end
   delete(cleanup)
end

function [run_name, run_dir, cleanup] = resolveRunDirectory(kwargs)
   %RESOLVERUNDIRECTORY Create a persistent or transient run directory.
   [~, ~, run_name] = ...
      icemodel.test.helpers.resolveRunStamp(kwargs.run_name);
   if kwargs.write_artifacts
      if kwargs.artifact_root == ""
         artifact_root = fullfile(string(icemodel.getpath('test')), ...
            'artifacts', 'promice-ablation-evaluation');
      else
         artifact_root = kwargs.artifact_root;
      end
      run_dir = fullfile(artifact_root, run_name);
      icemodel.helpers.ensureDirExists(run_dir);
      cleanup = onCleanup(@() []);
   else
      % Readiness normally writes files, so a no-write evaluation uses
      % a private temporary run directory and removes it before returning.
      run_dir = string(tempname);
      mkdir(run_dir)
      cleanup = onCleanup(@() removeTransientDirectory(run_dir));
   end
end

function removeTransientDirectory(run_dir)
   %REMOVETRANSIENTDIRECTORY Remove only the runner-owned temporary directory.
   if isfolder(run_dir)
      rmdir(run_dir, 's')
   end
end

function requested = requestedRows(rows, case_ids, years)
   %REQUESTEDROWS Apply explicit case and optional year selection.
   requested = false(height(rows), 1);
   case_ids = lower(case_ids);

   % An omitted selection is readiness-only mode: report which site-years are
   % admissible without launching the model.
   if isempty(case_ids)
      return
   end

   % The literal "all" is an explicit full-cohort request; mixing it with case
   % ids is rejected so command provenance remains unambiguous.
   if any(case_ids == "all")
      if numel(case_ids) ~= 1
         error('icemodel:verification:promiceAblationEvaluation:caseSelection', ...
            'case_ids="all" cannot be combined with individual case ids')
      end
      requested(:) = true;
   else
      unknown = setdiff(case_ids, unique(rows.case_id), 'stable');
      if ~isempty(unknown)
         error('icemodel:verification:promiceAblationEvaluation:unknownCase', ...
            'unknown PROMICE case id(s): %s', strjoin(unknown, ', '))
      end
      requested = ismember(rows.case_id, case_ids);
   end

   % An omitted year filter retains every annual row for the selected cases.
   if ~isempty(years)
      requested = requested & ismember(rows.year, years);
   end
   % A named selection that matches nothing is an error, not an empty run.
   % The year filter is not the only way to reach an empty mask.
   if ~any(requested)
      error('icemodel:verification:promiceAblationEvaluation:emptySelection', ...
         ['the requested case/year selection matches no readiness rows; ' ...
         'check the available case-year ledger before execution'])
   end
end

function provenance = plannedRunProvenance(row, policy)
   %PLANNEDRUNPROVENANCE Separate model initialization from evaluation support.
   evaluation_end = parseTime(row.snow_free_window_end);
   [display_start, display_end] = ...
      icemodel.verification.helpers.evaluationSeason( ...
      double(row.year), policy);
   provenance = struct( ...
      'initialization_start', parseTime(row.requested_window_start), ...
      'initialization_policy', "readiness requested_window_start", ...
      'evaluation_start', parseTime(row.snow_free_window_start), ...
      'evaluation_end', evaluation_end, ...
      'display_start', display_start, ...
      'display_end', display_end, ...
      'run_end_inclusive', ...
      display_end - seconds(policy.model_substep_seconds));
end

function provenance = validateRunProvenance(row, policy)
   %VALIDATERUNPROVENANCE Enforce ordered windows and forcing coverage.
   provenance = plannedRunProvenance(row, policy);
   requested_end = parseTime(row.requested_window_end);
   forcing_start = parseTime(row.forcing_acceptance_start);
   forcing_end = parseTime(row.forcing_acceptance_end);
   times = [provenance.initialization_start, provenance.evaluation_start, ...
      provenance.evaluation_end, provenance.run_end_inclusive, ...
      provenance.display_start, provenance.display_end, requested_end, ...
      forcing_start, forcing_end];
   ordered = ~any(isnat(times)) ...
      && provenance.initialization_start <= provenance.display_start ...
      && provenance.display_start <= provenance.evaluation_start ...
      && provenance.evaluation_start < provenance.evaluation_end ...
      && provenance.evaluation_end <= provenance.display_end ...
      && provenance.run_end_inclusive >= provenance.initialization_start ...
      && provenance.run_end_inclusive < provenance.display_end ...
      && requested_end >= provenance.display_end;
   if ~ordered
      error('icemodel:verification:promiceAblationEvaluation:windowOrdering', ...
         'readiness initialization, evaluation, and run windows are not ordered')
   end

   % Admission covers the annual initialization interval, not merely the later
   % snow-free evaluation endpoints.
   if forcing_start > provenance.initialization_start ...
         || forcing_end < provenance.display_end
      error('icemodel:verification:promiceAblationEvaluation:forcingCoverage', ...
         'promice_filled does not cover initialization through October 1')
   end
end

function [result, nested, perturbations] = evaluateRow( ...
      row, evaluation_data_root, input_data_root, model_provider, policy)
   %EVALUATEROW Run and compare one admitted case-year.
   result = emptySiteYearResult(row);
   provenance = validateRunProvenance(row, policy);
   initialization_start = provenance.initialization_start;
   t0 = provenance.evaluation_start;
   t1 = provenance.evaluation_end;
   display_start = provenance.display_start;
   display_end = provenance.display_end;

   % Resolve the canonical PROMICE case, then replace its raw PROMICE leg with
   % exactly the producer-pinned promice_filled artifact from readiness.
   manifest = icemodel.verification.loadmanifest(row.case_id, ...
      evaluation_data_root=evaluation_data_root, ...
      input_data_root=input_data_root, dataset_family="promice");
   manifest = overlayFilledForcing(manifest, row);
   observations = loadObservation( ...
      evaluation_data_root, row, policy, display_start, display_end);

   % Inclusive loadmet bounds would execute the future October 1 interval. Stop
   % one 15-minute substep earlier and add a state-only October 1 checkpoint.
   run_end = provenance.run_end_inclusive;

   % Rehash every producer link and scientific input immediately before the
   % run. A readiness check does not keep the identity of mutable staged
   % files valid over time.
   execution_artifacts = revalidateExecutionArtifacts( ...
      row, evaluation_data_root, input_data_root);
   manifest.report_inputs_file = execution_artifacts.report_inputs_file;
   manifest.readiness_file = execution_artifacts.readiness_file;
   [model, opts, source] = runModel( ...
      manifest, row, initialization_start, run_end, model_provider);
   [model, boundary_metadata] = appendBoundaryCheckpoint( ...
      model, display_end, provenance, source, policy);

   % The primary comparison must retain the readiness-selected direct endpoints;
   % moving either endpoint would change the scientific window.
   [comparison, aligned, diagnostics] = ...
      icemodel.verification.compareAblation( ...
      observations, model, window_start=t0, window_end=t1);
   requireExactEndpoints(comparison, t0, t1);
   nested = nestedDiagnostics( ...
      row, observations, model, comparison, policy, t0, t1);
   perturbations = endpointPerturbationDiagnostics( ...
      row, observations, model, policy);
   seasonal = seasonalDiagnostics( ...
      observations, model, provenance, policy);

   % Retain the actual inputs and outputs needed for saved-artifact reporting.
   result.status = "completed";
   result.reason = "";
   result.manifest = manifest;
   result.observations = observations;
   % The saved policy's channel list selects the saved columns, so that list
   % is the cohort's recorded schema. The report's schema gate compares it
   % against the channels the report reads and against the current namelist.
   keep_model = model.Time >= display_start & model.Time <= display_end;
   % Save the model schema recorded in the policy.
   saved_fields = policy.required_model_fields;
   result.model = model(keep_model, cellstr(saved_fields));
   result.model_options = opts;
   result.model_metadata = boundary_metadata;
   result.comparison = comparison;
   result.aligned = aligned;
   result.seasonal = seasonal;
   result.diagnostics = diagnostics;
   result.nested_windows = nested;
   result.endpoint_perturbations = perturbations;
end

function manifest = overlayFilledForcing(manifest, row)
   %OVERLAYFILLEDFORCING Pin one resolved manifest to its readiness artifact.
   artifact = replace(string(row.forcing_artifact), "\\", "/");
   prefix = "input/met/";
   if ~startsWith(artifact, prefix) || contains(artifact, "..")
      error('icemodel:verification:promiceAblationEvaluation:forcingPath', ...
         'readiness forcing_artifact is not an input/met-relative path')
   end
   met_file = extractAfter(artifact, prefix);
   if met_file == "" || ~startsWith(met_file, "promice_filled/")
      error('icemodel:verification:promiceAblationEvaluation:forcingIdentity', ...
         'readiness forcing_artifact is not a promice_filled product')
   end

   % setModelOptsForCase resolves this exact recorded path under input/met and
   % cannot fall back to filename discovery while the record is present.
   manifest.forcing_sources = "promice_filled";
   leg = struct( ...
      'kind', 'station_met_and_eval', ...
      'staged', true, ...
      'source', 'promice_filled', ...
      'source_id', 'promice_filled', ...
      'met_files', char(met_file), ...
      'window', struct( ...
      'start', char(row.forcing_acceptance_start), ...
      'end', char(row.forcing_acceptance_end)));
   manifest.colocation.promice_filled = leg;
end

function observations = loadObservation( ...
      evaluation_data_root, row, policy, display_start, display_end)
   %LOADOBSERVATION Load the compact June-through-October observation slice.
   relative_path = replace(string(row.observation_artifact), "\\", "/");
   pathname = ...
      icemodel.verification.helpers.assertRootRelativeArtifactSha256( ...
      evaluation_data_root, relative_path, string(row.observation_sha256));
   saved = load(pathname, 'targets');
   if ~isfield(saved, 'targets') || ~isstruct(saved.targets) ...
         || ~isfield(saved.targets, 'data') ...
         || ~istimetable(saved.targets.data)
      error('icemodel:verification:promiceAblationEvaluation:observation', ...
         'observation artifact lacks targets.data timetable')
   end
   data = saved.targets.data;
   data.Properties.RowTimes = ...
      icemodel.verification.setup.ensureUtc(data.Properties.RowTimes);
   required = policy.required_observation_fields;
   missing = setdiff(required, string(data.Properties.VariableNames), 'stable');
   if ~isempty(missing)
      error('icemodel:verification:promiceAblationEvaluation:observationSchema', ...
         'observation artifact is missing required field(s): %s', ...
         strjoin(missing, ', '))
   end

   % Persist no off-season rows or unused channels. Source path and full-file
   % hash remain in readiness, while this compact bundle is sufficient for the
   % comparator and saved-artifact report.
   keep = data.Time >= display_start & data.Time <= display_end;
   data = data(keep, cellstr(required));
   metadata = struct( ...
      'selection_start', display_start, 'selection_end', display_end, ...
      'retained_fields', required);
   source_metadata = metadataFields(saved.targets);
   names = fieldnames(source_metadata);
   for k = 1:numel(names)
      metadata.(names{k}) = source_metadata.(names{k});
   end
   observations = struct( ...
      'format', "timeseries", 'data', data, 'metadata', metadata);
end

function artifacts = revalidateExecutionArtifacts( ...
      row, evaluation_data_root, input_data_root)
   %REVALIDATEEXECUTIONARTIFACTS Check every producer hash still matches.
   required = ["forcing_producer_manifest", ...
      "forcing_producer_manifest_sha256", "forcing_readiness_artifact", ...
      "forcing_readiness_sha256", "forcing_artifact", "forcing_sha256", ...
      "observation_artifact", "observation_sha256"];
   missing = setdiff(required, string(row.Properties.VariableNames), 'stable');
   if ~isempty(missing)
      error('icemodel:verification:promiceAblationEvaluation:provenanceSchema', ...
         'readiness row lacks provenance field(s): %s', strjoin(missing, ', '))
   end

   % All producer paths are selected-data-root-relative. The observation path
   % is contained the same way under the paired evaluation-data root.
   selected_data_root = string(fileparts(input_data_root));
   report_inputs_file = ...
      icemodel.verification.helpers.assertRootRelativeArtifactSha256( ...
      selected_data_root, string(row.forcing_producer_manifest), ...
      string(row.forcing_producer_manifest_sha256));
   readiness_file = ...
      icemodel.verification.helpers.assertRootRelativeArtifactSha256( ...
      selected_data_root, string(row.forcing_readiness_artifact), ...
      string(row.forcing_readiness_sha256));
   icemodel.verification.helpers.assertRootRelativeArtifactSha256( ...
      selected_data_root, string(row.forcing_artifact), ...
      string(row.forcing_sha256));
   icemodel.verification.helpers.assertRootRelativeArtifactSha256( ...
      evaluation_data_root, string(row.observation_artifact), ...
      string(row.observation_sha256));

   % Pass only the root-confined paths returned by the identity verifier into
   % loadmet; raw readiness-row strings must never become runtime file options.
   artifacts = struct( ...
      'report_inputs_file', report_inputs_file, ...
      'readiness_file', readiness_file);
end

function metadata = metadataFields(targets)
   %METADATAFIELDS Retain only small source-identity observation metadata.
   metadata = struct();
   if ~isfield(targets, 'metadata') || ~isstruct(targets.metadata)
      return
   end
   allowed = {'source', 'source_family', 'station', 'site_id'};
   for k = 1:numel(allowed)
      name = allowed{k};
      if isfield(targets.metadata, name)
         metadata.(name) = targets.metadata.(name);
      end
   end
end

function [model, opts, source] = runModel( ...
      manifest, row, run_start, run_end, model_provider)
   %RUNMODEL Dispatch production through the shared case path or a test seam.
   if isempty(model_provider)
      [model, ~, opts] = icemodel.test.helpers.runModelCase( ...
         manifest, startdate=run_start, enddate=run_end, ...
         output_profile="diagnostic");
      source = "icemodel.test.helpers.runModelCase";
   else
      model = model_provider(manifest, row, run_start, run_end);
      opts = struct();
      source = "model_provider";
   end
   if isstruct(model) && isfield(model, 'data')
      model = model.data;
   end
   if ~istimetable(model)
      error('icemodel:verification:promiceAblationEvaluation:model', ...
         'model output must be a diagnostic timetable or data bundle')
   end
end

function [model, metadata] = appendBoundaryCheckpoint( ...
      model, boundary_time, provenance, source, policy)
   %APPENDBOUNDARYCHECKPOINT Add a state-only row at the display endpoint.

   % Require the ledger channels used by the comparison and the storage
   % channels needed to add the boundary row.
   required = unique([ ...
      icemodel.verification.namelists.ablationReportChannels('ledger'), ...
      string(icemodel.namelists.budgetoutputs('first')), ...
      string(icemodel.namelists.budgetoutputs('last'))], 'stable');
   missing = setdiff(required, string(model.Properties.VariableNames), 'stable');
   if ~isempty(missing)
      error('icemodel:verification:promiceAblationEvaluation:modelSchema', ...
         'diagnostic model is missing required field(s): %s', ...
         strjoin(missing, ', '))
   end
   model.Properties.RowTimes = ...
      icemodel.verification.setup.ensureUtc(model.Properties.RowTimes);
   if isempty(model) || any(model.Time >= boundary_time)
      error('icemodel:verification:promiceAblationEvaluation:futureInterval', ...
         'model output must end before the October 1 boundary checkpoint')
   end

   % The seasonal figure is a fixed June-to-October product. Refuse to turn a
   % truncated or internally gapped run into an apparently complete season by
   % appending a synthetic endpoint to incomplete hourly support.
   output_step = seconds(policy.model_output_cadence_seconds);
   expected_time = (provenance.display_start:output_step: ...
      boundary_time - output_step).';
   display_time = model.Time(model.Time >= provenance.display_start ...
      & model.Time < boundary_time);
   if ~isequal(display_time, expected_time)
      error( ...
         'icemodel:verification:promiceAblationEvaluation:noncontiguousDisplaySupport', ...
         ['model diagnostics must contain every %g-second interval-start ' ...
         'row on the fixed June-to-October display season'], ...
         policy.model_output_cadence_seconds)
   end

   % Clone one schema-compatible context row, zero every saved interval ledger,
   % and carry only the preceding endpoint storage into the synthetic checkpoint.
   boundary = model(end, :);
   boundary.Properties.RowTimes = boundary_time;
   sum_fields = intersect(string(icemodel.namelists.budgetoutputs('sum')), ...
      string(model.Properties.VariableNames), 'stable');
   for name = sum_fields
      boundary.(name) = 0;
   end
   boundary.mass_budget_solid_start_mwe = ...
      model.mass_budget_solid_end_mwe(end);
   boundary.mass_budget_solid_end_mwe = ...
      model.mass_budget_solid_end_mwe(end);
   boundary.mass_budget_liquid_start_mwe = ...
      model.mass_budget_liquid_end_mwe(end);
   boundary.mass_budget_liquid_end_mwe = ...
      model.mass_budget_liquid_end_mwe(end);
   model = [model; boundary];

   % Label the copied non-ledger values so they cannot be interpreted as a
   % simulated future forcing step in downstream reports.
   metadata = struct( ...
      'source', source, ...
      'initialization_start', provenance.initialization_start, ...
      'initialization_policy', provenance.initialization_policy, ...
      'evaluation_start', provenance.evaluation_start, ...
      'evaluation_end', provenance.evaluation_end, ...
      'display_start', provenance.display_start, ...
      'display_end', provenance.display_end, ...
      'run_end_inclusive', provenance.run_end_inclusive, ...
      'comparison_boundary', boundary_time, ...
      'future_interval_executed', false, ...
      'boundary_sum_fields_zero', true, ...
      'boundary_role', ...
      "state checkpoint; non-ledger values copy preceding row for schema context");
end

function seasonal = seasonalDiagnostics(observations, model, provenance, policy)
   %SEASONALDIAGNOSTICS Build the snow-aware June-to-October report series.
   data = observations.data;
   inside = model.Time >= provenance.display_start ...
      & model.Time <= provenance.display_end;
   time = model.Time(inside);
   model_index = find(inside);
   if isempty(time)
      error('icemodel:verification:promiceAblationEvaluation:seasonalSupport', ...
         'June-through-October observations and model output do not overlap')
   end

   % Direct exposed-ice postings define visible support. Trace snow may remain
   % inside the selected season, but both observed and modeled curves are
   % censored where the height record may not represent ice ablation. Use the
   % complete model grid so an omitted observation remains a visible NaN gap.
   [has_observation, observation_index] = ismember(time, data.Time);
   n_time = numel(time);
   observation = NaN(n_time, 1);
   snow = NaN(n_time, 1);
   direct = false(n_time, 1);
   if any(has_observation)
      rows = observation_index(has_observation);
      observation(has_observation) = ...
         data.(policy.observation_field)(rows);
      snow(has_observation) = data.(policy.snow_variable)(rows);
      % Apply the PROMICE flag rules that decide whether an observation row
      % is supported.
      flag_fields = ...
         icemodel.verification.helpers.observationSupportFields( ...
         policy.observation_field, policy);
      support = icemodel.verification.helpers.classifyObservationSupport( ...
         double(data{rows, cellstr(flag_fields)}), flag_fields, ...
         policy.observation_field, policy);
      direct(has_observation) = support.flag_clean;
   end
   [ice_exposed, ~, ~] = ...
      icemodel.verification.helpers.classifySnowDepth( ...
      snow, policy.ice_exposure_threshold_m);
   cumulative_fields = string( ...
      icemodel.namelists.cumulativeoutputs());
   cumulative_index = model_index - 1;
   output_step = seconds(policy.model_output_cadence_seconds);
   if any(cumulative_index < 1) || any( ...
         model.Time(cumulative_index) + output_step ~= time)
      error('icemodel:verification:promiceAblationEvaluation:seasonalSupport', ...
         ['seasonal cumulative diagnostics require the preceding hourly ' ...
         'endpoint state at every observation time'])
   end
   cumulative_values = model{cumulative_index, cellstr(cumulative_fields)};
   ledger_fields = ["mass_budget_phase_solid_mwe", ...
      "mass_budget_vapor_solid_mwe", ...
      "mass_budget_top_export_solid_mwe", ...
      "mass_budget_top_export_liquid_mwe"];
   model_values = [cumulative_values, ...
      model{model_index, cellstr(ledger_fields)}];
   eligible = direct & ice_exposed & all(isfinite(model_values), 2);
   reference = find(eligible, 1);
   if isempty(reference)
      error('icemodel:verification:promiceAblationEvaluation:seasonalSupport', ...
         'June-through-October support has no direct exposed-ice posting')
   end

   % Rebase every simple cumulative model diagnostic at the same observed
   % exposed-ice posting. Diagnosed runoff remains signed and is never clamped.
   observation_lowering = observation - observation(reference);
   density = policy.effective_density_kg_m3;
   ro_liq = icemodel.physicalConstant('ro_liq');
   observation_density = observation_lowering .* density(:).' ./ ro_liq;
   % Signed negative lowering reverses the density-endpoint order, so every
   % stored lower/upper field must use the pointwise numeric extrema.
   observation_lower = min(observation_density, [], 2);
   observation_upper = max(observation_density, [], 2);
   observation_reference = observation_lowering ...
      * policy.effective_density_reference_kg_m3 / ro_liq;
   model_melt = cumulative_values(:, cumulative_fields == "melt");
   model_runoff = cumulative_values(:, cumulative_fields == "runoff");
   model_freeze = cumulative_values(:, cumulative_fields == "freeze");
   model_layer_change = cumulative_values(:, cumulative_fields == "dlayer");
   model_melt = model_melt - model_melt(reference);
   model_runoff = model_runoff - model_runoff(reference);
   model_freeze = model_freeze - model_freeze(reference);
   model_layer_change = model_layer_change - model_layer_change(reference);

   % Interval ledgers become state-at-time prefixes. The signed net-solid curve
   % is not made monotonic: decreases expose refreezing/storage behavior but
   % must not be read as geometric surface rise. Cumulative surface
   % mass loss is the separate monotonic lowering analogue, built from the mass
   % that top-cell removal actually exported rather than from cell geometry.
   increments = ...
      icemodel.verification.helpers.ablationLedgerIncrements(model);
   balance_prefix = [0; cumsum(increments.solid_balance(1:end - 1))];
   surface_prefix = [0; cumsum(increments.surface_loss(1:end - 1))];
   model_net_solid = balance_prefix(model_index) ...
      - balance_prefix(model_index(reference));
   model_surface_loss = surface_prefix(model_index) ...
      - surface_prefix(model_index(reference));

   % Sublimation removes ice without ever becoming liquid, so runoff cannot see
   % it. Runoff plus the SOLID vapor loss is the continuous ablation proxy:
   % melt that actually left under the runoff residence-time limit, plus
   % sublimation, minus deposition. Liquid vapor exchange is excluded because
   % evaporation removes pore water already counted as melt that left, and
   % condensation adds pore liquid rather than ice.
   d_vapor_loss = increments.solid_vapor_loss;
   vapor_prefix = [0; cumsum(d_vapor_loss(1:end - 1))];
   model_vapor_loss = vapor_prefix(model_index) ...
      - vapor_prefix(model_index(reference));
   model_ablation_proxy = model_runoff + model_vapor_loss;

   visible = eligible & ((1:numel(time))' >= reference);
   curves = {observation_lowering, observation_lower, observation_upper, ...
      observation_reference, ...
      model_melt, model_runoff, model_freeze, model_net_solid, ...
      model_layer_change, model_surface_loss, model_ablation_proxy};
   for k = 1:numel(curves)
      values = curves{k};
      values(~visible) = NaN;
      curves{k} = values;
   end
   [observation_lowering, observation_lower, observation_upper, ...
      observation_reference, ...
      model_melt, model_runoff, model_freeze, model_net_solid, ...
      model_layer_change, model_surface_loss, model_ablation_proxy] = ...
      curves{:};
   evaluation_window = time >= provenance.evaluation_start ...
      & time <= provenance.evaluation_end;
   seasonal = timetable(observation_lowering, observation_lower, ...
      observation_upper, observation_reference, ...
      model_melt, model_runoff, model_freeze, ...
      model_net_solid, model_layer_change, model_surface_loss, ...
      model_ablation_proxy, snow, ...
      ice_exposed, direct, evaluation_window, 'RowTimes', time, ...
      'VariableNames', {'observation_lowering_m', ...
      'observation_lower_mwe', 'observation_upper_mwe', ...
      'observation_reference_mwe', ...
      'model_melt_mwe', 'model_runoff_mwe', 'model_freeze_mwe', ...
      'model_net_solid_loss_mwe', 'model_layer_change_mwe', ...
      'model_surface_mass_loss_mwe', 'model_ablation_proxy_mwe', ...
      'snow_depth_m', 'ice_exposed', ...
      'direct_observation', 'evaluation_window'});
end

function requireExactEndpoints(comparison, t0, t1)
   %REQUIREEXACTENDPOINTS Reject implicit comparator endpoint movement.
   if comparison.window_start ~= t0 || comparison.window_end ~= t1
      error('icemodel:verification:promiceAblationEvaluation:endpointShift', ...
         'comparison could not retain both readiness-selected direct endpoints')
   end
end

function rows = nestedDiagnostics( ...
      row, observations, model, primary, policy, t0, t1)
   %NESTEDDIAGNOSTICS Evaluate fixed cumulative windows plus the longest window.
   rows = table2struct(predeclaredNestedDiagnostics(row, policy));
   for k = 1:numel(rows)
      if rows(k).window_label == "longest"
         rows(k) = completedNestedRecord( ...
            rows(k), primary, t0, t1);
         continue
      end

      % Fixed windows are unavailable rather than shortened when their exact
      % predeclared endpoint lies outside the primary interval or is flagged.
      target_end = rows(k).window_end;
      if target_end > t1
         rows(k).reason = "fixed endpoint exceeds the primary window";
         continue
      end
      try
         comparison = icemodel.verification.compareAblation( ...
            observations, model, window_start=t0, window_end=target_end);
         if comparison.window_start ~= t0 ...
               || comparison.window_end ~= target_end
            rows(k).reason = ...
               "exact endpoint is unavailable or observation-flagged";
            continue
         end
         rows(k) = completedNestedRecord( ...
            rows(k), comparison, t0, target_end);
      catch err
         rows(k).reason = string(err.message);
         rows(k).error_identifier = string(err.identifier);
      end
   end
   rows = struct2table(rows);
end

function rows = predeclaredNestedDiagnostics(row, policy)
   %PREDECLAREDNESTEDDIAGNOSTICS Define every selected row before execution.
   fixed_days = reshape(policy.nested_window_days, [], 1);
   labels = [compose('%dd', fixed_days); "longest"];
   t0 = parseTime(row.snow_free_window_start);
   t1 = parseTime(row.snow_free_window_end);
   target_days = [fixed_days; days(t1 - t0)];
   records = repmat(emptyNestedRecord(), numel(labels), 1);
   for k = 1:numel(labels)
      records(k) = emptyNestedRecord(row, labels(k), target_days(k));
      records(k).window_start = t0;
      if labels(k) == "longest"
         records(k).window_end = t1;
      else
         records(k).window_end = t0 + days(target_days(k));
      end
   end
   rows = struct2table(records);
end

function rows = endpointPerturbationDiagnostics( ...
      row, observations, model, policy)
   %ENDPOINTPERTURBATIONDIAGNOSTICS Recompare one run at planned endpoints.
   records = table2struct(predeclaredEndpointPerturbations(row, policy));
   for k = 1:numel(records)
      requested_start = records(k).requested_window_start;
      requested_end = records(k).requested_window_end;

      % A planned perturbation remains visible even if it collapses a short
      % primary interval; it is never shortened or omitted.
      if isnat(requested_start) || isnat(requested_end) ...
            || requested_start >= requested_end
         records(k).reason = ...
            "perturbed window must satisfy requested start < requested end";
         records(k).error_identifier = ...
            "icemodel:verification:promiceAblationEvaluation:invalidPerturbedWindow";
         continue
      end

      % Reuse the completed observation/model pair. These sensitivity rows are
      % comparator calls only and must never launch another model simulation.
      try
         [comparison, ~, diagnostics] = ...
            icemodel.verification.compareAblation( ...
            observations, model, window_start=requested_start, ...
            window_end=requested_end);
         records(k).actual_window_start = comparison.window_start;
         records(k).actual_window_end = comparison.window_end;
         if comparison.window_start ~= requested_start ...
               || comparison.window_end ~= requested_end
            records(k).reason = ...
               "exact perturbed endpoints are unavailable or observation-flagged";
            records(k).error_identifier = ...
               "icemodel:verification:promiceAblationEvaluation:perturbationEndpointShift";
            continue
         end
         records(k) = completedEndpointPerturbationRecord( ...
            records(k), comparison, diagnostics);
      catch err
         [records(k).reason, records(k).error_identifier] = ...
            failureDetails(err);
      end
   end
   rows = struct2table(records);
end

function rows = predeclaredEndpointPerturbations(row, policy)
   %PREDECLAREDENDPOINTPERTURBATIONS Define sensitivity rows before execution.
   perturbation_days = reshape(double( ...
      policy.endpoint_perturbation_days), [], 1);
   records = repmat(emptyEndpointPerturbationRecord(), ...
      2 * numel(perturbation_days), 1);
   t0 = parseTime(row.snow_free_window_start);
   t1 = parseTime(row.snow_free_window_end);
   for k = 1:numel(perturbation_days)
      offset = perturbation_days(k);

      % Positive start offsets move the start later while retaining the primary
      % end; negative end offsets move the end earlier while retaining the start.
      records(k) = emptyEndpointPerturbationRecord(row, ...
         compose('start_plus_%gd', offset), "start", offset);
      records(k).requested_window_start = t0 + days(offset);
      records(k).requested_window_end = t1;
      end_index = numel(perturbation_days) + k;
      records(end_index) = emptyEndpointPerturbationRecord(row, ...
         compose('end_minus_%gd', offset), "end", -offset);
      records(end_index).requested_window_start = t0;
      records(end_index).requested_window_end = t1 - days(offset);
   end
   rows = struct2table(records);
end

function [reason, error_identifier] = failureDetails(err)
   %FAILUREDETAILS Preserve one stable failure identity across every artifact.
   reason = string(err.message);
   error_identifier = string(err.identifier);
   if error_identifier == ""
      error_identifier = ...
         "icemodel:verification:promiceAblationEvaluation:executionFailed";
   end
end

function record = completedNestedRecord(record, comparison, t0, t1)
   %COMPLETEDNESTEDRECORD Copy one exact-endpoint comparison into a flat row.
   record.available = true;
   record.reason = "";
   record.window_start = t0;
   record.window_end = t1;
   record.classification = comparison.classification;
   record.physical_comparable = comparison.physical_comparable;
   record.observation_intact_mwe = comparison.observation_intact_mwe;
   record.model_solid_loss_mwe = comparison.model_solid_loss_mwe;
   record.model_minus_observation_mwe = ...
      comparison.model_minus_observation_mwe;
   record.relative_difference = comparison.relative_difference;
end

function record = completedEndpointPerturbationRecord( ...
      record, comparison, diagnostics)
   %COMPLETEDENDPOINTPERTURBATIONRECORD Flatten one exact sensitivity result.
   record.available = true;
   record.reason = "";
   record.error_identifier = "";
   record.actual_window_start = comparison.window_start;
   record.actual_window_end = comparison.window_end;
   record.classification = comparison.classification;
   record.physical_comparable = comparison.physical_comparable;
   record.observation_lowering_m = comparison.observation_lowering_m;
   record.observation_intact_mwe = comparison.observation_intact_mwe;

   % The effective-density values are labeled sensitivity scenarios in the
   % comparator; min/max flatten them without promoting them to uncertainty.
   sensitivity = ...
      diagnostics.effective_density.observation_sensitivity_mwe;
   record.observation_sensitivity_min_mwe = min(sensitivity);
   record.observation_sensitivity_max_mwe = max(sensitivity);
   record.model_solid_loss_mwe = comparison.model_solid_loss_mwe;
   record.model_minus_observation_mwe = ...
      comparison.model_minus_observation_mwe;
   record.relative_difference = comparison.relative_difference;
end

function time = parseTime(value)
   %PARSETIME Convert one readiness timestamp to UTC.
   if strlength(string(value)) == 0
      time = NaT(1, 1, 'TimeZone', 'UTC');
   else
      time = icemodel.verification.setup.ensureUtc(string(value));
   end
end


function result = emptySiteYearResult(row)
   %EMPTYSITEYEARRESULT Return one stable per-site-year result schema.
   if nargin == 0
      row = table();
      case_id = "";
      site_id = "";
      year_value = 0;
   else
      case_id = row.case_id;
      site_id = row.site_id;
      year_value = row.year;
   end
   result = struct( ...
      'case_id', case_id, 'site_id', site_id, 'year', year_value, ...
      'status', "", 'reason', "", 'error_identifier', "", ...
      'readiness_row', row, 'manifest', struct(), ...
      'observations', struct(), 'model', timetable(), ...
      'model_options', struct(), 'model_metadata', struct(), ...
      'comparison', struct(), 'aligned', timetable(), ...
      'seasonal', timetable(), ...
      'diagnostics', struct(), 'nested_windows', emptyNestedTable(), ...
      'endpoint_perturbations', emptyEndpointPerturbationTable());
end

function record = emptySummaryRecord(row, selected, readiness_policy)
   %EMPTYSUMMARYRECORD Return one flat machine-readable cohort row.
   if nargin == 0
      row = table("", "", 0, false, "", "", "", "", 0, 0, 0, ...
         "", "", ...
         'VariableNames', {'case_id', 'site_id', 'year', 'admitted', ...
         'forcing_artifact', 'forcing_sha256', 'observation_artifact', ...
         'observation_sha256', 'possible_support_count', ...
         'direct_target_count', 'snow_free_direct_count', ...
         'snow_free_window_start', 'snow_free_window_end'});
      selected = false;
      readiness_policy = ...
         icemodel.verification.namelists.promiceAblationReadiness();
   end
   window_possible = windowPossibleSampleCount(row, readiness_policy);
   window_direct = double(row.snow_free_direct_count);
   record = struct( ...
      'case_id', string(row.case_id), ...
      'site_id', string(row.site_id), ...
      'year', double(row.year), ...
      'selected', logical(selected), ...
      'admitted', logical(row.admitted), ...
      'status', "", 'reason', "", 'error_identifier', "", ...
      'forcing_artifact', string(row.forcing_artifact), ...
      'forcing_sha256', string(row.forcing_sha256), ...
      'observation_artifact', string(row.observation_artifact), ...
      'observation_sha256', string(row.observation_sha256), ...
      'possible_support_count', double(row.possible_support_count), ...
      'direct_target_count', double(row.direct_target_count), ...
      'window_valid_sample_count', NaN, ...
      'window_direct_sample_count', window_direct, ...
      'window_possible_sample_count', window_possible, ...
      'coverage_fraction', NaN, ...
      'coverage_denominator', ...
      "inclusive [t0,t1] timestamps at the readiness observation cadence", ...
      'window_start', NaT(1, 1, 'TimeZone', 'UTC'), ...
      'window_end', NaT(1, 1, 'TimeZone', 'UTC'), ...
      'classification', "", 'physical_comparable', false, ...
      'eligible_sample_count', NaN, ...
      'observation_lowering_m', NaN, ...
      'observation_intact_mwe', NaN, ...
      'model_solid_loss_mwe', NaN, ...
      'model_minus_observation_mwe', NaN, ...
      'relative_difference', NaN);
end

function record = completedSummaryRecord( ...
      row, result, selected, readiness_policy)
   %COMPLETEDSUMMARYRECORD Flatten one completed primary comparison.
   record = emptySummaryRecord(row, selected, readiness_policy);
   comparison = result.comparison;
   record.status = result.status;
   record.window_start = comparison.window_start;
   record.window_end = comparison.window_end;
   record.classification = comparison.classification;
   record.physical_comparable = comparison.physical_comparable;
   record.eligible_sample_count = comparison.eligible_sample_count;
   record.window_valid_sample_count = comparison.eligible_sample_count;
   record.coverage_fraction = coverageFraction( ...
      record.window_valid_sample_count, record.window_possible_sample_count);
   record.observation_lowering_m = comparison.observation_lowering_m;
   record.observation_intact_mwe = comparison.observation_intact_mwe;
   record.model_solid_loss_mwe = comparison.model_solid_loss_mwe;
   record.model_minus_observation_mwe = ...
      comparison.model_minus_observation_mwe;
   record.relative_difference = comparison.relative_difference;
end

function count = windowPossibleSampleCount(row, readiness_policy)
   %WINDOWPOSSIBLESAMPLECOUNT Count hourly opportunities inside [t0,t1].
   t0 = parseTime(row.snow_free_window_start);
   t1 = parseTime(row.snow_free_window_end);
   cadence = double(readiness_policy.observation_cadence_seconds);
   if isnat(t0) || isnat(t1) || t1 < t0 ...
         || ~isfinite(cadence) || cadence <= 0
      count = 0;
      return
   end
   count = floor(seconds(t1 - t0) / cadence) + 1;
end

function fraction = coverageFraction(valid_count, possible_count)
   %COVERAGEFRACTION Return an explicit undefined fraction at zero denominator.
   if possible_count == 0
      fraction = NaN;
   else
      fraction = valid_count / possible_count;
   end
end

function record = emptyNestedRecord(row, label, target_days)
   %EMPTYNESTEDRECORD Return one stable cumulative-window summary row.
   if nargin == 0
      case_id = "";
      site_id = "";
      year_value = 0;
      label = "";
      target_days = NaN;
   else
      case_id = row.case_id;
      site_id = row.site_id;
      year_value = row.year;
   end
   record = struct( ...
      'case_id', string(case_id), 'site_id', string(site_id), ...
      'year', double(year_value), 'window_label', string(label), ...
      'target_days', double(target_days), 'available', false, ...
      'reason', "", 'error_identifier', "", ...
      'window_start', NaT(1, 1, 'TimeZone', 'UTC'), ...
      'window_end', NaT(1, 1, 'TimeZone', 'UTC'), ...
      'classification', "", 'physical_comparable', false, ...
      'observation_intact_mwe', NaN, 'model_solid_loss_mwe', NaN, ...
      'model_minus_observation_mwe', NaN, 'relative_difference', NaN);
end

function value = emptyNestedTable()
   %EMPTYNESTEDTABLE Return the zero-row nested-window public schema.
   value = struct2table(emptyNestedRecord());
   value = value([], :);
end

function record = emptyEndpointPerturbationRecord( ...
      row, label, axis_name, signed_days)
   %EMPTYENDPOINTPERTURBATIONRECORD Return one stable sensitivity row.
   if nargin == 0
      case_id = "";
      site_id = "";
      year_value = 0;
      label = "";
      axis_name = "";
      signed_days = NaN;
   else
      case_id = row.case_id;
      site_id = row.site_id;
      year_value = row.year;
   end
   record = struct( ...
      'case_id', string(case_id), 'site_id', string(site_id), ...
      'year', double(year_value), 'perturbation_label', string(label), ...
      'perturbation_axis', string(axis_name), ...
      'perturbation_days_signed', double(signed_days), ...
      'requested_window_start', NaT(1, 1, 'TimeZone', 'UTC'), ...
      'requested_window_end', NaT(1, 1, 'TimeZone', 'UTC'), ...
      'actual_window_start', NaT(1, 1, 'TimeZone', 'UTC'), ...
      'actual_window_end', NaT(1, 1, 'TimeZone', 'UTC'), ...
      'available', false, 'reason', "", 'error_identifier', "", ...
      'classification', "", 'physical_comparable', false, ...
      'observation_lowering_m', NaN, 'observation_intact_mwe', NaN, ...
      'observation_sensitivity_min_mwe', NaN, ...
      'observation_sensitivity_max_mwe', NaN, ...
      'model_solid_loss_mwe', NaN, ...
      'model_minus_observation_mwe', NaN, 'relative_difference', NaN);
end

function value = emptyEndpointPerturbationTable()
   %EMPTYENDPOINTPERTURBATIONTABLE Return the zero-row sensitivity schema.
   value = struct2table(emptyEndpointPerturbationRecord());
   value = value([], :);
end

function paths = resultPaths(run_dir, readiness_files, write_artifacts)
   %RESULTPATHS Return persisted output paths or the empty path schema.
   if ~write_artifacts
      paths = emptyPaths();
      return
   end
   paths = struct( ...
      'run_dir', string(run_dir), ...
      'readiness_csv', string(readiness_files.csv), ...
      'readiness_json', string(readiness_files.json), ...
      'summary_csv', string(fullfile(run_dir, 'summary.csv')), ...
      'nested_windows_csv', ...
      string(fullfile(run_dir, 'nested-windows.csv')), ...
      'endpoint_perturbations_csv', ...
      string(fullfile(run_dir, 'endpoint-perturbations.csv')), ...
      'results_mat', string(fullfile(run_dir, 'results.mat')));
end

function paths = emptyPaths()
   %EMPTYPATHS Return the stable no-artifact path schema.
   paths = struct( ...
      'run_dir', "", 'readiness_csv', "", 'readiness_json', "", ...
      'summary_csv', "", 'nested_windows_csv', "", ...
      'endpoint_perturbations_csv', "", 'results_mat', "");
end
