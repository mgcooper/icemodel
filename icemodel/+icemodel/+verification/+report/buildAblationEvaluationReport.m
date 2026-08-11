function report_file = buildAblationEvaluationReport(results_file, kwargs)
   %BUILDABLATIONEVALUATIONREPORT Render saved PROMICE ablation results.
   %
   %  report_file = ...
   %     icemodel.verification.report.buildAblationEvaluationReport( ...
   %     "test/artifacts/<run>/results.mat")
   %  report_file = ...
   %     icemodel.verification.report.buildAblationEvaluationReport( ...
   %     results_file, render=false, output_dir=tempdir)
   %
   % The MAT file's RESULTS struct is the only scientific input. This function
   % writes compact CSVs, MATLAB-exported PNGs, inspectable QMD, and optionally a
   % self-contained Quarto HTML report without loading canonical data or running
   % IceModel.

   arguments
      results_file (1, 1) string
      kwargs.render (1, 1) logical = true
      kwargs.output_dir (1, 1) string = ""
   end

   % Load one saved results file. Paths stored inside it are printed as-is;
   % this function never opens them.
   if ~isfile(results_file)
      error('icemodel:verification:report:missingAblationResults', ...
         'Saved ablation results are unavailable: %s', results_file)
   end
   source_sha256 = ...
      icemodel.verification.setup.fileSha256(results_file);
   saved_variables = whos('-file', results_file);
   if ~any(string({saved_variables.name}) == "results")
      error('icemodel:verification:report:invalidAblationResults', ...
         'The MAT file must contain one scalar RESULTS struct')
   end
   saved = load(results_file, 'results');
   if icemodel.verification.setup.fileSha256(results_file) ~= source_sha256
      error('icemodel:verification:report:changedAblationResults', ...
         'Saved ablation results changed while the report was loading: %s', ...
         results_file)
   end
   if ~isstruct(saved.results) ...
         || ~isscalar(saved.results)
      error('icemodel:verification:report:invalidAblationResults', ...
         'The MAT file must contain one scalar RESULTS struct')
   end
   results = validateResults(saved.results);

   % All derived report artifacts land beside the results file unless the
   % caller selects another output directory for focused inspection.
   output_dir = kwargs.output_dir;
   if output_dir == ""
      output_dir = string(fileparts(results_file));
      if output_dir == ""
         output_dir = ".";
      end
   end
   icemodel.helpers.ensureDirExists(output_dir)
   asset_dir = fullfile(output_dir, 'report-assets');
   icemodel.helpers.ensureDirExists(asset_dir)

   % Write report-specific copies of the four public ledgers and derive only
   % flat presentation tables from the saved per-site diagnostics.
   tables = reportTables(results);
   files = writeReportTables(tables, output_dir);
   assets = buildFigures(results, tables, asset_dir);

   % Plain generated QMD remains inspectable and contains no executable code.
   qmd_file = fullfile(output_dir, 'promice-ablation-evaluation-report.qmd');
   report_file = fullfile(output_dir, ...
      'promice-ablation-evaluation-report.html');
   manifest_file = fullfile(output_dir, ...
      'report-artifact-sha256.csv');
   lines = reportMarkdown(results, results_file, tables, files, assets, ...
      report_file, source_sha256, manifest_file);
   writelines(lines, qmd_file)

   % Rendering is optional for focused tests and source review only.
   if kwargs.render
      command = "quarto render " ...
         + icemodel.shellQuote(qmd_file);
      [status, output] = system(command);
      if status ~= 0
         error('icemodel:verification:report:quartoFailed', ...
            'Quarto failed to render %s:\n%s', qmd_file, output)
      end
      assert(isfile(report_file), ...
         'icemodel:verification:report:missingHtml', ...
         'Quarto did not create the expected report: %s', report_file)
   end

   % Hash the finalized report package after optional rendering. The manifest
   % excludes itself so its identity is nonrecursive and reproducible.
   writeReportArtifactManifest(results_file, qmd_file, report_file, ...
      kwargs.render, files, assets, manifest_file, source_sha256);
end

function results = validateResults(results)
   %VALIDATERESULTS Check the saved results have the fields the report needs.

   required = ["run_name", "readiness", "site_year_results", ...
      "summary", "nested_windows", "endpoint_perturbations", ...
      "policy", "paths"];
   missing = required(~isfield(results, required));
   if ~isempty(missing) || ~isstruct(results.readiness) ...
         || ~isfield(results.readiness, 'rows') ...
         || ~istable(results.readiness.rows) ...
         || ~istable(results.summary) || ~istable(results.nested_windows) ...
         || ~istable(results.endpoint_perturbations) ...
         || ~isstruct(results.site_year_results)
      error('icemodel:verification:report:invalidAblationResults', ...
         'RESULTS does not match the saved PROMICE evaluation contract')
   end

   % Flat summary fields are required because they define cohort state and the
   % cross-site comparison without reopening per-case artifacts.
   summary_fields = ["case_id", "site_id", "year", "selected", ...
      "admitted", "status", "classification", "physical_comparable", ...
      "observation_intact_mwe", "model_solid_loss_mwe", ...
      "window_start", "window_end"];
   if ~all(ismember(summary_fields, ...
         string(results.summary.Properties.VariableNames)))
      error('icemodel:verification:report:invalidAblationSummary', ...
         'RESULTS.summary is missing required cohort fields')
   end

   % Validate and normalize the saved observation registry even when no case
   % completed, so readiness-only reports cannot bypass the policy check.
   results.policy = validateObservationPolicy(results.policy);
   validateCompletedSeasonalPayloads( ...
      results.site_year_results, results.policy)

   % Endpoint perturbations are a planned six-row policy experiment per
   % selected site-year; unavailable rows must retain their full schema.
   endpoint_fields = ["case_id", "site_id", "year", ...
      "perturbation_label", "perturbation_axis", ...
      "perturbation_days_signed", "requested_window_start", ...
      "requested_window_end", "actual_window_start", "actual_window_end", ...
      "available", "reason", "error_identifier", "classification", ...
      "physical_comparable", "observation_lowering_m", ...
      "observation_intact_mwe", "observation_sensitivity_min_mwe", ...
      "observation_sensitivity_max_mwe", "model_solid_loss_mwe", ...
      "model_minus_observation_mwe", "relative_difference"];
   if ~all(ismember(endpoint_fields, ...
         string(results.endpoint_perturbations.Properties.VariableNames)))
      error('icemodel:verification:report:invalidEndpointPerturbations', ...
         'RESULTS.endpoint_perturbations is missing required fields')
   end
end

function policy = validateObservationPolicy(policy)
   %VALIDATEOBSERVATIONPOLICY Check the saved flag lists do not overlap or
   % miss a field.

   partition_fields = ["observation_field", "required_observation_fields", ...
      "support_flag_fields", "direct_zero_flag_fields", ...
      "datum_break_flag_fields", "ordinary_gap_flag_fields", ...
      "metadata_only_flag_fields"];
   scientific_fields = ["version", "snow_variable", ...
      "evaluation_season_start_month_day", ...
      "evaluation_season_end_month_day", "effective_density_kg_m3", ...
      "effective_density_reference_kg_m3", ...
      "snow_continuity_threshold_m", "ice_exposure_threshold_m"];
   if ~isstruct(policy) || ~isscalar(policy) ...
         || ~all(isfield(policy, [partition_fields, scientific_fields]))
      error('icemodel:verification:report:invalidAblationPolicy', ...
         'Saved policy lacks required observation-flag partitions')
   end
   if ~isnumeric(policy.evaluation_season_start_month_day) ...
         || ~isnumeric(policy.evaluation_season_end_month_day) ...
         || ~isnumeric(policy.effective_density_kg_m3) ...
         || ~isnumeric(policy.snow_continuity_threshold_m) ...
         || ~isnumeric(policy.ice_exposure_threshold_m)
      error('icemodel:verification:report:invalidAblationPolicy', ...
         'Saved policy has invalid season or density values')
   end
   observation_field = string(policy.observation_field);
   required_schema = string(policy.required_observation_fields(:)');
   support_fields = string(policy.support_flag_fields(:)');
   zero_fields = string(policy.direct_zero_flag_fields(:)');
   datum_fields = string(policy.datum_break_flag_fields(:)');
   gap_fields = string(policy.ordinary_gap_flag_fields(:)');
   metadata_fields = string(policy.metadata_only_flag_fields(:)');
   season_start = double(policy.evaluation_season_start_month_day(:)');
   season_end = double(policy.evaluation_season_end_month_day(:)');
   density = double(policy.effective_density_kg_m3(:)');
   snow_continuity = double(policy.snow_continuity_threshold_m);
   ice_exposure = double(policy.ice_exposure_threshold_m);

   % The registry is closed: support splits into direct-zero and metadata-only,
   % while direct-zero splits into disjoint datum-break and ordinary-gap roles.
   role_sets = {required_schema, support_fields, zero_fields, datum_fields, ...
      gap_fields, metadata_fields};
   has_duplicates = any(cellfun(@(value) ...
      numel(unique(value, 'stable')) ~= numel(value), role_sets));
   all_fields = [required_schema, support_fields, zero_fields, ...
      datum_fields, gap_fields, metadata_fields];
   invalid_names = any(ismissing(all_fields) | strlength(all_fields) == 0);
   expected_support = unique([zero_fields, metadata_fields], 'stable');
   expected_zero = unique([datum_fields, gap_fields], 'stable');
   snow_field = string(policy.snow_variable);
   invalid_snow_field = numel(snow_field) ~= 1 ...
      || ismissing(snow_field) || strlength(snow_field) == 0 ...
      || ismember(snow_field, support_fields);
   required_data_fields = [observation_field, snow_field];
   expected_required = unique( ...
      [required_data_fields, support_fields], 'stable');
   invalid_partition = numel(observation_field) ~= 1 ...
      || ismember(observation_field, support_fields) ...
      || invalid_snow_field ...
      || ~isempty(intersect(metadata_fields, zero_fields)) ...
      || ~isempty(intersect(datum_fields, gap_fields)) ...
      || ~isempty(setxor(support_fields, expected_support)) ...
      || ~isempty(setxor(zero_fields, expected_zero)) ...
      || ~isempty(setxor(required_schema, expected_required));
   invalid_scientific = ~validMonthDay(season_start) ...
      || ~validMonthDay(season_end) ...
      || numel(density) < 2 || any(~isfinite(density) | density <= 0) ...
      || any(diff(density) <= 0) ...
      || ~isscalar(snow_continuity) || ~isfinite(snow_continuity) ...
      || snow_continuity < 0 ...
      || ~isscalar(ice_exposure) || ~isfinite(ice_exposure) ...
      || ice_exposure < 0;
   canonical = ...
      icemodel.verification.namelists.promiceAblationPolicy();

   % Compare the VALUES that govern results, not the prose that explains them,
   % so rewording a citation does not invalidate a saved cohort. A report
   % cannot describe a run whose policy values differ from the current
   % namelist.
   invalid_fixed_policy = ~isequaln( ...
      strippedPolicyValues(policy, canonical), ...
      strippedPolicyValues(canonical, canonical));
   if has_duplicates || invalid_names || invalid_partition ...
         || invalid_scientific || invalid_fixed_policy
      error('icemodel:verification:report:invalidAblationPolicy', ...
         'Saved policy does not match the fixed PROMICE ablation policy')
   end

   % Hand the report one validated policy struct with row-shaped fields.
   policy.observation_field = observation_field;
   policy.required_observation_fields = required_schema;
   policy.support_flag_fields = support_fields;
   policy.direct_zero_flag_fields = zero_fields;
   policy.datum_break_flag_fields = datum_fields;
   policy.ordinary_gap_flag_fields = gap_fields;
   policy.metadata_only_flag_fields = metadata_fields;
   policy.evaluation_season_start_month_day = season_start;
   policy.evaluation_season_end_month_day = season_end;
   policy.effective_density_kg_m3 = density;
   policy.snow_variable = snow_field;
   policy.snow_continuity_threshold_m = snow_continuity;
   policy.ice_exposure_threshold_m = ice_exposure;
end

function values = strippedPolicyValues(policy, canonical)
   %STRIPPEDPOLICYVALUES Drop explanatory fields before comparing policies.
   %
   % The namelist names its own documentation fields. A saved policy that
   % lacks one of them is valid, so removal is guarded by isfield.

   values = policy;
   documentation_fields = canonical.documentation_fields;
   for k = 1:numel(documentation_fields)
      field = char(documentation_fields(k));
      if isfield(values, field)
         values = rmfield(values, field);
      end
   end

   % documentation_fields is itself documentation bookkeeping.
   if isfield(values, 'documentation_fields')
      values = rmfield(values, 'documentation_fields');
   end
end

function tf = validMonthDay(value)
   %VALIDMONTHDAY Return true for one calendar-valid [month, day] pair.

   tf = isnumeric(value) && numel(value) == 2 ...
      && all(isfinite(value)) && all(value == fix(value));
   if tf
      try
         datetime(2000, value(1), value(2));
      catch
         tf = false;
      end
   end
end

function validateCompletedSeasonalPayloads(site_results, policy)
   %VALIDATECOMPLETEDSEASONALPAYLOADS Require one hourly season per
   % completed site-year.

   required = ["observation_lowering_m", "observation_lower_mwe", ...
      "observation_upper_mwe", "observation_reference_mwe", ...
      "model_melt_mwe", "model_runoff_mwe", ...
      "model_freeze_mwe", "model_net_solid_loss_mwe", ...
      "model_layer_change_mwe", "model_surface_mass_loss_mwe", ...
      "model_ablation_proxy_mwe", "snow_depth_m", ...
      "ice_exposed", "direct_observation", "evaluation_window"];
   for k = 1:numel(site_results)
      result = site_results(k);
      if ~isfield(result, 'status') || string(result.status) ~= "completed"
         continue
      end
      valid_header = isfield(result, 'year') ...
         && isnumeric(result.year) && isscalar(result.year) ...
         && isfinite(result.year) && result.year == fix(result.year) ...
         && isfield(result, 'seasonal') && istimetable(result.seasonal) ...
         && all(ismember(required, string( ...
         result.seasonal.Properties.VariableNames))) ...
         && isfield(result, 'comparison') ...
         && isstruct(result.comparison) && isscalar(result.comparison) ...
         && all(isfield(result.comparison, ...
         ["window_start", "window_end"]));
      if ~valid_header
         error('icemodel:verification:report:invalidAblationSeasonal', ...
            'Completed result %d lacks the required seasonal contract', k)
      end

      % The persisted display payload must be the full inclusive policy grid;
      % exact equality also rejects duplicates, reordering, gaps, and truncation.
      [season_start, season_end] = ...
      icemodel.verification.helpers.evaluationSeason( ...
      double(result.year), policy);
      expected_time = (season_start:hours(1):season_end)';
      if ~isequal(result.seasonal.Time, expected_time)
         error('icemodel:verification:report:invalidAblationSeasonal', ...
            ['Completed result %d must retain the exact hourly policy grid ' ...
            'from 1 June through 1 October inclusive'], k)
      end

      % Evaluation labels are binary and nonempty, and their true rows exactly
      % reproduce the saved comparison endpoints. An all-true mask is valid
      % only when those bounds span the complete display season.
      labels = result.seasonal.evaluation_window;
      valid_labels = islogical(labels) || (isnumeric(labels) ...
         && all(isfinite(labels)) && all(labels == 0 | labels == 1));
      t0 = result.comparison.window_start;
      t1 = result.comparison.window_end;
      valid_bounds = isdatetime(t0) && isscalar(t0) && ~ismissing(t0) ...
         && isdatetime(t1) && isscalar(t1) && ~ismissing(t1) && t0 < t1 ...
         && any(expected_time == t0) && any(expected_time == t1);
      if ~valid_labels || ~valid_bounds
         error('icemodel:verification:report:invalidAblationSeasonal', ...
            'Completed result %d has invalid evaluation-window metadata', k)
      end
      labels = logical(labels);
      expected_labels = expected_time >= t0 & expected_time <= t1;
      full_season = t0 == season_start && t1 == season_end;
      invalid_full_mask = all(labels) && ~full_season;
      if ~isequal(labels, expected_labels) || ~any(labels) ...
            || invalid_full_mask
         error('icemodel:verification:report:invalidAblationSeasonal', ...
            ['Completed result %d evaluation_window must exactly match the ' ...
            'saved comparison bounds'], k)
      end
   end
end

function tables = reportTables(results)
   %REPORTTABLES Collect exact public ledgers and flat diagnostic tables.

   tables.readiness = results.readiness.rows;
   tables.summary = results.summary;
   tables.nested = results.nested_windows;
   tables.endpoint = kanFirstTable(results.endpoint_perturbations);
   tables.synthesis = synthesisTable(results.summary);
   tables.status = statusTable(tables.synthesis);
   tables.support = supportTable(results.site_year_results, results.policy);
   tables.components = componentTable(results.site_year_results);
   tables.grid_translation = gridTranslationTable(results.site_year_results);
   tables.initialization = initializationTable(results.site_year_results);
   % Scoring every diagnostic against the aligned observations is the primary
   % purpose of the evaluation, so it is a first-class report table.
   [tables.performance, tables.performance_summary] = ...
      icemodel.verification.ablationPerformanceMetrics( ...
      results.site_year_results);
   [tables.identities, tables.materiality, tables.scenarios, ...
      tables.effective_density, tables.exclusions] = ...
      diagnosticTables(results.site_year_results);
   % Observation-rate plausibility is a caveat on the observations themselves,
   % so it is reported beside the model diagnostics rather than inside them.
   tables.observation_rates = icemodel.verification.observationRateOutliers( ...
      results.summary, results.policy);
end

function files = writeReportTables(tables, output_dir)
   %WRITEREPORTTABLES Persist machine-readable tables beside the report.

   names = ["readiness", "summary", "nested", "endpoint", ...
      "synthesis", "status", "support", "components", "grid_translation", ...
      "initialization", "performance", "performance_summary", ...
      "identities", "materiality", "scenarios", ...
      "effective_density", "exclusions", "observation_rates"];
   basenames = ["report-readiness.csv", "report-summary.csv", ...
      "report-nested-windows.csv", "report-endpoint-perturbations.csv", ...
      "site-year-synthesis.csv", "synthesis-status.csv", ...
      "observation-support.csv", "signed-solid-components.csv", ...
      "grid-translation.csv", "model-initialization.csv", ...
      "performance-metrics.csv", "performance-summary.csv", ...
      "closure-identities.csv", ...
      "materiality.csv", "accounting-scenarios.csv", ...
      "effective-density-sensitivity.csv", "support-exclusions.csv", ...
      "observation-rate-outliers.csv"];
   files = struct();
   for k = 1:numel(names)
      % Writetable preserves saved numeric precision in a compact open format.
      files.(names(k)) = string(fullfile(output_dir, basenames(k)));
      writetable(tables.(names(k)), files.(names(k)))
   end
end

function manifest = writeReportArtifactManifest(results_file, qmd_file, ...
      report_file, include_html, files, assets, manifest_file, source_sha256)
   %WRITEREPORTARTIFACTMANIFEST Hash every finalized report input and output.

   % Table fields are inserted in the same declared order as WRITEREPORTTABLES.
   table_fields = string(fieldnames(files));
   table_files = strings(numel(table_fields), 1);
   for k = 1:numel(table_fields)
      table_files(k) = string(files.(table_fields(k)));
   end
   table_roles = "table_" + table_fields;
   table_paths = arrayfun(@relativeName, table_files);

   % Figure paths retain their output-local report-assets prefix. Caption
   % fields hold prose rather than paths and are already embedded in the
   % QMD/HTML, so they are the only excluded asset fields; every other field
   % of the struct is hashed.
   asset_fields = string(fieldnames(assets))';
   asset_fields = asset_fields(~endsWith(asset_fields, "_caption"));
   figure_file_parts = cell(numel(asset_fields), 1);
   figure_role_parts = cell(numel(asset_fields), 1);
   for k = 1:numel(asset_fields)
      values = string(assets.(asset_fields(k)));
      values = values(strlength(values) > 0);
      figure_file_parts{k} = values(:);
      figure_role_parts{k} = repmat( ...
         "figure_" + asset_fields(k), numel(values), 1);
   end
   figure_files = vertcat(figure_file_parts{:});
   figure_roles = vertcat(figure_role_parts{:});
   figure_paths = "report-assets/" ...
      + arrayfun(@relativeName, figure_files);

   report_files = string(qmd_file);
   report_roles = "report_qmd";
   report_paths = relativeName(qmd_file);
   if include_html
      report_files = [report_files; string(report_file)];
      report_roles = [report_roles; "report_html"];
      report_paths = [report_paths; relativeName(report_file)];
   end
   artifact_files = [string(results_file); report_files; ...
      table_files; figure_files];
   artifact_role = ["source_results_mat"; report_roles; ...
      table_roles; figure_roles];
   artifact_path = [string(results_file); report_paths; ...
      table_paths; figure_paths];
   bytes = zeros(numel(artifact_files), 1);
   sha256 = strings(numel(artifact_files), 1);
   for k = 1:numel(artifact_files)
      if ~isfile(artifact_files(k))
         error('icemodel:verification:report:missingReportArtifact', ...
            'Cannot hash missing report artifact: %s', artifact_files(k))
      end
      info = dir(artifact_files(k));
      bytes(k) = info.bytes;
      sha256(k) = ...
         icemodel.verification.setup.fileSha256(artifact_files(k));
   end
   if sha256(1) ~= source_sha256
      error('icemodel:verification:report:changedAblationResults', ...
         'Saved ablation results changed while the report was generated: %s', ...
         results_file)
   end
   manifest = table(artifact_role, artifact_path, bytes, sha256);
   writetable(manifest, manifest_file)
end

function synthesis = synthesisTable(summary)
   %SYNTHESISTABLE Normalize coverage and preserve every outcome category.

   synthesis = summary;
   vars = string(summary.Properties.VariableNames);
   possible = numericColumn(summary, ...
      ["window_possible_sample_count", "possible_count"]);
   valid = numericColumn(summary, ...
      ["window_valid_sample_count", "valid_count", "eligible_sample_count"]);
   direct = numericColumn(summary, ...
      ["window_direct_sample_count", "direct_target_count"]);
   coverage = numericColumn(summary, "coverage_fraction");
   completed = string(summary.status) == "completed";
   derive = ~isfinite(coverage) & completed & possible > 0 ...
      & isfinite(valid);
   coverage(derive) = valid(derive) ./ possible(derive);
   synthesis.window_possible_sample_count = possible;
   synthesis.window_valid_sample_count = valid;
   synthesis.window_direct_sample_count = direct;
   synthesis.coverage_fraction = coverage;
   if ismember("coverage_denominator", vars)
      synthesis.coverage_denominator = string(summary.coverage_denominator);
   else
      synthesis.coverage_denominator = repmat( ...
         "valid exact-common samples / possible timestamps in selected window", ...
         height(summary), 1);
   end

   % Interpretation state supports the scientific synthesis, while report
   % category keeps readiness, selection, execution, and scientific failures
   % distinct for operations and data-quality diagnosis.
   state = repmat("unavailable", height(summary), 1);
   status = string(summary.status);
   classification = string(summary.classification);
   non_identifiable = completed & (~logical(summary.physical_comparable) ...
      | ismember(classification, ...
      ["non_identifiable", "not_physically_comparable"]));
   state(completed & ~non_identifiable) = "interpretable";
   state(non_identifiable) = "non_identifiable";
   state(status == "excluded") = "excluded";
   synthesis.interpretation_state = state;

   category = repmat("unknown_status", height(summary), 1);
   category(completed & ~non_identifiable) = "interpretable";
   category(non_identifiable) = "scientifically_unavailable";
   category(status == "excluded") = "readiness_excluded";
   category(status == "not_selected") = "not_selected";
   category(status == "unavailable") = "execution_unavailable";
   synthesis.report_category = category;

   % Preserve the broad report category for operational roll-ups while adding
   % an exact outcome category for machine-readable summaries and figures.
   outcome_category = category;
   unavailable = category == "scientifically_unavailable";
   outcome_category(unavailable & classification == "non_identifiable") = ...
      "non_identifiable";
   outcome_category(unavailable ...
      & classification == "not_physically_comparable") = ...
      "not_physically_comparable";
   outcome_category(unavailable ...
      & ~ismember(classification, ...
      ["non_identifiable", "not_physically_comparable"])) = ...
      "scientifically_unavailable_other";
   synthesis.outcome_category = outcome_category;
end

function values = numericColumn(table_value, candidates)
   %NUMERICCOLUMN Return the first saved candidate or a NaN column.

   vars = string(table_value.Properties.VariableNames);
   selected = candidates(find(ismember(candidates, vars), 1));
   if isempty(selected)
      values = NaN(height(table_value), 1);
   else
      values = double(table_value.(selected));
   end
end

function summary = statusTable(synthesis)
   %STATUSTABLE Aggregate outcome, status, reason, and error id.

   if height(synthesis) == 0
      summary = table(strings(0, 1), strings(0, 1), strings(0, 1), ...
         strings(0, 1), strings(0, 1), strings(0, 1), zeros(0, 1), ...
         'VariableNames', {'report_category', 'outcome_category', ...
         'classification', 'status', 'reason', 'error_identifier', 'count'});
      return
   end
   vars = string(synthesis.Properties.VariableNames);
   reason = repmat("", height(synthesis), 1);
   error_identifier = repmat("", height(synthesis), 1);
   if ismember("reason", vars)
      reason = string(synthesis.reason);
   end
   if ismember("error_identifier", vars)
      error_identifier = string(synthesis.error_identifier);
   end
   keys = [string(synthesis.report_category), ...
      string(synthesis.outcome_category), string(synthesis.classification), ...
      string(synthesis.status), reason, error_identifier];
   groups = unique(keys, 'rows', 'stable');
   count = zeros(size(groups, 1), 1);
   for k = 1:size(groups, 1)
      count(k) = nnz(all(keys == groups(k, :), 2));
   end
   summary = table(groups(:, 1), groups(:, 2), groups(:, 3), groups(:, 4), ...
      groups(:, 5), groups(:, 6), count, 'VariableNames', ...
      {'report_category', 'outcome_category', 'classification', 'status', ...
      'reason', 'error_identifier', 'count'});
end

function support = supportTable(site_results, policy)
   %SUPPORTTABLE Classify every saved observation posting used by the report.

   parts = cell(numel(site_results), 1);
   n_parts = 0;
   for k = 1:numel(site_results)
      result = site_results(k);
      if ~isfield(result, 'status') || string(result.status) ~= "completed"
         continue
      end
      n_parts = n_parts + 1;
      parts{n_parts} = observationSupport(result, policy);
   end
   support = concatenateTableParts(emptySupportTable(), parts, n_parts);
end

function value = emptySupportTable()
   %EMPTYSUPPORTTABLE Return the zero-row observation-support schema.

   value = table(strings(0, 1), strings(0, 1), zeros(0, 1), ...
      NaT(0, 1, 'TimeZone', 'UTC'), false(0, 1), zeros(0, 1), zeros(0, 1), ...
      zeros(0, 1), strings(0, 1), strings(0, 1), 'VariableNames', ...
      {'case_id', 'site_id', 'year', 'Time', 'in_comparison_window', ...
      'observation_ablation_m', 'observation_lowering_m', 'snow_depth_m', ...
      'support_class', 'support_reason'});
end

function components = componentTable(site_results)
   %COMPONENTTABLE Sum signed physical solid-mass components from saved ledgers.

   empty_components = table(strings(0, 1), strings(0, 1), zeros(0, 1), ...
      zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
      'VariableNames', {'case_id', 'site_id', 'year', ...
      'phase_melt_solid_change_mwe', ...
      'phase_refreezing_solid_change_mwe', 'solid_vapor_change_mwe', ...
      'net_physical_solid_change_mwe'});
   parts = cell(numel(site_results), 1);
   n_parts = 0;
   for k = 1:numel(site_results)
      result = site_results(k);
      if ~isfield(result, 'status') || string(result.status) ~= "completed"
         continue
      end
      t0 = result.comparison.window_start;
      t1 = result.comparison.window_end;
      ledger = result.model(result.model.Time >= t0 ...
         & result.model.Time < t1, :);
      required = ["mass_budget_phase_solid_mwe", ...
         "mass_budget_vapor_solid_mwe"];
      if isempty(ledger) || ~all(ismember(required, ...
            string(ledger.Properties.VariableNames)))
         error('icemodel:verification:report:missingAblationComponents', ...
            'Completed case %s lacks saved physical component ledgers', ...
            result.case_id)
      end

      % Negative phase change is melt, positive phase change is refreezing, and
      % vapor exchange retains its signed deposition/sublimation convention.
      phase = ledger.mass_budget_phase_solid_mwe;
      vapor = ledger.mass_budget_vapor_solid_mwe;
      row = table(string(result.case_id), string(result.site_id), ...
         double(result.year), sum(min(phase, 0)), sum(max(phase, 0)), ...
         sum(vapor), sum(phase + vapor), 'VariableNames', ...
         empty_components.Properties.VariableNames);
      n_parts = n_parts + 1;
      parts{n_parts} = row;
   end
   components = concatenateTableParts(empty_components, parts, n_parts);
end

function values = gridTranslationTable(site_results)
   %GRIDTRANSLATIONTABLE Export interval-resolved H_grid event provenance.

   empty_values = table(strings(0, 1), strings(0, 1), zeros(0, 1), ...
      NaT(0, 1, 'TimeZone', 'UTC'), NaT(0, 1, 'TimeZone', 'UTC'), ...
      zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
      zeros(0, 1), zeros(0, 1), 'VariableNames', ...
      {'case_id', 'site_id', 'year', 'interval_start', 'interval_end', ...
      'top_deletion_count', 'top_deletion_height_m', ...
      'interior_merge_count', 'cumulative_top_deletion_count', ...
      'cumulative_top_deletion_height_m', ...
      'cumulative_interior_merge_count'});
   parts = cell(numel(site_results), 1);
   n_parts = 0;
   for k = 1:numel(site_results)
      result = site_results(k);
      if ~isfield(result, 'status') || string(result.status) ~= "completed"
         continue
      end
      n_parts = n_parts + 1;
      parts{n_parts} = gridTranslationRows(result, ...
         empty_values.Properties.VariableNames);
   end
   values = concatenateTableParts(empty_values, parts, n_parts);
end

function rows = gridTranslationRows(result, variable_names)
   %GRIDTRANSLATIONROWS Reconstruct one saved case's quantized geometry ledger.

   t0 = result.comparison.window_start;
   t1 = result.comparison.window_end;
   ledger = result.model(result.model.Time >= t0 ...
      & result.model.Time < t1, :);
   required = ["mass_budget_top_deletion_count", ...
      "mass_budget_top_deletion_height_m", ...
      "mass_budget_interior_merge_count"];
   if isempty(ledger) || ~all(ismember(required, ...
         string(ledger.Properties.VariableNames)))
      error('icemodel:verification:report:missingGridTranslationLedger', ...
         'Completed case %s lacks saved grid-translation event ledgers', ...
         result.case_id)
   end

   % Each saved model row starts one interval. The final end is the exclusive
   % comparison boundary, so cumulative H_grid can be reproduced without model
   % reruns or inference from the unrelated remeshing mass exchange.
   interval_start = ledger.Time;
   interval_end = [ledger.Time(2:end); t1];
   top_count = ledger.mass_budget_top_deletion_count;
   top_height = ledger.mass_budget_top_deletion_height_m;
   interior_count = ledger.mass_budget_interior_merge_count;
   n_rows = height(ledger);
   rows = table(repmat(string(result.case_id), n_rows, 1), ...
      repmat(string(result.site_id), n_rows, 1), ...
      repmat(double(result.year), n_rows, 1), interval_start, interval_end, ...
      top_count, top_height, interior_count, cumsum(top_count), ...
      cumsum(top_height), cumsum(interior_count), ...
      'VariableNames', variable_names);
end

function values = initializationTable(site_results)
   %INITIALIZATIONTABLE Preserve executed model-boundary provenance.

   empty_values = table(strings(0, 1), strings(0, 1), zeros(0, 1), ...
      strings(0, 1), NaT(0, 1, 'TimeZone', 'UTC'), strings(0, 1), ...
      NaT(0, 1, 'TimeZone', 'UTC'), NaT(0, 1, 'TimeZone', 'UTC'), ...
      NaT(0, 1, 'TimeZone', 'UTC'), 'VariableNames', ...
      {'case_id', 'site_id', 'year', 'status', 'initialization_start', ...
      'initialization_policy', 'evaluation_start', 'evaluation_end', ...
      'run_end_inclusive'});
   parts = cell(numel(site_results), 1);
   n_parts = 0;
   required = ["initialization_start", "initialization_policy", ...
      "evaluation_start", "evaluation_end", "run_end_inclusive"];
   for k = 1:numel(site_results)
      result = site_results(k);
      if ~isfield(result, 'model_metadata') ...
            || ~isstruct(result.model_metadata) ...
            || ~all(isfield(result.model_metadata, required))
         continue
      end
      metadata = result.model_metadata;
      row = table(string(result.case_id), string(result.site_id), ...
         double(result.year), string(result.status), ...
         metadata.initialization_start, string(metadata.initialization_policy), ...
         metadata.evaluation_start, metadata.evaluation_end, ...
         metadata.run_end_inclusive, 'VariableNames', ...
         empty_values.Properties.VariableNames);
      n_parts = n_parts + 1;
      parts{n_parts} = row;
   end
   values = concatenateTableParts(empty_values, parts, n_parts);
end

function values = kanFirstTable(values)
   %KANFIRSTTABLE Put KAN rows first without filtering any saved site.

   if height(values) == 0 || ~ismember("site_id", ...
         string(values.Properties.VariableNames))
      return
   end
   is_kan = startsWith(upper(string(values.site_id)), "KAN");
   [~, order] = sortrows([double(~is_kan), (1:height(values))'], [1, 2]);
   values = values(order, :);
end

function [identities, materiality, scenarios, density, exclusions] = ...
      diagnosticTables(site_results)
   %DIAGNOSTICTABLES Flatten completed per-site diagnostics with stable keys.

   identities = emptyIdentityTable();
   materiality = emptyMaterialityTable();
   scenarios = emptyScenarioTable();
   density = emptyDensityTable();
   exclusions = emptyExclusionTable();
   identity_parts = cell(numel(site_results), 1);
   materiality_parts = cell(numel(site_results), 1);
   scenario_parts = cell(numel(site_results), 1);
   density_parts = cell(numel(site_results), 1);
   exclusion_parts = cell(numel(site_results), 1);
   n_identity = 0;
   n_materiality = 0;
   n_scenario = 0;
   n_density = 0;
   n_exclusion = 0;
   for k = 1:numel(site_results)
      result = site_results(k);
      if ~isfield(result, 'status') || string(result.status) ~= "completed" ...
            || ~isfield(result, 'diagnostics') ...
            || ~isstruct(result.diagnostics)
         continue
      end

      % Each table retains case, site, and year so rows remain attributable
      % after all completed cases are concatenated.
      key = {string(result.case_id), string(result.site_id), ...
         double(result.year)};
      diagnostics = result.diagnostics;
      if isfield(diagnostics, 'identities') && istable(diagnostics.identities)
         n_identity = n_identity + 1;
         identity_parts{n_identity} = addKeys(diagnostics.identities, key);
      end
      if isfield(diagnostics, 'materiality') && istable(diagnostics.materiality)
         n_materiality = n_materiality + 1;
         materiality_parts{n_materiality} = addKeys( ...
            diagnostics.materiality, key);
      end
      if isfield(diagnostics, 'scenarios') && istable(diagnostics.scenarios)
         % A saved endpoint row cannot confer observation provenance on itself.
         % Keep its numeric sensitivity while enforcing the current comparator
         % naming in exported tables and report aggregation.
         scenario_values = diagnostics.scenarios;
         endpoint_rows = string(scenario_values.role) == "endpoint";
         scenario_values.credible(endpoint_rows) = false;
         n_scenario = n_scenario + 1;
         scenario_parts{n_scenario} = addKeys(scenario_values, key);
      end
      if isfield(diagnostics, 'effective_density') ...
            && istable(diagnostics.effective_density)
         n_density = n_density + 1;
         density_parts{n_density} = addKeys( ...
            diagnostics.effective_density, key);
      end
      if isfield(diagnostics, 'excluded') && isstruct(diagnostics.excluded)
         values = struct2table(diagnostics.excluded, 'AsArray', true);
         n_exclusion = n_exclusion + 1;
         exclusion_parts{n_exclusion} = addKeys(values, key);
      end
   end
   identities = concatenateTableParts(identities, identity_parts, n_identity);
   materiality = concatenateTableParts( ...
      materiality, materiality_parts, n_materiality);
   scenarios = concatenateTableParts(scenarios, scenario_parts, n_scenario);
   density = concatenateTableParts(density, density_parts, n_density);
   exclusions = concatenateTableParts( ...
      exclusions, exclusion_parts, n_exclusion);
end

function values = addKeys(values, key)
   %ADDKEYS Prepend stable case-year provenance to one diagnostic table.

   n_rows = height(values);
   keys = table(repmat(key{1}, n_rows, 1), repmat(key{2}, n_rows, 1), ...
      repmat(key{3}, n_rows, 1), ...
      'VariableNames', {'case_id', 'site_id', 'year'});
   values = [keys, values];
end

function values = concatenateTableParts(empty_value, parts, n_parts)
   %CONCATENATETABLEPARTS Join preallocated table fragments in stable order.

   values = empty_value;
   if n_parts > 0
      values = vertcat(values, parts{1:n_parts});
   end
end

function value = emptyIdentityTable()
   %EMPTYIDENTITYTABLE Return the zero-row closure schema.

   value = table(strings(0, 1), strings(0, 1), zeros(0, 1), ...
      strings(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), ...
      zeros(0, 1), false(0, 1), false(0, 1), zeros(0, 1), ...
      false(0, 1), 'VariableNames', ...
      {'case_id', 'site_id', 'year', 'identity', 'units', 'residual', ...
      'normalization', 'tolerance', 'window_passed', 'step_passed', ...
      'failed_step_count', 'passed'});
end

function value = emptyMaterialityTable()
   %EMPTYMATERIALITYTABLE Return the zero-row materiality schema.

   value = table(strings(0, 1), strings(0, 1), zeros(0, 1), ...
      strings(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
      zeros(0, 1), false(0, 1), 'VariableNames', ...
      {'case_id', 'site_id', 'year', 'channel', 'signed_net_mwe', ...
      'gross_mwe', 'signed_ratio', 'gross_ratio', 'material'});
end

function value = emptyScenarioTable()
   %EMPTYSCENARIOTABLE Return the zero-row accounting-scenario schema.

   value = table(strings(0, 1), strings(0, 1), zeros(0, 1), ...
      strings(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), ...
      zeros(0, 1), zeros(0, 1), strings(0, 1), false(0, 1), ...
      false(0, 1), false(0, 1), false(0, 1), 'VariableNames', ...
      {'case_id', 'site_id', 'year', 'scenario', 'role', 'model_mwe', ...
      'observation_mwe', 'difference_mwe', 'relative_difference', ...
      'classification', 'changes_sign', 'changes_classification', ...
      'material', 'credible'});
end

function value = emptyDensityTable()
   %EMPTYDENSITYTABLE Return the zero-row effective-density schema.

   value = table(strings(0, 1), strings(0, 1), zeros(0, 1), ...
      zeros(0, 1), zeros(0, 1), false(0, 1), 'VariableNames', ...
      {'case_id', 'site_id', 'year', 'effective_density_kg_m3', ...
      'observation_sensitivity_mwe', 'rigorous_bound'});
end

function value = emptyExclusionTable()
   %EMPTYEXCLUSIONTABLE Return the zero-row support-exclusion schema.

   value = table(strings(0, 1), strings(0, 1), zeros(0, 1), ...
      zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
      zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
      zeros(0, 1), zeros(0, 1), ...
      'VariableNames', ...
      {'case_id', 'site_id', 'year', 'gap_bridged', ...
      'station_transition', 'unresolved_step', ...
      'step_correctable_but_unresolved', 'nonfinite_observation', ...
      'unknown_quality_flag', 'unknown_snow_depth', 'snow_censored', ...
      'nonfinite_model', ...
      'total_unique_excluded'});
end

function assets = buildFigures(results, tables, asset_dir)
   %BUILDFIGURES Export only evidence available in the saved result.

   assets = struct('nested', "", 'endpoint', "", ...
      'site', strings(0, 1), 'site_caption', strings(0, 1), ...
      'components', "", 'closure', "", 'synthesis', "", ...
      'performance_scatter', "", 'performance_summary', "");
   completed = arrayfun(@(value) isfield(value, 'status') ...
      && string(value.status) == "completed", results.site_year_results);
   n_site = nnz(completed);
   assets.site = strings(n_site, 1);
   assets.site_caption = strings(n_site, 1);
   used_stems = strings(n_site, 1);
   site_index = 0;
   if height(tables.nested) > 0
      assets.nested = nestedFigure(tables.nested, asset_dir);
   end
   if height(tables.endpoint) > 0
      assets.endpoint = endpointFigure(tables.endpoint, asset_dir);
   end

   % Site figures precede synthesis and reconstruct the model only from saved
   % interval-start ledgers on the exact comparison window.
   for k = 1:numel(results.site_year_results)
      result = results.site_year_results(k);
      if isfield(result, 'status') && string(result.status) == "completed"
         site_index = site_index + 1;
         base_stem = safeFilename(string(result.case_id) + "-" ...
            + string(result.year));
         file_stem = base_stem;
         suffix = 1;
         while any(used_stems(1:site_index - 1) == file_stem)
            suffix = suffix + 1;
            file_stem = base_stem + "-" + string(suffix);
         end
         used_stems(site_index) = file_stem;
         [file, caption] = siteFigure( ...
            result, asset_dir, file_stem, results.policy);
         assets.site(site_index) = file;
         assets.site_caption(site_index) = caption;
      end
   end
   if height(tables.components) > 0
      assets.components = componentFigure(tables.components, asset_dir);
   end
   if height(tables.identities) > 0 || height(tables.materiality) > 0
      assets.closure = closureFigure(tables, asset_dir);
   end
   if height(tables.synthesis) > 0
      assets.synthesis = synthesisFigure(tables.synthesis, asset_dir);
   end

   % Model-versus-observation scoring is the primary deliverable, so both
   % performance views render whenever any case was scored.
   if any(tables.performance.scored)
      % Plot one density so the panels stay readable. The policy reference
      % density is the primary conversion; the other band values are a
      % sensitivity the verdict text quantifies.
      policy = icemodel.verification.namelists.promiceAblationPolicy();
      primary = policy.effective_density_reference_kg_m3;
      if ~any(tables.performance.density_kg_m3 == primary)
         primary = max(tables.performance.density_kg_m3);
      end
      performance = tables.performance( ...
         tables.performance.density_kg_m3 == primary, :);
      summary = tables.performance_summary( ...
         tables.performance_summary.density_kg_m3 == primary, :);
      assets.performance_scatter = performanceScatterFigure( ...
         performance, asset_dir);
      assets.performance_summary = performanceSummaryFigure( ...
         performance, summary, asset_dir);
   end
end

function file = performanceScatterFigure(performance, asset_dir)
   %PERFORMANCESCATTERFIGURE Plot endpoint model values against observations.

   % One panel per diagnostic keeps the shared 1:1 reference readable while the
   % panels stay directly comparable on identical axes.
   diagnostics = unique(performance.diagnostic, 'stable');
   n = numel(diagnostics);
   n_columns = min(2, n);
   n_rows = ceil(n / n_columns);
   file = string(fullfile(asset_dir, "performance-endpoint-scatter.png"));
   fig = icemodel.plot.newFigure(width=640 * n_columns, ...
      height=460 * n_rows);
   layout = tiledlayout(fig, n_rows, n_columns, ...
      TileSpacing='compact', Padding='compact');

   scored = performance(performance.scored, :);
   sites = unique(scored.site_id, 'stable');
   colors = stationColors(numel(sites));
   limit = max([scored.observation_endpoint_mwe; ...
      scored.model_endpoint_mwe; 0]);
   if ~isfinite(limit) || limit <= 0
      limit = 1;
   end
   limit = limit * 1.05;
   for k = 1:n
      ax = nexttile(layout);
      hold(ax, 'on')
      rows = scored.diagnostic == diagnostics(k);
      plot(ax, [0 limit], [0 limit], '--', Color=[0.4 0.4 0.4], ...
         LineWidth=1.2, DisplayName='1:1')
      for s = 1:numel(sites)
         in = rows & scored.site_id == sites(s);
         if ~any(in)
            continue
         end
         scatter(ax, scored.observation_endpoint_mwe(in), ...
            scored.model_endpoint_mwe(in), 44, colors(s, :), 'filled', ...
            MarkerFaceAlpha=0.75, DisplayName=sites(s));
      end
      icemodel.verification.report.formatReportAxes(ax)
      grid(ax, 'on')
      axis(ax, 'square')
      xlim(ax, [0 limit])
      ylim(ax, [0 limit])
      xlabel(ax, "Observed lowering (m w.e.)")
      ylabel(ax, "Modeled (m w.e.)")
      label = scored.label(find(rows, 1));
      if isempty(label)
         label = diagnostics(k);
      end
      title(ax, label, Interpreter='none')
      if k == 1
         legend(ax, Location='eastoutside')
      end
   end
   title(layout, ...
      "End-of-window cumulative model value versus measured lowering")
   icemodel.verification.report.exportAndClose(fig, file)
end

function file = performanceSummaryFigure(performance, summary, asset_dir)
   %PERFORMANCESUMMARYFIGURE Compare aggregate skill across all diagnostics.

   file = string(fullfile(asset_dir, "performance-summary.png"));
   fig = icemodel.plot.newFigure(width=1400, height=520);
   layout = tiledlayout(fig, 1, 3, ...
      TileSpacing='compact', Padding='compact');
   labels = summary.label;

   % Error magnitudes first: pooled RMSE weights long site-years correctly and
   % mean absolute error states the typical miss in the same units.
   ax = nexttile(layout);
   bar(ax, [summary.pooled_rmse_mwe, summary.mean_mae_mwe])
   configureDiagnosticAxis(ax, labels)
   grid(ax, 'on')
   ylabel(ax, "Error (m w.e.)")
   title(ax, "Pooled RMSE and mean absolute error")
   legend(ax, {'pooled RMSE', 'mean MAE'}, Location='best')

   % Signed endpoint error separates over-prediction from under-prediction,
   % which an absolute error cannot show.
   ax = nexttile(layout);
   bar(ax, summary.median_endpoint_error_mwe, FaceColor=[0.42 0.18 0.55])
   configureDiagnosticAxis(ax, labels)
   grid(ax, 'on')
   yline(ax, 0, '-', Color=[0.3 0.3 0.3], LineWidth=1.0)
   ylabel(ax, "Median endpoint error (m w.e.)")
   title(ax, "Signed endpoint error, positive is over-prediction")

   % Distribution of per-case endpoint error exposes spread that any single
   % aggregate statistic hides.
   ax = nexttile(layout);
   hold(ax, 'on')
   scored = performance(performance.scored, :);
   for k = 1:height(summary)
      rows = scored.diagnostic == summary.diagnostic(k);
      values = scored.endpoint_error_mwe(rows);
      if isempty(values)
         continue
      end
      % Spread the points deterministically. A random jitter would change the
      % exported bytes on every render and break the report hash manifest.
      offset = linspace(-0.28, 0.28, numel(values))';
      scatter(ax, k + offset, values, 26, [0.15 0.35 0.6], 'filled', ...
         MarkerFaceAlpha=0.45)
      plot(ax, k + [-0.3 0.3], median(values) * [1 1], '-', ...
         Color=[0.85 0.33 0.10], LineWidth=2.2)
   end
   configureDiagnosticAxis(ax, labels)
   grid(ax, 'on')
   yline(ax, 0, '-', Color=[0.3 0.3 0.3], LineWidth=1.0)
   ylabel(ax, "Per-site-year endpoint error (m w.e.)")
   title(ax, "Per-case spread; the orange bar is the median")
   icemodel.verification.report.exportAndClose(fig, file)
end

function file = nestedFigure(nested, asset_dir)
   %NESTEDFIGURE Plot cumulative-window stability for every saved site-year.

   use = logical(nested.available) ...
      & isfinite(nested.model_minus_observation_mwe);
   values = nested(use, :);
   file = string(fullfile(asset_dir, ...
      'ablation-nested-window-stability.png'));
   % Fixed dimensions keep the figure renderable for the full readiness-scale
   % nested cohort; the complete row inventory remains in the linked table.
   fig = icemodel.plot.newFigure(width=1100, height=650);
   ax = axes(fig);
   icemodel.verification.report.formatReportAxes(ax)
   hold(ax, 'on')
   group_key = string(nested.case_id) + " " + string(nested.year);
   groups = unique(group_key, 'stable');
   group_labels = strings(size(groups));
   plotted_groups = 0;
   for k = 1:numel(groups)
      first = find(group_key == groups(k), 1);
      group_labels(k) = string(nested.site_id(first)) + " " ...
         + string(nested.year(first));
   end
   colors = lines(max(numel(groups), 1));
   value_key = string(values.case_id) + " " + string(values.year);
   for k = 1:numel(groups)
      % A site-year enters the legend only through its saved available values;
      % unavailable windows remain explicit in the annotation and report table.
      select = value_key == groups(k);
      if ~any(select)
         continue
      end
      plotted_groups = plotted_groups + 1;
      x = double(values.target_days(select));
      y = double(values.model_minus_observation_mwe(select));
      [x, order] = sort(x);
      handle_visibility = 'off';
      if plotted_groups <= 12
         handle_visibility = 'on';
      end
      plot(ax, x, y(order), '-o', LineWidth=1.6, MarkerSize=6, ...
         Color=colors(k, :), ...
         DisplayName=icemodel.verification.report.safeLabel(group_labels(k)), ...
         HandleVisibility=handle_visibility)
   end

   % Empty axes and annotations distinguish unavailable evidence from a genuine
   % zero model-observation difference without inventing plot coordinates.
   if isempty(values)
      xlim(ax, [0, 1])
      ylim(ax, [-1, 1])
      text(ax, 0.5, 0, "No available nested-window values", ...
         HorizontalAlignment='center', FontWeight='bold', Color='k', ...
         FontSize=10)
   else
      yline(ax, 0, 'k-', HandleVisibility='off')
      if plotted_groups <= 12
         legend(ax, Location='best', Interpreter='none')
      else
         text(ax, 0.99, 0.99, string(plotted_groups) ...
            + " plotted site-years; see table for labels", ...
            Units='normalized', HorizontalAlignment='right', ...
            VerticalAlignment='top', Interpreter='none', Color='k', ...
            FontSize=10)
      end
   end
   unavailable = ~logical(nested.available) ...
      | ~isfinite(nested.model_minus_observation_mwe);
   unavailable_groups = group_labels(~ismember(groups, unique(value_key, 'stable')));
   note = "Saved unavailable windows: " + string(nnz(unavailable));
   if ~isempty(unavailable_groups)
      note = note + "; all-unavailable site-years: " ...
         + string(numel(unavailable_groups)) + " (see table)";
   end
   text(ax, 0.01, 0.01, icemodel.verification.report.safeLabel(note), ...
      Units='normalized', ...
      VerticalAlignment='bottom', Interpreter='none', Color='k', ...
      FontSize=10)
   grid(ax, 'on')
   xlabel(ax, "Cumulative window duration (days)")
   ylabel(ax, "Model minus intact-ice observation (m w.e.)")
   title(ax, "Nested-window stability across saved site-years")
   icemodel.verification.report.exportAndClose(fig, file)
end

function file = endpointFigure(endpoint, asset_dir)
   %ENDPOINTFIGURE Plot planned start and end perturbations for every site.

   file = string(fullfile(asset_dir, ...
      'ablation-endpoint-perturbation-stability.png'));
   fig = icemodel.plot.newFigure(width=1200, height=620);
   layout = tiledlayout(fig, 1, 2, ...
      TileSpacing='compact', Padding='compact');
   plotEndpointPanel(nexttile(layout), endpoint, ...
      "start", "Later comparison start", ...
      "Start offset (days)")
   plotEndpointPanel(nexttile(layout), endpoint, ...
      "end", "Earlier comparison end", ...
      "End offset (days)")
   title(layout, "Endpoint-perturbation stability; KAN first, all saved sites", ...
      Color='k')
   icemodel.verification.report.exportAndClose(fig, file)
end

function plotEndpointPanel(ax, endpoint, axis_name, panel_title, x_label)
   %PLOTENDPOINTPANEL Draw one set of window offsets, leaving gaps as gaps.

   offset = double(endpoint.perturbation_days_signed);
   represented = string(endpoint.perturbation_axis) == axis_name;
   available = represented & logical(endpoint.available) ...
      & isfinite(endpoint.model_minus_observation_mwe);
   group_key = string(endpoint.case_id) + " " + string(endpoint.year);
   groups = unique(group_key, 'stable');
   colors = lines(max(numel(groups), 1));
   hold(ax, 'on')
   plotted = 0;
   for k = 1:numel(groups)
      use = available & group_key == groups(k);
      if ~any(use)
         continue
      end
      plotted = plotted + 1;
      x = offset(use);
      y = endpoint.model_minus_observation_mwe(use);
      [x, order] = sort(x);
      visibility = 'off';
      if plotted <= 12
         visibility = 'on';
      end
      first = find(group_key == groups(k), 1);
      label = string(endpoint.site_id(first)) + " " ...
         + string(endpoint.year(first));
      plot(ax, x, y(order), '-o', Color=colors(k, :), ...
         LineWidth=1.5, MarkerSize=5, ...
         DisplayName=icemodel.verification.report.safeLabel(label), ...
         HandleVisibility=visibility)
   end
   if plotted == 0
      text(ax, 0.5, 0.5, "No available perturbation values", ...
         Units='normalized', HorizontalAlignment='center', Color='k', ...
         FontSize=10)
   else
      yline(ax, 0, 'k-', HandleVisibility='off')
      if plotted <= 12
         legend(ax, Location='best', Interpreter='none')
      end
   end
   unavailable = represented & (~logical(endpoint.available) ...
      | ~isfinite(endpoint.model_minus_observation_mwe));
   text(ax, 0.01, 0.01, "Saved unavailable rows: " ...
      + string(nnz(unavailable)) + "; see table", Units='normalized', ...
      VerticalAlignment='bottom', Interpreter='none', Color='k', ...
      FontSize=10)
   icemodel.verification.report.formatReportAxes(ax)
   grid(ax, 'on')
   xlabel(ax, x_label)
   ylabel(ax, "Model minus intact-ice observation (m w.e.)")
   title(ax, panel_title)
end

function [file, caption] = siteFigure(result, asset_dir, file_stem, policy)
   %SITEFIGURE Compare cumulative ablation diagnostics for one season.
   %
   % The figure carries two panels. The upper panel keeps the full June-to-
   % October season and marks the selected evaluation window. The lower panel
   % crops to that window and rebases every series at its first row. Both panels
   % use metres water equivalent on a single axis; a geometric metre axis in the
   % same frame would put grid translation and mass on one scale.

   required = ["observation_lowering_m", "observation_lower_mwe", ...
      "observation_upper_mwe", "observation_reference_mwe", ...
      "model_melt_mwe", "model_runoff_mwe", ...
      "model_freeze_mwe", "model_net_solid_loss_mwe", ...
      "model_layer_change_mwe", "model_surface_mass_loss_mwe", ...
      "model_ablation_proxy_mwe", ...
      "snow_depth_m", "ice_exposed", "direct_observation", ...
      "evaluation_window"];
   if ~isfield(result, 'seasonal') || ~istimetable(result.seasonal) ...
         || isempty(result.seasonal) || ~all(ismember(required, ...
         string(result.seasonal.Properties.VariableNames)))
      error('icemodel:verification:report:missingAblationSeasonal', ...
         'Completed case %s lacks the saved seasonal comparison fields', ...
         result.case_id)
   end

   % The scientific display is fixed to the common melt-season interval. Every
   % comparison curve is broken unless the saved row has exposed ice and direct
   % finite observation support.
   [season_start, season_end] = ...
      icemodel.verification.helpers.evaluationSeason( ...
      double(result.year), policy);
   seasonal = sortrows(result.seasonal);
   season_rows = seasonal.Time >= season_start ...
      & seasonal.Time <= season_end;
   if ~any(season_rows)
      error('icemodel:verification:report:missingAblationSeasonal', ...
         'Completed case %s has no saved rows in the policy season', ...
         result.case_id)
   end
   evaluation_highlight = season_rows ...
      & logical(seasonal.evaluation_window);
   direct_posting = logical(seasonal.direct_observation);
   [season_ice_exposed, season_snow_censored, season_unknown_snow] = ...
      icemodel.verification.helpers.classifySnowDepth( ...
      seasonal.snow_depth_m, policy.ice_exposure_threshold_m);
   ice_support = season_rows & direct_posting & season_ice_exposed;
   direct_support = ice_support ...
      & isfinite(seasonal.observation_lowering_m);
   snow_censored = season_rows & direct_posting & season_snow_censored;
   unknown_snow = season_rows & direct_posting & season_unknown_snow;
   missing_or_quality = season_rows & ~direct_posting;
   unknown_or_missing = unknown_snow | missing_or_quality;

   % Insert NaNs at unsupported rows so neither the observation band nor any
   % modeled curve visually bridges snow-covered intervals. The signed physical
   % phase-plus-vapor diagnostic is never clamped or made monotonic.
   series = seasonalSeriesSet(seasonal, direct_support);

   file = string(fullfile(asset_dir, ...
      file_stem + "-cumulative-comparison.png"));
   fig = icemodel.plot.newFigure(width=1300, height=900);
   layout = tiledlayout(fig, 2, 1, ...
      TileSpacing='compact', Padding='compact');

   % Full-season panel: rebased at each series' first plotted row, so the whole
   % June-to-October record is visible and the selected window is marked in it.
   ax = nexttile(layout);
   full_series = rebaseSeasonalSeries(series, direct_support);
   plotSeasonalPanel(ax, seasonal.Time, full_series, ...
      evaluation_highlight, snow_censored, unknown_or_missing, ...
      season_start, season_end, direct_support, policy)
   title(ax, icemodel.verification.report.safeLabel(result.site_id) ...
      + " " + string(result.year) ...
      + ": full season, pale blue marks the selected window", ...
      Interpreter='none')
   legend(ax, Location='best', Interpreter='tex')

   % Selected-window panel: the same diagnostics cropped and rebased at the
   % first row of the window, so a season that begins under snow does not carry
   % a pre-window offset into the comparison the metrics actually score.
   ax = nexttile(layout);
   window_support = direct_support & evaluation_highlight;
   window_series = rebaseSeasonalSeries(series, window_support);
   window_series = maskSeriesSet(window_series, evaluation_highlight);
   [window_start, window_end] = seasonalWindowLimits( ...
      seasonal.Time, evaluation_highlight, season_start, season_end);
   plotSeasonalPanel(ax, seasonal.Time, window_series, ...
      false(size(evaluation_highlight)), ...
      snow_censored & evaluation_highlight, ...
      unknown_or_missing & evaluation_highlight, ...
      window_start, window_end, window_support, policy)
   title(ax, "Selected snow-free evaluation window, rebased at its first row")
   icemodel.verification.report.exportAndClose(fig, file)

   % The caption names each quantity and states the central runoff limitation.
   caption = icemodel.verification.report.safeLabel(result.site_id) ...
      + " " + string(result.year) ...
      + ". Upper panel: the full " + seasonRangeText(policy) ...
      + "with the selected evaluation window shaded pale blue. Lower panel: " ...
      + "the same diagnostics cropped to that window and rebased at its " ...
      + "first row. Both panels use metres water equivalent only. Green " ...
      + "shading converts measured surface lowering with " ...
      + densityRangeText(policy) + " kg m^-3. Lines show the legacy " ...
      + "cumulative melt and refreezing diagnostics, the " ...
      + residenceWindowText(result.model_options) + " " ...
      + "runoff diagnostic, runoff plus vapor loss, and the signed " ...
      + "net solid mass balance from modeled phase and vapor terms " ...
      + "(excluding remeshing and domain exchange). The signed balance can " ...
      + "decrease when internal refreezing exceeds melt, so the runoff " ...
      + "diagnostics are the surface-lowering comparators. Cumulative merge " ...
      + "export is not plotted: merges assign the joined cell the mean of " ...
      + "the pair, so it over-counts what the top cell held. " ...
      + "Unsupported rows are left as gaps. " ...
      + "Finite snow-covered observation postings are shaded gray; rows with " ...
      + "unknown snow or no direct-quality observation are shaded amber and " ...
      + "are not classified as snow-censored.";
end

function series = seasonalSeriesSet(seasonal, support)
   %SEASONALSERIESSET Collect the plotted metres-water-equivalent diagnostics.

   % Producers store ordered bounds; normalize again so legacy saved artifacts
   % cannot invert the plotted polygon at negative signed-lowering rows.
   bound_a = maskSeasonalSeries(seasonal.observation_lower_mwe, support);
   bound_b = maskSeasonalSeries(seasonal.observation_upper_mwe, support);
   series = struct( ...
      'observation_lower', min(bound_a, bound_b), ...
      'observation_upper', max(bound_a, bound_b), ...
      'observation_reference', ...
      maskSeasonalSeries(seasonal.observation_reference_mwe, support), ...
      'model_melt', maskSeasonalSeries(seasonal.model_melt_mwe, support), ...
      'model_runoff', maskSeasonalSeries(seasonal.model_runoff_mwe, support), ...
      'model_freeze', maskSeasonalSeries(seasonal.model_freeze_mwe, support), ...
      'model_balance', ...
      maskSeasonalSeries(seasonal.model_net_solid_loss_mwe, support), ...
      'model_ablation_proxy', ...
      maskSeasonalSeries(seasonal.model_ablation_proxy_mwe, support));
end

function series = rebaseSeasonalSeries(series, support)
   %REBASESEASONALSERIES Zero every cumulative series at its first plotted row.

   % Each series subtracts its own value at the reference row. Both density
   % edges of the observation band are the same measured lowering scaled by a
   % different density, so subtracting each edge's own reference leaves both
   % edges at zero there and keeps every later value equal to the rebased
   % lowering times that edge's density. A shared reference would instead
   % offset one edge by the other's density and distort the band.
   names = fieldnames(series);
   reference = find(support, 1);
   if isempty(reference)
      return
   end
   for n = 1:numel(names)
      values = series.(names{n});
      if ~isfinite(values(reference))
         continue
      end
      series.(names{n}) = values - values(reference);
   end
end

function series = maskSeriesSet(series, keep)
   %MASKSERIESSET Restrict every plotted series to one contiguous interval.

   names = fieldnames(series);
   for n = 1:numel(names)
      values = series.(names{n});
      values(~keep) = NaN;
      series.(names{n}) = values;
   end
end

function [lower_limit, upper_limit] = seasonalWindowLimits( ...
      time, highlight, season_start, season_end)
   %SEASONALWINDOWLIMITS Return the plotted x-limits of the selected window.

   % Fall back to the full season when no window row survives, so the panel
   % still renders with a valid axis instead of an empty frame.
   if ~any(highlight)
      lower_limit = season_start;
      upper_limit = season_end;
      return
   end
   lower_limit = min(time(highlight));
   upper_limit = max(time(highlight));
   if upper_limit <= lower_limit
      upper_limit = lower_limit + hours(1);
   end
end

function plotSeasonalPanel(ax, time, series, evaluation_highlight, ...
      snow_censored, unknown_or_missing, lower_limit, upper_limit, ...
      support, policy)
   %PLOTSEASONALPANEL Draw one metres-water-equivalent cumulative panel.

   hold(ax, 'on')
   shadeSeasonalMask(ax, time, evaluation_highlight, ...
      lower_limit, upper_limit, [0.65 0.78 0.95], 0.12)
   shadeSeasonalMask(ax, time, snow_censored, ...
      lower_limit, upper_limit, [0.75 0.78 0.82], 0.25)
   shadeSeasonalMask(ax, time, unknown_or_missing, ...
      lower_limit, upper_limit, [0.95 0.72 0.30], 0.18)
   plotObservationBand(ax, time, series.observation_lower, ...
      series.observation_upper, support, densityRangeText(policy))
   plot(ax, time, series.observation_reference, '-', LineWidth=2.4, ...
      Color=[0.05 0.35 0.22], ...
      DisplayName="observed lowering at " ...
      + compose('%g', policy.effective_density_reference_kg_m3) ...
      + " kg m^{-3}")
   plot(ax, time, series.model_melt, '-', LineWidth=1.8, ...
      Color=[0.85 0.33 0.10], ...
      DisplayName='legacy cumulative melt diagnostic')
   plot(ax, time, series.model_runoff, '--', LineWidth=1.8, ...
      Color=[0.10 0.45 0.75], ...
      DisplayName='runoff (refreezing limited to recent melt)')
   plot(ax, time, series.model_ablation_proxy, '-', LineWidth=2.0, ...
      Color=[0.00 0.30 0.50], ...
      DisplayName='runoff + vapor loss (ablation proxy)')
   plot(ax, time, series.model_freeze, ':', LineWidth=1.7, ...
      Color=[0.12 0.55 0.32], ...
      DisplayName='legacy cumulative refreezing diagnostic')
   plot(ax, time, series.model_balance, '-.', LineWidth=2.0, ...
      Color=[0.42 0.18 0.55], ...
      DisplayName='signed net solid balance (phase + vapor, can decrease)')
   icemodel.verification.report.formatReportAxes(ax)
   grid(ax, 'on')
   xlim(ax, [lower_limit upper_limit])
   xlabel(ax, "UTC time")
   ylabel(ax, "Cumulative water equivalent (m w.e.)")

   % Zero is the common rebase origin, so it must stay visible even when every
   % plotted series is strictly positive.
   limits = ylim(ax);
   ylim(ax, [min(0, limits(1)), max(0, limits(2))])
end

function values = maskSeasonalSeries(values, support)
   %MASKSEASONALSERIES Break a plotted series at unsupported seasonal rows.

   values = double(values);
   values(~support | ~isfinite(values)) = NaN;
end

function plotObservationBand(ax, time, lower, upper, support, density_range)
   %PLOTOBSERVATIONBAND Shade each contiguous supported density-band segment.

   valid = support & isfinite(lower) & isfinite(upper);
   starts = find(valid & [true; ~valid(1:end - 1)]);
   ends = find(valid & [~valid(2:end); true]);
   for k = 1:numel(starts)
      index = starts(k):ends(k);
      x = [time(index); flipud(time(index))];
      y = [lower(index); flipud(upper(index))];
      band = fill(ax, x, y, [0.26 0.67 0.45], ...
         FaceAlpha=0.24, EdgeColor='none');
      if k == 1
         band.DisplayName = char( ...
            "observed " + density_range + " kg m^{-3} band");
      else
         band.HandleVisibility = 'off';
      end
   end
end

function shadeSeasonalMask(ax, time, mask, lower_limit, upper_limit, ...
      face_color, face_alpha)
   %SHADESEASONALMASK Shade contiguous seasonal rows with a supplied style.

   if ~any(mask)
      return
   end
   if isscalar(time)
      edges = [lower_limit; upper_limit];
   else
      midpoint = time(1:end - 1) + diff(time) / 2;
      edges = [lower_limit; midpoint; upper_limit];
   end
   starts = find(mask & [true; ~mask(1:end - 1)]);
   ends = find(mask & [~mask(2:end); true]);
   for k = 1:numel(starts)
      icemodel.plot.markTimeSpan(ax, edges(starts(k)), edges(ends(k) + 1), ...
         style="fill", color=face_color, face_alpha=face_alpha);
   end
end

function support = observationSupport(result, policy)
   %OBSERVATIONSUPPORT Classify saved source points under fixed policy rules.

   observations = result.observations;
   if istimetable(observations)
      data = observations;
   elseif isstruct(observations) && isfield(observations, 'data') ...
         && istimetable(observations.data)
      data = observations.data;
   else
      error('icemodel:verification:report:missingAblationObservations', ...
         'Completed results must retain their saved observation timetable')
   end
   % VALIDATERESULTS has already normalized these saved fields, including for
   % readiness-only runs; support classification only consumes it here.
   observation_field = string(policy.observation_field);
   snow_field = string(policy.snow_variable);
   ice_threshold = double(policy.ice_exposure_threshold_m);
   zero_fields = string(policy.direct_zero_flag_fields);
   datum_fields = string(policy.datum_break_flag_fields);
   gap_fields = string(policy.ordinary_gap_flag_fields);
   names = string(data.Properties.VariableNames);
   required = unique([ ...
      icemodel.verification.helpers.observationSupportFields( ...
      observation_field, policy), snow_field], 'stable');
   if ~all(ismember(required, names))
      error('icemodel:verification:report:missingAblationObservations', ...
         'Saved observations do not contain the policy-required fields')
   end

   % Classify every retained observation posting in the fixed display season;
   % the comparison-window label preserves the narrower evaluation selection.
   [season_start, season_end] = ...
      icemodel.verification.helpers.evaluationSeason( ...
      double(result.year), policy);
   use = data.Time >= season_start & data.Time <= season_end;
   data = data(use, :);
   t0 = result.comparison.window_start;
   t1 = result.comparison.window_end;
   in_comparison_window = data.Time >= t0 & data.Time <= t1;
   raw = data.(observation_field);
   snow_depth = data.(snow_field);
   if ~isnumeric(snow_depth) && ~islogical(snow_depth)
      error('icemodel:verification:report:invalidAblationSnowDepth', ...
         'Saved snow field %s must be numeric or logical', snow_field)
   end
   snow_depth = double(snow_depth);
   % The support rules read more than policy.support_flag_fields, so take the
   % set from observationSupportFields.
   quality_fields = setdiff( ...
      icemodel.verification.helpers.observationSupportFields( ...
      observation_field, policy), observation_field, 'stable');
   quality = zeros(height(data), numel(quality_fields));
   for k = 1:numel(quality_fields)
      values = data.(quality_fields(k));
      if ~isnumeric(values) && ~islogical(values)
         error('icemodel:verification:report:invalidAblationQualityFlag', ...
            'Saved quality field %s must be numeric or logical', ...
            quality_fields(k))
      end
      quality(:, k) = double(values);
   end
   % classifyObservationSupport applies every flag rule. The quality matrix is
   % still needed below to name which flags made a given posting unknown or
   % flagged, so it is built here and passed in.
   [~, zero_index] = ismember(zero_fields, quality_fields);
   support = icemodel.verification.helpers.classifyObservationSupport( ...
      [double(raw), quality], [observation_field, quality_fields], ...
      observation_field, policy);
   quality_finite = support.quality_finite;
   target_finite = support.target_finite;

   % Direct support requires finite zero-valued quality flags and exposed ice.
   % Finite flag-clean observations above the exposure threshold are retained
   % as snow-censored support rather than treated as direct ice lowering.
   nonzero_flag = ~support.direct_flags_zero;
   flag_clean = support.flag_clean;
   [ice_exposed, censored_snow, invalid_snow] = ...
      icemodel.verification.helpers.classifySnowDepth( ...
      snow_depth, ice_threshold);
   direct = flag_clean & ice_exposed;
   if ~any(direct)
      error('icemodel:verification:report:missingAblationObservations', ...
         'Saved completed observations contain no direct finite support')
   end
   snow_censored = flag_clean & censored_snow;
   unknown_snow = flag_clean & invalid_snow;
   flagged = target_finite & quality_finite & nonzero_flag;
   unknown = target_finite & ~quality_finite;

   % Unknown quality takes precedence over flag comparisons so NaN and Inf
   % never become direct support. Per-row reasons retain the fields that made a
   % posting unknown or explicitly flagged.
   support_class = repmat("nonfinite_target", height(data), 1);
   support_class(direct) = "direct";
   support_class(snow_censored) = "snow_censored";
   support_class(unknown_snow) = "unknown_snow";
   support_class(flagged) = "flagged";
   support_class(unknown) = "unknown_quality";
   support_reason = repmat("nonfinite observation ablation", height(data), 1);
   support_reason(direct) = "direct finite exposed-ice support";
   for row = find(snow_censored)'
      support_reason(row) = sprintf( ...
         "snow depth %.6g m exceeds ice-exposure threshold %.6g m", ...
         snow_depth(row), ice_threshold);
   end
   support_reason(unknown_snow) = ...
      "nonfinite required snow field: " + snow_field;
   negative_snow = unknown_snow & isfinite(snow_depth) & snow_depth < 0;
   support_reason(negative_snow) = ...
      "negative invalid snow field: " + snow_field;
   for row = find(flagged)'
      active = zero_fields(quality(row, zero_index) ~= 0);
      datum = active(ismember(active, datum_fields));
      gap = active(ismember(active, gap_fields));
      other = setdiff(active, [datum_fields, gap_fields], 'stable');
      parts = strings(3, 1);
      n_parts = 0;
      if ~isempty(datum)
         n_parts = n_parts + 1;
         parts(n_parts) = "datum break: " + join(datum, ", ");
      end
      if ~isempty(gap)
         n_parts = n_parts + 1;
         parts(n_parts) = "ordinary gap: " + join(gap, ", ");
      end
      if ~isempty(other)
         n_parts = n_parts + 1;
         parts(n_parts) = "other direct-zero flag: " ...
            + join(other, ", ");
      end
      support_reason(row) = "nonzero exclusion flag(s) [" ...
         + join(parts(1:n_parts), "; ") + "]";
   end
   for row = find(unknown)'
      nonfinite = quality_fields(~isfinite(quality(row, :)));
      support_reason(row) = "nonfinite required flag(s): " ...
         + join(nonfinite, ", ");
   end

   % Rebase geometric lowering to the first direct finite posting while retaining
   % excluded values in the machine-readable support ledger.
   lowering = raw - raw(find(direct, 1));
   n_rows = height(data);
   support = table(repmat(string(result.case_id), n_rows, 1), ...
      repmat(string(result.site_id), n_rows, 1), ...
      repmat(double(result.year), n_rows, 1), data.Time, ...
      in_comparison_window, raw, lowering, snow_depth, support_class, ...
      support_reason, 'VariableNames', ...
      {'case_id', 'site_id', 'year', 'Time', 'in_comparison_window', ...
      'observation_ablation_m', 'observation_lowering_m', 'snow_depth_m', ...
      'support_class', 'support_reason'});
end

function file = componentFigure(components, asset_dir)
   %COMPONENTFIGURE Plot signed physical solid-mass changes at fixed size.

   file = string(fullfile(asset_dir, 'ablation-signed-components.png'));
   fig = icemodel.plot.newFigure(width=1250, height=680);
   ax = axes(fig);
   hold(ax, 'on')
   x = (1:height(components))';
   plot(ax, x, components.phase_melt_solid_change_mwe, 'v', ...
      LineStyle='none', MarkerSize=7, Color=[0.78 0.19 0.19], ...
      DisplayName='phase melt (solid change)')
   plot(ax, x, components.phase_refreezing_solid_change_mwe, '^', ...
      LineStyle='none', MarkerSize=7, Color=[0.08 0.52 0.34], ...
      DisplayName='phase refreezing (solid change)')
   plot(ax, x, components.solid_vapor_change_mwe, 'd', ...
      LineStyle='none', MarkerSize=6, Color=[0.55 0.23 0.67], ...
      DisplayName='solid-vapor exchange')
   plot(ax, x, components.net_physical_solid_change_mwe, 'o', ...
      LineStyle='none', MarkerSize=6, Color=[0.10 0.32 0.64], ...
      DisplayName='net physical solid change')
   yline(ax, 0, 'k-', HandleVisibility='off')
   icemodel.verification.report.formatReportAxes(ax)
   grid(ax, 'on')
   xlabel(ax, "Completed site-year")
   ylabel(ax, "Signed solid-mass change (m w.e.)")
   title(ax, "Signed phase, refreezing, and solid-vapor components")
   if height(components) <= 24
      xticks(ax, x)
      xticklabels(ax, string(components.site_id) + " " ...
         + string(components.year))
      ax.TickLabelInterpreter = 'none';
   end
   legend(ax, Location='southoutside', Orientation='horizontal')
   icemodel.verification.report.exportAndClose(fig, file)
end

function file = closureFigure(tables, asset_dir)
   %CLOSUREFIGURE Plot scalable site closure and channel materiality summaries.

   file = string(fullfile(asset_dir, ...
      'ablation-closure-materiality.png'));
   fig = icemodel.plot.newFigure(width=1200, height=650);
   layout = tiledlayout(fig, 1, 2, ...
      TileSpacing='compact', Padding='compact');

   % Reduce arbitrarily many identity channels to one maximum window-closure
   % ratio per site-year. A categorical sentinel places any per-step failure
   % visibly above the limit because the saved diagnostics do not retain its
   % exact worst-step ratio; exact verdicts and failed-step counts remain in CSV.
   ax = nexttile(layout);
   hold(ax, 'on')
   identities = tables.identities;
   if height(identities) > 0
      window_ratio = abs(identities.residual) ./ identities.tolerance;
      step_failure = ~logical(identities.step_passed);
      step_failure_sentinel = 1.1;
      key = string(identities.case_id) + " " + string(identities.year);
      groups = unique(key, 'stable');
      maximum = NaN(numel(groups), 1);
      group_step_failure = false(numel(groups), 1);
      labels = strings(numel(groups), 1);
      for k = 1:numel(groups)
         use_group = key == groups(k);
         use = use_group & isfinite(window_ratio);
         if any(use)
            maximum(k) = max(window_ratio(use));
         end
         group_step_failure(k) = any(use_group & step_failure);
         if group_step_failure(k)
            if isfinite(maximum(k))
               maximum(k) = max(maximum(k), step_failure_sentinel);
            else
               maximum(k) = step_failure_sentinel;
            end
         end
         first = find(use_group, 1);
         labels(k) = string(identities.site_id(first)) + " " ...
            + string(identities.year(first));
      end
      finite = isfinite(maximum);
      window_only = finite & ~group_step_failure;
      if any(window_only)
         plot(ax, find(window_only), ...
            log10(max(maximum(window_only), eps)), 'o', ...
            LineStyle='none', MarkerFaceColor=[0.25 0.54 0.73], ...
            MarkerEdgeColor=[0.25 0.54 0.73], ...
            DisplayName='window ratio; all saved steps pass')
      end
      failed_steps = finite & group_step_failure;
      if any(failed_steps)
         plot(ax, find(failed_steps), ...
            log10(maximum(failed_steps)), 'x', LineStyle='none', ...
            MarkerSize=8, LineWidth=1.5, Color=[0.78 0.19 0.19], ...
            DisplayName='one or more saved forcing steps fail')
      end
      yline(ax, 0, 'r--', DisplayName='acceptance limit')
      icemodel.verification.report.formatReportAxes(ax)
      grid(ax, 'on')
       xlabel(ax, "Completed site-year")
      ylabel(ax, "log10(max window ratio; step-failure sentinel)")
      title(ax, "Window closure and per-step acceptance")
      if any(window_only) || any(failed_steps)
         legend(ax, Location='best')
      end
      if numel(groups) <= 24
         xticks(ax, 1:numel(groups))
         xticklabels(ax, icemodel.verification.report.safeLabel(labels))
         ax.TickLabelInterpreter = 'none';
      end
   else
      axis(ax, 'off')
      text(ax, 0.5, 0.5, 'No saved closure identities', ...
         HorizontalAlignment='center', Color='k', FontSize=10)
   end

   % Summarize each accounting channel across site-years so the display scales
   % with the fixed diagnostic vocabulary instead of the cohort size.
   ax = nexttile(layout);
   materiality = tables.materiality;
   if height(materiality) > 0
      channels = unique(string(materiality.channel), 'stable');
      maximum_signed = NaN(numel(channels), 1);
      maximum_gross = NaN(numel(channels), 1);
      for k = 1:numel(channels)
         use = string(materiality.channel) == channels(k);
         signed = abs(materiality.signed_ratio(use));
         gross = materiality.gross_ratio(use);
         if any(isfinite(signed))
            maximum_signed(k) = max(signed(isfinite(signed)));
         end
         if any(isfinite(gross))
            maximum_gross(k) = max(gross(isfinite(gross)));
         end
      end
      barh(ax, [maximum_signed, maximum_gross])
      channel_labels = replace(channels, "_", " ");
      channel_labels = replace(channel_labels, ...
         ["remesh", "merge delete", "cloned bottom", ...
         "condensation overflow", "unapplied vapor solid equivalent"], ...
         ["remeshing", "layer merge/delete", "cloned-bottom", ...
         "excess condensation", "unapplied vapor (solid equivalent)"]);
      icemodel.verification.report.configureCategoryAxis(ax, channel_labels)
      xlabel(ax, "Maximum ratio to comparison signal G (dimensionless)")
      title(ax, "Channel-level accounting materiality")
      legend(ax, ["absolute signed ratio", "gross ratio"], ...
         Location='best')
   else
      axis(ax, 'off')
      text(ax, 0.5, 0.5, 'No saved materiality diagnostics', ...
         HorizontalAlignment='center', Color='k', FontSize=10)
   end
   icemodel.verification.report.exportAndClose(fig, file)
end

function file = synthesisFigure(synthesis, asset_dir)
   %SYNTHESISFIGURE Compare completed selected site-years at fixed size.

   file = string(fullfile(asset_dir, ...
      'ablation-cross-site-synthesis.png'));
   vars = string(synthesis.Properties.VariableNames);
   use_completed = string(synthesis.status) == "completed";
   if ismember("selected", vars)
      use_completed = use_completed & logical(synthesis.selected);
   end
   synthesis = synthesis(use_completed, :);
   if height(synthesis) == 0
      file = "";
      return
   end

   % Operational rows remain in CSV evidence; the scientific graphic contains
   % only completed selected comparisons with real model and observation values.
   categories = ["interpretable", "non_identifiable", ...
      "not_physically_comparable", "scientifically_unavailable_other"];
   palette = categoryColors(categories);
   markers = categoryMarkers(categories);
   category = string(synthesis.outcome_category);
   observation = synthesis.observation_intact_mwe;
   model = synthesis.model_solid_loss_mwe;
   if ismember("model_minus_observation_mwe", ...
         string(synthesis.Properties.VariableNames))
      difference = synthesis.model_minus_observation_mwe;
   else
      difference = model - observation;
   end
   rows = (1:height(synthesis))';
   labels = icemodel.verification.report.safeLabel(string(synthesis.site_id) + " " ...
      + string(synthesis.year));

   % The two panels show only physical comparisons and their signed differences;
   % operational status inventories stay in the linked machine-readable tables.
   fig = icemodel.plot.newFigure(width=1200, height=620);
   layout = tiledlayout(fig, 1, 2, ...
      TileSpacing='compact', Padding='compact');
   ax = nexttile(layout);
   hold(ax, 'on')
   finite_comparison = isfinite(observation) & isfinite(model);
   for k = 1:numel(categories)
      use = category == categories(k) & finite_comparison;
      if any(use)
         plot(ax, observation(use), model(use), LineStyle='none', ...
             Marker=markers(k), MarkerSize=7, Color=palette(k, :), ...
             MarkerFaceColor=palette(k, :), ...
             DisplayName=replace(categories(k), "_", " "))
      end
   end
   if any(finite_comparison)
      limits = [observation(finite_comparison); model(finite_comparison)];
      lower = min(limits);
      upper = max(limits);
      padding = max(0.025, 0.05 * max(upper - lower, eps));
      plot(ax, [lower - padding, upper + padding], ...
         [lower - padding, upper + padding], 'k--', ...
         HandleVisibility='off')
      legend(ax, Location='best', Interpreter='none')
   else
      text(ax, 0.5, 0.5, "No saved finite model-observation pairs", ...
         Units='normalized', HorizontalAlignment='center', Color='k', ...
         FontSize=10)
   end
   icemodel.verification.report.formatReportAxes(ax)
   grid(ax, 'on')
   xlabel(ax, "Observed intact-ice conversion (m w.e.)")
   ylabel(ax, "Model phase + vapor solid loss (m w.e.)")
   title(ax, "Completed site-years: model versus observation")

   ax = nexttile(layout);
   hold(ax, 'on')
   for k = 1:numel(categories)
      use = category == categories(k) & isfinite(difference);
      if any(use)
         plot(ax, rows(use), difference(use), LineStyle='none', ...
            Marker=markers(k), MarkerSize=6, Color=palette(k, :), ...
            MarkerFaceColor=palette(k, :), HandleVisibility='off')
      end
   end
   if any(isfinite(difference))
      yline(ax, 0, 'k-', HandleVisibility='off')
   else
      text(ax, 0.5, 0.5, "No saved finite signed differences", ...
         Units='normalized', HorizontalAlignment='center', Color='k', ...
         FontSize=10)
   end
   icemodel.verification.report.formatReportAxes(ax)
   grid(ax, 'on')
   xlabel(ax, "Completed selected site-year")
   ylabel(ax, "Model minus observation (m w.e.)")
   title(ax, "Signed site-year difference")
   if height(synthesis) <= 24
      xticks(ax, rows)
      xticklabels(ax, labels)
      ax.TickLabelInterpreter = 'none';
   end
   icemodel.verification.report.exportAndClose(fig, file)
end

function colors = stationColors(n_stations)
   %STATIONCOLORS Return distinguishable colors for an arbitrary station list.

   % The cohort spans 17 stations, so the fixed outcome-category palette does
   % not apply. Walk the hue circle at fixed saturation and value instead, which
   % stays deterministic and separates neighbouring stations at any count.
   if n_stations < 1
      colors = zeros(0, 3);
      return
   end
   hue = ((0:n_stations - 1)' / n_stations);
   colors = hsv2rgb([hue, repmat(0.62, n_stations, 1), ...
      repmat(0.78, n_stations, 1)]);
end

function colors = categoryColors(categories)
   %CATEGORYCOLORS Return fixed colors for distinct report outcome categories.

   colors = zeros(numel(categories), 3);
   map = struct( ...
      'interpretable', [0.08 0.52 0.34], ...
      'non_identifiable', [0.83 0.52 0.10], ...
      'not_physically_comparable', [0.78 0.19 0.19], ...
      'scientifically_unavailable_other', [0.65 0.43 0.26]);
   for k = 1:numel(categories)
      colors(k, :) = map.(categories(k));
   end
end

function markers = categoryMarkers(categories)
   %CATEGORYMARKERS Return fixed symbols for distinct outcome categories.

   map = struct('interpretable', "o", ...
      'non_identifiable', "s", ...
      'not_physically_comparable', "x", ...
      'scientifically_unavailable_other', "d");
   markers = strings(size(categories));
   for k = 1:numel(categories)
      markers(k) = map.(categories(k));
   end
end

function lines = reportMarkdown(results, results_file, tables, files, ...
      assets, report_file, source_sha256, manifest_file)
   %REPORTMARKDOWN Build the scientific report from saved values only.

   generated = string(datetime('now', TimeZone='UTC', ...
      Format="yyyy-MM-dd HH:mm:ss 'UTC'"));
   [~, output_name, output_ext] = fileparts(report_file);
   counts = categoryCounts(tables.synthesis);
   n_completed = nnz(string(tables.summary.status) == "completed");
   closure_counts = closureSiteYearCounts(tables.identities);
   directional = directionalClassificationTable(tables.synthesis);
   drivers = accountingDriverTable(tables.scenarios);
   has_endpoint_deficit = hasEndpointDeficitRows(tables.scenarios);
   key_findings = keyFindingsText( ...
      directional, drivers, has_endpoint_deficit);
   scenario_interpretation = scenarioInterpretationText(has_endpoint_deficit);
   endpoint_diagram_line = endpointDeficitDiagramLine(has_endpoint_deficit);
   season_range = seasonRangeText(results.policy);
   density_range = densityRangeText(results.policy);
   density_values = results.policy.effective_density_kg_m3;
   density_lower = string(sprintf('%g', density_values(1)));
   density_upper = string(sprintf('%g', density_values(end)));

   % Static scientific prose is combined with sanitized saved identifiers and
   % exact evidence counts; no report text depends on live repository state.
   lines = [ ...
      "---"
      "title: ""PROMICE ablation-zone evaluation"""
      "date: """ + generated + """"
      "format:"
      "  html:"
      "    embed-resources: true"
      "    toc: true"
      "    toc-depth: 3"
      "output-file: """ + output_name + output_ext + """"
      "---"
      ""
      "## Structured Abstract"
      ""
       "**Background.** PROMICE instruments measure surface lowering, while " ...
         + "IceModel calculates melt, refreezing, vapor exchange, and changes " ...
         + "in stored mass. These quantities are related, but they are not " ...
         + "interchangeable over short periods when a weathering crust can grow " ...
         + "or decay."
       ""
       "**Methods.** The primary figures cover " + season_range + ". " ...
         + "Finite snow-covered postings are shaded gray and omitted from the " ...
         + "plotted curves; unknown snow or missing observations are shaded " ...
         + "amber instead. " ...
         + "Measured lowering is shown both geometrically and as a " ...
         + density_range + " kg m^-3 water-equivalent band. The legacy " ...
         + "cumulative melt and refreezing diagnostics, " ...
         + residenceWindowText( ...
         firstCompletedModelOptions(results)) + " " ...
         + "runoff diagnostic, runoff plus vapor loss, and net physical " ...
         + "solid loss from modeled phase and vapor terms (excluding " ...
         + "remeshing and domain exchange) are shown separately. " ...
         + "Cumulative merge export is NOT plotted " ...
         + "because it is a regridding quantity rather than a surface mass " ...
         + "flux. Every " ...
         + "plotted quantity is in metres water equivalent. Quantized " ...
         + "top-cell deletion height is grid geometry rather than mass and is " ...
         + "not a continuous surface prediction, so it is reported only in " ...
         + "the grid-translation ledger."
       ""
       "**Results.** " + string(n_completed) + " site-year(s) completed. " ...
          + string(closure_counts.passed) + " of " ...
          + string(closure_counts.evaluated) + " site-year(s) with saved " ...
          + "closure evidence passed every saved closure identity. " ...
          + string(counts.interpretable) + " yielded an identifiable " ...
          + "directional conclusion, while " ...
         + string(counts.scientifically_unavailable) + " did not (" ...
         + string(counts.non_identifiable) + " non-identifiable, " ...
         + string(counts.not_physically_comparable) ...
         + " not physically comparable, and " ...
         + string(counts.scientifically_unavailable_other) ...
         + " with another saved classification). " + key_findings
       ""
       "**Conclusions.** Measured lowering should be evaluated against the " ...
         + "collection of cumulative diagnostics, not against gross melt alone. " ...
         + "The runoff series is a postprocessed water budget for the pore " ...
         + "reservoir with a trailing refreezing residence limit, not a " ...
         + "modeled boundary flux: the prognostic column never drains. " ...
         + "Because the prognostic column retains that water, the signed " ...
         + "solid balance can fall late in a season, which is why the runoff " ...
         + "diagnostics rather than the signed balance are compared against " ...
         + "observed lowering."
      ""
       "## Executive Summary"
       ""
       "- Completed seasonal comparisons: **" + string(n_completed) + "**."
       "- Closure: **" + string(closure_counts.passed) + "/" ...
          + string(closure_counts.evaluated) + "** evaluated site-year(s) " ...
          + "passed every saved closure identity."
       "- Identifiability: **" + string(counts.interpretable) + "/" ...
          + string(n_completed) + "** completed site-year(s) yielded an " ...
          + "identifiable directional conclusion."
       "- Observed surface lowering is shown as a " + density_range ...
          + " kg m^-3 water-equivalent sensitivity band. The " + density_lower ...
          + " kg m^-3 density endpoint represents porous weathering-crust " ...
          + "material, not intact glacier ice; the two band edges are ordered " ...
          + "pointwise for signed lowering."
       "- Legacy melt, cumulative refreezing, the residence-limited " ...
         + "runoff diagnostic, and net physical solid loss from modeled phase " ...
         + "and vapor terms are kept separate because they answer different " ...
         + "physical questions; the last excludes remeshing and domain exchange."
       "- The runoff diagnostic is postprocessed, not a modeled boundary flux; " ...
         + "the prognostic column never removes that water during the " ...
         + "simulation. The budget credits refreezing only up to the liquid " ...
         + "supplied within a trailing residence window, and it accounts for " ...
         + "condensation, evaporation, and condensation that exceeded the top " ...
         + "cell's pore capacity."
       "- Cumulative merge export is the mass the remeshing step removed, in " ...
         + "metres water equivalent. It is NOT a surface mass flux: the merge " ...
         + "rule assigns the joined cell the mean of the pair, so the export " ...
         + "over-counts what the top cell actually held. The quantized cell " ...
         + "height is not convertible to mass and is never plotted as such; " ...
         + "it remains a separate grid-geometry counter."
      ""
      "## Data and Methods"
      ""
       "### Observations and comparison window"
      ""
       "Each seasonal figure uses the full saved " + season_range ...
         + " interval. " ...
         + "Pale-blue shading marks the selected longest evaluation window but " ...
         + "does not crop the seasonal record. Finite snow-covered observation " ...
         + "postings are shaded gray; unknown snow or missing observations are " ...
         + "shaded amber. " ...
         + "All comparison curves are shown only where exposed ice coincides " ...
         + "with direct finite observation support."
      ""
      "The observation-support CSV retains the full quality classification " ...
         + "for every posting in the fixed display season and labels whether " ...
         + "each row belongs to the narrower comparison window. Coverage and " ...
         + "excluded-row counts remain available there and in the site-year " ...
         + "synthesis CSV."
      ""
      "### Control-volume definition"
      ""
      "```{mermaid}"
      "flowchart LR"
      "  I[""PROMICE sonic ranger or pressure transducer""] --> D[""Instrument reference datum""]"
      "  D --> H[""Geometric surface lowering, m positive down""]"
      "  H --> C[""Intact-ice or effective-density conversion""]"
       endpoint_diagram_line
      "  P[""Phase exchange: melt or refreezing""] --> S[""Modeled solid phase term""]"
      "  V[""Solid-vapor exchange: sublimation or deposition""] --> Q[""Modeled solid-vapor term""]"
      "  P --> L[""Liquid storage, reported separately""]"
      "  S --> A[""Net physical solid loss: phase + vapor only, m w.e.""]"
      "  Q --> A"
      "  C --> O[""Observed lowering conversion, m w.e.""]"
      "  A --> X[""Model-observation comparison""]"
      "  O --> X"
      "  R[""Remeshing: clone minus merge export""] --> G[""Secondary grid geometry""]"
      "  R -. ""closure only; excluded from plotted loss"" .-> A"
      "```"
      ""
       "The legacy cumulative melt diagnostic includes all modeled phase " ...
         + "change. Saved cumulative refreezing is plotted explicitly so its " ...
         + "relationship to melt and runoff remains visible. The runoff " ...
         + "series is postprocessed rather than a modeled boundary flux, but " ...
         + "it is a genuine water budget for the pore reservoir. Per step the " ...
         + "reservoir gains melt and condensation and loses refreezing, " ...
         + "evaporation, and runoff, while condensation the top cell could " ...
         + "not store runs off directly. Refreezing is credited only up to " ...
         + "the liquid supplied within a trailing residence window, so water " ...
         + "older than that window is treated as having already left. That " ...
         + "residence limit is what " ...
         + "prevents the physically impossible case in which the column's " ...
         + "entire accumulated meltwater refreezes when the melt season ends, " ...
         + "and it also limits how much of a day's melt can refreeze " ...
         + "overnight. Runoff therefore tracks the signed solid balance " ...
         + "closely until refreezing becomes large enough to exceed that " ...
         + "limit, typically in late August or September, after which the " ...
         + "signed balance credits refreezing that runoff does not. " ...
         + "Runoff plus vapor loss adds only the SOLID vapor exchange, that " ...
         + "is sublimation minus deposition, to give the continuous ablation " ...
         + "proxy. Liquid vapor exchange is excluded from that proxy because " ...
         + "runoff already accounts for it: evaporated pore water is water " ...
         + "runoff would otherwise have carried, so adding it again would " ...
         + "count the same mass twice. The signed net solid balance is calculated " ...
         + "from the modeled phase and solid-vapor terms only; it excludes " ...
         + "remeshing and domain exchange. It is a mass-balance term and not " ...
         + "a surface-lowering proxy: it credits ALL internal refreezing, " ...
         + "including refreezing of meltwater that would already have drained " ...
         + "at a bare-ice site, so refreezing subtracts from it and the curve " ...
         + "falls. Refreezing at depth does not raise the surface, which is " ...
         + "why a falling balance must never be read as surface rise. It can " ...
         + "increase or decrease as melt, " ...
         + "refreezing, and vapor exchange act, and the report plots that saved " ...
         + "series without clipping or enforcing monotonicity. Because a " ...
         + "signed balance that decreases cannot be read as surface lowering, " ...
         + "the runoff diagnostics are the surface-lowering comparators. " ...
         + "Cumulative merge export is not a surface mass flux and is " ...
         + "neither scored nor plotted. A merge gives the joined cell the " ...
         + "MEAN of the two cells it replaces, so removing a nearly empty " ...
         + "top cell still exports about half the pair's mass." ...
         + mergeExportScaleText(results) ...
         + " It is retained in the mass ledger for closure checking. The " ...
         + "quantized cell height is likewise never converted to mass."
      ""
      "### Model initialization and evaluation boundary"
      ""
      "The runner initializes each selected case at readiness " ...
         + "`requested_window_start`--the requested year's January 1--with " ...
         + "zero earlier spin-up, then begins scientific comparison at the " ...
         + "saved snow-free `evaluation_start`. This does not provide " ...
         + "production-snow-physics winter/spring preconditioning. The saved " ...
         + "initialization table is provenance, not an initialization-" ...
         + "sensitivity experiment."
      ""
       "[Download model initialization provenance](" ...
         + relativeName(files.initialization) + ")"
      ""
      "### Geometric-to-water-equivalent conversion"
      ""
       "Measured surface lowering is multiplied by densities of " ...
         + density_lower + " and " + density_upper ...
         + " kg m^-3 and divided by water density. The resulting shaded band " ...
         + "is a broad sensitivity envelope from porous weathering-crust " ...
         + "material toward intact ice. The " + density_lower ...
          + " kg m^-3 density endpoint is not an intact glacier-ice density. " ...
          + "The two numeric band edges are ordered pointwise, so their " ...
         + "density identity reverses when signed lowering is negative; the " ...
         + "band does not resolve the weathering crust as a modeled state."
       ""
       "## Results"
       ""
       "### Model performance against the measurements"
       ""
       performanceVerdictText(tables.performance_summary)
       ""];

   % Scoring comes before the curve gallery so the ranking of which diagnostic
   % reproduces the measurements appears ahead of the per-site figures.
   lines = [lines; figureOrEmpty(assets.performance_scatter, ...
      "End-of-window cumulative modeled value against measured lowering " ...
      + "converted at the " ...
      + compose('%g', results.policy.effective_density_reference_kg_m3) ...
      + " kg m^-3 reference density, one panel per diagnostic and one " ...
      + "point per site-year, colored by station. Points above the dashed " ...
      + "1:1 line are over-predictions.", ...
      "No site-year was scored, so no endpoint scatter is available.")];
   lines = [lines; figureOrEmpty(assets.performance_summary, ...
      "Aggregate skill across every scored site-year. Left: pooled RMSE " ...
      + "(sample-count weighted) and mean absolute error. Center: median " ...
      + "signed endpoint error, where positive is over-prediction. Right: " ...
      + "the per-site-year endpoint errors behind those aggregates.", ...
      "No site-year was scored, so no aggregate skill figure is available.")];
   lines = [lines; ""; "Aggregate metrics for every scored diagnostic:"; ""];
   lines = [lines; ...
      icemodel.verification.report.markdownTable(tables.performance_summary); ""];
   lines = [lines; "Per-site-year metrics, including every case that could " ...
      + "not be scored and the reason, are in " ...
      + icemodel.verification.report.markdownCode( ...
      relativeName(files.performance)) + "; the aggregate " ...
      + "rows are in " + icemodel.verification.report.markdownCode( ...
      relativeName(files.performance_summary)) ...
      + "."; ""];
   lines = [lines;
       "### Cumulative ablation by site and year"
       ""];

   % Each result subsection states an explicit empty state when saved evidence
   % is absent, instead of emitting a plot or conclusion.
   if isempty(assets.site)
      lines(end + 1, 1) = ...
         "No site-year completed, so no cumulative scientific comparison is available.";
   else
      site_blocks = cell(numel(assets.site), 1);
      for k = 1:numel(assets.site)
         site_blocks{k} = imageMarkdown(assets.site(k), ...
            assets.site_caption(k));
      end
      lines = [lines; vertcat(site_blocks{:})];
   end
   lines = [lines
      ""
      supportSummary(tables.support, files.support)
      ""
      "### Endpoint-perturbation stability"
      ""
      figureOrEmpty(assets.endpoint, ...
         "Start +" + strjoin(compose('%g', ...
         results.policy.endpoint_perturbation_days), "/+") ...
         + "-day and end -" + strjoin(compose('%g', ...
         results.policy.endpoint_perturbation_days), "/-") ...
         + "-day perturbations of the " ...
         + "readiness-selected comparison window. KAN site-years are shown " ...
         + "first, followed by every other saved site-year. Points are saved " ...
         + "model-minus-observation differences; unavailable planned rows are " ...
         + "not plotted as zeros and remain explicit below.", ...
         "No endpoint-perturbation rows were saved.")
      ""
      endpointScope(tables.endpoint)
      ""
       "[Download every endpoint perturbation](" ...
         + relativeName(files.endpoint) + ")"
      ""
      "### Nested-window stability"
      ""
      figureOrEmpty(assets.nested, ...
         "Cumulative-window stability for every represented saved site-year. " ...
         + "The ordinate is net physical solid loss from modeled phase and " ...
         + "vapor terms minus the intact-ice observation conversion in m " ...
         + "w.e.; each point uses its labeled duration in days, and unavailable " ...
         + "fixed endpoints are omitted rather than shortened.", ...
         "No completed nested windows were saved.")
      ""
      nestedScope(tables.nested)
      ""
       "[Download all saved nested windows](" ...
         + relativeName(files.nested) + ")"
      ""
      "### Signed physical solid-mass components"
      ""
      figureOrEmpty(assets.components, ...
         "Signed physical solid-mass components for each completed site-year " ...
         + "in m w.e. Negative phase values are melt-driven solid loss; positive " ...
         + "phase values are refreezing-driven solid gain; solid-vapor exchange " ...
         + "is positive for deposition and negative for sublimation. Net physical " ...
         + "solid change is their signed sum over each row's own [t0,t1) window.", ...
         "No completed signed solid-mass component ledgers were saved.")
      ""
       "[Download signed solid-mass components](" ...
         + relativeName(files.components) + ")"
      ""
       "### Conservation closure and materiality"
       ""
       materialityDenominatorText(results.policy)
       ""
       closureAcceptanceText(tables.identities)
      ""
       figureOrEmpty(assets.closure, ...
         "Scalable closure and materiality summary. Blue circles in the first " ...
         + "panel are each site-year's maximum window identity residual divided " ...
         + "by tolerance (log10; acceptance <= 0). A red cross denotes a " ...
         + "per-step failure; when its window ratio is below the limit, it is " ...
         + "placed at a categorical sentinel just above the limit. It is not a " ...
         + "measured step ratio because the worst step ratio is not saved. The " ...
         + "second panel is " ...
         + "the maximum absolute signed and non-cancelling gross ratio, " ...
         + "normalized by the saved comparison signal G, for each accounting " ...
         + "channel across site-years. Exact rows remain in CSV.", ...
          "No completed closure or materiality diagnostics were saved.")
       ""
       "### Directional outcomes and identifiability drivers"
       ""
       key_findings
       ""
       "Interpretable final classifications:"
       ""
       icemodel.verification.report.markdownTable(directional)
       ""
       "Governing model-accounting sensitivity by case-year:"
       ""
       icemodel.verification.report.markdownTable(drivers)
       ""
       scenario_interpretation
       ""
       "[Download all saved scenario sensitivities](" ...
          + relativeName(files.scenarios) + ")"
       ""
       "### Cross-site site-year synthesis"
      ""
       figureOrEmpty(assets.synthesis, ...
          "Cross-site scientific synthesis for completed selected site-years " ...
          + "only. The left panel compares net physical solid loss from modeled " ...
          + "phase and vapor terms with the observed intact-ice conversion; the " ...
          + "right panel shows signed model-minus-observation differences by " ...
          + "site-year. Operational exclusions and unavailable rows remain in " ...
          + "the linked CSV evidence and are not plotted as scientific values.", ...
          "No completed selected site-year comparisons were saved.")
      ""
       "[Download per-site-year synthesis](" ...
         + relativeName(files.synthesis) + ")"
      ""
      "[Download aggregated synthesis status](" ...
         + relativeName(files.status) + ")"
      ""
      "## Discussion"
      ""
      "Cumulative windows target the timescale on which weathering-crust " ...
         + "storage and refreezing can partly reconcile, while the saved ledger " ...
         + "keeps physical phase and vapor exchanges separate from numerical " ...
         + "domain exchange. Agreement or disagreement is interpreted only " ...
         + "after closure, endpoint storage, and non-cancelling materiality are " ...
         + "examined. Nested windows show whether the sign or magnitude depends " ...
         + "strongly on the chosen cumulative duration. Endpoint perturbations " ...
         + "test dependence on nearby direct start and end postings; they do " ...
         + "not test the initialization policy."
      ""
      "## Limitations"
      ""
      "- Surface-height instruments observe geometry, not whole-column phase " ...
         + "change; weathering-crust density and internal melt are not directly " ...
         + "resolved by the comparison."
      "- IceModel does not explicitly drain subsurface liquid water to runoff, " ...
         + "so retained water can affect later refreezing and thermal state."
       "- The " + density_range ...
         + " kg m^-3 band is a broad porous/weathering-crust sensitivity " ...
         + "envelope, not an intact-ice density range or complete uncertainty " ...
         + "model."
      "- Direct finite endpoint rules avoid interpolation but can shorten or " ...
         + "exclude otherwise plausible windows."
      "- Top-cell deletion height is a quantized remeshing diagnostic; it is " ...
         + "not a continuous surface-displacement state and is not a mass. " ...
         + "Cumulative merge export carries the removed mass, but that " ...
         + "quantity over-counts because merges average the joined pair, so " ...
         + "neither is a surface mass flux and neither is scored."
      "- The prognostic column retains liquid water, so late-season " ...
         + "refreezing can exceed melt and pull the signed solid balance down. " ...
         + "At a bare-ice site that water would have drained, so the signed " ...
         + "balance understates real surface lowering in those intervals. The " ...
         + "runoff diagnostic exists precisely to bound this: it credits " ...
         + "refreezing only up to the melt produced within its trailing " ...
         + "residence window. Runoff and runoff-plus-vapor-loss are therefore " ...
         + "the better surface-lowering comparators, and the performance " ...
         + "section scores that claim rather than asserting it."
      observationRateCaveatText(tables.observation_rates, results.policy)
      "- Requested-year January 1 initialization has zero earlier spin-up and " ...
         + "no production-snow-physics winter/spring preconditioning. Thermal " ...
         + "state, retained liquid water, and later refreezing can therefore " ...
         + "depend on initialization; no initialization sensitivity was run."
      ""
       "## Conclusions"
       ""
       conclusionText(n_completed, counts, closure_counts)
       ""
       key_findings
      ""
       "## Operational Accounting Appendix"
       ""
       "This appendix records which saved site-years reached the model run. " ...
         + "It is kept separate from the scientific comparison above."
       ""
       cohortAccountingText(tables)
       ""
       "- [Readiness and admission ledger](" ...
          + relativeName(files.readiness) + ")"
       "- [Model-run summary](" + relativeName(files.summary) + ")"
       ""
       "Exact row-level selections, exclusions, and reasons are retained in " ...
         + "the downloadable CSV files rather than repeated as large tables."
       ""
       "## Reproducibility Appendix"
      ""
      "This report was generated solely from the saved MAT artifact below. " ...
         + "The builder did not reopen readiness inputs, canonical forcing, " ...
         + "observations, or model configuration paths."
      ""
      "- Run name: " ...
         + icemodel.verification.report.markdownCode(string(results.run_name))
      "- Results MAT: " + icemodel.verification.report.markdownCode(results_file)
      "- Results MAT SHA-256: " ...
         + icemodel.verification.report.markdownCode(source_sha256)
      "- Saved policy version: " + policyVersion(results.policy)
      "- Generated: " + icemodel.verification.report.markdownCode(generated)
      ""
      "### Machine-readable evidence"
      ""
      "- [Readiness ledger](" + relativeName(files.readiness) + ")"
      "- [Site-year summary](" + relativeName(files.summary) + ")"
      "- [Nested windows](" + relativeName(files.nested) + ")"
      "- [Endpoint perturbations](" + relativeName(files.endpoint) + ")"
      "- [Per-site-year synthesis](" + relativeName(files.synthesis) + ")"
      "- [Synthesis status and reasons](" ...
         + relativeName(files.status) + ")"
      "- [Observation support](" + relativeName(files.support) + ")"
      "- [Signed solid-mass components](" ...
         + relativeName(files.components) + ")"
      "- [Grid-translation event ledger](" ...
         + relativeName(files.grid_translation) + ")"
      "- [Model initialization provenance](" ...
         + relativeName(files.initialization) + ")"
      "- [Closure identities](" + relativeName(files.identities) + ")"
      "- [Materiality diagnostics](" + relativeName(files.materiality) + ")"
      "- [Saved scenario sensitivities](" + relativeName(files.scenarios) + ")"
      "- [Effective-density sensitivity](" ...
         + relativeName(files.effective_density) + ")"
      "- [Support exclusions](" + relativeName(files.exclusions) + ")"
      "- [Observation-rate outliers](" ...
         + relativeName(files.observation_rates) + ")"
      "- [Report artifact SHA-256 manifest](" ...
         + relativeName(manifest_file) + ")"
      ""
      "### Saved artifact provenance"
      ""
      provenanceLines(results.paths)
      ""
      "The summary CSV retains forcing and observation artifact names and " ...
         + "SHA-256 fingerprints for each runner row. Report tables preserve " ...
         + "stored numeric precision; displayed Markdown values are compact. " ...
         + "The report artifact manifest hashes the source MAT, generated QMD, " ...
         + "every CSV and figure, and the rendered HTML when present; it omits " ...
         + "itself to avoid recursive identity. " ...
         + "The grid-translation CSV retains interval-level top-deletion and " ...
         + "interior-merge provenance used by the H_grid report panel. " ...
         + "The initialization CSV preserves `initialization_start`, " ...
         + "`initialization_policy`, evaluation boundaries, and the inclusive " ...
         + "run endpoint exactly as saved by the runner."
      ""];
end

function text = performanceVerdictText(summary)
   %PERFORMANCEVERDICTTEXT State which diagnostic best matches the observations.

   % The tolerance and ranking density come from the policy;
   % ablationPerformanceMetrics scores against the same values.
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   tolerance_pct = compose('%g', ...
      100 * policy.scientific.endpoint_tolerance_fraction);

   % The ranking is stated in prose together with the number that produced it.
   all_scored = summary(summary.n_scored > 0, :);
   if isempty(all_scored)
      text = "No site-year could be scored against its observations, so " ...
         + "this report states no performance ranking.";
      return
   end

   % Rank at the reference density, the representative bubbly near-surface
   % ice value, rather than at a band endpoint. Converting observed lowering
   % scales the observation only, so the ranking is density dependent.
   primary = policy.effective_density_reference_kg_m3;
   if ~any(all_scored.density_kg_m3 == primary)
      primary = max(all_scored.density_kg_m3);
   end
   scored = all_scored(all_scored.density_kg_m3 == primary, :);
   if ~any(isfinite(scored.pooled_rmse_mwe))
      % Every diagnostic was scored but none produced a finite pooled RMSE,
      % which happens when no two eligible observations sit one output step
      % apart. Ranking on NaN would name the same diagnostic best and worst.
      text = "**Which diagnostic reproduces the measurements.** No " ...
         + "diagnostic produced a finite pooled RMSE, so no ranking is " ...
         + "available.";
      return
   end
   [~, best] = min(scored.pooled_rmse_mwe);
   [~, worst] = max(scored.pooled_rmse_mwe);
   text = "**Which diagnostic reproduces the measurements.** Ranked by " ...
      + "pooled RMSE against measured lowering converted at the " ...
      + compose('%g', policy.effective_density_reference_kg_m3) ...
      + " kg m^-3 reference density, **" + scored.label(best) ...
      + "** agrees best (" ...
      + compose('%.3f', scored.pooled_rmse_mwe(best)) + " m w.e. pooled " ...
      + "RMSE, median endpoint error " ...
      + compose('%+.3f', scored.median_endpoint_error_mwe(best)) ...
      + " m w.e., " + string(scored.n_within_tolerance(best)) + "/" ...
      + string(scored.n_scored(best)) + " site-years within " ...
      + tolerance_pct + "% of the observed endpoint or of the " ...
      + compose('%g', policy.scientific.signal_floor_mwe) ...
      + " m w.e. signal floor, whichever is larger), and **" ...
      + scored.label(worst) + "** agrees " ...
      + "worst (" + compose('%.3f', scored.pooled_rmse_mwe(worst)) ...
      + " m w.e. pooled RMSE, median endpoint error " ...
      + compose('%+.3f', scored.median_endpoint_error_mwe(worst)) ...
      + " m w.e.). A positive endpoint error is an over-prediction of " ...
      + "ablation. These are agreement statistics over the selected " ...
      + "snow-free windows only; they do not by themselves establish that " ...
      + "the best-agreeing diagnostic is the physically correct comparator, " ...
      + "because the closure and materiality evidence below can still leave " ...
      + "a site-year non-identifiable." ...
      + " Ranking uses the " + compose('%g', primary) + " kg m^-3 " ...
      + "conversion." ...
      + directionalConsensusText(scored) + rateSkillText(scored) ...
      + densitySensitivityText(all_scored, primary);
end

function text = densitySensitivityText(all_scored, primary)
   %DENSITYSENSITIVITYTEXT Say how much the conversion density moves the score.

   % Observed lowering is geometric, so comparing it against a modeled water
   % equivalent needs a density, and that scales the observation while leaving
   % the model untouched. Over a season the material lost is mostly ice, but a
   % single hour of lowering can be porous weathering crust, so the per-step
   % metrics are more sensitive to this choice than the endpoints are.
   others = unique(all_scored.density_kg_m3(all_scored.density_kg_m3 ~= primary));
   if isempty(others)
      text = "";
      return
   end
   text = " **Density sensitivity, and why it governs the conclusion.** " ...
      + "Observed lowering is geometric, so comparing it with a modeled " ...
      + "water-equivalent series needs a density, and that density scales " ...
      + "the observation while leaving the model untouched.";
   for k = 1:numel(others)
      rows = all_scored(all_scored.density_kg_m3 == others(k), :);
      if ~any(isfinite(rows.pooled_rmse_mwe))
         % No finite pooled RMSE at this density, so there is nothing to rank.
         continue
      end
      [~, best_other] = min(rows.pooled_rmse_mwe);
      text = text + " At " + compose('%g', others(k)) + " kg m^-3 the " ...
         + "best-agreeing diagnostic is **" + rows.label(best_other) ...
         + "** with median endpoint error " ...
         + compose('%+.3f', rows.median_endpoint_error_mwe(best_other)) ...
         + " m w.e.";
   end
   text = text + signFlipText(all_scored) ...
      + " Every scored diagnostic and density is in the per-site-year " ...
      + "metrics table, so the full sensitivity is auditable.";
end

function text = signFlipText(all_scored)
   %SIGNFLIPTEXT Warn when the density choice reverses the direction of bias.

   % If the model reads low at one admissible density and high at another, no
   % directional claim survives the conversion uncertainty, and saying the
   % model under-predicts would be an artifact of picking an endpoint.
   densities = unique(all_scored.density_kg_m3);
   low = false;
   high = false;
   for k = 1:numel(densities)
      rows = all_scored(all_scored.density_kg_m3 == densities(k), :);
      low = low || any(rows.median_endpoint_error_mwe < 0);
      high = high || any(rows.median_endpoint_error_mwe > 0);
   end
   if ~(low && high)
      text = "";
      return
   end
   text = " **The direction of the model-observation difference is not " ...
      + "determined by these data.** The sign of the median endpoint error " ...
      + "reverses across the admissible density range: the model reads low " ...
      + "against a dense-ice conversion and high against a porous-crust " ...
      + "conversion. Any statement that the model under-predicts or " ...
      + "over-predicts ablation is therefore a statement about the assumed " ...
      + "density of the material actually removed, not an independent " ...
      + "result. Resolving the direction requires constraining that density, " ...
      + "not more model runs.";
end

function text = directionalConsensusText(scored)
   %DIRECTIONALCONSENSUSTEXT Warn when the ranking looks like compensating error.

   % When every diagnostic misses the same way, the ranking measures how much
   % each one subtracts rather than which is physically right.
   n_low = nnz(scored.median_endpoint_error_mwe < 0);
   n_high = nnz(scored.median_endpoint_error_mwe > 0);
   if n_low == 0 || n_high == 0 || min(n_low, n_high) > 1
      text = "";
      return
   end
   if n_low > n_high
      majority = "under-predict";
      minority = "over-predicts";
      n_majority = n_low;
   else
      majority = "over-predict";
      minority = "under-predicts";
      n_majority = n_high;
   end
   text = " **Read the ranking with care.** " + string(n_majority) + " of " ...
      + string(height(scored)) + " scored diagnostics " + majority ...
      + " and only one " + minority + ", and the diagnostics order by how " ...
      + "much each one subtracts from gross melt. A ranking produced that " ...
      + "way measures how far each definition happens to sit from a common " ...
      + "bias, not which definition is physically correct: subtracting less " ...
      + "from an already low estimate moves it toward the observation " ...
      + "without making it a better description of the surface. Treat the " ...
      + "shared sign as the finding and the ordering as a consequence of it.";
end

function text = observationRateCaveatText(rates, policy)
   %OBSERVATIONRATECAVEATTEXT State which site-years look measurement-limited.
   %
   % Written as a Limitations bullet because a compressed observation record
   % is a limitation of the comparison, not a model result.

   flagged = rates(rates.rate_outlier_flag, :);
   if isempty(flagged)
      text = "- No completed site-year has an observed ablation rate far " ...
         + "below its own station's family, so no record shows the " ...
         + "signature of a compressed measurement scale.";
      return
   end

   % Name the site-years so the reader can act on the caveat.
   labels = strings(height(flagged), 1);
   for k = 1:height(flagged)
      labels(k) = flagged.site_id(k) + " " + compose('%d', flagged.year(k)) ...
         + " (" + compose('%.2f', flagged.rate_ratio_to_station_median(k)) ...
         + " of its station median)";
   end
   text = "- " + compose('%d', height(flagged)) + " completed site-year(s) " ...
      + "have an observed ablation rate below " ...
      + compose('%.2f', policy.observation_rate_outlier_ratio) ...
      + " of the median rate for the same station " ...
      + "across its own completed years: " + strjoin(labels, "; ") + ". A " ...
      + "systematically compressed observation record reads as a model " ...
      + "error, and readiness gates only reject flagged transitions " ...
      + "and unresolved steps, so an unflagged sensor or datum problem " ...
      + "passes admission. These site-years are reported and scored like " ...
      + "any other; excluding them is a human decision, and doing so " ...
      + "silently would improve apparent model skill by deleting " ...
      + "inconvenient observations.";
end

function text = mergeExportScaleText(results)
   %MERGEEXPORTSCALETEXT Measure how far merge export exceeds melt.
   %
   % The scale of the regridding artifact is measured from the run being
   % reported, so the quoted ratio always belongs to this results file.

   completed = results.site_year_results( ...
      string({results.site_year_results.status}) == "completed");
   ratio = NaN(numel(completed), 1);
   for k = 1:numel(completed)
      seasonal = completed(k).seasonal;
      window = seasonal.evaluation_window;
      if ~any(window)
         continue
      end
      first_row = find(window, 1, 'first');
      last_row = find(window, 1, 'last');
      export_mwe = seasonal.model_surface_mass_loss_mwe(last_row) ...
         - seasonal.model_surface_mass_loss_mwe(first_row);
      melt_mwe = seasonal.model_melt_mwe(last_row) ...
         - seasonal.model_melt_mwe(first_row);
      if isfinite(export_mwe) && isfinite(melt_mwe) && melt_mwe > 0
         ratio(k) = export_mwe / melt_mwe;
      end
   end

   usable = ratio(isfinite(ratio));
   if isempty(usable)
      text = "";
      return
   end
   text = " Across this cohort the series runs about " ...
      + compose('%.1f', median(usable)) + " times melt and exceeds melt in " ...
      + compose('%.0f', 100 * mean(usable > 1)) + " percent of the " ...
      + compose('%d', numel(usable)) + " site-years with a usable window " ...
      + "and positive melt, which is the " ...
      + "regridding rule showing through rather than mass leaving the " ...
      + "surface.";
end

function text = residenceWindowText(options)
   %RESIDENCEWINDOWTEXT Describe the runoff diagnostic's trailing window.
   %
   % The window is opts.tlag * opts.dt, so it is read from the model options
   % of the run being reported.

   if ~isstruct(options) || ~isfield(options, 'tlag') ...
         || ~isfield(options, 'dt')
      text = "legacy trailing-window";
      return
   end
   residence_hours = options.tlag * options.dt / 3600;
   text = "legacy " + compose('%g', residence_hours) + "-hour";
end

function options = firstCompletedModelOptions(results)
   %FIRSTCOMPLETEDMODELOPTIONS Model options of the first completed case.
   %
   % Every case in a cohort runs with the same options, so the first
   % completed case describes the run for prose that names a model setting.

   completed = results.site_year_results( ...
      string({results.site_year_results.status}) == "completed");
   if isempty(completed)
      options = struct();
      return
   end
   options = completed(1).model_options;
end

function text = rateSkillText(scored)
   %RATESKILLTEXT State whether any diagnostic reproduces the ablation rate.

   % Endpoint agreement and rate agreement are different claims, and a
   % quantized series can win the first while failing the second badly. The
   % verdict below follows the computed NSE rather than a fixed conclusion.
   best_nse = max(scored.median_nse);
   if ~isfinite(best_nse)
      text = "";
      return
   end

   % NSE > 0 means the diagnostic beats the observed mean as a predictor of
   % per-step increments. At or below zero it does not.
   if best_nse > 0.5
      verdict = " Rate agreement is good as well as cumulative: the best " ...
         + "median Nash-Sutcliffe efficiency on per-step increments is " ...
         + compose('%.3f', best_nse) + ", so at least one diagnostic " ...
         + "reproduces the hourly ablation rate and not only the seasonal " ...
         + "total.";
   elseif best_nse > 0
      verdict = " Rate agreement is modest: the best median " ...
         + "Nash-Sutcliffe efficiency on per-step increments is " ...
         + compose('%.3f', best_nse) + ", so the best diagnostic beats the " ...
         + "observed mean as a predictor of the hourly rate, but not by " ...
         + "much.";
   else
      verdict = " Rate agreement is weaker than endpoint agreement " ...
         + "throughout: the best median Nash-Sutcliffe efficiency on " ...
         + "per-step increments is " + compose('%.3f', best_nse) ...
         + ", so no diagnostic reproduces the hourly ablation rate well " ...
         + "even where cumulative totals are close.";
   end

   text = verdict + " A diagnostic that advances in discrete remeshing " ...
      + "events can score a competitive endpoint while scoring poorly on " ...
      + "rate, because quantization misplaces the loss in time without " ...
      + "changing its total.";
end

function text = closureAcceptanceText(identities)
   %CLOSUREACCEPTANCETEXT Report window and forcing-step verdicts separately.

   if height(identities) == 0
      text = "No saved closure identities were available for acceptance.";
      return
   end

   % The combined verdict requires both levels. Failed-step counts remain exact
   % even though the saved identity table does not retain the worst step ratio.
   window_failure = ~logical(identities.window_passed);
   step_failure = ~logical(identities.step_passed);
   combined_pass = logical(identities.passed);
   failed_step_count = sum(double(identities.failed_step_count));
   text = "Saved closure acceptance: **" + string(nnz(combined_pass)) + "/" ...
      + string(height(identities)) + "** identity rows passed both the window " ...
      + "and per-forcing-step checks. Window acceptance failed for **" ...
      + string(nnz(window_failure)) + "** identity row(s); per-step acceptance " ...
      + "failed for **" + string(nnz(step_failure)) + ...
      "** identity row(s) across **" + string(failed_step_count) ...
      + "** forcing row(s).";
end

function text = cohortAccountingText(tables)
   %COHORTACCOUNTINGTEXT State model-run and supplemental-data counts plainly.

   rows = tables.readiness;
   if height(rows) == 0
      text = "No saved readiness rows were available for cohort accounting.";
      return
   end

   % Admission, selection, and completion are saved model-run decisions. The
   % supplemental snow-data count is reported separately without naming it an
   % independent observation source or an intervening gate.
   admitted = nnz(logical(rows.admitted));
   selected = nnz(logical(tables.summary.admitted) ...
      & logical(tables.summary.selected));
   completed = nnz(string(tables.summary.status) == "completed");
   text = "The saved inventory contains **" + string(height(rows)) ...
      + "** candidate site-year(s). **" + string(admitted) ...
      + "** met the saved admission rules, **" + string(selected) ...
      + "** were selected for a run, and **" + string(completed) ...
      + "** completed.";

   vars = string(rows.Properties.VariableNames);
   if ismember("snow_input_ready", vars)
      snow_ready = nnz(logical(rows.snow_input_ready));
      text = text + " Supplemental snow data were available for **" ...
         + string(snow_ready) + "** and unavailable for **" ...
         + string(height(rows) - snow_ready) ...
         + "**; this inventory is reported separately from admission.";
   else
      text = text + " No snow-input readiness inventory was saved.";
   end
end

function text = materialityDenominatorText(policy)
   %MATERIALITYDENOMINATORTEXT Define the saved comparison-signal denominator.

   formula = "`G = max(abs(A_model), abs(A_obs), signal_floor)`";
   text = "Saved materiality ratios use comparison signal " + formula + ".";
   if isstruct(policy) && isfield(policy, 'scientific') ...
         && isstruct(policy.scientific) ...
         && isfield(policy.scientific, 'signal_floor_mwe')
      signal_floor = policy.scientific.signal_floor_mwe;
      if isnumeric(signal_floor) && isscalar(signal_floor) ...
            && isfinite(signal_floor)
         text = text + " The saved `signal_floor` is **" ...
            + string(sprintf('%.5g', signal_floor)) + " m w.e.**";
      end
   end
end

function lines = endpointScope(endpoint)
   %ENDPOINTSCOPE Summarize perturbation availability without row-level prose.

   if height(endpoint) == 0
      lines = "Saved endpoint-perturbation site-years: none.";
      return
   end
   keys = string(endpoint.site_id) + "|" + string(endpoint.year);
   n_available = nnz(logical(endpoint.available) ...
      & isfinite(endpoint.model_minus_observation_mwe));
   lines = "Saved endpoint perturbations: **" + string(n_available) + "/" ...
      + string(height(endpoint)) + " available** across **" ...
      + string(numel(unique(keys, 'stable'))) + "** site-year(s). " ...
      + "The CSV retains every unavailable row and reason.";
end

function lines = nestedScope(nested)
   %NESTEDSCOPE Summarize nested windows without one bullet per site-year.

   if height(nested) == 0
      lines = "Saved nested-window site-years: none.";
      return
   end

   keys = string(nested.site_id) + "|" + string(nested.year);
   n_available = nnz(logical(nested.available) ...
      & isfinite(nested.model_minus_observation_mwe));
   lines = "Saved nested windows: **" + string(n_available) + "/" ...
      + string(height(nested)) + " available** across **" ...
      + string(numel(unique(keys, 'stable'))) + "** site-year(s). " ...
      + "The CSV retains every unavailable window and reason.";
end

function lines = supportSummary(support, support_file)
   %SUPPORTSUMMARY Count every saved observation-support class.

   if height(support) == 0
      lines = "No completed observation-support postings were saved.";
      return
   end
   classes = string(support.support_class);
   lines = [ ...
      "Saved full-season observation support: " ...
         + string(nnz(classes == "direct")) + " direct, " ...
         + string(nnz(classes == "snow_censored")) ...
         + " snow-censored, " ...
         + string(nnz(classes == "unknown_snow")) ...
         + " excluded for unknown snow depth, " ...
         + string(nnz(classes == "flagged")) + " flagged, " ...
         + string(nnz(classes == "unknown_quality")) ...
         + " excluded for unknown quality, and " ...
         + string(nnz(classes == "nonfinite_target")) ...
         + " excluded for a nonfinite target. **" ...
         + string(nnz(logical(support.in_comparison_window))) + "/" ...
         + string(height(support)) + "** postings are in the saved " ...
         + "comparison window."
      ""
      "[Download posting-level observation support](" ...
         + relativeName(support_file) + ")"];
end

function counts = categoryCounts(synthesis)
   %CATEGORYCOUNTS Count the mutually exclusive saved outcome categories.

   category = string(synthesis.report_category);
   classification = string(synthesis.classification);
   scientifically_unavailable = category == "scientifically_unavailable";
   non_identifiable = scientifically_unavailable ...
      & classification == "non_identifiable";
   not_physically_comparable = scientifically_unavailable ...
      & classification == "not_physically_comparable";
   counts = struct( ...
      'interpretable', nnz(category == "interpretable"), ...
      'scientifically_unavailable', ...
         nnz(scientifically_unavailable), ...
      'non_identifiable', nnz(non_identifiable), ...
      'not_physically_comparable', nnz(not_physically_comparable), ...
      'scientifically_unavailable_other', nnz(scientifically_unavailable ...
         & ~non_identifiable & ~not_physically_comparable), ...
      'readiness_excluded', nnz(category == "readiness_excluded"), ...
      'not_selected', nnz(category == "not_selected"), ...
      'execution_unavailable', ...
         nnz(category == "execution_unavailable"), ...
      'unknown_status', nnz(category == "unknown_status"));
end

function summary = directionalClassificationTable(synthesis)
   %DIRECTIONALCLASSIFICATIONTABLE Count interpretable final directions.

   use = string(synthesis.report_category) == "interpretable";
   values = string(synthesis.classification(use));
   groups = unique(values, 'stable');
   count = zeros(numel(groups), 1);
   for k = 1:numel(groups)
      count(k) = nnz(values == groups(k));
   end
   summary = table(groups, count, ...
      'VariableNames', {'classification', 'count'});
   summary = sortrows(summary, {'count', 'classification'}, ...
      {'descend', 'ascend'});
end

function counts = closureSiteYearCounts(identities)
   %CLOSURESITEYEARCOUNTS Count site-years that pass every saved identity.

   counts = struct('evaluated', 0, 'passed', 0);
   if height(identities) == 0
      return
   end

   % Closure is a site-year verdict only when every saved identity for that
   % site-year passes; row counts remain available in the detailed table.
   keys = string(identities.case_id) + "|" + string(identities.site_id) ...
      + "|" + string(identities.year);
   groups = unique(keys, 'stable');
   counts.evaluated = numel(groups);
   for k = 1:numel(groups)
      counts.passed = counts.passed ...
         + all(logical(identities.passed(keys == groups(k))));
   end
end

function summary = accountingDriverTable(scenarios)
   %ACCOUNTINGDRIVERTABLE Count credible material accounting sensitivities.

   % This reproduces the comparator gate exactly and also protects reports made
   % from older or synthetic saved rows: only accounting alternatives can govern,
   % and only after the non-cancelling materiality gate fires.
   role_value = string(scenarios.role);
   use = logical(scenarios.credible) ...
      & role_value == "accounting" & logical(scenarios.material);
   values = scenarios(use, :);
   group_key = string(values.role) + "|" + string(values.scenario);
   groups = unique(group_key, 'stable');
   role = strings(numel(groups), 1);
   scenario = strings(numel(groups), 1);
   evaluated = zeros(numel(groups), 1);
   changes_sign = zeros(numel(groups), 1);
   changes_classification = zeros(numel(groups), 1);
   destabilized = zeros(numel(groups), 1);
   for k = 1:numel(groups)
      rows = values(group_key == groups(k), :);
      role(k) = string(rows.role(1));
      scenario(k) = string(rows.scenario(1));
      keys = string(rows.case_id) + "|" + string(rows.site_id) + "|" ...
         + string(rows.year);
      evaluated(k) = numel(unique(keys));
      changes_sign(k) = numel(unique(keys(logical(rows.changes_sign))));
      changes_classification(k) = numel(unique( ...
         keys(logical(rows.changes_classification))));
      destabilized(k) = numel(unique(keys(logical(rows.changes_sign) ...
         | logical(rows.changes_classification))));
   end
   summary = table(role, scenario, evaluated, changes_sign, ...
      changes_classification, destabilized, 'VariableNames', {'role', 'scenario', ...
      'evaluated_case_years', 'changes_sign_case_years', ...
      'changes_classification_case_years', 'destabilized_case_years'});
   summary = sortrows(summary, ...
      {'destabilized_case_years', 'role', 'scenario'}, ...
      {'descend', 'ascend', 'ascend'});
end

function text = keyFindingsText(directional, drivers, has_endpoint_deficit)
   %KEYFINDINGSTEXT Summarize directional and identifiability evidence.

   if isempty(directional)
      direction_text = "No completed case-year supports a directional " ...
         + "model--observation interpretation.";
   else
      labels = replace(string(directional.classification), "_", "-") ...
         + " n=" + string(directional.count);
      direction_text = "Interpretable directional classifications: " ...
         + strjoin(labels, ", ") + ".";
   end

   destabilized = drivers(drivers.destabilized_case_years > 0, :);
   if isempty(destabilized)
      driver_text = "No credible material model-accounting scenario changes " ...
         + "sign or classification.";
   else
      % The report describes the numerical sensitivity in plain language while
      % the linked CSV preserves exact machine-readable scenario identifiers.
      merge_delete = string(destabilized.scenario) == "merge_delete_solid";
      if any(merge_delete)
         n_changed = sum( ...
            destabilized.destabilized_case_years(merge_delete));
         case_year_label = "case-year";
         if n_changed ~= 1
            case_year_label = "case-years";
         end
         driver_text = "Material numerical merge/delete solid exchange changed " ...
            + "the comparison in " + string(n_changed) + " " ...
            + case_year_label + ". " ...
            + "This accounting sensitivity does not identify a physical cause.";
         if any(~merge_delete)
            driver_text = driver_text + " Other material numerical accounting " ...
               + "alternatives also changed sign or classification; exact " ...
               + "non-exclusive counts remain in the scenario table and CSV.";
         end
      else
         driver_text = "Material numerical accounting alternatives changed " ...
            + "sign or classification in saved case-year comparisons; exact " ...
            + "non-exclusive counts remain in the scenario table and CSV. " ...
            + "These accounting sensitivities do not identify a physical cause.";
      end
   end

   % Endpoint-deficit interpretation is run-specific and must not appear when
   % no saved endpoint-deficit scenario exists.
   endpoint_text = "";
   if has_endpoint_deficit
      endpoint_text = " Caller-supplied endpoint-deficit scenarios are " ...
         + "unverified and never gate classification.";
   end
   text = direction_text + " " + driver_text + endpoint_text ...
      + " Temporal endpoint perturbations are separate, non-gating " ...
      + "sensitivity diagnostics. Scenario counts are non-exclusive and do " ...
      + "not establish physical causation.";
end

function present = hasEndpointDeficitRows(scenarios)
   %HASENDPOINTDEFICITROWS Detect saved endpoint-deficit scenario evidence.

   present = false;
   if height(scenarios) == 0
      return
   end

   % Role and scenario are both required so similarly named accounting rows do
   % not trigger observation-endpoint interpretation.
   present = any(string(scenarios.role) == "endpoint" ...
      & startsWith(string(scenarios.scenario), "endpoint_deficit"));
end

function text = endpointDeficitDiagramLine(present)
   %ENDPOINTDEFICITDIAGRAMLINE Return the optional conceptual-diagram node.

   text = "";
   if present
      text = "  W[""Unverified endpoint weathering-crust deficit " ...
         + "sensitivity""] -.-> C";
   end
end

function text = scenarioInterpretationText(has_endpoint_deficit)
   %SCENARIOINTERPRETATIONTEXT Explain saved scenario evidence without causation.

   text = "Driver counts are non-exclusive because one site-year can be " ...
      + "destabilized by multiple accounting scenarios. Only credible, " ...
      + "material model-accounting alternatives gate classification. ";
   if has_endpoint_deficit
      text = text + "Caller-supplied endpoint-deficit values remain visible " ...
         + "sensitivity cases but cannot gate classification because the " ...
         + "current runner has no validated, hashed observation-derived " ...
         + "endpoint-deficit payload. ";
   end
   text = text + "Temporal endpoint perturbations are reported separately " ...
      + "and do not gate classification. These counts diagnose sensitivity " ...
      + "to saved scenarios, not physical causation.";
end

function text = conclusionText(n_completed, counts, closure_counts)
   %CONCLUSIONTEXT Return an evidence-scaled conclusion.

   if n_completed == 0
      text = "No site-year completed. The saved artifact supports readiness " ...
         + "and exclusion reporting only; it does not support a model--observation " ...
         + "ablation conclusion.";
   else
      text = string(n_completed) + " case-year(s) completed. " ...
         + string(closure_counts.passed) + " of " ...
         + string(closure_counts.evaluated) + " site-year(s) with saved " ...
         + "closure evidence passed every saved closure identity. " ...
         + string(counts.interpretable) + " yielded an identifiable " ...
         + "directional conclusion; " ...
         + string(counts.scientifically_unavailable) ...
         + " did not support a unique directional interpretation. Site-level " ...
         + "curves and nested-window evidence must remain primary to any " ...
         + "cross-site generalization.";
   end
end

function lines = provenanceLines(paths)
   %PROVENANCELINES Format every saved runner path without opening it.

   if ~isstruct(paths)
      lines = "No saved path registry was present.";
      return
   end
   names = string(fieldnames(paths));
   lines = strings(numel(names), 1);
   n_lines = 0;
   for k = 1:numel(names)
      value = string(paths.(names(k)));
      if isscalar(value) && strlength(value) > 0
         n_lines = n_lines + 1;
         lines(n_lines) = "- " ...
            + icemodel.verification.report.escapeMarkdownText(names(k)) ...
            + ": " + icemodel.verification.report.markdownCode(value);
      end
   end
   lines = lines(1:n_lines);
   if n_lines == 0
      lines = "No persisted artifact paths were recorded.";
   end
end

function text = policyVersion(policy)
   %POLICYVERSION Format the saved policy version, or state that it is
   % unavailable.

   if isstruct(policy) && isfield(policy, 'version')
      text = icemodel.verification.report.markdownCode(string(policy.version));
   else
      text = "unavailable";
   end
end

function text = seasonRangeText(policy)
   %SEASONRANGETEXT Format saved month-day bounds for scientific prose.

   first = policy.evaluation_season_start_month_day;
   last = policy.evaluation_season_end_month_day;
   first_time = datetime(2000, first(1), first(2));
   last_time = datetime(2000, last(1), last(2));
   text = string(first_time, "d MMMM") + " through " ...
      + string(last_time, "d MMMM");
end

function text = densityRangeText(policy)
   %DENSITYRANGETEXT Format the saved density conversion envelope.
   %
   % The band may carry interior points such as the reference density, so the
   % envelope is the first and last value rather than the first two.

   values = policy.effective_density_kg_m3;
   text = string(sprintf('%g--%g', values(1), values(end)));
end

function lines = figureOrEmpty(asset, caption, empty_text)
   %FIGUREOREMPTY Return one figure or an explicit absence statement.

   if strlength(asset) == 0
      lines = string(empty_text);
   else
      lines = imageMarkdown(asset, caption);
   end
end

function lines = imageMarkdown(asset, caption)
   %IMAGEMARKDOWN Emit a parse-safe image and separate sanitized caption.

   [~, name, ext] = fileparts(asset);
   alt_text = replace(string(name), ["-", "_"], " ");
   lines = [""; "![" + alt_text + "](" ...
      + "report-assets/" + name + ext + ")"; ...
      ""; "*" + icemodel.verification.report.escapeMarkdownText(caption) + "*"];
end

function name = safeFilename(value)
   %SAFEFILENAME Restrict saved identifiers to one portable basename.

   name = lower(regexprep(string(value), '[^A-Za-z0-9._-]+', '_'));
   name = strip(name, '_');
   if strlength(name) == 0
      name = "case";
   end
end


function name = relativeName(filename)
   %RELATIVENAME Return the output-local file name for a report link.

   [~, stem, ext] = fileparts(filename);
   name = stem + ext;
end

function configureDiagnosticAxis(ax, labels)
   %CONFIGUREDIAGNOSTICAXIS Label vertical diagnostic categories on the x axis.

   % The performance panels use vertical bars and scatter columns, so the
   % categories belong on x. configureCategoryAxis labels a reversed y axis for
   % the horizontal evidence rows elsewhere in this report.
   xticks(ax, 1:numel(labels))
   xticklabels(ax, icemodel.verification.report.safeLabel(labels))
   xlim(ax, [0.4, numel(labels) + 0.6])
   ax.TickLabelInterpreter = 'none';
   ax.XAxis.FontSize = 9;
   icemodel.verification.report.formatReportAxes(ax)
end
