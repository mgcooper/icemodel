function report = writePromiceAblationReadiness(kwargs)
   %WRITEPROMICEABLATIONREADINESS Write the PROMICE ablation readiness ledger.
   %
   %  report = ...
   %     icemodel.verification.setup.writePromiceAblationReadiness( ...
   %     evaluation_data_root=eval_root, input_data_root=input_root, ...
   %     output_dir=output_dir)
   %
   % The CSV contains one row for every canonical PROMICE case and calendar
   % year intersecting its manifest period. Verdicts inspect the actual staged
   % observations and promice_filled artifacts. The manifest and final producer
   % ledger provide inventory and provenance, never substitute for payloads.

   arguments
      kwargs.evaluation_data_root (1, 1) string
      kwargs.input_data_root (1, 1) string
      kwargs.output_dir (1, 1) string
   end

   policy = ...
      icemodel.verification.namelists.promiceAblationReadiness();
   cases = icemodel.verification.listcases( ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      input_data_root=kwargs.input_data_root, dataset_family="promice");
   if isempty(cases)
      error('icemodel:verification:promiceAblationReadiness:manifestMissing', ...
         'no canonical PROMICE cases found under %s', ...
         kwargs.evaluation_data_root)
   end
   case_ids = lower(string({cases.case_id}));
   if numel(unique(case_ids)) ~= numel(case_ids)
      error('icemodel:verification:promiceAblationReadiness:duplicateCase', ...
         'PROMICE manifest contains duplicate case ids')
   end

   % Load one station at a time so the large staged tables are promptly freed.
   % Each station's annual rows are collected into a preallocated slot and
   % joined once, so the ledger never grows an element at a time.
   per_case = cell(numel(cases), 1);
   for k = 1:numel(cases)
      observation = observationPayload(cases(k), ...
         kwargs.evaluation_data_root);
      forcing = forcingPayload(cases(k), kwargs.input_data_root, policy);
      first_year = year(toUtc(cases(k).period.start));
      last_year = year(toUtc(cases(k).period.end));
      years = first_year:last_year;
      case_rows = repmat(rowTemplate(policy), numel(years), 1);
      for y = 1:numel(years)
         case_rows(y) = annualRow( ...
            cases(k), years(y), observation, forcing, policy);
      end
      per_case{k} = case_rows;
      clear observation forcing
   end
   rows = vertcat(per_case{:});
   ledger = sortrows(struct2table(rows), {'case_id', 'year'});
   summary = struct( ...
      'row_count', height(ledger), ...
      'canonical_case_count', numel(unique(ledger.case_id)), ...
      'ice_model_forcing_ready_count', nnz(ledger.ice_model_forcing_ready), ...
      'evaluation_target_ready_count', nnz(ledger.evaluation_target_ready), ...
      'snow_input_ready_count', nnz(ledger.snow_input_ready), ...
      'production_snow_candidate_ready_count', ...
         nnz(ledger.production_snow_candidate_ready), ...
      'admitted_count', nnz(ledger.admitted), ...
      'inventory', caseInventory(cases, kwargs.evaluation_data_root, ...
         kwargs.input_data_root), ...
      'reconstruction_source_audit', reconstructionWriteAudit());

   if ~isfolder(kwargs.output_dir)
      mkdir(kwargs.output_dir)
   end
   csv_file = fullfile(kwargs.output_dir, ...
      "promice_ablation_readiness.csv");
   json_file = fullfile(kwargs.output_dir, ...
      "promice_ablation_readiness.json");
   writetable(ledger, csv_file)
   payload = struct( ...
      'schema_version', "1.1", 'path_base', "selected_data_root", ...
      'policy', policy, 'summary', summary, ...
      'observation_hash_scope', ...
         "current hashes are a forward baseline; historical byte identity " ...
         + "cannot be proven without an earlier ledger", ...
      'csv_sha256', ...
         icemodel.verification.setup.fileSha256(csv_file), ...
      'rows', table2struct(ledger));
   icemodel.verification.setup.writeJson(json_file, payload)
   report = struct('rows', ledger, 'policy', policy, 'summary', summary, ...
      'files', struct('csv', string(csv_file), 'json', string(json_file)));
end

function observation = observationPayload(c, evaluation_data_root)
   %OBSERVATIONPAYLOAD Hash and load the actual staged targets.data timetable.
   pathname = string(c.evaluation_path);
   observation = struct('path', "", 'sha256', "", 'ok', false, ...
      'reason', "observation artifact is missing", 'data', timetable());
   if pathname == ""
      return
   end
   if ~icemodel.isPathInside(pathname, evaluation_data_root)
      observation.reason = ...
         "observation artifact escapes selected evaluation root";
      return
   end
   observation.path = ...
      icemodel.verification.setup.fixtureRelativePosix( ...
      evaluation_data_root, pathname);
   if ~isfile(pathname)
      return
   end
   observation.sha256 = ...
      icemodel.verification.setup.fileSha256(pathname);
   try
      saved = load(pathname, 'targets');
   catch err
      observation.reason = "observation artifact is unreadable: " ...
         + string(err.message);
      return
   end
   if ~isfield(saved, 'targets') || ~isstruct(saved.targets) ...
         || ~isfield(saved.targets, 'data') ...
         || ~istimetable(saved.targets.data)
      observation.reason = ...
         "observation artifact lacks targets.data timetable";
      return
   end
   expected_case = ...
      icemodel.forcing.helpers.normalizedFileToken(c.case_id);
   expected_site = ...
      icemodel.forcing.helpers.normalizedFileToken(c.site_id);
   valid_identity = isfield(saved.targets, 'metadata') ...
      && isstruct(saved.targets.metadata) ...
      && isscalar(saved.targets.metadata) ...
      && all(isfield(saved.targets.metadata, {'station', 'site_id'})) ...
      && isscalar(string(saved.targets.metadata.station)) ...
      && isscalar(string(saved.targets.metadata.site_id)) ...
      && expected_site == expected_case ...
      && icemodel.forcing.helpers.normalizedFileToken( ...
      saved.targets.metadata.station) == expected_case ...
      && icemodel.forcing.helpers.normalizedFileToken( ...
      saved.targets.metadata.site_id) == expected_case;
   if ~valid_identity
      observation.reason = ...
         "observation artifact station/site identity does not match case";
      return
   end
   observation.ok = true;
   observation.reason = "";
   observation.data = saved.targets.data;
end

function forcing = forcingPayload(c, input_data_root, policy)
   %FORCINGPAYLOAD Verify producer hashes and inspect complete payload windows.
   data_root = string(fileparts(input_data_root));
   case_id = lower(string(c.case_id));
   manifest_file = fullfile(data_root, 'preview', 'qa', 'gapfill', ...
      'plans', case_id + "-report-inputs.json");
   forcing = struct('manifest', "", 'manifest_sha256', "", 'artifact', "", ...
      'artifact_sha256', "", 'readiness_artifact', "", ...
      'readiness_sha256', "", 'acceptance_start', "", ...
      'acceptance_end', "", 'ok', false, ...
      'reason', "promice_filled producer manifest is missing", ...
      'windows', struct([]), 'ledger', table());
   forcing.manifest = ...
      icemodel.verification.setup.fixtureRelativePosix( ...
      data_root, manifest_file);
   if ~isfile(manifest_file)
      return
   end
   forcing.manifest_sha256 = ...
      icemodel.verification.setup.fileSha256(manifest_file);
   try
      manifest = jsondecode(fileread(manifest_file));
   catch err
      forcing.reason = "promice_filled producer manifest is unreadable: " ...
         + string(err.message);
      return
   end
   valid_manifest = isstruct(manifest) && isscalar(manifest) ...
      && isfield(manifest, 'site') && isfield(manifest, 'path_base') ...
      && isfield(manifest, 'artifacts') ...
      && isfield(manifest, 'acceptance_window') ...
      && strcmpi(string(manifest.site), case_id) ...
      && string(manifest.path_base) == "selected_data_root";
   if ~valid_manifest
      forcing.reason = "promice_filled producer manifest identity is invalid";
      return
   end
   acceptance = manifest.acceptance_window;
   valid_acceptance = isstruct(acceptance) && isscalar(acceptance) ...
      && all(isfield(acceptance, {'start', 'end'})) ...
      && isscalar(string(acceptance.start)) ...
      && isscalar(string(acceptance.end));
   if ~valid_acceptance
      forcing.reason = "producer acceptance-window identity is invalid";
      return
   end
   try
      acceptance_start = toUtc(acceptance.start);
      acceptance_end = toUtc(acceptance.end);
   catch
      forcing.reason = "producer acceptance-window identity is invalid";
      return
   end
   if isnat(acceptance_start) || isnat(acceptance_end) ...
         || acceptance_start > acceptance_end
      forcing.reason = "producer acceptance-window identity is invalid";
      return
   end
   artifacts = reshape(manifest.artifacts, [], 1);
   if ~isstruct(artifacts) ...
         || ~all(isfield(artifacts, {'role', 'path', 'sha256'}))
      forcing.reason = "producer manifest artifact schema is incomplete";
      return
   end
   roles = string({artifacts.role});
   filled_index = find(roles == "filled");
   readiness_index = find(roles == "readiness");
   if numel(filled_index) ~= 1 || numel(readiness_index) ~= 1
      forcing.reason = ...
         "producer manifest must name one filled and one readiness artifact";
      return
   end
   filled = artifacts(filled_index);
   readiness = artifacts(readiness_index);
   forcing.artifact = string(filled.path);
   forcing.artifact_sha256 = string(filled.sha256);
   forcing.readiness_artifact = string(readiness.path);
   forcing.readiness_sha256 = string(readiness.sha256);
   try
      filled_file = ...
         icemodel.verification.helpers.assertRootRelativeArtifactSha256( ...
         data_root, string(filled.path), string(filled.sha256));
      readiness_file = ...
         icemodel.verification.helpers.assertRootRelativeArtifactSha256( ...
         data_root, string(readiness.path), string(readiness.sha256));
   catch err
      forcing.reason = "producer-pinned artifact identity is invalid: " ...
         + string(err.message);
      return
   end
   try
      loaded = load(filled_file, 'met');
      % Reuse the runtime's strongest product/provenance gate so scientific
      % readiness cannot admit a stale policy, engine, registry, or channel
      % ledger that the model itself would reject.
      icemodel.forcing.reconstruct.assertPromiceFilledArtifact( ...
         filled_file, loaded.met, case_id)
      % Whole-artifact readiness is intentionally not the cohort gate: a
      % requested year is admissible only when forcingCoverage below finds one
      % complete contiguous window that encloses that exact annual request.
      [~, payload_reason, forcing.windows, payload_cadence_seconds] = ...
         icemodel.verification.setup.metArtifactReadiness(filled_file);
      forcing.ledger = readtable(readiness_file, TextType='string');
   catch err
      forcing.reason = "promice_filled payload is invalid: " ...
         + string(err.message);
      return
   end
   if ~isscalar(payload_cadence_seconds) ...
         || ~isfinite(payload_cadence_seconds) ...
         || payload_cadence_seconds ~= policy.forcing_cadence_seconds
      forcing.reason = "promice_filled payload cadence must equal " ...
         + string(policy.forcing_cadence_seconds) + " seconds; found " ...
         + string(payload_cadence_seconds);
      return
   end
   if isempty(forcing.windows)
      forcing.reason = "promice_filled has no complete forcing window: " ...
         + payload_reason;
      return
   end
   payload_start = toUtc(forcing.windows(1).start_time);
   payload_end = toUtc(forcing.windows(end).end_time);
   if acceptance_start ~= payload_start || acceptance_end ~= payload_end
      forcing.reason = ...
         "producer acceptance window does not match filled payload support";
      return
   end
   ledger_names = string(forcing.ledger.Properties.VariableNames);
   if ~all(ismember(["site", "year"], ledger_names))
      forcing.reason = ...
         "producer readiness ledger identity schema is incomplete";
      return
   end
   ledger_sites = string(forcing.ledger.site);
   valid_ledger_site = height(forcing.ledger) > 0 ...
      && all(strlength(ledger_sites) > 0) ...
      && all(icemodel.forcing.helpers.normalizedFileToken( ...
      ledger_sites) == ...
      icemodel.forcing.helpers.normalizedFileToken(case_id));
   if ~valid_ledger_site
      forcing.reason = ...
         "producer readiness ledger site identity does not match case";
      return
   end
   ledger_years = forcing.ledger.year;
   expected_years = (year(acceptance_start):year(acceptance_end)).';
   valid_ledger_years = isnumeric(ledger_years) ...
      && all(isfinite(ledger_years)) ...
      && all(ledger_years == fix(ledger_years)) ...
      && isequal(sort(double(ledger_years(:))), expected_years);
   if ~valid_ledger_years
      forcing.reason = ...
         "producer readiness ledger years do not match acceptance window";
      return
   end
   forcing.ok = true;
   forcing.reason = "";
   forcing.acceptance_start = formatTime(acceptance_start);
   forcing.acceptance_end = formatTime(acceptance_end);
end

function row = annualRow(c, y, observation, forcing, policy)
   %ANNUALROW Derive one case-year from actual payload support and flags.
   row = rowTemplate(policy);
   row.case_id = lower(string(c.case_id));
   row.site_id = string(c.site_id);
   row.surface_zone = string(c.surface_zone);
   row.year = y;
   [window_start, window_end] = annualWindow(y, policy);
   row.requested_window_start = formatTime(window_start);
   row.requested_window_end = formatTime(window_end);
   row.observation_artifact = observation.path;
   row.observation_sha256 = observation.sha256;
   row.forcing_producer_manifest = forcing.manifest;
   row.forcing_producer_manifest_sha256 = forcing.manifest_sha256;
   row.forcing_artifact = forcing.artifact;
   row.forcing_sha256 = forcing.artifact_sha256;
   row.forcing_readiness_artifact = forcing.readiness_artifact;
   row.forcing_readiness_sha256 = forcing.readiness_sha256;
   row.forcing_acceptance_start = forcing.acceptance_start;
   row.forcing_acceptance_end = forcing.acceptance_end;
   row.target_variable = policy.target_variable;
   % target_units records payload metadata and remains empty when the source
   % omits it; policy.target_units is the comparison requirement below.
   row.target_positive_direction = policy.target_positive_direction;
   row.timestamp_convention = policy.timestamp_convention;
   row.interval_convention = policy.interval_convention;
   row.reference_semantics = policy.reference_semantics;
   row.adjustment_semantics = policy.adjustment_semantics;
   row.possible_support_count = floor( ...
      seconds(window_end - window_start) ...
      / policy.observation_cadence_seconds) + 1;

   [row, row.evaluation_target_reason] = observationStatus( ...
      row, observation, window_start, window_end, policy);
   period_start = toUtc(c.period.start);
   period_end = toUtc(c.period.end);
   if period_start > window_start || period_end < window_end
      period_reason = ...
         "canonical case period does not contain requested annual window";
      row.evaluation_target_reason = strjoin( ...
         [row.evaluation_target_reason( ...
         row.evaluation_target_reason ~= ""); period_reason], "; ");
   end
   row.evaluation_target_ready = row.evaluation_target_reason == "";
   [coverage_ready, coverage_reason] = ...
      forcingCoverage(forcing, window_start, window_end);
   [row.ice_model_forcing_ready, row.ice_model_forcing_reason] = ...
      consumerStatus(forcing, y, "icemodel", coverage_ready, coverage_reason);
   [row.snow_input_ready, row.snow_input_reason] = ...
      consumerStatus(forcing, y, "snowmodel", coverage_ready, coverage_reason);
   row.production_snow_candidate_ready = false;
   row.production_snow_candidate_reason = ...
      "no production snow-model candidate is owned by this readiness audit";
   row.admitted = row.ice_model_forcing_ready ...
      && row.evaluation_target_ready;
   reasons = strings(0, 1);
   if ~row.ice_model_forcing_ready
      reasons(end + 1, 1) = "forcing: " + row.ice_model_forcing_reason;
   end
   if ~row.evaluation_target_ready
      reasons(end + 1, 1) = "evaluation target: " ...
         + row.evaluation_target_reason;
   end
   row.exclusion_reason = strjoin(reasons, "; ");
end

function [row, reason] = observationStatus( ...
      row, observation, window_start, window_end, policy)
   %OBSERVATIONSTATUS Audit the target's measure, support, flags, and
   % sensitivity.
   if ~observation.ok
      reason = observation.reason;
      return
   end
   all_data = observation.data;
   [all_data, time_reason] = normalizeObservationRowTimes(all_data);
   if time_reason ~= ""
      reason = time_reason;
      return
   end
   data = all_data(all_data.Time >= window_start ...
      & all_data.Time <= window_end, :);
   row.observation_sample_count = height(data);
   names = string(data.Properties.VariableNames);
   target = policy.target_field;
   snow = policy.snow_variable;
   required = [target, snow, policy.support_flag_fields];
   row.target_available = ismember(target, names);
   missing = setdiff(required, names, 'stable');
   if ~isempty(missing)
      reason = "actual observation payload is missing required variable(s): " ...
         + strjoin(missing, ", ");
      return
   end

   % A malformed target or flag is a case-level scientific exclusion. Detect
   % it before isfinite/table concatenation so one bad payload cannot abort
   % writing the complete cohort CSV and JSON audit.
   numeric_variable = arrayfun(@(name) ...
      isnumeric(data.(name)) || islogical(data.(name)), required);
   nonnumeric = required(~numeric_variable);
   if ~isempty(nonnumeric)
      reason = ...
         "actual observation payload has nonnumeric required variable(s): " ...
         + strjoin(nonnumeric, ", ");
      return
   end

   if row.target_available
      index = find(names == target, 1);
      % Older observation timetables may omit these optional table metadata
      % arrays. Leave the public fields empty so this case receives the normal
      % semantic exclusion below instead of aborting the complete cohort audit.
      units = string(data.Properties.VariableUnits);
      if numel(units) >= index
         row.target_units = units(index);
      end
      descriptions = string(data.Properties.VariableDescriptions);
      if numel(descriptions) >= index
         row.target_description = descriptions(index);
      end
      row.finite_target_count = nnz(isfinite(data.(target)));
   end
   if ismember(snow, names)
      values = data.(snow);
      values = sort(values(isfinite(values)));
      row.finite_snow_count = numel(values);
      if ~isempty(values)
         row.snow_depth_min_m = values(1);
         row.snow_depth_p01_m = ...
            icemodel.verification.helpers.sampleQuantile(values, 0.01);
         row.snow_depth_p05_m = ...
            icemodel.verification.helpers.sampleQuantile(values, 0.05);
         row.snow_depth_p50_m = ...
            icemodel.verification.helpers.sampleQuantile(values, 0.50);
      end
   end

   % Primary support excludes every detected step because this audit applies no
   % de-step. A correctable classification is evidence, not a correction.
   % One owner applies every flag rule, so this cannot drift from the
   % comparator, the runner, or the report builder.
   flag_fields = icemodel.verification.helpers.observationSupportFields( ...
      target, policy);
   support = icemodel.verification.helpers.classifyObservationSupport( ...
      double(data{:, cellstr(flag_fields)}), flag_fields, target, policy);
   detected = support.unresolved_step;
   correctable = support.metadata_flagged;
   direct = support.flag_clean;
   datum_intact = support.datum_intact;
   row.gap_flag_count = nnz(support.gap_flagged);
   row.station_transition_flag_count = nnz(support.station_transition);
   row.step_detected_flag_count = nnz(detected);
   row.step_correctable_flag_count = nnz(correctable);
   row.unresolved_step_flag_count = nnz(detected);
   row.direct_target_count = nnz(direct);
   [season_start, season_end] = ...
      icemodel.verification.helpers.evaluationSeason( ...
      year(window_start), policy);
   in_season = data.Time >= season_start & data.Time <= season_end;
   season_data = data(in_season, :);
   season_direct = direct(in_season);
   season_datum_intact = datum_intact(in_season);
   thresholds = policy.snow_free_sensitivity_thresholds_m;
   windows = repmat(emptyWindow(), numel(thresholds), 1);
   for k = 1:numel(thresholds)
      windows(k) = longestWindow( ...
         season_data, season_direct, season_datum_intact, ...
         thresholds(k), policy);
      row.(thresholdField(thresholds(k), "window_days")) = windows(k).days;
      row.(thresholdField(thresholds(k), "direct_count")) = ...
         windows(k).direct_count;
   end
   primary_index = find(abs( ...
      thresholds - policy.snow_continuity_threshold_m) ...
      < eps(policy.snow_continuity_threshold_m), 1);
   if isempty(primary_index)
      error('icemodel:verification:promiceAblationReadiness:policyMismatch', ...
         'primary snow-free threshold is absent from sensitivity thresholds')
   end
   primary = windows(primary_index);
   row.snow_free_window_start = formatTime(primary.start);
   row.snow_free_window_end = formatTime(primary.end);
   row.snow_free_window_days = primary.days;
   row.snow_free_direct_count = primary.direct_count;
   row.snow_free_snow_count = primary.snow_count;

   % Require actual units, UTC hourly postings, a verified zero reference, and
   % the predeclared minimum long-window support.
   reasons = strings(0, 1);
   if row.surface_zone ~= "ablation"
      reasons(end + 1, 1) = "manifest surface zone is not ablation";
   end
   if row.target_units ~= policy.target_units
      reasons(end + 1, 1) = "target units are not " + policy.target_units;
   end
   if ~contains(lower(row.target_description), "lowering")
      reasons(end + 1, 1) = ...
         "target description does not identify geometric lowering";
   end
   coordinate_ok = height(data) > 0 ...
      && string(data.Time.TimeZone) == "UTC";
   if coordinate_ok
      day_start = dateshift(data.Time, 'start', 'day');
      elapsed_seconds = seconds(data.Time - day_start);
      coordinate_ok = all(mod(elapsed_seconds, ...
         policy.observation_cadence_seconds) == 0);
   end
   if height(data) > 1
      spacing_seconds = seconds(diff(data.Time));
      coordinate_ok = coordinate_ok && all(spacing_seconds > 0) ...
         && all(mod(spacing_seconds, ...
         policy.observation_cadence_seconds) == 0);
   end
   if ~coordinate_ok
      reasons(end + 1, 1) = ...
         "observation timestamps are not a strictly increasing native whole-hour UTC coordinate";
   end
   first_finite = find(isfinite(all_data.(target)), 1);
   if isempty(first_finite) || abs(all_data.(target)(first_finite)) ...
         > policy.zero_reference_tolerance_m
      reasons(end + 1, 1) = ...
         "target is not zero-referenced at its first finite sample";
   end
   if primary.days < policy.minimum_window_days
      reasons(end + 1, 1) = "longest data-based snow-free window is " ...
         + compose('%.2f', primary.days) + " days; at least " ...
         + policy.minimum_window_days + " days are required";
   end
   if primary.direct_count < policy.minimum_direct_samples
      reasons(end + 1, 1) = "snow-free window has fewer than " ...
         + policy.minimum_direct_samples + " direct target samples";
   end
   reason = strjoin(reasons, "; ");
end

function [data, reason] = normalizeObservationRowTimes(data)
   %NORMALIZEOBSERVATIONROWTIMES Convert one staged coordinate safely to UTC.
   reason = "";

   % Tag unzoned wall times as UTC and preserve zoned instants through the
   % shared staging normalizer before converting their display zone to UTC.
   try
      row_times = icemodel.verification.setup.ensureUtc( ...
         data.Properties.RowTimes);
      if ~isdatetime(row_times) || any(isnat(row_times), 'all')
         reason = ...
            "observation timestamps cannot be converted to a complete UTC datetime coordinate";
         return
      end
      row_times.TimeZone = 'UTC';
      data.Properties.RowTimes = row_times;
   catch
      % One malformed artifact excludes its readiness row without aborting the
      % cohort-level CSV/JSON audit for every other staged case.
      reason = ...
         "observation timestamps cannot be converted to a complete UTC datetime coordinate";
   end
end

function window = longestWindow( ...
      data, direct, datum_intact, threshold, policy)
   %LONGESTWINDOW Find the longest summer run bounded by exposed-ice data.
   window = emptyWindow();
   snow = data.(policy.snow_variable);
   [ice_exposed, ~, ~] = ...
      icemodel.verification.helpers.classifySnowDepth( ...
      snow, policy.ice_exposure_threshold_m);
   [~, snow_present, ~] = ...
      icemodel.verification.helpers.classifySnowDepth(snow, threshold);
   endpoint_ok = ice_exposed & direct & datum_intact;
   % Missing postings and unknown snow support are valid inside a cumulative
   % interval. Split only at a row carrying finite snow presence or failed
   % datum continuity; an omitted hourly row therefore cannot create a break.
   % Finite exposed-ice direct rows still define both endpoints.
   breaks = find(snow_present | ~datum_intact);
   bounds = [0; breaks; height(data) + 1];
   for k = 1:numel(bounds) - 1
      run = (bounds(k) + 1:bounds(k + 1) - 1).';
      supported = run(endpoint_ok(run));
      if numel(supported) < policy.minimum_direct_samples
         continue
      end
      candidate = emptyWindow();
      candidate.start = data.Time(supported(1));
      candidate.end = data.Time(supported(end));
      candidate.days = days(candidate.end - candidate.start);
      inside = data.Time >= candidate.start & data.Time <= candidate.end;
      candidate.direct_count = nnz(endpoint_ok & inside);
      candidate.snow_count = nnz(ice_exposed & inside);
      if candidate.days > window.days ...
            || (candidate.days == window.days ...
            && candidate.direct_count > window.direct_count)
         window = candidate;
      end
   end
end


function [tf, reason] = forcingCoverage(forcing, window_start, window_end)
   %FORCINGCOVERAGE Test whether one actual complete window encloses a request.
   if ~forcing.ok
      tf = false;
      reason = forcing.reason;
      return
   end
   tf = false;
   for k = 1:numel(forcing.windows)
      first = toUtc(forcing.windows(k).start_time);
      last = toUtc(forcing.windows(k).end_time);
      if first <= window_start && last >= window_end
         tf = true;
         break
      end
   end
   if tf
      reason = "";
   else
      reason = "actual promice_filled payload does not cover the requested window";
   end
end

function [tf, reason] = consumerStatus( ...
      forcing, y, consumer, coverage_ready, coverage_reason)
   %CONSUMERSTATUS Combine actual coverage with one producer consumer verdict.
   if ~coverage_ready
      tf = false;
      reason = coverage_reason;
      return
   end
   verdict_name = "verdict_" + consumer;
   reason_name = "reason_" + consumer;
   names = string(forcing.ledger.Properties.VariableNames);
   if ~all(ismember(["year", verdict_name, reason_name], names))
      tf = false;
      reason = "producer readiness ledger schema is incomplete";
      return
   end
   match = forcing.ledger.year == y;
   if nnz(match) ~= 1
      tf = false;
      reason = "producer readiness ledger has no unique row for year " + y;
      return
   end
   tf = lower(string(forcing.ledger.(verdict_name)(match))) == "ready";
   reason = string(forcing.ledger.(reason_name)(match));
   if ismissing(reason) || strlength(reason) == 0
      reason = "";
      if ~tf
         reason = "producer ledger supplied no reason";
      end
   end
end

function inventory = caseInventory(cases, evaluation_data_root, input_data_root)
   %CASEINVENTORY Detect missing and extra staged cases around the manifest.
   canonical = sort(lower(string({cases.case_id}))).';
   files = dir(fullfile(evaluation_data_root, 'promice', '*', ...
      'observations.mat'));
   observed = strings(numel(files), 1);
   for k = 1:numel(files)
      [~, observed(k)] = fileparts(files(k).folder);
   end
   data_root = string(fileparts(input_data_root));
   files = dir(fullfile(data_root, 'preview', 'qa', 'gapfill', 'plans', ...
      '*-report-inputs.json'));
   forced = erase(lower(string({files.name})).', ...
      "-report-inputs.json");
   inventory = struct( ...
      'canonical_case_count', numel(canonical), ...
      'missing_observation_cases', setdiff(canonical, observed), ...
      'extra_observation_cases', setdiff(observed, canonical), ...
      'missing_forcing_cases', setdiff(canonical, forced), ...
      'extra_forcing_cases', setdiff(forced, canonical));
end

function audit = reconstructionWriteAudit()
   %RECONSTRUCTIONWRITEAUDIT Inventory writes behind runtime path guards.
   root = fullfile(icemodel.internal.fullpath, 'icemodel', '+icemodel', ...
      '+forcing', '+reconstruct');
   files = dir(fullfile(root, '**', '*.m'));
   per_file = repmat({strings(0, 1)}, numel(files), 1);
   n_calls = 0;
   n_guard_calls = 0;
   for k = 1:numel(files)
      code = regexprep(fileread(fullfile(files(k).folder, files(k).name)), ...
         '%[^\n]*', '');
      calls = regexp(code, ['(?m)^\s*[^\n]*' ...
         '(save|writetable|writematrix|writecell|fopen)\s*\([^\n]*'], ...
         'match');
      n_calls = n_calls + numel(calls);
      if string(files(k).name) ~= "assertNotEvaluationDestination.m"
         n_guard_calls = n_guard_calls + numel(regexp(code, ...
            'assertNotEvaluationDestination\s*\(', 'match'));
      end
      % Select the offending calls for this file in one pass and stash them
      % in a preallocated slot, so the match list is joined once at the end.
      is_match = ~cellfun(@isempty, regexpi(calls, ...
         'observations\.mat|data[\\/]eval|evaluation_data_root', 'once'));
      if any(is_match)
         per_file{k} = transpose(string(files(k).name) + ": " ...
            + strtrim(string(calls(is_match))));
      end
   end
   matches = vertcat(per_file{:});
   verified = numel(files) > 0 && n_calls > 0 ...
      && n_guard_calls >= 2 && isempty(matches);
   audit = struct('reconstruction_source_file_count', numel(files), ...
      'reconstruction_write_call_count', n_calls, ...
      'runtime_evaluation_guard_call_count', n_guard_calls, ...
      'observation_write_match_count', numel(matches), ...
      'observation_write_matches', matches, ...
      'observation_write_path_verified', verified, ...
      'evidence', ...
         "canonical runtime guards reject reconstruction output, QA, and " ...
         + "split-manifest destinations under evaluation roots; the " ...
         + "supplemental non-comment write inventory has no explicit " ...
         + "observation target; current hashes are the forward baseline");
end

function row = rowTemplate(policy)
   %ROWTEMPLATE Define the stable public CSV schema and column order.
   row = struct( ...
      'case_id', "", 'site_id', "", 'surface_zone', "", 'year', 0, ...
      'requested_window_start', "", 'requested_window_end', "", ...
      'forcing_producer_manifest', "", ...
      'forcing_producer_manifest_sha256', "", ...
      'forcing_readiness_artifact', "", ...
      'forcing_readiness_sha256', "", ...
      'forcing_artifact', "", 'forcing_sha256', "", ...
      'forcing_acceptance_start', "", 'forcing_acceptance_end', "", ...
      'observation_artifact', "", 'observation_sha256', "", ...
      'target_variable', "", 'target_available', false, ...
      'target_units', "", 'target_positive_direction', "", ...
      'target_description', "", 'timestamp_convention', "", ...
      'interval_convention', "", 'reference_semantics', "", ...
      'adjustment_semantics', "", 'possible_support_count', 0, ...
      'observation_sample_count', 0, 'finite_target_count', 0, ...
      'finite_snow_count', 0, 'direct_target_count', 0, ...
      'gap_flag_count', 0, 'station_transition_flag_count', 0, ...
      'step_detected_flag_count', 0, 'step_correctable_flag_count', 0, ...
      'unresolved_step_flag_count', 0, ...
      'snow_depth_min_m', NaN, 'snow_depth_p01_m', NaN, ...
      'snow_depth_p05_m', NaN, 'snow_depth_p50_m', NaN, ...
      'snow_free_window_start', "", 'snow_free_window_end', "", ...
      'snow_free_window_days', 0, 'snow_free_direct_count', 0, ...
      'snow_free_snow_count', 0, 'ice_model_forcing_ready', false, ...
      'ice_model_forcing_reason', "", 'evaluation_target_ready', false, ...
      'evaluation_target_reason', "", 'snow_input_ready', false, ...
      'snow_input_reason', "", ...
      'production_snow_candidate_ready', false, ...
      'production_snow_candidate_reason', "", 'admitted', false, ...
      'exclusion_reason', "");
   for threshold = policy.snow_free_sensitivity_thresholds_m
      row.(thresholdField(threshold, "window_days")) = 0;
      row.(thresholdField(threshold, "direct_count")) = 0;
   end
end

function window = emptyWindow()
   %EMPTYWINDOW Return one no-support window record.
   window = struct('start', NaT(1, 1, 'TimeZone', 'UTC'), ...
      'end', NaT(1, 1, 'TimeZone', 'UTC'), 'days', 0, ...
      'direct_count', 0, 'snow_count', 0);
end

function name = thresholdField(threshold, suffix)
   %THRESHOLDFIELD Encode one metre threshold as a stable millimetre field.
   name = char("snow_free_" + compose('%03d', round(1000 * threshold)) ...
      + "mm_" + suffix);
end


function [first, last] = annualWindow(y, policy)
   %ANNUALWINDOW Return the full initialization year for one readiness row.
   first = datetime(y, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
   last = datetime(y + 1, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      - seconds(policy.observation_cadence_seconds);
end

function value = toUtc(value)
   %TOUTC Reuse the staged-data UTC normalizer and convert zoned instants.
   if (isstring(value) || ischar(value)) ...
         && endsWith(string(value), "Z")
      value = datetime(string(value), 'InputFormat', ...
         "yyyy-MM-dd'T'HH:mm:ss'Z'", 'TimeZone', 'UTC');
   else
      value = icemodel.verification.setup.ensureUtc(value);
   end
   value.TimeZone = 'UTC';
end

function text = formatTime(value)
   %FORMATTIME Reuse the canonical manifest formatter for ledger timestamps.
   text = string(icemodel.verification.setup.formatManifestTime(value));
end
