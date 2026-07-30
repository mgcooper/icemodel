function report = auditArtifacts(kwargs)
   %AUDITARTIFACTS Read-only QA/QC for manifest-referenced verification artifacts.
   %
   %  report = icemodel.verification.auditArtifacts( ...
   %     evaluation_data_root=eval_root, input_data_root=input_root, ...
   %     families=["promice", "esm_snowmip", "laugh_tests"])
   %  report = ... auditArtifacts(..., report_dir=output_dir)
   %
   % The audit reads artifacts referenced by the selected family manifests plus
   % the exact runtime-resolved met artifact for atomic ESM-SnowMIP cases. It
   % checks manifest/schema consistency, artifact payload and metadata shape,
   % time axes and periods, canonical names and units, physical ranges,
   % placeholder semantics, contiguous missing runs, immutable file identity,
   % and source-specific MERRA, MAR, RACMO, MODIS, and PROMICE contracts. Returned
   % records are deterministic struct arrays suitable for plotting or report
   % generation. Supplying REPORT_DIR is the only write path; it emits
   % artifact_qa.json and artifact_qa.md without changing staged data.
   %
   % Inputs
   %  data_root              Whole data tree containing eval/ and input/.
   %  evaluation_data_root   Explicit eval/ root. Blank uses the verification
   %                         helper's repo-local default.
   %  input_data_root        Explicit input/ root paired with the eval tree.
   %  icemodel_config_casename   Optional nonmutating config-case selector.
   %  families               Canonical dataset-family ids or "all".
   %  report_dir             Optional generated-output directory.
   %
   % Outputs
   %  report   Struct with summary, artifact, channel, and finding records.
   %
   % See also: icemodel.verification.comparisonCompatibility,
   %  icemodel.verification.helpers.esmRuntimeMetFiles,
   %  icemodel.forcing.helpers.metchecks,
   %  icemodel.verification.setup.repairRcmArtifactMetadata

   arguments
      kwargs.data_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
      kwargs.families (1, :) string ...
         {icemodel.verification.validators.mustBeDatasetFamilySelection} = "all"
      kwargs.report_dir (1, 1) string = ""
   end

   % Resolve roots without mutating global configuration. An explicit input root
   % is carried through listcases/comparisonCompatibility for portable manifests.
   [eval_root, input_root] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      data_root=kwargs.data_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      input_data_root=kwargs.input_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);

   % Discover and validate manifests before using the operational readers. This
   % lets malformed or missing manifests become findings instead of aborting the
   % entire multi-family audit.
   manifest_files = icemodel.verification.helpers.familyManifestFiles( ...
      evaluation_data_root=eval_root);
   families = selectedFamilies(kwargs.families, manifest_files);
   artifacts = repmat(emptyArtifact(), 0, 1);
   channels = repmat(emptyChannel(), 0, 1);
   findings = repmat(emptyFinding(), 0, 1);
   case_count = 0;

   for family = reshape(families, 1, [])
      manifest_path = familyManifestPath(family, manifest_files);
      if manifest_path == ""
         findings(end + 1) = finding("error", "missing_manifest", family, ...
            "", "", "manifest", fullfile(eval_root, family, "manifest.json"), ...
            "", "selected family has no manifest.json"); %#ok<AGROW>
         continue
      end

      % Keep the manifest itself in the artifact ledger and validate its portable
      % schema before downstream path resolution.
      manifest_record = emptyArtifact();
      manifest_record.dataset_family = family;
      manifest_record.kind = "manifest";
      manifest_record.path = manifest_path;
      manifest_record.exists = true;
      manifest_record.status = "checked";
      artifacts(end + 1) = manifest_record; %#ok<AGROW>
      try
         % Record the exact manifest bytes before parsing so a later report can
         % reject a manifest that changed after this audit pass.
         manifest_info = dir(manifest_path);
         artifacts(end).artifact_size_bytes = manifest_info.bytes;
         artifacts(end).artifact_sha256 = ...
            icemodel.verification.setup.fileSha256(manifest_path);
         manifest = icemodel.verification.helpers.readFamilyManifest(manifest_path);
      catch err
         artifacts(end).status = "read_error";
         findings(end + 1) = finding("error", "manifest_read_error", family, ...
            "", "", "manifest", manifest_path, "", string(err.message)); %#ok<AGROW>
         continue
      end
      findings = appendStructs(findings, validateManifest(manifest));

      % Reuse the existing resolved-case and comparison discovery paths. The
      % latter owns portable path resolution and exact manifest artifact lists.
      try
         cases = icemodel.verification.listcases( ...
            evaluation_data_root=eval_root, input_data_root=input_root, ...
            dataset_family=family);
         compatibility = icemodel.verification.comparisonCompatibility( ...
            evaluation_data_root=eval_root, input_data_root=input_root, ...
            dataset_family=family);
      catch err
         findings(end + 1) = finding("error", "manifest_resolution_error", ...
            family, "", "", "manifest", manifest_path, "", ...
            string(err.message)); %#ok<AGROW>
         continue
      end
      case_count = case_count + numel(cases);

      % Audit source-list declarations separately from file content so a stale
      % manifest cannot be hidden by a valid artifact left elsewhere on disk.
      for k = 1:numel(cases)
         findings = appendStructs(findings, validateCaseSources(cases(k)));
      end

      % Inspect every concrete artifact discovered by comparisonCompatibility.
      % Append exact manifest met references excluded from comparison pairing and
      % the standard runtime met fallback for atomic ESM-SnowMIP cases.
      for k = 1:numel(compatibility.cases)
         case_report = compatibility.cases(k);
         case_manifest = caseForReport(cases, case_report);
         discovered = case_report.artifacts(:);
         [omitted_met, met_resolution_findings] = ...
            omittedMetArtifacts(case_manifest, discovered, input_root);
         discovered = appendStructs(discovered, omitted_met);
         findings = appendStructs(findings, met_resolution_findings);
         for n = 1:numel(discovered)
            spec = discovered(n);

            % Compatibility discovery uses an empty metadata record for a
            % staged status-only association. It is evidence about a source,
            % not a file that can be opened or reported missing.
            if string(spec.kind) == "metadata" ...
                  && strlength(string(spec.path)) == 0
               continue
            end
            spec.path = canonicalMissingPath(spec, case_manifest, input_root);
            period = artifactPeriod(case_manifest, spec.source, spec.kind);
            [artifact_record, channel_records, artifact_findings] = ...
               inspectArtifact(case_manifest, spec, period);
            artifacts(end + 1) = artifact_record; %#ok<AGROW>
            channels = appendStructs(channels, channel_records);
            findings = appendStructs(findings, artifact_findings);
         end
      end
   end

   % Stable ordering makes JSON diffs and downstream report ingestion repeatable.
   artifacts = sortRecords(artifacts, ...
      ["dataset_family", "case_id", "source", "kind", "path"]);
   channels = sortRecords(channels, ...
      ["dataset_family", "case_id", "source", "kind", "path", ...
      "table_path", "channel"]);
   findings = sortFindings(findings);

   % Keep the public result self-contained and machine serializable.
   generated = datetime("now", TimeZone="UTC", ...
      Format="yyyy-MM-dd'T'HH:mm:ss'Z'");
   report = struct();
   report.schema_version = "1.0";
   report.generated_at = string(generated);
   report.evaluation_data_root = eval_root;
   report.input_data_root = input_root;
   report.families = families(:);
   report.artifacts = artifacts;
   report.channels = channels;
   report.findings = findings;
   report.summary = summarizeReport(families, case_count, artifacts, ...
      channels, findings);
   report.report_files = struct('json', "", 'markdown', "");

   % Writing is explicit and confined to the caller-selected generated-output
   % directory. Staged manifests and artifacts remain read-only.
   if ~isblanktext(kwargs.report_dir)
      report = writeReports(report, kwargs.report_dir);
   end
end

%% Manifest discovery and schema checks
function families = selectedFamilies(requested, manifest_files)
   %SELECTEDFAMILIES Expand "all" to the families actually present on disk.
   requested = unique(reshape(string(requested), 1, []), 'stable');
   if isempty(requested) || all(isblanktext(requested)) || ismember("all", requested)
      families = strings(numel(manifest_files), 1);
      for k = 1:numel(manifest_files)
         [~, families(k)] = fileparts(fileparts(manifest_files(k)));
      end
      families = unique(families, 'stable');
   else
      families = requested(:);
   end
end

function pathname = familyManifestPath(family, manifest_files)
   %FAMILYMANIFESTPATH Return the discovered manifest for one family.
   pathname = "";
   for candidate = reshape(manifest_files, 1, [])
      [~, candidate_family] = fileparts(fileparts(candidate));
      if string(candidate_family) == family
         pathname = candidate;
         return
      end
   end
end

function findings = validateManifest(manifest)
   %VALIDATEMANIFEST Check required family/case fields and period values.
   findings = repmat(emptyFinding(), 0, 1);
   family = string(fieldOr(manifest, 'dataset_family', ""));

   % Family fields are centralized in the setup schema helper.
   required = string(icemodel.verification.setup.familyManifestFieldNames());
   missing = setdiff(required, string(fieldnames(manifest)), 'stable');
   for name = reshape(missing, 1, [])
      findings(end + 1) = finding("error", "missing_manifest_field", ...
         family, "", "", "manifest", string(fieldOr(manifest, ...
         'manifest_path', "")), name, "missing required family field"); %#ok<AGROW>
   end
   if ~isfield(manifest, 'cases') || ~isstruct(manifest.cases)
      findings(end + 1) = finding("error", "invalid_manifest_cases", ...
         family, "", "", "manifest", string(fieldOr(manifest, ...
         'manifest_path', "")), "cases", "cases must be a struct array");
      return
   end

   % Validate each case against the appropriate existing schema source.
   ids = strings(numel(manifest.cases), 1);
   for k = 1:numel(manifest.cases)
      c = manifest.cases(k);
      case_id = string(fieldOr(c, 'case_id', ""));
      ids(k) = case_id;
      if string(fieldOr(c, 'case_type', "")) == "firn_observational"
         case_required = icemodel.verification.setup.firnCaseManifestFieldNames();
      else
         case_required = icemodel.verification.setup.caseManifestFieldNames();
      end
      case_missing = setdiff(case_required, string(fieldnames(c)), 'stable');
      for name = reshape(case_missing, 1, [])
         findings(end + 1) = finding("error", "missing_case_field", ...
            family, case_id, "", "manifest", string(manifest.manifest_path), ...
            name, "missing required case field"); %#ok<AGROW>
      end
      findings = appendStructs(findings, validateCasePeriod(c, family, ...
         string(manifest.manifest_path)));
   end

   % Duplicate ids make loadmanifest's first-match behavior ambiguous.
   unique_ids = unique(ids(ids ~= ""));
   for id = reshape(unique_ids, 1, [])
      if sum(ids == id) > 1
         findings(end + 1) = finding("error", "duplicate_case_id", ...
            family, id, "", "manifest", string(manifest.manifest_path), ...
            "case_id", "case_id appears more than once in this manifest"); %#ok<AGROW>
      end
   end
end

function findings = validateCasePeriod(c, family, manifest_path)
   %VALIDATECASEPERIOD Require a parseable, ordered case period.
   findings = repmat(emptyFinding(), 0, 1);
   case_id = string(fieldOr(c, 'case_id', ""));
   if ~isfield(c, 'period') || ~isstruct(c.period) ...
         || ~all(isfield(c.period, ["start", "end"]))
      findings(end + 1) = finding("error", "missing_period", family, ...
         case_id, "", "manifest", manifest_path, "period", ...
         "case period must contain start and end");
      return
   end
   try
      [t1, t2] = icemodel.verification.setup.periodBounds(c.period);
      valid = ~isnat(t1) && ~isnat(t2) && t1 <= t2;
   catch
      valid = false;
   end
   if ~valid
      findings(end + 1) = finding("error", "invalid_period", family, ...
         case_id, "", "manifest", manifest_path, "period", ...
         "case period is not a parseable ordered UTC interval");
   end
end

function findings = validateCaseSources(c)
   %VALIDATECASESOURCES Check source lists against staged colocation declarations.
   findings = repmat(emptyFinding(), 0, 1);
   family = string(c.dataset_family);
   case_id = string(c.case_id);
   manifest_path = string(c.manifest_path);
   colocation = struct();
   if isfield(c, 'colocation') && isstruct(c.colocation)
      colocation = c.colocation;
   end

   % Every declared forcing/evaluation source must resolve to a concrete leg or
   % the case-level evaluation artifact.
   forcing_sources = stringArray(fieldOr(c, 'forcing_sources', strings(0, 1)));
   for source = reshape(forcing_sources, 1, [])
      [found, leg] = sourceLeg(colocation, source);
      has_met = found && hasFiles(leg, "met_files");
      if ~has_met
         findings(end + 1) = finding("error", "forcing_source_without_met", ...
            family, case_id, source, "manifest", manifest_path, ...
            "forcing_sources", "declared forcing source has no staged met file"); %#ok<AGROW>
      end
   end

   eval_sources = stringArray(fieldOr(c, 'eval_sources', strings(0, 1)));
   for source = reshape(eval_sources, 1, [])
      base_source = erase(erase(source, "_obs"), "_protocol");
      [found, leg] = sourceLeg(colocation, source);
      if ~found
         [found, leg] = sourceLeg(colocation, base_source);
      end
      case_eval = endsWith(source, "_obs") ...
         && strlength(string(fieldOr(c, 'evaluation_path', ""))) > 0;
      leg_eval = found && any([hasFiles(leg, "obs_file"), ...
         hasFiles(leg, "evaluation_file"), hasFiles(leg, "data_files"), ...
         hasFiles(leg, "model_output_files"), hasFiles(leg, "reference_file")]);
      if ~case_eval && ~leg_eval
         findings(end + 1) = finding("error", "eval_source_without_artifact", ...
            family, case_id, source, "manifest", manifest_path, ...
            "eval_sources", "declared evaluation source has no artifact reference"); %#ok<AGROW>
      end
   end

   % A staged leg must carry at least one file reference; an unstaged leg must not
   % claim files. This catches stale source-list/file drift directly.
   names = string(fieldnames(colocation));
   for name = reshape(names, 1, [])
      leg = colocation.(char(name));
      if ~isstruct(leg)
         continue
      end
      staged = logicalScalar(fieldOr(leg, 'staged', false));
      has_any = any([hasFiles(leg, "obs_file"), hasFiles(leg, "evaluation_file"), ...
         hasFiles(leg, "met_files"), hasFiles(leg, "data_files"), ...
         hasFiles(leg, "model_output_files"), hasFiles(leg, "reference_file")]);

      % RETMIP's native_met leg records readiness for the native source while
      % the owning retmip leg carries the actual met/data file paths. A status-
      % only leg has no artifact-ownership contract to violate.
      status_only = ~isfield(leg, 'kind') && isfield(leg, 'status');
      if staged && ~has_any && ~status_only
         findings(end + 1) = finding("error", "staged_leg_without_artifact", ...
            family, case_id, name, "manifest", manifest_path, ...
            "colocation", "staged colocation leg has no file reference"); %#ok<AGROW>
      elseif ~staged && has_any
         findings(end + 1) = finding("error", "unstaged_leg_with_artifact", ...
            family, case_id, name, "manifest", manifest_path, ...
            "colocation", "unstaged colocation leg still references files"); %#ok<AGROW>
      end
   end
end

%% Artifact discovery and inspection
function c = caseForReport(cases, case_report)
   %CASEFORREPORT Match a compatibility row to its resolved case manifest.
   family = string({cases.dataset_family});
   ids = string({cases.case_id});
   idx = find(family == string(case_report.dataset_family) ...
      & ids == string(case_report.case_id), 1);
   c = cases(idx);
end

function [artifacts, findings] = omittedMetArtifacts(c, discovered, input_root)
   %OMITTEDMETARTIFACTS Add exact met artifacts omitted from comparison pairing.
   artifacts = repmat(discovered(1), 0, 1);
   findings = repmat(emptyFinding(), 0, 1);

   % Comparison discovery intentionally omits RCM met files. Preserve manifest
   % authority by adding only each staged leg's exact recorded paths.
   if isfield(c, 'colocation') && isstruct(c.colocation)
      names = string(fieldnames(c.colocation));
      for name = reshape(names, 1, [])
         leg = c.colocation.(char(name));
         if ~isstruct(leg) || ~logicalScalar(fieldOr(leg, 'staged', false)) ...
               || ~hasFiles(leg, "met_files")
            continue
         end
         files = stringArray(leg.met_files);
         source = sourceLabel(name, leg);
         for file = reshape(files, 1, [])
            known = [string({discovered.path}), string({artifacts.path})];
            same_path = any(known == file) ...
               || any(endsWith(known, filesep + file)) ...
               || any(endsWith(file, filesep + known));
            if same_path
               continue
            end
            a = discovered(1);
            a.source = source;
            a.role = "data";
            a.kind = "met";
            a.path = file;
            a.exists = false;
            a.artifact_variables = strings(1, 0);
            a.declared_variables = strings(1, 0);
            a.variables = strings(1, 0);
            a.evidence = "manifest";
            artifacts = appendStructs(artifacts, a);
         end
      end
   end

   % Atomic ESM manifests deliberately have no forcing leg. Resolve their one
   % runnable met artifact through the same standard chain as a model run, but
   % never override a future explicit met declaration.
   found_kinds = [string({discovered.kind}), string({artifacts.kind})];
   if string(c.dataset_family) ~= "esm_snowmip" ...
         || any(found_kinds == "met")
      return
   end
   try
      files = icemodel.verification.helpers.esmRuntimeMetFiles(c, input_root);
   catch err
      findings(end + 1) = finding("error", ...
         "runtime_met_resolution_error", string(c.dataset_family), ...
         string(c.case_id), "esm_snowmip", "met", ...
         fullfile(input_root, "met"), "", ...
         string(err.identifier) + ": " + string(err.message));
      return
   end
   for file = reshape(files, 1, [])
      known = [string({discovered.path}), string({artifacts.path})];
      same_path = any(known == file) ...
         || any(endsWith(known, filesep + file)) ...
         || any(endsWith(file, filesep + known));
      if same_path
         continue
      end
      a = discovered(1);
      a.source = "esm_snowmip";
      a.role = "data";
      a.kind = "met";
      a.path = file;
      a.exists = false;
      a.artifact_variables = strings(1, 0);
      a.declared_variables = strings(1, 0);
      a.variables = strings(1, 0);
      a.evidence = "runtime";
      artifacts = appendStructs(artifacts, a);
   end
end

function pathname = canonicalMissingPath(spec, c, input_root)
   %CANONICALMISSINGPATH Resolve a still-relative missing artifact canonically.
   pathname = string(spec.path);
   if pathname == "" || isfile(pathname) || startsWith(pathname, filesep)
      return
   end
   if spec.kind == "met"
      pathname = fullfile(input_root, "met", pathname);
   elseif spec.kind == "userdata"
      pathname = fullfile(input_root, "userdata", pathname);
   elseif isfield(c, 'family_root')
      pathname = fullfile(string(c.family_root), pathname);
   end
end

function period = artifactPeriod(c, source, kind)
   %ARTIFACTPERIOD Prefer a source-leg window over the broader case period.
   period = fieldOr(c, 'period', struct('start', "", 'end', ""));

   % Exact-date model profiles and protocol outputs are evaluation products.
   % Their dates follow the observation/case interval, not the separately
   % staged forcing leg's reusable coverage window.
   if string(kind) == "model_output"
      return
   end
   if ~isfield(c, 'colocation') || ~isstruct(c.colocation)
      return
   end
   [found, leg] = sourceLeg(c.colocation, source);
   if found && isfield(leg, 'window') && isstruct(leg.window)
      period = leg.window;
   end
end

function [record, channels, findings] = inspectArtifact(c, spec, period)
   %INSPECTARTIFACT Read one artifact using its native file format.
   record = emptyArtifact();
   record.dataset_family = string(c.dataset_family);
   record.case_id = string(c.case_id);
   record.source = string(spec.source);
   record.kind = string(spec.kind);
   record.path = string(spec.path);
   record.expected_start = string(fieldOr(period, 'start', ""));
   record.expected_end = string(fieldOr(period, 'end', ""));
   record.exists = isfile(record.path);
   channels = repmat(emptyChannel(), 0, 1);
   findings = repmat(emptyFinding(), 0, 1);

   if ~record.exists
      record.status = "missing";
      findings(end + 1) = finding("error", "missing_artifact", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "", "required staged artifact does not exist");
      return
   end

   % Hash the exact staged bytes before format-specific inspection. The report
   % can therefore detect any replacement or edit made after this audit pass.
   try
      file_info = dir(record.path);
      record.artifact_size_bytes = file_info.bytes;
      record.artifact_sha256 = ...
         icemodel.verification.setup.fileSha256(record.path);
   catch err
      record.status = "read_error";
      findings(end + 1) = finding("error", "artifact_read_error", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "", string(err.message));
      return
   end

   % RetMIP protocol outputs are native NetCDF products, not MAT containers.
   [~, ~, extension] = fileparts(record.path);
   if strcmpi(extension, '.nc')
      [record, findings] = inspectNetcdfArtifact(record);
      return
   end

   % Restrict MAT loading to the recognized primary payload plus metadata.
   try
      info = whos('-file', record.path);
      names = string({info.name});
      payload = primaryPayload(record.kind, names);
      if payload == ""
         record.status = "invalid_payload";
         findings(end + 1) = finding("error", "missing_payload", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, "", "MAT file has no recognized verification payload");
         return
      end
      load_names = payload;
      if ismember("artifact_metadata", names)
         load_names(end + 1) = "artifact_metadata";
      end
      loaded = load(record.path, load_names{:});
      value = loaded.(char(payload));
   catch err
      record.status = "read_error";
      findings(end + 1) = finding("error", "artifact_read_error", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "", string(err.message));
      return
   end
   record.payload = payload;
   bundle_format = artifactFormat(value);

   % Writers duplicate table UserData at top level. Require exact parity for met
   % and Data payloads so lightweight readers see the same provenance.
   if ismember(payload, ["met", "Data"])
      [metadata, metadata_findings] = artifactMetadata(loaded, value, record);
      findings = appendStructs(findings, metadata_findings);
   else
      metadata = struct();
   end

   % Reuse the canonical met validator so required channels, fixed cadence, and
   % precipitation units stay aligned with the runtime contract.
   if payload == "met"
      try
         icemodel.forcing.helpers.validatemet(value)
      catch err
         findings(end + 1) = finding("error", "met_contract", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, "", string(err.identifier) + ": " + string(err.message));
      end
   end

   % Recursively inspect every table/timetable inside target/reference bundles.
   nodes = tableNodes(value, payload);
   if isempty(nodes)
      findings(end + 1) = finding("error", "payload_without_tables", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "", "recognized payload contains no table or timetable");
      record.status = "invalid_payload";
      return
   end
   all_times = NaT(0, 1, TimeZone="UTC");
   first_cadence = NaN;
   for k = 1:numel(nodes)
      [node_channels, node_findings, times, cadence] = ...
         inspectTable(nodes(k).value, nodes(k).path, record, metadata, ...
         bundle_format);
      channels = appendStructs(channels, node_channels);
      findings = appendStructs(findings, node_findings);
      if ~isempty(times)
         utc_times = icemodel.verification.setup.ensureUtc(times(:));
         all_times = [all_times; utc_times(~isnat(utc_times))]; %#ok<AGROW>
      end
      if isnan(first_cadence) && ~isnan(cadence)
         first_cadence = cadence;
      end
   end

   % Optional MAR profile sidecars have a stricter scientific contract than a
   % generic nested table: exact profile/date identity, units, depth ordering,
   % and source provenance are all required.
   if record.kind == "model_output" && isstruct(value) ...
         && isfield(value, 'format') ...
         && string(value.format) == "subsurface_profile_bundle"
      findings = appendStructs(findings, ...
         modelOutputProfileChecks(value, record));
   end
   record.n_tables = numel(nodes);
   record.n_channels = numel(channels);
   record.n_samples = sum([nodes.n_rows]);
   record.cadence_seconds = first_cadence;
   if ~isempty(all_times)
      record.time_start = timeString(min(all_times));
      record.time_end = timeString(max(all_times));
      findings = appendStructs(findings, ...
         validateArtifactPeriod(record, all_times, period));
   end

   % Source policies operate on the primary timetable and its synchronized
   % metadata. Target/reference bundles are covered by generic table checks.
   if istimetable(value)
      findings = appendStructs(findings, sourceChecks(value, metadata, record));
   end
   record.status = artifactStatus(findings);
end

function [record, findings] = inspectNetcdfArtifact(record)
   %INSPECTNETCDFARTIFACT Validate a native protocol NetCDF header.
   findings = repmat(emptyFinding(), 0, 1);
   try
      info = ncinfo(record.path);
   catch err
      record.status = "read_error";
      findings(end + 1) = finding("error", "artifact_read_error", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "", string(err.message));
      return
   end

   record.payload = "netcdf";
   record.n_channels = numel(info.Variables);
   names = string({info.Variables.Name});
   sizes = [info.Variables.Size];
   if ~isempty(sizes)
      record.n_samples = max(sizes, [], 'all');
   end

   % Participant files preserve their native time conventions: some use t,
   % some expose only a time dimension, and some omit coordinate attributes.
   % The portable protocol contract is therefore the declared column/value
   % payload, not a canonicalized time coordinate.
   if record.dataset_family == "retmip" && record.kind == "model_output"
      [~, basename] = fileparts(record.path);
      if contains(lower(basename), "_columns")
         expected = ["rho", "depth", "temp", "lwc", "icecon", "rfrz", "dz"];
      elseif contains(lower(basename), "_values")
         expected = ["trunoff", "trfrz", "tfac", "tlwc", ...
            "T10m", "z550", "z830", "FAC20m"];
      else
         expected = ["rho", "depth", "temp", "lwc", "trunoff", "tfac"];
      end
      if ~any(ismember(names, expected))
         findings(end + 1) = finding("error", "netcdf_protocol_variables", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, "", ...
            "NetCDF lacks the expected RetMIP column/value variables");
      end
   end
   record.status = artifactStatus(findings);
end

function format = artifactFormat(value)
   %ARTIFACTFORMAT Return a nested bundle's declared public format.
   format = "";
   if isstruct(value) && isscalar(value) && isfield(value, 'format')
      format = string(value.format);
   elseif isstruct(value) && isscalar(value) && isfield(value, 'data') ...
         && isstruct(value.data) && isscalar(value.data) ...
         && isfield(value.data, 'format')
      format = string(value.data.format);
   end
end

function findings = modelOutputProfileChecks(bundle, record)
   %MODELOUTPUTPROFILECHECKS Validate MAR fixed-depth density sidecars.
   findings = repmat(emptyFinding(), 0, 1);
   if ~isfield(bundle, 'data') || ~isstruct(bundle.data) ...
         || ~isfield(bundle.data, 'density') ...
         || ~istable(bundle.data.density)
      findings(end + 1) = finding("error", "profile_payload", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "density", "profile bundle lacks a density table");
      return
   end

   density = bundle.data.density;
   names = string(density.Properties.VariableNames);
   required = ["profile_id", "datetime", "depth", "density", ...
      "source_id", "source_variable", "source_file", "sample_method"];
   missing = required(~ismember(required, names));
   if ~isempty(missing)
      findings(end + 1) = finding("error", "profile_schema", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "density", "missing profile fields: " + strjoin(missing, ','));
      return
   end

   depth_index = find(names == "depth", 1);
   density_index = find(names == "density", 1);
   units = string(density.Properties.VariableUnits);
   if numel(units) < width(density) || units(depth_index) ~= "m" ...
         || units(density_index) ~= "kg m-3"
      findings(end + 1) = finding("error", "profile_units", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "density", "profile depth/density units are not m and kg m-3");
   end
   if ~isdatetime(density.datetime) ...
         || any(string(density.datetime.TimeZone) ~= "UTC")
      findings(end + 1) = finding("error", "profile_time", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "density", "profile datetime must be UTC");
   end
   if any(string(density.source_id) ~= "mar3.11") ...
         || any(string(density.source_variable) ~= "RO1") ...
         || any(string(density.sample_method) ~= "nearest") ...
         || any(strlength(string(density.source_file)) == 0)
      findings(end + 1) = finding("error", "profile_provenance", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "density", "MAR RO1 nearest-cell provenance is incomplete");
   end

   groups = icemodel.verification.helpers.profileGroups(density);
   for group = reshape(groups, 1, [])
      depth = double(group.data.depth);
      values = double(group.data.density);
      if any(diff(depth) <= 0) || any(depth < 0 | depth > 20) ...
            || any(~isfinite(values) | values < 250 | values > 1000)
         findings(end + 1) = finding("error", "profile_values", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, "density", ...
            "profile depths or densities violate the MAR RO1 contract"); %#ok<AGROW>
         return
      end
   end
end

function payload = primaryPayload(kind, names)
   %PRIMARYPAYLOAD Choose the one public payload appropriate for an artifact kind.
   if kind == "met"
      preferred = "met";
   elseif kind == "userdata"
      preferred = "Data";
   elseif ismember(kind, ["reference", "model_output"])
      preferred = ["reference", "targets", "observations", "Data"];
   else
      preferred = ["targets", "observations", "evaluation", "Data", "reference"];
   end
   hit = intersect(preferred, names, 'stable');
   if isempty(hit)
      payload = "";
   else
      payload = hit(1);
   end
end

function [metadata, findings] = artifactMetadata(loaded, value, record)
   %ARTIFACTMETADATA Require synchronized top-level and UserData provenance.
   findings = repmat(emptyFinding(), 0, 1);
   metadata = struct();
   if isprop(value.Properties, 'UserData') && isstruct(value.Properties.UserData)
      metadata = value.Properties.UserData;
   end
   if ~isfield(loaded, 'artifact_metadata') || ~isstruct(loaded.artifact_metadata)
      findings(end + 1) = finding("error", "missing_artifact_metadata", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "", "met/Data payload lacks top-level artifact_metadata");
      return
   end
   if isempty(metadata)
      findings(end + 1) = finding("error", "missing_userdata_metadata", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "", "payload Properties.UserData lacks provenance metadata");
      return
   end
   if ~metadataEqual(metadata, loaded.artifact_metadata)
      findings(end + 1) = finding("error", "artifact_metadata_mismatch", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "", ...
         "artifact_metadata and payload Properties.UserData differ");
   end
end

function nodes = tableNodes(value, label)
   %TABLENODES Recursively collect table payloads without descending into arrays.
   prototype = struct('path', "", 'value', [], 'n_rows', 0);
   nodes = repmat(prototype, 0, 1);
   if istable(value) || istimetable(value)
      node = prototype;
      node.path = string(label);
      node.value = value;
      node.n_rows = height(value);
      nodes = node;
      return
   end
   if isstruct(value)
      for n = 1:numel(value)
         fields = string(fieldnames(value(n)));
         for field = reshape(fields, 1, [])
            child_label = string(label) + "." + field;
            nodes = appendStructs(nodes, ...
               tableNodes(value(n).(char(field)), child_label));
         end
      end
   elseif iscell(value)
      for n = 1:numel(value)
         nodes = appendStructs(nodes, ...
            tableNodes(value{n}, string(label) + "{" + n + "}"));
      end
   end
end

function [channels, findings, times, cadence] = inspectTable(T, table_path, ...
      record, metadata, bundle_format)
   %INSPECTTABLE Audit one table's time axis, channel metadata, and values.
   channels = repmat(emptyChannel(), 0, 1);
   findings = repmat(emptyFinding(), 0, 1);
   times = tableTimes(T);
   cadence = NaN;
   table_names = string(T.Properties.VariableNames);
   is_event_bundle = bundle_format == "subsurface_profile_bundle" ...
      && ~ismember(record.kind, ["met", "userdata"]);

   % Profile/event observations legitimately repeat survey dates across depth
   % rows and need not arrive in chronological record order. Forcing and true
   % timeseries axes retain the strict monotonic/cadence contract.
   if ~isempty(times)
      if is_event_bundle
         time_findings = validateEventTimes(T, times, record, table_path);
      else
         [time_findings, cadence] = validateTimeAxis(times, record, table_path);
      end
      findings = appendStructs(findings, time_findings);
   end

   names = table_names;
   units = tableUnits(T, numel(names));
   numeric_mask = false(size(names));
   for k = 1:numel(names)
      numeric_mask(k) = isnumeric(T.(names(k))) || islogical(T.(names(k)));
   end

   % Reuse metchecks in read-only mode for the standard missing/complex counts.
   checks = struct();
   numeric_names = names(numeric_mask);
   if istimetable(T) && ~isempty(numeric_names)
      [~, checks] = icemodel.forcing.helpers.metchecks(T(:, numeric_names), ...
         fillgaps=false, clamp=false);
   end

   for k = 1:numel(names)
      name = names(k);
      data = T.(name);
      if ~numeric_mask(k)
         continue
      end
      channel = emptyChannel();
      channel.dataset_family = record.dataset_family;
      channel.case_id = record.case_id;
      channel.source = record.source;
      channel.kind = record.kind;
      channel.path = record.path;
      channel.table_path = table_path;
      channel.channel = name;
      channel.unit = units(k);
      channel.sample_count = numel(data);
      if istimetable(T)
         channel.missing_count = checks.numnan.(name);
         channel.complex_count = checks.numcomplex.(name);
      else
         channel.missing_count = sum(isnan(data), 'all');
         channel.complex_count = sum(imag(data) ~= 0, 'all');
      end

      % A time-varying channel reports contiguous gaps independently for every
      % trailing component. Tables without one timestamp per row remain not
      % applicable rather than receiving misleading flattened run statistics.
      if ~isempty(times) && numel(times) == height(T)
         [channel.missing_run_count, ...
            channel.longest_missing_run_samples] = missingRuns(data);
      end
      finite = isfinite(data) & imag(data) == 0;
      channel.finite_count = nnz(finite);
      if any(finite, 'all')
         real_values = real(data(finite));
         channel.minimum = min(real_values, [], 'all');
         channel.maximum = max(real_values, [], 'all');
      end

      % Canonical metadata is checked through the single existing variable map.
      [canonical, expected_unit] = canonicalUnit(name, T, units);
      compatible_units = expected_unit;
      if is_event_bundle && name == "subsurface_temperature"
         compatible_units = ["K", "degC"];
      end
      if ~canonical
         severity = forcingSeverity(record.kind, "warning");
         findings(end + 1) = finding(severity, "noncanonical_channel", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, name, "numeric channel is absent from the canonical map"); %#ok<AGROW>
      elseif units(k) == ""
         findings(end + 1) = finding("error", "missing_unit", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, name, "expected unit '" + expected_unit ...
            + "', found no unit metadata"); %#ok<AGROW>
      elseif ~ismember(units(k), compatible_units)
         findings(end + 1) = finding("error", "wrong_unit", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, name, "expected unit '" + strjoin(compatible_units, ...
            "' or '") ...
            + "', found '" + units(k) + "'"); %#ok<AGROW>
      end

      % Complex values are never valid observations or forcing samples.
      if channel.complex_count > 0
         findings(end + 1) = finding("error", "complex_values", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, name, "channel contains complex-valued samples"); %#ok<AGROW>
      end

      % Preserve all-missing channels as explicit placeholders when documented;
      % do not count them as observed data.
      if channel.finite_count == 0
         documented = documentedPlaceholder(metadata, name);
         if documented || ~ismember(record.kind, ["met", "userdata"])
            findings(end + 1) = finding("placeholder", "intentional_placeholder", ...
               record.dataset_family, record.case_id, record.source, record.kind, ...
               record.path, name, "all-missing channel is not counted as observations"); %#ok<AGROW>
         else
            if requiredMetChannel(name, record.kind, T)
               severity = "error";
            else
               severity = "warning";
            end
            findings(end + 1) = finding(severity, "undocumented_all_missing", ...
               record.dataset_family, record.case_id, record.source, record.kind, ...
               record.path, name, ...
               "all-missing forcing channel lacks placeholder policy"); %#ok<AGROW>
         end
      elseif channel.missing_count > 0 ...
            && requiredMetChannel(name, record.kind, T)
         % A partially missing required channel is structurally valid but is
         % not direct-run forcing. Keep the audit nonblocking while making the
         % runtime repair requirement first-class rather than table-only.
         coverage = 100 * channel.finite_count / channel.sample_count;
         message = sprintf([ ...
            '%d of %d required met samples are missing (%.6g%% finite); ' ...
            'artifact is not directly runnable without an explicit repair policy'], ...
            channel.missing_count, channel.sample_count, coverage);
         findings(end + 1) = finding("warning", "required_met_gap", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, name, message); %#ok<AGROW>
      end

      % Apply conservative, source-independent physical bounds only.
      [bounded, minimum, maximum] = physicalRange(name, units(k));
      if bounded && any(real(data(finite)) < minimum ...
            | real(data(finite)) > maximum)
         findings(end + 1) = finding("error", "physical_range", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, name, sprintf('finite values must lie in [%g, %g]', ...
            minimum, maximum)); %#ok<AGROW>
      end

      % MAR quantiles and UTC-day reconstruction diagnostics do not impose a
      % magnitude threshold.
      if isMarSource(record.source) && record.kind == "userdata" ...
            && ismember(name, ["runoff", "smb"])
         [channel, mar_findings] = marChannelDiagnostics(channel, data, ...
            times, metadata, record);
         findings = appendStructs(findings, mar_findings);
      end
      channels(end + 1) = channel; %#ok<AGROW>
   end

   % The split precipitation channels must sum to total precipitation wherever
   % all three are actual observations.
   findings = appendStructs(findings, precipitationRelationship(T, record));

   % PROMICE observations are wrapped inside targets while userdata is a direct
   % timetable. Prefer each table's own provenance so the same thermistor QC
   % contract is enforced at both persistence boundaries without special-casing
   % the outer MAT payload shape.
   table_metadata = metadata;
   if isprop(T.Properties, 'UserData') ...
         && isstruct(T.Properties.UserData) ...
         && ~isempty(fieldnames(T.Properties.UserData))
      table_metadata = T.Properties.UserData;
   end
   if record.dataset_family == "promice"
      findings = appendStructs(findings, ...
         promiceThermistorChecks(T, table_metadata, record));
   end
end

function findings = validateEventTimes(T, times, record, table_path)
   %VALIDATEEVENTTIMES Validate dates without imposing timeseries row ordering.
   findings = repmat(emptyFinding(), 0, 1);
   names = string(T.Properties.VariableNames);
   if all(ismember(["start_date", "end_date"], names))
      start_missing = isnat(T.start_date);
      end_missing = isnat(T.end_date);
      if ismember("start_year", names)
         start_missing = start_missing & ~isfinite(double(T.start_year));
      end
      if ismember("end_year", names)
         end_missing = end_missing & ~isfinite(double(T.end_year));
      end
      invalid_interval = (~isnat(T.start_date) & ~isnat(T.end_date) ...
         & T.start_date > T.end_date);
      if ismember("start_year", names) && ismember("end_year", names)
         invalid_interval = invalid_interval ...
            | (isfinite(double(T.start_year)) & isfinite(double(T.end_year)) ...
            & double(T.start_year) > double(T.end_year));
      end
      invalid = any(start_missing | end_missing | invalid_interval);
   else
      invalid = any(isnat(times));
   end
   if invalid
      findings(end + 1) = finding("error", "missing_time", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, table_path, ...
         "event rows require a finite date or year and ordered intervals");
   end
end

function times = tableTimes(T)
   %TABLETIMES Return a table's canonical datetime axis when one exists.
   times = datetime.empty(0, 1);
   if istimetable(T)
      times = T.Properties.RowTimes;
      return
   end
   names = string(T.Properties.VariableNames);
   candidates = ["Time", "time", "datetime", "date"];
   hit = intersect(candidates, names, 'stable');
   if ~isempty(hit) && isdatetime(T.(hit(1)))
      times = T.(hit(1));
   elseif all(ismember(["start_date", "end_date"], names)) ...
         && isdatetime(T.start_date) && isdatetime(T.end_date)
      times = [T.start_date(:); T.end_date(:)];
   end
end

function [run_count, longest_run] = missingRuns(data)
   %MISSINGRUNS Count contiguous NaN runs along each time-series component.

   % Table variables may carry multiple trailing components. Reshape only those
   % trailing dimensions so the end of one component never joins the beginning
   % of another component into one artificial run.
   missing = reshape(isnan(data), size(data, 1), []);
   edges = diff([false(1, size(missing, 2)); missing; ...
      false(1, size(missing, 2))], 1, 1);
   starts = find(edges == 1);
   stops = find(edges == -1);
   run_count = numel(starts);

   % A complete channel has no missing runs and therefore a zero longest run.
   longest_run = 0;
   if ~isempty(starts)
      longest_run = max(stops - starts, [], 'all');
   end
end

function [findings, cadence] = validateTimeAxis(times, record, table_path)
   %VALIDATETIMEAXIS Check finite ordering and fixed forcing cadence.
   findings = repmat(emptyFinding(), 0, 1);
   cadence = NaN;
   if any(isnat(times))
      findings(end + 1) = finding("error", "missing_time", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, table_path, "time axis contains NaT");
      return
   end
   if numel(times) < 2
      return
   end
   steps = seconds(diff(times));
   if any(steps <= 0)
      findings(end + 1) = finding("error", "nonmonotonic_time", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, table_path, ...
         "time axis must be strictly increasing and duplicate-free");
      return
   end
   cadence = median(steps);
   if ismember(record.kind, ["met", "userdata"]) ...
         && any(abs(steps - steps(1)) > 1e-6)
      findings(end + 1) = finding("error", "irregular_cadence", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, table_path, "forcing time axis has nonuniform cadence");
   end
end

function findings = validateArtifactPeriod(record, times, period)
   %VALIDATEARTIFACTPERIOD Compare artifact bounds with its manifest window.
   findings = repmat(emptyFinding(), 0, 1);
   if ~isstruct(period) || ~all(isfield(period, ["start", "end"]))
      return
   end
   try
      [expected_start, expected_end] = ...
         icemodel.verification.setup.periodBounds(period);
   catch
      return
   end
   if isnat(expected_start) || isnat(expected_end)
      return
   end
   actual_start = min(icemodel.verification.setup.ensureUtc(times));
   actual_end = max(icemodel.verification.setup.ensureUtc(times));
   tolerance = 1;
   if isfinite(record.cadence_seconds)
      tolerance = max(tolerance, record.cadence_seconds);
   end
   if ismember(record.kind, ["met", "userdata"])
      % Cached forcing artifacts may intentionally enclose a narrower manifest
      % leg. They are valid when they cover the full leg; downstream loading
      % subsets them to the declared window.
      if seconds(actual_start - expected_start) > tolerance ...
            || seconds(expected_end - actual_end) > tolerance
         findings(end + 1) = finding("error", "period_coverage_drift", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, "Time", ...
            "forcing artifact does not cover its manifest window");
      end
   elseif seconds(expected_start - actual_start) > tolerance ...
         || seconds(actual_end - expected_end) > tolerance
      findings(end + 1) = finding("error", "period_outside_manifest", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "Time", "artifact time bounds extend outside manifest period");
   end
end

%% Source-specific policies
function findings = sourceChecks(T, metadata, record)
   %SOURCECHECKS Apply durable source and derived-met contracts.
   findings = repmat(emptyFinding(), 0, 1);
   if isMerraSource(record.source)
      findings = appendStructs(findings, merraChecks(T, metadata, record));
   end
   if isRacmoSource(record.source)
      findings = appendStructs(findings, racmoChecks(T, metadata, record));
   end
   if isMarSource(record.source)
      if record.kind == "userdata"
         % Native 24-posting RU/SMB ledgers and optional mass diagnostics belong
         % to Data/userdata. Derived met retains provenance but is validated by
         % its own generic and resampling contracts below.
         findings = appendStructs(findings, ...
            marMetadataChecks(T, metadata, record));
         findings = appendStructs(findings, ...
            marDiagnosticChecks(T, metadata, record));
      end
      findings = appendStructs(findings, ...
         marSnowDepthChecks(T, metadata, record));
   end
   findings = appendStructs(findings, modisChecks(T, metadata, record));
   if record.source == "promice" && record.kind == "met"
      findings = appendStructs(findings, ...
         promicePrecipitationChecks(T, metadata, record));
   end
   if record.source == "promice" && record.kind == "userdata"
      findings = appendStructs(findings, ...
         promiceShortwaveChecks(T, metadata, record));
   end
   findings = appendStructs(findings, ...
      albedoTransientChecks(T, metadata, record));
   findings = appendStructs(findings, metResampleChecks(T, metadata, record));
end

function findings = albedoTransientChecks(T, metadata, record)
   %ALBEDOTRANSIENTCHECKS Verify shared transient-albedo mask provenance.
   findings = repmat(emptyFinding(), 0, 1);
   names = string(T.Properties.VariableNames);
   source_family = lower(string(fieldOr(metadata, 'source_family', "")));
   applies = record.kind == "userdata" ...
      && all(ismember(["swd", "swu"], names)) ...
      && (any(ismember(source_family, ["gcnet_vandecrux", "imau"])) ...
      || isfield(metadata, 'albedo_transient_qc'));
   if ~applies
      return
   end

   % Recompute the compact contract from raw radiation so an artifact cannot
   % make its mask self-validating by stamping stale counts or dates.
   [flags, expected] = ...
      icemodel.forcing.helpers.dailyAlbedoAnomalyFlags( ...
      T.Time, double(T.swd), double(T.swu));
   expected = rmfield(expected, 'diagnostics');
   actual = fieldOr(metadata, 'albedo_transient_qc', struct());
   if ~validAlbedoTransientMetadata(actual, expected)
      findings(end + 1) = finding("error", ...
         "albedo_transient_qc_contract", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, ...
         "albedo", ...
         "transient-albedo counts, dates, or policy disagree with raw radiation");
   end

   % Raw swd/swu remain source-faithful; only derived radiation channels must
   % be missing on a confirmed episode.
   derived = intersect(["albedo", "swn", "netr"], names, 'stable');
   leaking = strings(0, 1);
   for name = reshape(derived, 1, [])
      values = T.(char(name));
      if any(isfinite(values(flags)))
         leaking(end + 1, 1) = name; %#ok<AGROW>
      end
   end
   if ~isempty(leaking)
      findings(end + 1) = finding("error", ...
         "albedo_transient_qc_mask", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, ...
         strjoin(leaking, ","), ...
         "confirmed transient-albedo rows retain finite derived radiation");
   end
end

function tf = validAlbedoTransientMetadata(actual, expected)
   %VALIDALBEDOTRANSIENTMETADATA Compare one compact artifact-safe contract.
   required = ["policy", "seed_day_count", "episode_day_count", ...
      "flagged_row_count", "seed_dates", "episode_dates"];
   tf = isstruct(actual) && isscalar(actual) ...
      && all(isfield(actual, required));
   if ~tf
      return
   end

   actual_policy = string(actual.policy);
   tf = isscalar(actual_policy) ...
      && actual_policy == string(expected.policy) ...
      && doubleScalar(actual.seed_day_count) == expected.seed_day_count ...
      && doubleScalar(actual.episode_day_count) == expected.episode_day_count ...
      && doubleScalar(actual.flagged_row_count) == expected.flagged_row_count ...
      && isequal(stringArray(actual.seed_dates), ...
      stringArray(expected.seed_dates)) ...
      && isequal(stringArray(actual.episode_dates), ...
      stringArray(expected.episode_dates));
end

function findings = racmoChecks(T, metadata, record)
   %RACMOCHECKS Require canonical positive-loss sublimation provenance.
   findings = repmat(emptyFinding(), 0, 1);
   if ~ismember("subl", string(T.Properties.VariableNames))
      return
   end

   native = string(fieldOr(metadata, ...
      'racmo_subl_native_sign_convention', ""));
   canonical = string(fieldOr(metadata, ...
      'racmo_subl_sign_convention', ""));
   if ~isscalar(native) || ~isscalar(canonical) || native == "" ...
         || canonical == ""
      findings(end + 1) = finding("error", ...
         "racmo_subl_sign_unverified", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "subl", ...
         "RACMO sublimation lacks scalar native and canonical sign markers");
   elseif native ~= "negative_loss_positive_deposition" ...
         || canonical ~= "positive_loss_negative_deposition"
      findings(end + 1) = finding("error", ...
         "racmo_subl_sign_not_canonical", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "subl", ...
         "RACMO sublimation sign markers do not match the canonical contract");
   end
end

function findings = marSnowDepthChecks(T, metadata, record)
   %MARSNOWDEPTHCHECKS Enforce SHSN2 semantics and the durable year mask.
   findings = repmat(emptyFinding(), 0, 1);
   if ~ismember("snowd", string(T.Properties.VariableNames))
      return
   end
   required = ["mar_snowd_qc_method", "mar_snowd_qc_status", ...
      "mar_snowd_source", "mar_snowd_semantics", ...
      "mar_snowd_shsn3_policy", "mar_snowd_jump_ratio", ...
      "mar_snowd_minimum_jump_m", ...
      "mar_snowd_discontinuous_boundary_years", ...
      "mar_snowd_boundary_jump_m", ...
      "mar_snowd_boundary_reference_p99_m", ...
      "mar_snowd_masked_years", "mar_snowd_unverified_years", ...
      "mar_snowd_masked_sample_count", "mar_snowd_qc_basis"];
   missing = required(~isfield(metadata, required));
   if ~isempty(missing)
      findings(end + 1) = finding("error", ...
         "mar_snowd_qc_metadata_missing", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "snowd", ...
         "MAR SHSN2 lacks required source-discontinuity provenance");
      return
   end

   boundary_years = double(metadata.mar_snowd_discontinuous_boundary_years(:));
   jumps = double(metadata.mar_snowd_boundary_jump_m(:));
   references = double(metadata.mar_snowd_boundary_reference_p99_m(:));
   masked_years = double(metadata.mar_snowd_masked_years(:));
   unverified_years = double(metadata.mar_snowd_unverified_years(:));
   status = string(metadata.mar_snowd_qc_status);
   valid = string(metadata.mar_snowd_qc_method) ...
      == "shsn2_archive_calibrated_year_boundary_screen" ...
      && string(metadata.mar_snowd_source) == "SHSN2" ...
      && string(metadata.mar_snowd_semantics) ...
      == "snow_pack_height_above_ice" ...
      && string(metadata.mar_snowd_shsn3_policy) ...
      == "not_used_total_multilayer_snow_firn_thickness" ...
      && ismember(status, ["applied", "insufficient_context"]) ...
      && doubleScalar(metadata.mar_snowd_jump_ratio) > 0 ...
      && doubleScalar(metadata.mar_snowd_minimum_jump_m) >= 0 ...
      && numel(boundary_years) == numel(jumps) ...
      && numel(boundary_years) == numel(references) ...
      && all(isfinite(boundary_years)) && all(fix(boundary_years) == boundary_years) ...
      && all(isfinite(jumps) & jumps >= 0) ...
      && all(isfinite(references) & references >= 0) ...
      && isempty(intersect(masked_years, unverified_years));
   if ~valid
      findings(end + 1) = finding("error", ...
         "mar_snowd_qc_metadata_invalid", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "snowd", ...
         "MAR SHSN2 definition, screen, or year classification is invalid");
      return
   end

   artifact_years = unique(year(T.Time));
   if any(~ismember([masked_years; unverified_years], artifact_years)) ...
         || any(isfinite(T.snowd(ismember(year(T.Time), masked_years))))
      findings(end + 1) = finding("error", ...
         "mar_snowd_mask_invalid", record.dataset_family, record.case_id, ...
         record.source, record.kind, record.path, "snowd", ...
         "MAR SHSN2 masked/unverified years disagree with the artifact");
   end
   if ~isempty(unverified_years)
      findings(end + 1) = finding("warning", ...
         "mar_snowd_unverified_years", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "snowd", ...
         "MAR SHSN2 retains one-sided or ambiguous discontinuity years");
   end
end

function findings = metResampleChecks(T, metadata, record)
   %METRESAMPLECHECKS Reject seasonal met artifacts that filled source gaps.
   findings = repmat(emptyFinding(), 0, 1);
   if record.kind ~= "met" || ~istimetable(T) ...
         || height(T) < 2
      return
   end

   cadence_s = median(seconds(diff(T.Time)), 'omitnan');
   if abs(cadence_s - 900) > 1e-6
      return
   end

   required = ["met_resample_policy", ...
      "met_resample_expected_missing_counts", ...
      "met_resample_time_semantics", ...
      "met_resample_support_end_exclusive"];
   if ~isstruct(metadata) || ~all(isfield(metadata, required))
      findings(end + 1) = finding("error", ...
         "met_resample_gap_metadata_missing", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "Time", ...
         "15-minute seasonal met lacks source-gap resampling provenance");
      return
   end

   allowed = ["interval_start_zero_order_hold", "native_15m_unchanged"];
   expected = metadata.met_resample_expected_missing_counts;
   if ~ismember(string(metadata.met_resample_policy), allowed) ...
         || string(metadata.met_resample_time_semantics) ~= "interval_start" ...
         || ~isstruct(expected)
      findings(end + 1) = finding("error", ...
         "met_resample_gap_metadata_invalid", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "Time", ...
         "15-minute seasonal met has invalid source-gap resampling provenance");
      return
   end

   % The output must include the complete final source interval without crossing
   % its exclusive support boundary. A linear-era artifact ends 45 minutes early.
   support_end = metadata.met_resample_support_end_exclusive;
   if ~isdatetime(support_end) || ~isscalar(support_end) ...
         || isnat(support_end) ...
         || T.Time(end) ~= support_end - minutes(15)
      findings(end + 1) = finding("error", ...
         "met_resample_support_end_invalid", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "Time", ...
         "15-minute met does not end one output step before exclusive source support");
   end

   % The writer records a source-derived lower bound independently of artifact
   % finiteness. Any channel with fewer missing values has bridged a native NaN
   % or omitted-time interval and is unsafe for verification.
   names = intersect(string(fieldnames(expected)), ...
      string(T.Properties.VariableNames), 'stable');
   for name = reshape(names, 1, [])
      values = T.(char(name));
      expected_count = doubleScalar(expected.(char(name)));
      if ~isnumeric(values) || ~isfinite(expected_count)
         continue
      end
      observed_count = nnz(~isfinite(values));
      if observed_count < expected_count
         message = sprintf(['derived met has %d missing values but source ' ...
            'support requires at least %d'], observed_count, expected_count);
         findings(end + 1) = finding("error", "met_gap_interpolation", ...
            record.dataset_family, record.case_id, record.source, ...
            record.kind, record.path, name, message); %#ok<AGROW>
      end
   end
end

function findings = merraChecks(T, metadata, record)
   %MERRACHECKS Require durable proof of flux and timestamp conventions.
   findings = repmat(emptyFinding(), 0, 1);
   channels = intersect(["shf", "lhf"], string(T.Properties.VariableNames));
   if ~isempty(channels)
      marker = string(fieldOr(metadata, 'merra_flux_sign_convention', ""));
      if numel(marker) ~= 1 || marker == ""
         findings(end + 1) = finding("blocker", "merra_sign_unverified", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, strjoin(channels, ","), ...
            "MERRA fluxes lack a scalar sign-orientation marker");
      elseif marker == "positive_upward"
         findings(end + 1) = finding("error", "merra_sign_not_canonical", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, strjoin(channels, ","), ...
            "MERRA shf/lhf remain positive upward instead of toward the surface");
      elseif marker ~= "positive_toward_surface"
         findings(end + 1) = finding("blocker", "merra_sign_unknown", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, strjoin(channels, ","), ...
            "MERRA flux-sign marker is not a recognized orientation");
      end
   end

   % Timing provenance applies to reduced glacier-channel artifacts even when no
   % turbulent flux is present; never couple its audit to the sign-check branch.
   required = ["merra_source_time_coordinate", ...
      "merra_time_relabel_policy", "merra_time_upsample_policy", ...
      "merra_collection_support_hours"];
   if ~isstruct(metadata) || ~all(isfield(metadata, required))
      findings(end + 1) = finding("error", "merra_time_semantics_missing", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "Time", ...
         "MERRA artifact lacks native/source-support timestamp provenance");
      return
   end

   if ~icemodel.forcing.helpers.hasCanonicalMerraTimeSupport(metadata)
      findings(end + 1) = finding("error", "merra_time_semantics_invalid", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "Time", ...
          "MERRA timestamp provenance is not the native-center to start/hold policy");
   end

   % A legacy value at a 00/03/... stamp is not necessarily native: an old
   % regularizer could have invented it across an omitted glc source interval.
   if ~icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(T, metadata)
      findings(end + 1) = finding("error", ...
         "merra_tavg3_source_grid_unproven", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "Time", ...
         "MERRA tavg3 channels lack an exact native glc timestamp inventory");
   end

   % Numeric block constancy is independent proof: stale canonical markers must
   % not allow an hourly/15-minute linear ramp to pass source-specific QA.
   if ~icemodel.forcing.helpers.hasConstantMerraTavg3Support(T)
      findings(end + 1) = finding("error", ...
         "merra_tavg3_support_invalid", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, ...
         "runoff,albedo,snowd,swe", ...
         "MERRA tavg3 channels are not constant over UTC three-hour support");
   end

   % Application artifacts begin on the canonical clock-hour grid. A :30 first
   % row is direct evidence that a native center escaped the builder boundary.
   if istimetable(T) && ~isempty(T.Time) ...
         && (minute(T.Time(1)) ~= 0 || second(T.Time(1)) ~= 0)
      findings(end + 1) = finding("error", ...
         "merra_time_axis_not_interval_start", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "Time", ...
         "MERRA application artifact does not begin on an interval-start hour");
   end
end

function findings = marMetadataChecks(T, metadata, record)
   %MARMETADATACHECKS Enforce the native-daily RU/SMB metadata contract.
   findings = repmat(emptyFinding(), 0, 1);
   names = string(T.Properties.VariableNames);
   present = intersect(["runoff", "smb"], names, 'stable');
   required = ["mar_qc_method", "mar_qc_status", "mar_qc_fallback", ...
      "mar_qc_channels", ...
      "mar_qc_sector", "mar_qc_sector_name", "mar_qc_runoff_source", ...
      "mar_qc_smb_source", "mar_qc_hourly_distribution", ...
      "mar_qc_partial_day_policy", "mar_qc_abs_tolerance_mwe_day", ...
      "mar_qc_rel_tolerance", "mar_qc_day_status_codes", ...
      "mar_qc_daily_reference_units", ...
      "mar_qc_runoff_day_status", "mar_qc_smb_day_status", ...
      "mar_qc_runoff_daily_reference_mwe", ...
      "mar_qc_smb_daily_reference_mwe", ...
      "mar_qc_preserved_runoff_day_count", ...
      "mar_qc_replaced_runoff_day_count", ...
      "mar_qc_unverified_runoff_day_count", ...
      "mar_qc_preserved_smb_day_count", ...
      "mar_qc_replaced_smb_day_count", ...
      "mar_qc_unverified_smb_day_count", ...
      "mar_qc_complete_utc_day_count", ...
      "mar_qc_partial_utc_day_count", "mar_qc_replaced_runoff_count", ...
      "mar_qc_replaced_smb_count", "mar_qc_basis"];
   missing = required(~isfield(metadata, required));
   if ~isempty(missing)
      findings(end + 1) = finding("error", "mar_qc_metadata_missing", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, strjoin(missing, ","), ...
         "MAR artifact lacks required native-daily QC metadata");
      return
   end

   % Source ids and distribution policy are frozen by the evidence-backed MAR
   % builder. Sector names must agree with the helper's two allowed sectors.
   sector = doubleScalar(metadata.mar_qc_sector);
   sector_name = "";
   if sector == 1
      sector_name = "permanent_ice";
   elseif sector == 2
      sector_name = "tundra";
   end
   applied_basis = ...
      "MAR hourly RUH/SMBH preserved where complete UTC-day sums agree " ...
      + "with native daily RU/SMB; missing, partial, or inconsistent days " ...
      + "use the native daily rate";
   fallback_basis = ...
      "MAR native daily RU/SMB unavailable; retained hourly RUH/SMBH";
   common = isequal(string(metadata.mar_qc_method), ...
      "daily_constrained_hourly") ...
      && ismember(sector, [1, 2]) ...
      && isequal(string(metadata.mar_qc_sector_name), sector_name) ...
      && isequal(string(metadata.mar_qc_runoff_source), "RU") ...
      && isequal(string(metadata.mar_qc_smb_source), "SMB") ...
      && isequal(string(metadata.mar_qc_hourly_distribution), ...
      "preserve_matching_hourly_else_daily_divided_by_24") ...
      && isequal(string(metadata.mar_qc_partial_day_policy), ...
      "native_daily_rate_applied_to_available_rows_marked_replaced") ...
      && isequal(string(metadata.mar_qc_day_status_codes), ...
      "1=preserved_hourly;2=replaced_daily;3=unverified") ...
      && isequal(string(metadata.mar_qc_daily_reference_units), "mWE/day");
   status = string(metadata.mar_qc_status);
   fallback = string(metadata.mar_qc_fallback);
   declared = stringArray(metadata.mar_qc_channels);
   is_applied = isequal(status, "applied");
   is_fallback = isequal(status, "not_applicable");
   if is_applied
      expected_channels = present(:);
      mode_consistent = ~isempty(present) && isequal(fallback, "none") ...
         && isequal(string(metadata.mar_qc_basis), applied_basis);
   elseif is_fallback
      expected_channels = strings(0, 1);
      mode_consistent = isequal(fallback, ...
         "hourly_RUH_SMBH_retained_native_daily_unavailable") ...
         && isequal(string(metadata.mar_qc_basis), fallback_basis);
   else
      expected_channels = strings(0, 1);
      mode_consistent = false;
   end
   consistent = common && mode_consistent;
   if ~consistent
      findings(end + 1) = finding("error", "mar_qc_metadata_invalid", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "mar_qc", ...
         "MAR native-daily method/status/fallback/source/sector/policy/basis disagree");
   end
   if (is_applied || is_fallback) ...
         && ~isequal(sort(declared), sort(expected_channels))
      findings(end + 1) = finding("error", "mar_qc_channels_mismatch", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "mar_qc_channels", ...
         "MAR QC channel list does not match the applied/fallback contract");
   end
   if is_fallback && consistent
      findings(end + 1) = finding("warning", "mar_native_daily_unavailable", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, strjoin(present, ","), ...
         "MAR retained hourly RUH/SMBH because native daily RU/SMB were unavailable");
   end
   [~, complete] = marUtcDayGroups(T.Time);
   if doubleScalar(metadata.mar_qc_complete_utc_day_count) ~= nnz(complete) ...
         || doubleScalar(metadata.mar_qc_partial_utc_day_count) ~= nnz(~complete)
      findings(end + 1) = finding("error", "mar_qc_count_mismatch", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "mar_qc_counts", ...
         "MAR QC complete/partial UTC-day counts disagree with the time axis");
   end
   % Replacement counts are cumulative idempotent provenance. Validate both
   % fields even when a reduced-source fallback has no declared replacements.
   for channel = ["runoff", "smb"]
      field = "mar_qc_replaced_" + channel + "_count";
      value = doubleScalar(metadata.(field));
      if ~isfinite(value) || value < 0 || fix(value) ~= value
         findings(end + 1) = finding("error", ...
            "mar_qc_replacement_count_invalid", record.dataset_family, ...
            record.case_id, record.source, record.kind, record.path, channel, ...
            "MAR replacement count must be a nonnegative integer"); %#ok<AGROW>
      end

      % The compact per-day ledger is the durable source-light proof. Its
      % length, status vocabulary, reference availability, support class, and
      % summary counts must agree exactly.
      day_status = uint8(metadata.("mar_qc_" + channel + "_day_status"));
      day_status = day_status(:);
      daily_reference = double(metadata.( ...
         "mar_qc_" + channel + "_daily_reference_mwe"));
      daily_reference = daily_reference(:);
      ndays = numel(unique(dateshift(T.Time, 'start', 'day')));
      ledger_valid = numel(day_status) == ndays ...
         && numel(daily_reference) == ndays ...
         && all(ismember(day_status, uint8([1 2 3]))) ...
         && all(isnan(daily_reference(day_status == 3))) ...
         && all(isfinite(daily_reference(day_status ~= 3))) ...
         && ~any(day_status == 1 & ~complete) ...
         && doubleScalar(metadata.("mar_qc_preserved_" + channel ...
         + "_day_count")) == nnz(day_status == 1) ...
         && doubleScalar(metadata.("mar_qc_replaced_" + channel ...
         + "_day_count")) == nnz(day_status == 2) ...
         && doubleScalar(metadata.("mar_qc_unverified_" + channel ...
         + "_day_count")) == nnz(day_status == 3);
      if ~ledger_valid
         findings(end + 1) = finding("error", ...
            "mar_qc_day_ledger_invalid", record.dataset_family, ...
            record.case_id, record.source, record.kind, record.path, channel, ...
            "MAR per-day status/reference ledger disagrees with its summary"); %#ok<AGROW>
      elseif is_applied && any(day_status == 3)
         findings(end + 1) = finding("warning", ...
            "mar_qc_unverified_days", record.dataset_family, ...
            record.case_id, record.source, record.kind, record.path, channel, ...
            sprintf('%d MAR UTC days lack a finite native daily reference', ...
            nnz(day_status == 3))); %#ok<AGROW>
      end
   end
end

function findings = marDiagnosticChecks(T, metadata, record)
   %MARDIAGNOSTICCHECKS Enforce optional SUH/SU/RZ/ME source semantics.
   findings = repmat(emptyFinding(), 0, 1);
   required = ["mar_diagnostic_method", "mar_diagnostic_status", ...
      "mar_diagnostic_channels", "mar_diagnostic_native_variables", ...
      "mar_diagnostic_units", "mar_diagnostic_subl_source", ...
      "mar_diagnostic_subl_evap_source", ...
      "mar_diagnostic_refreeze_deposition_source", ...
      "mar_diagnostic_melt_hourly_source", ...
      "mar_diagnostic_melt_daily_source", "mar_diagnostic_su_sector", ...
      "mar_diagnostic_su_sector_name", "mar_diagnostic_rz_sector", ...
      "mar_diagnostic_me_sector", "mar_diagnostic_daily_distribution", ...
      "mar_diagnostic_subl_sign", "mar_diagnostic_subl_evap_sign", ...
      "mar_diagnostic_refreeze_deposition_sign", ...
      "mar_diagnostic_suh_su_relationship", ...
       "mar_diagnostic_rz_relationship", ...
       "mar_diagnostic_abs_limit_mwe_h", ...
       "mar_diagnostic_refreeze_negative_day_count", ...
       "mar_diagnostic_refreeze_negative_minimum_mwe_h", ...
       "mar_diagnostic_refreeze_negative_statistics_basis", ...
       "mar_diagnostic_refreeze_material_negative_threshold_mwe_h", ...
       "mar_diagnostic_refreeze_material_negative_day_count", ...
       "mar_diagnostic_refreeze_material_negative_minimum_mwe_h", ...
       "mar_diagnostic_refreeze_material_negative_statistics_basis", ...
      "mar_diagnostic_melt_validation_status", ...
      "mar_diagnostic_melt_abs_tolerance_mwe_day", ...
      "mar_diagnostic_melt_rel_tolerance", ...
      "mar_diagnostic_melt_day_status_codes", ...
      "mar_diagnostic_melt_day_status", ...
      "mar_diagnostic_melt_daily_reference_mwe", ...
      "mar_diagnostic_melt_residual_mwe_day", ...
      "mar_diagnostic_melt_validated_day_count", ...
      "mar_diagnostic_melt_mismatch_day_count", ...
      "mar_diagnostic_melt_unverified_day_count", ...
      "mar_diagnostic_melt_max_abs_error_mwe_day", ...
      "mar_diagnostic_basis"];
   missing = required(~isfield(metadata, required));
   if ~isempty(missing)
      findings(end + 1) = finding("error", ...
         "mar_diagnostic_metadata_missing", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, ...
         strjoin(missing, ","), ...
         "MAR artifact lacks optional-diagnostic provenance");
      return
   end

   % Public names retain only scientifically defensible identities: SUH is
   % pure hourly sublimation, whereas SU and RZ remain explicitly combined.
   names = string(T.Properties.VariableNames);
   public_names = ["subl", "subl_evap", "refreeze_deposition"];
   source_names = ["SUH", "SU", "RZ"];
   present = intersect(public_names, names, 'stable')';
   expected_native = source_names(ismember(public_names, present))';
   melt_status = string(metadata.mar_diagnostic_melt_validation_status);
   if melt_status ~= "not_available"
      expected_native(end + 1, 1) = "ME";
   end
   if isempty(expected_native)
      expected_status = "not_available";
   elseif isequal(sort(expected_native), sort(["SUH"; "SU"; "RZ"; "ME"]))
      expected_status = "applied";
   else
      expected_status = "partial";
   end
   declared = stringArray(metadata.mar_diagnostic_channels);
   declared_native = stringArray(metadata.mar_diagnostic_native_variables);

   sector = doubleScalar(metadata.mar_diagnostic_su_sector);
   sector_name = "";
   if sector == 1
      sector_name = "permanent_ice";
   elseif sector == 2
      sector_name = "tundra";
   end
   basis = "MAR SUH is hourly sublimation; SU and RZ are native daily rates " ...
      + "divided by 24 and previous-held; daily ME validates but does not " ...
      + "replace hourly MEH";
   abs_limit = doubleScalar(metadata.mar_diagnostic_abs_limit_mwe_h);
   expected_refreeze = icemodel.forcing.helpers.marRefreezeMetadata(T);
   expected_refreeze_threshold = doubleScalar( ...
      expected_refreeze. ...
      mar_diagnostic_refreeze_material_negative_threshold_mwe_h);
   strict_statistics_match = doubleScalar( ...
      metadata.mar_diagnostic_refreeze_negative_day_count) ...
      == expected_refreeze.mar_diagnostic_refreeze_negative_day_count ...
      && isequaln(doubleScalar( ...
      metadata.mar_diagnostic_refreeze_negative_minimum_mwe_h), ...
      expected_refreeze.mar_diagnostic_refreeze_negative_minimum_mwe_h);
   material_statistics_match = doubleScalar( ...
      metadata.mar_diagnostic_refreeze_material_negative_day_count) ...
      == expected_refreeze.mar_diagnostic_refreeze_material_negative_day_count ...
      && isequaln(doubleScalar( ...
      metadata.mar_diagnostic_refreeze_material_negative_minimum_mwe_h), ...
      expected_refreeze.mar_diagnostic_refreeze_material_negative_minimum_mwe_h);
   refreeze_statistics_match = strict_statistics_match ...
      && material_statistics_match;
   common = string(metadata.mar_diagnostic_method) ...
      == "native_hourly_and_daily_previous_hold" ...
      && string(metadata.mar_diagnostic_status) == expected_status ...
      && string(metadata.mar_diagnostic_units) == "mWE/h" ...
      && string(metadata.mar_diagnostic_subl_source) == "SUH" ...
      && string(metadata.mar_diagnostic_subl_evap_source) == "SU" ...
      && string(metadata.mar_diagnostic_refreeze_deposition_source) == "RZ" ...
      && string(metadata.mar_diagnostic_melt_hourly_source) == "MEH" ...
      && string(metadata.mar_diagnostic_melt_daily_source) == "ME" ...
      && ismember(sector, [1 2]) ...
      && string(metadata.mar_diagnostic_su_sector_name) == sector_name ...
      && doubleScalar(metadata.mar_diagnostic_rz_sector) == 1 ...
      && doubleScalar(metadata.mar_diagnostic_me_sector) == 1 ...
      && string(metadata.mar_diagnostic_daily_distribution) ...
      == "native_daily_divided_by_24_previous_hold" ...
      && string(metadata.mar_diagnostic_subl_sign) ...
      == "positive_loss_negative_deposition" ...
      && string(metadata.mar_diagnostic_subl_evap_sign) ...
      == "positive_loss_negative_deposition" ...
      && string(metadata.mar_diagnostic_refreeze_deposition_sign) ...
      == "native_signed_combined_term_preserved_not_pure_refreeze" ...
      && string(metadata.mar_diagnostic_suh_su_relationship) ...
      == "distinct_native_products_not_interchangeable" ...
      && string(metadata.mar_diagnostic_rz_relationship) ...
      == "combined_refreezing_and_deposition_not_pure_refreeze" ...
      && isfinite(abs_limit) && abs_limit > 0 ...
      && string(metadata.mar_diagnostic_refreeze_negative_statistics_basis) ...
      == "distinct_utc_days_and_minimum_finite_negative_artifact_value" ...
      && doubleScalar( ...
      metadata.mar_diagnostic_refreeze_material_negative_threshold_mwe_h) ...
      == expected_refreeze_threshold ...
      && string( ...
      metadata.mar_diagnostic_refreeze_material_negative_statistics_basis) ...
      == "distinct_utc_days_and_minimum_below_declared_reporting_threshold" ...
      && string(metadata.mar_diagnostic_melt_day_status_codes) ...
      == "1=validated;2=mismatch;3=unverified" ...
      && string(metadata.mar_diagnostic_basis) == basis;
   if ~common
      findings(end + 1) = finding("error", ...
         "mar_diagnostic_metadata_invalid", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, ...
         "mar_diagnostic", ...
         "MAR diagnostic source, sector, support, sign, or basis is invalid");
   end
   if ~refreeze_statistics_match
      findings(end + 1) = finding("error", ...
         "mar_refreeze_metadata_mismatch", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, ...
         "refreeze_deposition", ...
         "MAR signed RZ strict/material negative statistics disagree with artifact");
   end
   if ~isequal(sort(declared), sort(present)) ...
         || ~isequal(sort(declared_native), sort(expected_native))
      findings(end + 1) = finding("error", ...
         "mar_diagnostic_channels_mismatch", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, ...
         "mar_diagnostic_channels", ...
         "MAR diagnostic channel/native-variable inventory is inconsistent");
   end

   % A broad symmetric source-informed bound catches unconverted mmWE values
   % without rejecting rare source-native signed RZ or ordinary signed SUH/SU.
   for channel = reshape(present, 1, [])
      values = real(T.(channel));
      finite = isfinite(values) & imag(T.(channel)) == 0;
      if any(abs(values(finite)) > abs_limit)
         findings(end + 1) = finding("error", ...
            "mar_diagnostic_range", record.dataset_family, record.case_id, ...
            record.source, record.kind, record.path, channel, ...
            sprintf('absolute MAR diagnostic values must not exceed %g mWE/h', ...
            abs_limit)); %#ok<AGROW>
      end
    end

   % Daily SU/RZ support must remain a constant /24 rate within each UTC day.
   [groups, ~] = marUtcDayGroups(T.Time);
   for channel = ["subl_evap", "refreeze_deposition"]
      if ~ismember(channel, present)
         continue
      end
      values = real(T.(channel));
      invalid = false;
      for day = 1:max(groups)
         sample = values(groups == day & isfinite(values));
         if numel(sample) > 1
            scale = max(1, max(abs(sample)));
            invalid = invalid || max(sample) - min(sample) > 16 * eps(scale);
         end
      end
      if invalid
         findings(end + 1) = finding("error", ...
            "mar_daily_diagnostic_support", record.dataset_family, ...
            record.case_id, record.source, record.kind, record.path, channel, ...
            "MAR daily diagnostic is not constant within its UTC-day support"); %#ok<AGROW>
      end
   end

   % Recompute the daily ME/MEH ledger from the staged artifact. This is the
   % only cross-product identity enforced: SUH/SU and RZ/pure refreeze remain
   % deliberately non-equivalent.
   day_status = uint8(metadata.mar_diagnostic_melt_day_status(:));
   reference = double(metadata.mar_diagnostic_melt_daily_reference_mwe(:));
   residual = double(metadata.mar_diagnostic_melt_residual_mwe_day(:));
   [groups, complete] = marUtcDayGroups(T.Time);
   ndays = numel(complete);
   expected_day_status = repmat(uint8(3), ndays, 1);
   expected_residual = nan(ndays, 1);
   abs_tolerance = doubleScalar( ...
      metadata.mar_diagnostic_melt_abs_tolerance_mwe_day);
   rel_tolerance = doubleScalar(metadata.mar_diagnostic_melt_rel_tolerance);
   ledger_valid = numel(day_status) == ndays ...
      && numel(reference) == ndays && numel(residual) == ndays ...
      && all(ismember(day_status, uint8([1 2 3]))) ...
      && isfinite(abs_tolerance) && abs_tolerance >= 0 ...
      && isfinite(rel_tolerance) && rel_tolerance >= 0;
   if ledger_valid && melt_status ~= "not_available"
      for day = 1:ndays
         rows = groups == day;
         if ~complete(day) || ~isfinite(reference(day)) ...
               || ~ismember("melt", names) || ~all(isfinite(T.melt(rows)))
            continue
         end
         expected_residual(day) = sum(T.melt(rows)) - reference(day);
         tolerance = abs_tolerance + rel_tolerance ...
            * max(abs([sum(T.melt(rows)), reference(day)]));
         if abs(expected_residual(day)) <= tolerance
            expected_day_status(day) = uint8(1);
         else
            expected_day_status(day) = uint8(2);
         end
      end
   end
   if melt_status == "not_available"
      expected_melt_status = "not_available";
      expected_day_status = zeros(0, 1, 'uint8');
      expected_residual = zeros(0, 1);
      ledger_valid = isempty(day_status) && isempty(reference) ...
         && isempty(residual);
   elseif any(expected_day_status == 2)
      expected_melt_status = "mismatch";
   elseif any(expected_day_status == 1)
      expected_melt_status = "validated";
   else
      expected_melt_status = "unverified";
   end
   finite_residual = expected_residual(isfinite(expected_residual));
   if isempty(finite_residual)
      expected_maximum = NaN;
   else
      expected_maximum = max(abs(finite_residual));
   end
   summary_valid = doubleScalar( ...
      metadata.mar_diagnostic_melt_validated_day_count) ...
      == nnz(expected_day_status == 1) ...
      && doubleScalar(metadata.mar_diagnostic_melt_mismatch_day_count) ...
      == nnz(expected_day_status == 2) ...
      && doubleScalar(metadata.mar_diagnostic_melt_unverified_day_count) ...
      == nnz(expected_day_status == 3) ...
      && isequaln(doubleScalar( ...
      metadata.mar_diagnostic_melt_max_abs_error_mwe_day), expected_maximum);
   if ~ledger_valid || melt_status ~= expected_melt_status ...
         || ~isequal(day_status, expected_day_status) ...
         || ~isequaln(residual, expected_residual) || ~summary_valid
      findings(end + 1) = finding("error", ...
         "mar_melt_daily_validation", record.dataset_family, record.case_id, ...
         record.source, record.kind, record.path, "melt", ...
         "MAR daily ME provenance does not reproduce hourly MEH validation");
   end
end

function [channel, findings] = marChannelDiagnostics(channel, data, times, ...
      metadata, record)
   %MARCHANNELDIAGNOSTICS Verify native-daily rates within each UTC day.
   findings = repmat(emptyFinding(), 0, 1);
   if isempty(times) || numel(times) ~= size(data, 1)
      return
   end
   finite = isfinite(data) & imag(data) == 0;
   finite_values = real(data(finite));
   if ~isempty(finite_values)
      channel.p01 = sampleQuantile(finite_values, 0.01);
      channel.p99 = sampleQuantile(finite_values, 0.99);
   end

   % Time-axis counts remain useful for both contracts. Applied hybrid QC is
   % checked against its saved native daily reference; explicit reduced-source
   % fallback remains unverified and retains hourly RUH/SMBH.
   [groups, complete] = marUtcDayGroups(times);
   channel.complete_utc_day_count = nnz(complete);
   channel.partial_utc_day_count = nnz(~complete);
   if ~isfield(metadata, 'mar_qc_status') ...
         || ~isequal(string(metadata.mar_qc_status), "applied")
      return
   end

   status_field = "mar_qc_" + channel.channel + "_day_status";
   reference_field = "mar_qc_" + channel.channel ...
      + "_daily_reference_mwe";
   if ~isfield(metadata, status_field) || ~isfield(metadata, reference_field) ...
         || ~isfield(metadata, 'mar_qc_abs_tolerance_mwe_day') ...
         || ~isfield(metadata, 'mar_qc_rel_tolerance')
      return
   end
   status = uint8(metadata.(status_field));
   reference = double(metadata.(reference_field));
   if numel(status) ~= numel(complete) || numel(reference) ~= numel(complete)
      return
   end
   abs_tolerance = doubleScalar(metadata.mar_qc_abs_tolerance_mwe_day);
   rel_tolerance = doubleScalar(metadata.mar_qc_rel_tolerance);
   invalid = 0;
   for day = 1:numel(complete)
      if status(day) == 3
         continue
      end
      values = data(groups == day);
      valid = isfinite(values) & imag(values) == 0;
      invalid_day = ~all(valid) || ~isfinite(reference(day));
      if ~invalid_day
         values = real(values);
         if status(day) == 2
            scale = max(1, max(abs(values)));
            invalid_day = max(values) - min(values) > 16 * eps(scale);
         end
         if complete(day)
            tolerance = abs_tolerance + rel_tolerance ...
               * max(abs([sum(values), reference(day)]));
            numeric_slack = 16 * eps(max(1, ...
               max(abs([sum(values), reference(day)]))));
            invalid_day = invalid_day ...
               || abs(sum(values) - reference(day)) ...
               > tolerance + numeric_slack;
         end
      end
      if invalid_day
         invalid = invalid + 1;
      end
   end
   channel.nonconstant_utc_day_count = invalid;
   if invalid > 0
      findings(end + 1) = finding("error", ...
         "mar_daily_constraint_inconsistent", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, channel.channel, ...
         "MAR constrained day disagrees with its saved daily reference/status");
   end

end

function [groups, complete] = marUtcDayGroups(times)
   %MARUTCDAYGROUPS Group rows by UTC day and identify exact complete days.
   times = icemodel.verification.setup.ensureUtc(times(:));
   days = dateshift(times, 'start', 'day');
   [groups, unique_days] = findgroups(days);
   complete = false(numel(unique_days), 1);
   for day = 1:numel(unique_days)
      actual = times(groups == day);
      expected = unique_days(day) + hours((0:23)');
      complete(day) = isequal(actual, expected);
   end
end

function findings = modisChecks(T, metadata, record)
   %MODISCHECKS Match MODIS coverage metadata to physical-data years.
   findings = repmat(emptyFinding(), 0, 1);
   names = string(T.Properties.VariableNames);
   has_channel = ismember("modis", names);
   metadata_fields = ["modis_product", "modis_status", "modis_coverage_years"];
   has_metadata = all(isfield(metadata, metadata_fields));
   if ~has_channel && ~has_metadata
      return
   elseif ~has_metadata
      findings(end + 1) = finding("error", "modis_metadata_missing", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "modis", "MODIS channel lacks coverage metadata");
      return
   end
   status = string(metadata.modis_status);
   coverage = sort(double(metadata.modis_coverage_years(:)))';
   valid_status = ismember(status, ["no_source_coverage", "source_coverage", ...
      "partial_source_coverage"]);
   if string(metadata.modis_product) ~= "GEUS Greenland Reflectivity 5km C6" ...
         || ~valid_status
      findings(end + 1) = finding("error", "modis_metadata_invalid", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "modis", "MODIS product/status metadata is not canonical");
      return
   end
   artifact_years = unique(year(T.Time))';
   if any(~ismember(coverage, artifact_years))
      findings(end + 1) = finding("error", "modis_coverage_outside_artifact", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "modis", "MODIS coverage metadata names years outside artifact");
   end
   valid = false(height(T), 1);
   if has_channel
      valid = isfinite(T.modis) & imag(T.modis) == 0 ...
         & real(T.modis) >= 0 & real(T.modis) <= 1;
   end
   if status == "no_source_coverage"
      if ~isempty(coverage) || (has_channel && any(valid))
         findings(end + 1) = finding("error", "modis_no_coverage_has_data", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, "modis", "finite MODIS data contradicts no-source coverage");
      end
      return
   end
   if ~has_channel
      findings(end + 1) = finding("error", "modis_coverage_without_channel", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "modis", "MODIS source coverage lacks a staged channel");
      return
   end
   finite_years = unique(year(T.Time(valid)))';
   if any(~ismember(coverage, finite_years)) || any(~ismember(finite_years, coverage))
      findings(end + 1) = finding("error", "modis_coverage_data_mismatch", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "modis", "MODIS finite-data years differ from coverage metadata");
   end
   if status == "source_coverage" && ~isequal(coverage, artifact_years)
      findings(end + 1) = finding("error", "modis_full_coverage_incomplete", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "modis", "full MODIS coverage does not include every artifact year");
   elseif status == "partial_source_coverage" ...
         && (isempty(setdiff(artifact_years, coverage)) || isempty(coverage))
      findings(end + 1) = finding("error", "modis_partial_coverage_invalid", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "modis", "partial MODIS coverage requires covered and missing years");
   end
end

function findings = promicePrecipitationChecks(T, metadata, record)
   %PROMICEPRECIPITATIONCHECKS Keep liquid rain source-faithful and splits explicit.
   findings = repmat(emptyFinding(), 0, 1);
   names = string(T.Properties.VariableNames);
   required_metadata = ["precip_policy", "rainf_source_present", ...
      "rainf_observations_present", "rainf_policy"];
   if ~all(isfield(metadata, required_metadata))
      findings(end + 1) = finding("error", "promice_precip_metadata_missing", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "precip", "PROMICE precipitation policy metadata is incomplete");
      return
   end
   policy = lower(string(metadata.precip_policy));
   if ~contains(policy, "rainf") || ~contains(policy, "snowf") ...
         || ~contains(policy, "ppt") || ~contains(policy, "placeholder")
      findings(end + 1) = finding("error", "promice_precip_policy_invalid", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "precip", "PROMICE policy must identify rainf and placeholder snowf/ppt");
   end
   for name = ["snowf", "ppt"]
      if ~ismember(name, names) || any(isfinite(T.(name)))
         findings(end + 1) = finding("error", "promice_invented_precip", ...
            record.dataset_family, record.case_id, record.source, record.kind, ...
            record.path, name, "PROMICE snowf/ppt must remain all-NaN placeholders"); %#ok<AGROW>
      end
   end
   source_present = logicalScalar(metadata.rainf_source_present);
   observations_present = logicalScalar(metadata.rainf_observations_present);
   finite_rain = ismember("rainf", names) && any(isfinite(T.rainf));
   if observations_present ~= finite_rain || (observations_present && ~source_present)
      findings(end + 1) = finding("error", "promice_rain_metadata_mismatch", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "rainf", "rainf finiteness disagrees with PROMICE source metadata");
   end
end

function findings = promiceShortwaveChecks(T, metadata, record)
   %PROMICESHORTWAVECHECKS Reject missing physical zeros in hourly userdata.
   findings = repmat(emptyFinding(), 0, 1);
   names = string(T.Properties.VariableNames);
   channels = ["swd", "swu"];
   channels = channels(ismember(channels, names));
   if isempty(channels) || height(T) < 2
      return
   end

   % The shared selector defines one-hour interval darkness. Apply it only to
   % hourly native userdata; 15-minute met has a separate resampling contract.
   steps = seconds(diff(T.Properties.RowTimes));
   if any(abs(steps - 3600) > 1e-6)
      return
   end
   latitude = doubleScalar(fieldOr(metadata, 'lat', ...
      fieldOr(metadata, 'lat_wgs84', NaN)));
   longitude = doubleScalar(fieldOr(metadata, 'lon', ...
      fieldOr(metadata, 'lon_wgs84', NaN)));
   if ~isfinite(latitude) || ~isfinite(longitude)
      findings(end + 1) = finding("error", ...
         "promice_shortwave_location_missing", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, ...
         strjoin(channels, ","), ...
         "hourly PROMICE shortwave audit requires finite station coordinates");
      return
   end

   % Reuse the production selector with staged public channels as raw inputs.
   % This identifies only missing bins whose complete source hour is deep-dark.
   Time = T.Properties.RowTimes;
   aws = timetable(Time);
   for channel = channels
      aws.(channel) = T.(channel);
   end
   swd_file_support = doubleScalar(fieldOr(metadata, ...
      'swd_source_file_observations_present', NaN));
   swu_file_support = doubleScalar(fieldOr(metadata, ...
      'swu_source_file_observations_present', NaN));
   [expected_swd, expected_swu] = ...
      icemodel.forcing.helpers.promiceShortwave(aws, fill_darkness=true, ...
      latitude=latitude, longitude=longitude, ...
      swd_source_file_observations_present=swd_file_support, ...
      swu_source_file_observations_present=swu_file_support);
   for channel = channels
      if channel == "swd"
         expected = expected_swd;
      else
         expected = expected_swu;
      end
      missing_dark = ~isfinite(T.(channel)) & expected == 0;
      if any(missing_dark)
         bands = icemodel.forcing.reconstruct.solarElevationBands();
         findings(end + 1) = finding("error", ...
            "promice_shortwave_darkness_gap", record.dataset_family, ...
            record.case_id, record.source, record.kind, record.path, channel, ...
            sprintf(['%d hourly samples are missing although the complete ' ...
            'interval is below the %g degree civil-night threshold'], ...
            nnz(missing_dark), bands.civil_twilight_deg)); %#ok<AGROW>
      end
   end
end

function findings = promiceThermistorChecks(T, metadata, record)
   %PROMICETHERMISTORCHECKS Enforce the canonical masked 10 m target contract.
   findings = repmat(emptyFinding(), 0, 1);
   names = string(T.Properties.VariableNames);
   contract = ["tice10m", "tice10m_source", "tice10m_qc_flag"];
   if ~any(ismember(contract, names))
      return
   end

   % A pre-QC artifact with only tice10m is unsafe: the source value, explicit
   % flag, and durable method/count provenance must travel together.
   required = ["tice10m_qc_status", "tice10m_qc_method", ...
      "tice10m_qc_source_variable", "tice10m_qc_source_channel", ...
      "tice10m_qc_jump_threshold_K", ...
       "tice10m_qc_persistent_jump_threshold_K", ...
       "tice10m_qc_other_sensor_median_threshold_K", ...
       "tice10m_qc_target_depth_tolerance_m", ...
       "tice10m_qc_recovery_window_hours", ...
      "tice10m_qc_depth_reset_threshold_m", "tice10m_qc_flag_codes", ...
      "tice10m_qc_flagged_sample_count", ...
      "tice10m_qc_failed_sample_count", ...
      "tice10m_qc_unreviewed_sample_count", ...
      "tice10m_qc_persistent_sample_count", "tice10m_qc_basis"];
   missing_channels = contract(~ismember(contract, names));
   missing_metadata = required(~isfield(metadata, required));
   if ~isempty(missing_channels) || ~isempty(missing_metadata)
      findings(end + 1) = finding("error", ...
         "promice_tice10m_qc_contract_missing", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "tice10m", ...
         "PROMICE 10 m target lacks source/flag channels or QC provenance");
      return
   end

   codes = metadata.tice10m_qc_flag_codes;
   threshold = doubleScalar(metadata.tice10m_qc_jump_threshold_K);
   valid_metadata = string(metadata.tice10m_qc_status) == "applied" ...
      && string(metadata.tice10m_qc_method) ...
      == "mask_gt_1K_hourly_endpoints_and_large_isolated_sensor_epochs" ...
      && string(metadata.tice10m_qc_source_variable) == "t_i_10m" ...
      && string(metadata.tice10m_qc_source_channel) == "tice10m_source" ...
      && isequal(threshold, 1) ...
      && doubleScalar(metadata.tice10m_qc_persistent_jump_threshold_K) == 4 ...
      && doubleScalar( ...
       metadata.tice10m_qc_other_sensor_median_threshold_K) == 0.25 ...
       && doubleScalar(metadata.tice10m_qc_target_depth_tolerance_m) == 2 ...
       && doubleScalar(metadata.tice10m_qc_recovery_window_hours) == 24 ...
      && doubleScalar(metadata.tice10m_qc_depth_reset_threshold_m) == 0.5 ...
      && isstruct(codes) ...
      && all(isfield(codes, ["accepted", "failed", "unreviewed", ...
      "persistent_unreviewed"])) ...
      && isequal([double(codes.accepted), double(codes.failed), ...
      double(codes.unreviewed), double(codes.persistent_unreviewed)], ...
      [0, 1, 2, 3]);
   if ~valid_metadata
      findings(end + 1) = finding("error", ...
         "promice_tice10m_qc_metadata_invalid", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "tice10m", ...
         "PROMICE 10 m QC method, threshold, source, or flag codes are invalid");
      return
   end

   source = T.tice10m_source;
   target = T.tice10m;
   flag = T.tice10m_qc_flag;
   valid_flag = isnumeric(flag) && isreal(flag) ...
      && all(isfinite(flag)) && all(ismember(flag, [0, 1, 2, 3]));
   if ~valid_flag || ~isequal(size(source), size(target), size(flag))
      findings(end + 1) = finding("error", ...
         "promice_tice10m_qc_flag_invalid", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, ...
         "tice10m_qc_flag", "PROMICE 10 m QC flag shape or values are invalid");
      return
   end

   % The canonical target must be source-identical where accepted and missing
   % wherever QC is nonzero. This catches silent edits and unmasked failures.
   accepted = flag == 0;
   if any(isfinite(target(~accepted))) ...
         || ~isequaln(target(accepted), source(accepted))
      findings(end + 1) = finding("error", ...
         "promice_tice10m_qc_mask_invalid", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "tice10m", ...
         "PROMICE canonical 10 m target disagrees with its source/flag mask");
   end

   % Recompute the cadence-independent source envelope so even a one-row
   % artifact cannot pass by carrying a self-consistent but false accepted flag.
   Tf = icemodel.physicalConstant('Tf');
   real_source = real(source);
   out_of_range = isfinite(source) & imag(source) == 0 ...
      & (real_source < Tf - 80 | real_source > Tf + 1);
   if any(out_of_range & accepted)
      findings(end + 1) = finding("error", ...
         "promice_tice10m_source_range_unmasked", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "tice10m", ...
         "PROMICE source outside -80..1 degC is marked accepted");
   end

   % Independently recompute the source-light temporal invariant. Every >1 K
   % exactly-hourly jump endpoint must be flagged even if metadata claims QC.
   if istimetable(T) && height(T) > 1
      contiguous = abs(seconds(diff(T.Time)) - 3600) <= 1;
      jumps = abs(diff(source));
      events = find(contiguous & isfinite(source(1:end-1)) ...
         & isfinite(source(2:end)) & jumps > threshold);
      if any(flag(unique([events; events + 1])) == 0)
         findings(end + 1) = finding("error", ...
            "promice_tice10m_discontinuity_unmasked", record.dataset_family, ...
            record.case_id, record.source, record.kind, record.path, "tice10m", ...
            "PROMICE source has an unmasked >1 K consecutive-hour jump");
      end
   end

   failed = nnz(flag == 1);
   unreviewed = nnz(flag >= 2);
   persistent_count = nnz(flag == 3);
   flagged = failed + unreviewed;
   counts_match = doubleScalar( ...
      metadata.tice10m_qc_flagged_sample_count) == flagged ...
      && doubleScalar(metadata.tice10m_qc_failed_sample_count) == failed ...
      && doubleScalar(metadata.tice10m_qc_unreviewed_sample_count) ...
      == unreviewed ...
      && doubleScalar(metadata.tice10m_qc_persistent_sample_count) ...
      == persistent_count;
   if ~counts_match
      findings(end + 1) = finding("error", ...
         "promice_tice10m_qc_count_mismatch", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, ...
         "tice10m_qc_flag", "PROMICE 10 m QC counts disagree with the flag");
   end
   if failed > 0
      findings(end + 1) = finding("warning", ...
         "promice_tice10m_failed_samples_masked", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "tice10m", ...
         sprintf(['%d source-range or reviewed thermistor-discontinuity ' ...
         'samples are masked'], ...
         failed));
   end
   if unreviewed > 0
      findings(end + 1) = finding("warning", ...
         "promice_tice10m_unreviewed_samples_masked", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "tice10m", ...
         sprintf('%d unreviewed discontinuity or epoch samples are masked', ...
         unreviewed));
   end
   if persistent_count > 0
      findings(end + 1) = finding("warning", ...
         "promice_tice10m_persistent_epoch_masked", record.dataset_family, ...
         record.case_id, record.source, record.kind, record.path, "tice10m", ...
         sprintf('%d unresolved isolated-sensor epoch samples are masked', ...
         persistent_count));
   end
end

function findings = precipitationRelationship(T, record)
   %PRECIPITATIONRELATIONSHIP Check ppt = rainf + snowf where all are observed.
   findings = repmat(emptyFinding(), 0, 1);
   names = string(T.Properties.VariableNames);
   if ~all(ismember(["ppt", "rainf", "snowf"], names))
      return
   end
   valid = isfinite(T.ppt) & isfinite(T.rainf) & isfinite(T.snowf);
   if ~any(valid)
      return
   end
   if any(~icemodel.forcing.helpers.precipitationConsistency( ...
         T.ppt(valid), T.rainf(valid), T.snowf(valid)))
      findings(end + 1) = finding("error", "precip_partition_mismatch", ...
         record.dataset_family, record.case_id, record.source, record.kind, ...
         record.path, "ppt", "ppt differs from rainf + snowf on finite samples");
   end
end

%% Metadata, range, and report helpers
function [canonical, unit] = canonicalUnit(name, T, units)
   %CANONICALUNIT Return the canonical or parent-channel uncertainty unit.
   canonical = true;
   unit = "";
   names = string(T.Properties.VariableNames);
   if name == "error"
      parents = ["density", "subsurface_temperature", "smb"];
      hit = find(ismember(names, parents), 1);
      if ~isempty(hit)
         unit = units(hit);
         return
      end
   end
   if name == "smb" && all(ismember(["start_date", "end_date"], names))
      unit = "m w.e.";
      return
   end
   try
      info = icemodel.netcdf.defaults.variable(name);
      unit = string(info.unit);
   catch err
      if strcmp(err.identifier, 'icemodel:netcdf:variable:unknownChannel')
         canonical = false;
      else
         rethrow(err)
      end
   end
end

function units = tableUnits(T, n)
   %TABLEUNITS Normalize table unit metadata to one string per variable.
   units = strings(1, n);
   raw = string(T.Properties.VariableUnits);
   units(1:min(n, numel(raw))) = raw(1:min(n, numel(raw)));
end

function [bounded, minimum, maximum] = physicalRange(name, unit)
   %PHYSICALRANGE Conservative bounds shared across verification sources.
   bounded = true;
   minimum = -Inf;
   maximum = Inf;
   token = lower(string(name));
   switch token
      case {"albedo", "modis", "cfrac"}
         minimum = 0;
         maximum = 1;
      case "rh"
         minimum = 0;
         maximum = 100;
      case "wspd"
         minimum = 0;
      case "wdir"
         minimum = 0;
         maximum = 360;
      case "psfc"
         minimum = 20000;
         maximum = 120000;
      case {"tair"}
         minimum = 180;
         maximum = 330;
      case "tsfc"
         if unit == "K"
            minimum = 150;
            maximum = 273.16;
         else
            minimum = -123.16;
            maximum = 0;
         end
      case {"surface_temp_c", "subsurface_temperature"}
         % ESM-SnowMIP includes land/forest surface observations that can exceed
         % the snow melting point; source profile temperatures may use degC.
         if unit == "K"
            minimum = 150;
            maximum = 353.15;
         else
            minimum = -123.16;
            maximum = 80;
         end
      case {"swd", "lwd", "ppt", "rainf", "snowf", "snow_depth", ...
            "snowd", "swe", "swe_kg_m2", "density", "lwc", ...
            "bottom_outflow_mps"}
         minimum = 0;
      otherwise
         bounded = false;
   end
end

function tf = documentedPlaceholder(metadata, channel)
   %DOCUMENTEDPLACEHOLDER Find an explicit channel-specific placeholder policy.
   tf = false;
   if isempty(metadata) || ~isstruct(metadata)
      return
   end
   fields = string(fieldnames(metadata));
   policy_fields = fields(contains(lower(fields), "policy") ...
      | contains(lower(fields), "reason"));
   for field = reshape(policy_fields, 1, [])
      value = metadata.(char(field));
      if ischar(value) || isstring(value) || iscellstr(value)
         text = lower(strjoin(stringArray(value), " "));
         if contains(text, lower(channel)) && contains(text, "placeholder")
            tf = true;
            return
         end
      end
   end
end

function tf = requiredMetChannel(name, kind, T)
   %REQUIREDMETCHANNEL True for a required channel in a met artifact.
   required = icemodel.forcing.helpers.metvariables();
   tf = kind == "met" && ismember(name, required) ...
      && ismember(name, string(T.Properties.VariableNames));
end

function severity = forcingSeverity(kind, fallback)
   %FORCINGSEVERITY Escalate forcing schema defects while preserving obs metadata.
   if ismember(kind, ["met", "userdata"])
      severity = "error";
   else
      severity = fallback;
   end
end

function status = artifactStatus(findings)
   %ARTIFACTSTATUS Reduce finding severities to one compact artifact status.
   severities = string({findings.severity});
   if any(severities == "error")
      status = "error";
   elseif any(severities == "blocker")
      status = "blocked";
   elseif any(severities == "warning")
      status = "warning";
   elseif any(severities == "placeholder")
      status = "placeholder";
   else
      status = "passed";
   end
end

function tf = metadataEqual(a, b)
   %METADATAEQUAL Compare scalar metadata structs independent of field order.
   if ~isstruct(a) || ~isstruct(b) || ~isscalar(a) || ~isscalar(b)
      tf = isequaln(a, b);
      return
   end
   if ~isequal(sort(string(fieldnames(a))), sort(string(fieldnames(b))))
      tf = false;
      return
   end
   tf = isequaln(orderfields(a), orderfields(b));
end

function [found, leg] = sourceLeg(colocation, source)
   %SOURCELEG Match a source key, public source id, or builder source name.
   found = false;
   leg = struct();
   if ~isstruct(colocation)
      return
   end
   source = string(source);
   names = string(fieldnames(colocation));
   for name = reshape(names, 1, [])
      candidate = colocation.(char(name));
      if ~isstruct(candidate)
         continue
      end
      labels = [name, string(fieldOr(candidate, 'source', "")), ...
         string(fieldOr(candidate, 'source_id', ""))];
      if any(labels == source)
         found = true;
         leg = candidate;
         return
      end
   end
end

function label = sourceLabel(name, leg)
   %SOURCELABEL Prefer a public product id over an internal colocation key.
   label = string(fieldOr(leg, 'source_id', ""));
   if label == ""
      label = string(name);
   end
end

function tf = hasFiles(s, field)
   %HASFILES True when a manifest file field has a nonblank value.
   tf = isfield(s, field) && any(strlength(stringArray(s.(char(field)))) > 0);
end

function values = stringArray(value)
   %STRINGARRAY Normalize decoded scalar/cell text to a nonblank string column.
   if isempty(value)
      values = strings(0, 1);
      return
   end
   if ischar(value)
      values = string(value);
   else
      values = string(value(:));
   end
   values = values(strlength(values) > 0 & ~ismissing(values));
end

function value = fieldOr(s, field, fallback)
   %FIELDOR Return a struct field or a caller-supplied fallback.
   if isstruct(s) && isfield(s, field)
      value = s.(char(field));
   else
      value = fallback;
   end
end

function combined = appendStructs(existing, added)
   %APPENDSTRUCTS Concatenate homogeneous records independent of row orientation.
   combined = [existing(:); added(:)];
end

function value = logicalScalar(value)
   %LOGICALSCALAR Coerce scalar logical-like metadata safely.
   value = logical(value);
   value = isscalar(value) && value;
end

function value = doubleScalar(value)
   %DOUBLESCALAR Return a scalar double or NaN for malformed metadata.
   value = double(value);
   if ~isscalar(value)
      value = NaN;
   end
end

function tf = isMerraSource(source)
   %ISMERRASOURCE Recognize internal and public MERRA source ids.
   tf = ismember(lower(string(source)), ["merra", "merra2"]);
end

function tf = isMarSource(source)
   %ISMARSOURCE Recognize internal and public MAR source ids.
   tf = ismember(lower(string(source)), ["mar", "mar3.11"]);
end

function tf = isRacmoSource(source)
   %ISRACMOSOURCE Recognize internal and public RACMO source ids.
   tf = ismember(lower(string(source)), ["racmo", "racmo2.3p3"]);
end

function value = sampleQuantile(values, probability)
   %SAMPLEQUANTILE Linear-interpolated sample quantile without extra toolboxes.
   values = sort(double(values(:)));
   if isempty(values)
      value = NaN;
      return
   end
   position = 1 + (numel(values) - 1) * probability;
   lower_index = floor(position);
   upper_index = ceil(position);
   fraction = position - lower_index;
   value = values(lower_index) * (1 - fraction) ...
      + values(upper_index) * fraction;
end

function value = timeString(t)
   %TIMESTRING Render a datetime in the portable manifest timestamp format.
   t = icemodel.verification.setup.ensureUtc(t);
   value = string(t, 'yyyy-MM-dd HH:mm:ss');
end

%% Stable records and report writing
function record = emptyArtifact()
   %EMPTYARTIFACT Prototype for one audited manifest/MAT artifact.
   record = struct( ...
      'dataset_family', "", 'case_id', "", 'source', "", 'kind', "", ...
      'path', "", 'exists', false, 'payload', "", 'status', "", ...
      'artifact_sha256', "", 'artifact_size_bytes', NaN, ...
      'expected_start', "", 'expected_end', "", 'time_start', "", ...
      'time_end', "", 'cadence_seconds', NaN, 'n_tables', 0, ...
      'n_samples', 0, 'n_channels', 0);
end

function record = emptyChannel()
   %EMPTYCHANNEL Prototype for one numeric channel diagnostic record.
   record = struct( ...
      'dataset_family', "", 'case_id', "", 'source', "", 'kind', "", ...
      'path', "", 'table_path', "", 'channel', "", 'unit', "", ...
      'sample_count', 0, 'finite_count', 0, 'missing_count', 0, ...
      'missing_run_count', NaN, 'longest_missing_run_samples', NaN, ...
      'complex_count', 0, 'minimum', NaN, 'maximum', NaN, 'p01', NaN, ...
      'p99', NaN, 'complete_utc_day_count', 0, ...
      'partial_utc_day_count', 0, 'nonconstant_utc_day_count', 0);
end

function record = emptyFinding()
   %EMPTYFINDING Prototype for one stable QA/QC finding.
   record = struct( ...
      'severity', "", 'code', "", 'dataset_family', "", 'case_id', "", ...
      'source', "", 'kind', "", 'path', "", 'channel', "", 'message', "");
end

function record = finding(severity, code, family, case_id, source, kind, ...
      path, channel, message)
   %FINDING Construct one normalized QA/QC finding record.
   record = emptyFinding();
   record.severity = string(severity);
   record.code = string(code);
   record.dataset_family = string(family);
   record.case_id = string(case_id);
   record.source = string(source);
   record.kind = string(kind);
   record.path = string(path);
   record.channel = string(channel);
   record.message = string(message);
end

function records = sortRecords(records, fields)
   %SORTRECORDS Sort a homogeneous struct array by concatenated string fields.
   if isempty(records)
      return
   end
   keys = strings(numel(records), 1);
   for k = 1:numel(records)
      parts = strings(1, numel(fields));
      for n = 1:numel(fields)
         parts(n) = string(records(k).(char(fields(n))));
      end
      keys(k) = strjoin(parts, char(31));
   end
   [~, order] = sort(keys);
   records = records(order);
end

function findings = sortFindings(findings)
   %SORTFINDINGS Put blocking severities first, then use stable content keys.
   if isempty(findings)
      return
   end
   severity_order = ["error", "blocker", "warning", "placeholder"];
   keys = strings(numel(findings), 1);
   for k = 1:numel(findings)
      rank = find(severity_order == findings(k).severity, 1);
      if isempty(rank)
         rank = numel(severity_order) + 1;
      end
      content = strjoin([findings(k).dataset_family, findings(k).case_id, ...
         findings(k).source, findings(k).kind, findings(k).path, ...
         findings(k).channel, findings(k).code], char(31));
      keys(k) = string(sprintf('%02d%c%s', rank, char(31), char(content)));
   end
   [~, order] = sort(keys);
   findings = findings(order);
end

function summary = summarizeReport(families, case_count, artifacts, ...
      channels, findings)
   %SUMMARIZEREPORT Count report records and pass/fail severities.
   severities = string({findings.severity});
   summary = struct();
   summary.family_count = numel(families);
   summary.case_count = case_count;
   summary.artifact_count = numel(artifacts);
   summary.channel_count = numel(channels);
   summary.finding_count = numel(findings);
   summary.error_count = nnz(severities == "error");
   summary.blocker_count = nnz(severities == "blocker");
   summary.warning_count = nnz(severities == "warning");
   summary.placeholder_count = nnz(severities == "placeholder");
   summary.passed = summary.error_count == 0 && summary.blocker_count == 0;
end

function report = writeReports(report, report_dir)
   %WRITEREPORTS Emit machine-readable JSON and a compact Markdown summary.
   if ~isfolder(report_dir)
      mkdir(report_dir)
   end
   json_path = fullfile(report_dir, "artifact_qa.json");
   markdown_path = fullfile(report_dir, "artifact_qa.md");
   report.report_files = struct('json', json_path, 'markdown', markdown_path);

   % JSON carries every stable record; Markdown stays compact for humans and
   % downstream Quarto inclusion.
   writeText(json_path, jsonencode(report, PrettyPrint=true));
   writeText(markdown_path, markdownReport(report));
end

function text = markdownReport(report)
   %MARKDOWNREPORT Render the compact human-facing audit summary.
   s = report.summary;
   lines = [ ...
      "# Verification artifact QA/QC"
      ""
      "- Generated: " + report.generated_at
      "- Evaluation root: `" + report.evaluation_data_root + "`"
      "- Input root: `" + report.input_data_root + "`"
      "- Families: " + strjoin(report.families, ", ")
      "- Passed: **" + string(s.passed) + "**"
      ""
      "| artifacts | channels | errors | blockers | warnings | placeholders |"
      "|---:|---:|---:|---:|---:|---:|"
      sprintf('| %d | %d | %d | %d | %d | %d |', s.artifact_count, ...
      s.channel_count, s.error_count, s.blocker_count, s.warning_count, ...
      s.placeholder_count)
      ""
      "## Findings"
      ""];
   visible = report.findings(string({report.findings.severity}) ~= "placeholder");
   if isempty(visible)
      lines(end + 1) = "No errors, blockers, or warnings.";
   else
      lines = [lines; ...
         "| severity | code | family | source | channel | count | message |"; ...
         "|---|---|---|---|---|---:|---|"];

      % Group repeated per-file findings so the human report remains compact;
      % the JSON retains every exact case/path record for repair routing.
      keys = strings(numel(visible), 1);
      for k = 1:numel(visible)
         keys(k) = strjoin([visible(k).severity, visible(k).code, ...
            visible(k).dataset_family, visible(k).source, ...
            visible(k).channel], char(31));
      end
      unique_keys = unique(keys, 'stable');
      for key = reshape(unique_keys, 1, [])
         keep = keys == key;
         f = visible(find(keep, 1));
         row = [f.severity, f.code, f.dataset_family, f.source, ...
            f.channel, string(nnz(keep)), f.message];
         row = replace(row, "|", "\\|");
         row = replace(row, newline, " ");
         lines(end + 1) = "| " + strjoin(row, " | ") + " |"; %#ok<AGROW>
      end
   end
   if s.placeholder_count > 0
      note = sprintf(['%d intentional placeholder findings are retained in ' ...
         'the JSON record and excluded from observations.'], s.placeholder_count);
      lines = [lines; ""; string(note)];
   end
   text = strjoin(lines, newline) + newline;
end

function writeText(pathname, text)
   %WRITETEXT Write one complete UTF-8 report and close it on every path.
   fid = fopen(pathname, 'w', 'n', 'UTF-8');
   if fid < 0
      error('icemodel:verification:auditArtifacts:reportOpen', ...
         'could not open report for writing: %s', pathname)
   end
   cleanup = onCleanup(@() fclose(fid));
   fwrite(fid, char(text), 'char');
   clear cleanup
end
