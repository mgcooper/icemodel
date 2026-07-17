function report = repairRcmArtifactMetadata(input_root, kwargs)
   %REPAIRRCMARTIFACTMETADATA Repair reusable staged RCM MAT artifacts.
   %
   %  report = icemodel.verification.setup.repairRcmArtifactMetadata()
   %  report = ... repairRcmArtifactMetadata(input_root, eval_root=..., ...
   %     dataset_family="promice", source_id="racmo2.3p3", ...
   %     repair_racmo_subl=true, mar_dir=..., modis_dir=..., dry_run=false)
   %
   % Role
   %  Repair manifest-referenced staged RCM MAT artifacts without rebuilding
   %  their expensive source payloads. The legacy function name now covers both
   %  current metadata stamping and this explicit source-aware transform set:
   %  RACMO precip -> nonnegative canonical ppt plus canonical sublimation sign;
   %  MERRA turbulent-flux sign and tavg3 support normalization; MAR native
   %  daily RU/SMB diagnostics and signed RZ provenance; and optional GEUS MODIS
   %  augmentation. Current manifests supply requested-site coordinates, while
   %  existing payload CustomProperties supply sampled source-grid coordinates.
   %  Every supported transform is dry-runnable and idempotent, and unrelated
   %  numeric payloads plus time axes are preserved exactly.
   %
   % Boundaries
   %  dry_run defaults true. Unsupported, ambiguous, or unprovable migrations
   %  are reported as restage-required; they are never guessed. This helper does
   %  not change manifests or observations.mat and is not a generic staged-
   %  artifact coordinator. Use refreshManifestSourceLists for manifest-only
   %  source-list policy changes, repairMetTimeSupport for its narrow met-grid
   %  migration, and a canonical importer/restage for raw payload, observation,
   %  grid-cell selection, or point-colocation changes. Unmapped aliases are
   %  reported for explicit classification.
   %
   % Output
   %  report.records classifies each artifact; report.summary and report.actions
   %  aggregate status/actions. Source-read counters make the required second
   %  unchanged dry run auditable.
   %
   % See also: icemodel.forcing.helpers.writeuserdata,
   %  icemodel.forcing.helpers.writemet,
   %  icemodel.verification.setup.refreshManifestSourceLists,
   %  icemodel.verification.setup.repairMetTimeSupport

   arguments
      input_root (1, 1) string = string(fullfile( ...
         icemodel.internal.fullpath("data"), "input"))
      kwargs.eval_root (1, 1) string = string(fullfile( ...
         icemodel.internal.fullpath("data"), "eval"))
      kwargs.dataset_family (1, :) string = strings(1, 0)
      kwargs.source_id (1, :) string = strings(1, 0)
      kwargs.repair_racmo_subl (1, 1) logical = true
      kwargs.modis_dir (1, 1) string = ""
      kwargs.mar_dir (1, 1) string = ""
      kwargs.merra_dir (1, 1) string = ""
      kwargs.verify_merra_sign (1, 1) logical = false
      kwargs.dry_run (1, 1) logical = true
   end

   % Build current artifact/alias -> requested-location lookups from manifests.
   % Exact file references win; alias fallback is allowed only when all current
   % manifests using that alias agree on coordinates.
   families = unique(reshape(kwargs.dataset_family, 1, []), "stable");
   locations = manifestLocations(kwargs.eval_root, input_root, families);
   files = rcmArtifactFiles(input_root, locations, kwargs.source_id);
   state = repairState(kwargs.modis_dir, kwargs.mar_dir, kwargs.merra_dir, ...
      kwargs.verify_merra_sign);
   state.modis_plan = modisPointPlan(files, locations, state.modis_inventory);
   state.mar_plan = marPointPlan(files, locations, state.mar_dir);
   records = repmat(emptyRecord(), numel(files), 1);

   % Process each artifact independently so one malformed file does not hide the
   % full classification of the cache.
   for k = 1:numel(files)
      records(k) = repairOne(files(k), locations, state, kwargs.dry_run, ...
         kwargs.repair_racmo_subl);
   end

   report = struct();
   report.input_root = input_root;
   report.eval_root = kwargs.eval_root;
   report.dataset_family = families;
   report.source_id = kwargs.source_id;
   report.repair_racmo_subl = kwargs.repair_racmo_subl;
   report.modis_dir = kwargs.modis_dir;
   report.mar_dir = state.mar_dir;
   report.dry_run = kwargs.dry_run;
   report.records = records;
   report.summary = summarizeRecords(records);
   report.actions = summarizeActions(records);
   report.modis_source_reads = double(state.modis_reads.Count);
   report.modis_cached_series = double(state.modis_series.Count);
   report.mar_source_reads = double(state.mar_reads.Count);
   report.mar_cached_series = double(state.mar_series.Count);
   report.merra_source_reads = double( ...
      state.merra_flux.Count + state.merra_time.Count);
end

function state = repairState(modis_dir, mar_dir, merra_dir, verify_merra_sign)
   %REPAIRSTATE Create shared source inventories and per-run read caches.
   state = struct();
   state.modis_dir = modis_dir;
   state.modis_inventory = modisInventory(modis_dir);
   state.modis_series = containers.Map('KeyType', 'char', 'ValueType', 'any');
   state.modis_reads = containers.Map('KeyType', 'char', 'ValueType', 'logical');
   state.mar_dir = resolveMarDir(mar_dir);
   state.mar_series = containers.Map('KeyType', 'char', 'ValueType', 'any');
   state.mar_reads = containers.Map('KeyType', 'char', 'ValueType', 'logical');
   state.merra_dir = resolveMerraDir(merra_dir);
   state.merra_flux = containers.Map('KeyType', 'char', 'ValueType', 'any');
   state.merra_time = containers.Map('KeyType', 'char', 'ValueType', 'any');
   state.merra_grid = containers.Map('KeyType', 'char', 'ValueType', 'any');
   state.merra_glc_files = merraGlcInventory(state.merra_dir);
   state.verify_merra_sign = verify_merra_sign;
end

function locations = manifestLocations(eval_root, input_root, families)
   %MANIFESTLOCATIONS Build current artifact and alias location lookups.
   manifests = dir(fullfile(eval_root, "*", "manifest.json"));
   locations = struct();
   locations.by_file = containers.Map('KeyType', 'char', 'ValueType', 'any');
   locations.by_alias = containers.Map('KeyType', 'char', 'ValueType', 'any');
   locations.ambiguous_file = containers.Map('KeyType', 'char', ...
      'ValueType', 'logical');
   locations.ambiguous_alias = containers.Map('KeyType', 'char', ...
      'ValueType', 'logical');
   locations.restrict_to_files = ~isempty(families);
   for k = 1:numel(manifests)
      [~, family] = fileparts(manifests(k).folder);
      if ~isempty(families) && ~ismember(string(family), families)
         continue
      end
      manifest_file = fullfile(manifests(k).folder, manifests(k).name);
      manifest = jsondecode(fileread(manifest_file));
      cases = manifest.cases;
      if isempty(cases)
         continue
      end
      for n = 1:numel(cases)
         case_id = string(cases(n).case_id);
         if ~isfield(cases(n), "site_location")
            continue
         end
         loc = cases(n).site_location;
         if ~all(isfield(loc, ["lat_wgs84", "lon_wgs84"]))
            continue
         end
         locations = addAliasLocation(locations, case_id, loc);
         locations = addArtifactLocations(locations, input_root, cases(n), loc);
      end
   end
end

function locations = addAliasLocation(locations, alias, loc)
   %ADDALIASLOCATION Add one alias location or mark conflicting aliases.
   key = char(matlab.lang.makeValidName(alias));
   entry = manifestLocationEntry(loc, "nearest");
   if ~isKey(locations.by_alias, key)
      locations.by_alias(key) = entry;
      locations.ambiguous_alias(key) = false;
      return
   end
   prior = locations.by_alias(key);
   if ~sameManifestLocationEntry(prior, entry)
      locations.ambiguous_alias(key) = true;
   end
end

function locations = addArtifactLocations(locations, input_root, c, loc)
   %ADDARTIFACTLOCATIONS Map exact manifest artifact refs to case locations.
   if ~isfield(c, "colocation") || ~isstruct(c.colocation)
      return
   end
   names = string(fieldnames(c.colocation));
   for source = reshape(names, 1, [])
      leg = c.colocation.(char(source));
      if ~isstruct(leg)
         continue
      end
      locations = addArtifactField(locations, input_root, leg, ...
         "met_files", "met", loc);
      locations = addArtifactField(locations, input_root, leg, ...
         "data_files", "userdata", loc);
   end
end

function locations = addArtifactField(locations, input_root, leg, field, kind, loc)
   %ADDARTIFACTFIELD Add exact file refs from one manifest field.
   if ~isfield(leg, field) || isempty(leg.(field))
      return
   end
   method = "nearest";
   if isfield(leg, "sample_method") && strlength(string(leg.sample_method)) > 0
      method = string(leg.sample_method);
   end
   entry = manifestLocationEntry(loc, method);
   values = string(leg.(field));
   for value = reshape(values, 1, [])
      filename = absolutePath(fullfile(input_root, kind, value));
      key = char(filename);
      if ~isKey(locations.by_file, key)
         locations.by_file(key) = entry;
         locations.ambiguous_file(key) = false;
         continue
      end
      prior = locations.by_file(key);
      if ~sameManifestLocationEntry(prior, entry)
         locations.ambiguous_file(key) = true;
      end
   end
end

function entry = manifestLocationEntry(loc, method)
   %MANIFESTLOCATIONENTRY Bundle requested point and sampling method metadata.
   entry = struct('location', loc, 'sample_method', string(method));
end

function tf = sameManifestLocationEntry(a, b)
   %SAMEMANIFESTLOCATIONENTRY True when exact-file metadata agrees.
   tf = sameLocation(a.location, b.location) ...
      && string(a.sample_method) == string(b.sample_method);
end

function tf = sameLocation(a, b)
   %SAMELOCATION True when two manifest locations identify the same point.
   tf = abs(double(a.lat_wgs84) - double(b.lat_wgs84)) <= 1e-8 ...
      && abs(double(a.lon_wgs84) - double(b.lon_wgs84)) <= 1e-8;
end

function files = rcmArtifactFiles(input_root, locations, source_ids)
   %RCMARTIFACTFILES Return current product-token RCM artifact files.
   patterns = [
      fullfile(input_root, "met", "mar3.11", "*.mat")
      fullfile(input_root, "met", "merra2", "*.mat")
      fullfile(input_root, "userdata", "mar3.11", "*.mat")
      fullfile(input_root, "userdata", "merra2", "*.mat")
      fullfile(input_root, "userdata", "racmo2.3p3", "*.mat")];
   files = strings(0, 1);
   for pattern = reshape(patterns, 1, [])
      hits = dir(pattern);
      hits = hits(~[hits.isdir]);
      next = strings(numel(hits), 1);
      for k = 1:numel(hits)
         next(k) = string(fullfile(hits(k).folder, hits(k).name));
      end
      files = [files; next]; %#ok<AGROW>
   end
   files = sort(files);
   if ~isempty(source_ids)
      % A source-scoped bounded repair must not restamp unrelated RCM files.
      keep = false(size(files));
      for k = 1:numel(files)
         record = parseArtifactFilename(files(k));
         keep(k) = ismember(record.source_id, source_ids);
      end
      files = files(keep);
   end
   if locations.restrict_to_files
      % A family-scoped repair must touch exact manifest references only;
      % alias fallback is deliberately disabled outside those references.
      keep = false(size(files));
      for k = 1:numel(files)
         keep(k) = isKey(locations.by_file, char(absolutePath(files(k))));
      end
      files = files(keep);
   end
end

function record = repairOne(filename, locations, state, dry_run, repair_racmo_subl)
   %REPAIRONE Classify and optionally repair one staged MAT file.
   record = parseArtifactFilename(filename);
   if record.status ~= "parsed"
      return
   end
   record.hash_before = ...
      icemodel.verification.setup.fileSha256(filename);

   [found, location, sample_method, ambiguous] = artifactLocation( ...
      filename, record.alias, locations);
   if ambiguous
      record.status = "ambiguous";
      record.reason = "current manifests disagree for this artifact";
      return
   end
   if ~found
      record.status = "unmapped";
      record.reason = "alias not present in current manifests";
      return
   end

   try
      info = whos("-file", filename);
      names = string({info.name});
      payload = intersect(names, ["met", "Data"], "stable");
      if isempty(payload)
         record.status = "skipped";
         record.reason = "no met or Data payload";
         return
      end

      % Restamp the table/timetable exactly as current writers do, then attach
      % requested-site and sampled-grid metadata.
      variable = payload(1);
      if ismember("artifact_metadata", names)
         S = load(filename, variable, "artifact_metadata", "-mat");
      else
         S = load(filename, variable, "-mat");
      end
      T = S.(variable);
      before = T;
      existing_metadata = struct();
      if isfield(S, "artifact_metadata") && isstruct(S.artifact_metadata)
         existing_metadata = S.artifact_metadata;
      end
      metadata = mergeArtifactMetadata(T.Properties.UserData, ...
         existing_metadata);
      metadata_before = metadata;

      % Apply only current-builder transformations, then restamp metadata.
      [T, metadata, actions, record, ok, reason] = ...
         repairSourceSchema(T, metadata, record, location, sample_method, ...
         state, variable, repair_racmo_subl);
      if ~ok
         if record.modis_status == "repair_failed"
            record.status = "repair_failed";
         else
            record.status = "restage_required";
         end
         record.reason = reason;
         return
      end
      transformed = T;
      T = icemodel.forcing.helpers.stampMetadata(T, strict=false);
      metadata = artifactMetadata(T, location, sample_method, metadata);
      T.Properties.UserData = metadata;
      if ~isequaln(transformed, T)
         actions(end + 1) = "restamp_metadata";
      end
      record.unrelated_payload_preserved = preservesUnrelatedPayload(before, T);
      if ~record.unrelated_payload_preserved
         record.status = "error";
         record.reason = "unrelated payload changed during artifact repair";
         return
      end
      record.changed_variables = cellstr(changedVariables(before, T));
      record.changed_metadata_fields = cellstr( ...
         changedMetadataFields(metadata_before, metadata));

      record.variable = variable;
      record.has_artifact_metadata = ismember("artifact_metadata", names);
      top_level_changed = ~record.has_artifact_metadata ...
         || ~isequaln(existing_metadata, metadata);
      if top_level_changed
         actions(end + 1) = "sync_artifact_metadata";
      end
      changed = ~isequaln(before, T) || top_level_changed;
      record.actions = cellstr(unique(actions, "stable"));
      if ~changed
         record.status = "unchanged";
         record.reason = "";
         record.hash_after = record.hash_before;
      elseif dry_run
         record.status = "would_repair";
         record.reason = "";
      else
         artifact_metadata = T.Properties.UserData;
         saveRepairedArtifact(filename, variable, T, artifact_metadata);
         record.status = "repaired";
         record.reason = "";
         record.hash_after = ...
            icemodel.verification.setup.fileSha256(filename);
      end
   catch err
      record.status = "error";
      record.reason = string(err.message);
   end
end

function names = changedVariables(before, after)
   %CHANGEDVARIABLES List payload variables whose values or names changed.
   before_names = string(before.Properties.VariableNames);
   after_names = string(after.Properties.VariableNames);
   names = setxor(before_names, after_names, "stable");
   common = intersect(before_names, after_names, "stable");
   for name = reshape(common, 1, [])
      if ~isequaln(before.(name), after.(name))
         names(end + 1) = name; %#ok<AGROW>
      end
   end
end

function names = changedMetadataFields(before, after)
   %CHANGEDMETADATAFIELDS List added, removed, or value-changed metadata.
   before_names = string(fieldnames(before));
   after_names = string(fieldnames(after));
   names = setxor(before_names, after_names, "stable");
   common = intersect(before_names, after_names, "stable");
   for name = reshape(common, 1, [])
      if ~isequaln(before.(name), after.(name))
         names(end + 1) = name; %#ok<AGROW>
      end
   end
end

function saveRepairedArtifact(filename, variable, T, artifact_metadata)
   %SAVEREPAIREDARTIFACT Replace an artifact only after a full sibling save.
   temporary = string(tempname(fileparts(filename))) + ".mat";
   cleanup = onCleanup(@() deleteIfPresent(temporary));
   [copied, message] = copyfile(filename, temporary, 'f');
   if ~copied
      error('icemodel:verification:repairRcmArtifactMetadata:copyFailed', ...
         'could not copy %s before repair: %s', filename, message)
   end
   if variable == "met"
      met = T;
      save(temporary, "met", "artifact_metadata", '-append')
   else
      Data = T;
      save(temporary, "Data", "artifact_metadata", '-append')
   end
   [moved, message] = movefile(temporary, filename, 'f');
   if ~moved
      error('icemodel:verification:repairRcmArtifactMetadata:replaceFailed', ...
         'could not replace %s: %s', filename, message)
   end
   clear cleanup
end

function deleteIfPresent(filename)
   %DELETEIFPRESENT Remove an incomplete sibling repair file after failure.
   if isfile(filename)
      delete(filename)
   end
end

function metadata = mergeArtifactMetadata(userdata, artifact_metadata)
   %MERGEARTIFACTMETADATA Preserve top-level fields absent from UserData.
   if isempty(userdata) || ~isstruct(userdata)
      metadata = struct();
   else
      metadata = userdata;
   end
   if isempty(artifact_metadata) || ~isstruct(artifact_metadata)
      return
   end
   names = fieldnames(artifact_metadata);
   for k = 1:numel(names)
      name = names{k};
      if ~isfield(metadata, name) || isempty(metadata.(name))
         metadata.(name) = artifact_metadata.(name);
      end
   end
end

function [T, metadata, actions, record, ok, reason] = repairSourceSchema( ...
      T, metadata, record, location, sample_method, state, variable, ...
      repair_racmo_subl)
   %REPAIRSOURCESCHEMA Apply source-aware canonical transforms exactly once.
   actions = strings(1, 0);
   ok = true;
   reason = "";

   % Align saved MAR daily ledgers before any source-aware transform. Current
   % clipped artifacts need only metadata repair, including derived met; no raw
   % MAR directory is consulted for this bounded operation.
   if record.source_id == "mar3.11"
      [source_days, source_ok, source_reason] = ...
         inferSavedMarSourceDays(metadata, T.Properties.RowTimes);
      if ~source_ok
         ok = false;
         reason = source_reason;
         return
      end
      before_alignment = metadata;
      try
         metadata = icemodel.forcing.helpers.alignMarDailyMetadata( ...
            metadata, source_days, T.Properties.RowTimes);
      catch err
         ok = false;
         reason = "MAR daily metadata alignment failed: " + string(err.message);
         return
      end
      if ~isequaln(metadata, before_alignment)
         actions(end + 1) = "align_mar_daily_metadata";
      end

      % Signed native RZ provenance is derivable from the saved artifact. Keep
      % this repair source-light for both hourly Data and derived met while
      % leaving every value and time posting untouched.
      before_refreeze = metadata;
      metadata = icemodel.forcing.helpers.marRefreezeMetadata(T, metadata);
      if ~isequaln(metadata, before_refreeze)
         actions(end + 1) = "stamp_mar_refreeze_signed_metadata";
      end
   end

   % MAR public diagnostics come from native daily delayed RU and daily SMB,
   % expanded as constant hourly rates. The source-year cache reads each
   % annual daily product once for every exact referenced point.
   if record.source_id == "mar3.11" && variable == "Data"
      % Native daily RU/SMB constrains only hourly source Data. A derived met
      % artifact may retain that provenance, but applying the 24-posting repair
      % to 15-minute values would multiply its daily integral.
      [T, metadata, mar_actions, record, ok, reason] = ...
         repairMarDailyDiagnostics(T, metadata, record, state);
      actions = [actions, mar_actions];
      if ~ok
         return
      end
   end

   % RACMO native precipitation is one canonical total-precipitation channel.
   names = string(T.Properties.VariableNames);
   if record.source_id == "racmo2.3p3"
      if all(ismember(["precip", "ppt"], names))
         ok = false;
         reason = "RACMO artifact contains both precip and ppt";
         return
      end
      if ismember("precip", names)
         T = renamevars(T, "precip", "ppt");
         actions(end + 1) = "rename_racmo_precip_to_ppt";
      end

      if ismember("ppt", string(T.Properties.VariableNames))
         % Enforce the same physical precipitation invariant as the current
         % builder only when this source-specific artifact carries ppt. Passing
         % prior metadata retains the original minimum/count on pass two.
         ppt = T.ppt;
         replaced_this_pass = any(isfinite(ppt) & ppt < 0, 'all');
         [T, ppt_qc] = ...
            icemodel.forcing.helpers.applyRacmoPrecipitationQualityControl( ...
            T, metadata);
         if replaced_this_pass
            actions(end + 1) = "zero_negative_racmo_ppt";
         end
         qc_fields = fieldnames(ppt_qc);
         qc_changed = false;
         for k = 1:numel(qc_fields)
            field = qc_fields{k};
            qc_changed = qc_changed || ~isfield(metadata, field) ...
               || ~isequaln(metadata.(field), ppt_qc.(field));
            metadata.(field) = ppt_qc.(field);
         end
         if qc_changed
            actions(end + 1) = "stamp_racmo_ppt_qc";
         end
      end

      % Legacy RACMO artifacts retain native negative-loss sublimation. Flip
      % that channel only when the paired durable markers prove it has not
      % already crossed the current builder boundary.
      if repair_racmo_subl
         [T, metadata, subl_actions, ok, reason] = ...
            repairRacmoSublimationSign(T, metadata);
         actions = [actions, subl_actions];
         if ~ok
            return
         end
      end
   end

   % MERRA sign orientation is either marked by a current builder/repair or
   % proven against one native daily source file before any numeric change.
   if record.source_id == "merra2"
      [T, metadata, merra_actions, record.merra_flux_orientation, ...
         ok, reason] = repairMerraFluxSigns(T, metadata, sample_method, state);
      actions = [actions, merra_actions];
      if ~ok
         return
      end

      if variable == "Data"
         % Only hourly Data retains exact tavg3 source rows at 3-hour stamps;
         % its legacy linear ramps are repairable only after the native glc
         % timestamp inventory proves which saved rows are genuine source rows.
         if height(T) > 1 && any(diff(T.Time) ~= hours(1))
            ok = false;
            reason = "MERRA Data time-support repair requires hourly input";
            return
         end
         [metadata, proof_changed, ok, reason] = ...
            ensureMerraTavg3SourceGrid(T, metadata, state);
         if ~ok
            return
         end
         if proof_changed
            actions(end + 1) = "stamp_merra_tavg3_source_grid";
         end
         [T, metadata, time_diagnostics] = ...
            icemodel.forcing.helpers.applyMerraTimeSupport(T, metadata);
         if time_diagnostics.replaced_count > 0
            actions(end + 1) = "hold_merra_tavg3_support";
         end
         if time_diagnostics.metadata_changed
            actions(end + 1) = "stamp_merra_time_support";
         end
      else
         % A 15-minute met file has already crossed the resampling boundary. Do
         % not reinterpret it as hourly Data; validate its own durable support
         % metadata and leave every numeric value untouched.
         [ok, reason] = validMerraMetTimeSupport(T, metadata);
         if ~ok
            return
         end
      end
   end

   % Optional MODIS augmentation is cached by requested location, method, and
   % year, so paired MAR/MERRA/RACMO artifacts never reread the same GEUS data.
   if state.modis_dir ~= ""
      [T, metadata, modis_actions, record, ok, reason] = repairModis( ...
         T, metadata, record, location, sample_method, state);
      actions = [actions, modis_actions];
   end
end

function [T, metadata, actions, ok, reason] = ...
      repairRacmoSublimationSign(T, metadata)
   %REPAIRRACMOSUBLIMATIONSIGN Normalize native RACMO sublimation exactly once.
   actions = strings(1, 0);
   ok = true;
   reason = "";
   if ~ismember("subl", string(T.Properties.VariableNames))
      return
   end

   native = "negative_loss_positive_deposition";
   canonical = "positive_loss_negative_deposition";
   native_marker = "";
   canonical_marker = "";
   if isfield(metadata, 'racmo_subl_native_sign_convention')
      native_marker = string(metadata.racmo_subl_native_sign_convention);
   end
   if isfield(metadata, 'racmo_subl_sign_convention')
      canonical_marker = string(metadata.racmo_subl_sign_convention);
   end
   if ~isscalar(native_marker) || ~isscalar(canonical_marker)
      ok = false;
      reason = "RACMO sublimation sign markers must be scalar";
      return
   end

   if native_marker == native && canonical_marker == canonical
      return
   end
   if native_marker ~= "" || canonical_marker ~= ""
      ok = false;
      reason = "RACMO sublimation sign markers are incomplete or unknown";
      return
   end

   T.subl = -T.subl;
   metadata.racmo_subl_native_sign_convention = char(native);
   metadata.racmo_subl_sign_convention = char(canonical);
   actions(end + 1) = "flip_racmo_subl_sign";
   actions(end + 1) = "stamp_racmo_subl_sign";
end

function [source_days, ok, reason] = inferSavedMarSourceDays(metadata, times)
   %INFERSAVEDMARSOURCEDAYS Recover only an unambiguous saved ledger day axis.
   ok = true;
   reason = "";
   times = times(:);
   times.TimeZone = 'UTC';
   if any(isnat(times))
      source_days = NaT(0, 1, 'TimeZone', 'UTC');
      ok = false;
      reason = "MAR artifact time axis contains NaT";
      return
   end
   retained_days = unique(dateshift(times, 'start', 'day'), 'stable');

   % Discover exactly the RU/SMB and ME/MEH per-day vectors. Summary counts and
   % textual status-code fields are deliberately excluded from axis inference.
   names = string(fieldnames(metadata));
   qc = startsWith(names, "mar_qc_") ...
      & (endsWith(names, "_day_status") ...
      | endsWith(names, "_daily_reference_mwe"));
   melt = ismember(names, ["mar_diagnostic_melt_day_status", ...
      "mar_diagnostic_melt_daily_reference_mwe", ...
      "mar_diagnostic_melt_residual_mwe_day"]);
   ledger_names = names(qc | melt);
   lengths = zeros(0, 1);
   for name = reshape(ledger_names, 1, [])
      values = metadata.(name);
      if ~isempty(values)
         lengths(end + 1, 1) = numel(values); %#ok<AGROW>
      end
   end

   % A current-method marker without any per-day vectors is not repairable from
   % filenames alone. Truly ledger-free legacy metadata remains a safe no-op.
   if isempty(lengths)
      marked = isfield(metadata, 'mar_qc_method') ...
         && string(metadata.mar_qc_method) == "daily_constrained_hourly";
      if marked && ~isempty(retained_days)
         source_days = retained_days;
         ok = false;
         reason = "MAR daily metadata marker has no saved per-day ledger";
      else
         source_days = retained_days;
      end
      return
   end
   if numel(unique(lengths)) ~= 1
      source_days = retained_days;
      ok = false;
      reason = "MAR saved per-day ledger lengths disagree";
      return
   end

   % An equal retained-day count is the idempotent already-aligned case. Only an
   % exact whole-calendar count authorizes reconstruction of the pre-clip axis.
   source_count = lengths(1);
   if source_count == numel(retained_days)
      source_days = retained_days;
      return
   end
   if isempty(retained_days)
      source_days = retained_days;
      ok = false;
      reason = "MAR nonempty daily ledger has no retained time axis";
      return
   end
   first_day = datetime(year(retained_days(1)), 1, 1, 'TimeZone', 'UTC');
   last_day = datetime(year(retained_days(end)), 12, 31, 'TimeZone', 'UTC');
   calendar_days = (first_day:caldays(1):last_day)';
   if source_count ~= numel(calendar_days)
      source_days = retained_days;
      ok = false;
      reason = "MAR saved ledger is neither retained-axis nor whole-calendar length";
      return
   end
   source_days = calendar_days;
end

function [ok, reason] = validMerraMetTimeSupport(T, metadata)
   %VALIDMERRAMETTIMESUPPORT Fail closed on noncanonical derived MERRA met.
   ok = false;
   reason = "";
   if ~istimetable(T) || height(T) < 2 ...
         || any(diff(T.Time) ~= minutes(15)) ...
         || any(mod(minute(T.Time), 15) ~= 0 | second(T.Time) ~= 0)
      reason = "MERRA met metadata repair accepts only uniform canonical " ...
         + "15-minute artifacts; regenerate this met from corrected hourly Data";
      return
   end

   % The generic met-resample contract must be independently complete before
   % source-specific MERRA metadata can be synchronized.
   required = ["met_resample_policy", ...
      "met_resample_expected_missing_counts", ...
      "met_resample_time_semantics", ...
      "met_resample_support_end_exclusive"];
   if ~isstruct(metadata) || ~all(isfield(metadata, required))
      reason = "MERRA 15-minute met lacks complete resampling provenance";
      return
   end
   policy = string(metadata.met_resample_policy);
   if isscalar(policy) && policy == "linear_adjacent_finite_only"
      reason = "legacy linear MERRA met requires " ...
         + "repairMetTimeSupport or regeneration from corrected hourly Data";
      return
   end
   allowed = ["interval_start_zero_order_hold", "native_15m_unchanged"];
   support_end = metadata.met_resample_support_end_exclusive;
   valid_generic = isscalar(policy) && ismember(policy, allowed) ...
      && isscalar(string(metadata.met_resample_time_semantics)) ...
      && string(metadata.met_resample_time_semantics) == "interval_start" ...
      && isstruct(metadata.met_resample_expected_missing_counts) ...
      && isdatetime(support_end) && isscalar(support_end) ...
      && ~isnat(support_end) && T.Time(end) == support_end - minutes(15);
    if ~valid_generic ...
          || ~icemodel.forcing.helpers.hasCanonicalMerraTimeSupport(metadata) ...
          || ~icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(T, metadata) ...
          || ~icemodel.forcing.helpers.hasConstantMerraTavg3Support(T)
      reason = "MERRA 15-minute met carries invalid timing/support provenance";
       return
    end

    % Missing-count provenance is a source-derived lower bound. Fewer observed
    % nonfinite values prove that an omitted support interval was bridged.
    expected = metadata.met_resample_expected_missing_counts;
    names = intersect(string(fieldnames(expected)), ...
       string(T.Properties.VariableNames), 'stable');
    for name = reshape(names, 1, [])
       values = T.(char(name));
       expected_count = double(expected.(char(name)));
       if ~isnumeric(values) || ~isscalar(expected_count) ...
             || ~isfinite(expected_count) || expected_count < 0 ...
             || expected_count ~= fix(expected_count) ...
             || nnz(~isfinite(values)) < expected_count
          reason = "MERRA 15-minute met bridges required missing source support";
          return
       end
    end
    ok = true;
end

function [metadata, changed, ok, reason] = ...
      ensureMerraTavg3SourceGrid(T, metadata, state)
   %ENSUREMERRATAVG3SOURCEGRID Prove native rows before legacy reconstruction.
   changed = false;
   ok = true;
   reason = "";
   if icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(T, metadata)
      return
   end
   channels = intersect(["runoff", "albedo", "snowd", "swe"], ...
      string(T.Properties.VariableNames), 'stable');
   if isempty(channels)
      return
   end
   if state.merra_dir == "" || ~isfolder(state.merra_dir)
      ok = false;
      reason = "native MERRA glc source is unavailable to prove tavg3 rows";
      return
   end

   % Identical artifact axes share one native inventory proof across all sites.
   source_times = T.Time(mod(hour(T.Time), 3) == 0);
   if isempty(source_times)
      ok = false;
      reason = "MERRA artifact has no recoverable tavg3 source stamps";
      return
   end
   key = char(compose('%s|%s|%d', ...
      string(source_times(1), 'yyyyMMddHHmmss'), ...
      string(source_times(end), 'yyyyMMddHHmmss'), numel(source_times)));
   if isKey(state.merra_grid, key)
      cached = state.merra_grid(key);
   else
      cached = proveMerraTavg3SourceGrid(source_times, state);
      state.merra_grid(key) = cached;
   end
   if ~cached.ok
      ok = false;
      reason = cached.reason;
      return
   end

   % Preserve unrelated provenance and add only the independently derived proof.
   fields = fieldnames(cached.proof);
   for k = 1:numel(fields)
      field = fields{k};
      metadata.(field) = cached.proof.(field);
   end
   changed = true;
   if ~icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(T, metadata)
      ok = false;
      reason = "native MERRA glc proof does not match the artifact time grid";
   end
end

function result = proveMerraTavg3SourceGrid(source_times, state)
   %PROVEMERRATAVG3SOURCEGRID Verify every required native daily time coordinate.
   result = struct('ok', false, 'reason', "", 'proof', struct());
   present = false(size(source_times));
   days = unique(dateshift(source_times, 'start', 'day'), 'stable');
   for day = reshape(days, 1, [])
      day_key = char(string(day, 'yyyyMMdd'));
      if ~isKey(state.merra_glc_files, day_key)
         result.reason = "native MERRA glc archive has a missing daily file";
         return
      end
      files = state.merra_glc_files(day_key);
      if numel(files) ~= 1
         result.reason = "native MERRA glc archive has ambiguous daily files";
         return
      end
      verified = cachedMerraGlcStarts(files(1), state);
      if ~verified.ok
         result.reason = verified.reason;
         return
      end
      in_day = dateshift(source_times, 'start', 'day') == day;
      present(in_day) = ismember(source_times(in_day), verified.starts);
   end
   if ~all(present)
      result.reason = "native MERRA glc archive has an incomplete tavg3 time grid";
      return
   end

   % A gap-free inventory makes every saved 00/03/... row an exact native row.
   missing = NaT(0, 1, 'TimeZone', 'UTC');
   result.proof = struct( ...
      'merra_tavg3_source_grid_policy', ...
      'native_glc_timestamp_inventory', ...
      'merra_tavg3_expected_source_row_count', numel(source_times), ...
      'merra_tavg3_source_row_count', numel(source_times), ...
      'merra_tavg3_source_time_gap_count', 0, ...
      'merra_tavg3_missing_source_times', missing);
   result.ok = true;
end

function verified = cachedMerraGlcStarts(filename, state)
   %CACHEDMERRAGLCSTARTS Read and validate one native daily glc coordinate once.
   key = char(filename);
   if isKey(state.merra_time, key)
      verified = state.merra_time(key);
      return
   end
   try
      native = icemodel.forcing.helpers.readMerra2Time(filename);
      day = dateshift(native(1), 'start', 'day');
      expected = day + hours((1.5:3:22.5)');
      ok = isequal(native, expected);
      starts = native - hours(1.5);
      if ok
         reason = "";
      else
         reason = "native MERRA glc file has noncanonical tavg3 timestamps";
      end
   catch err
      ok = false;
      starts = NaT(0, 1, 'TimeZone', 'UTC');
      reason = "could not read native MERRA glc timestamps: " + string(err.message);
   end
   verified = struct('ok', ok, 'reason', reason, 'starts', starts);
   state.merra_time(key) = verified;
end

function inventory = merraGlcInventory(merra_dir)
   %MERRAGLCINVENTORY Index daily native glc files by their YYYYMMDD token.
   inventory = containers.Map('KeyType', 'char', 'ValueType', 'any');
   if merra_dir == "" || ~isfolder(merra_dir)
      return
   end
   files = dir(fullfile(merra_dir, 'glc', '*_Nx.*.nc4*'));
   for k = 1:numel(files)
      token = regexp(files(k).name, '_Nx\.(\d{8})\.', 'tokens', 'once');
      if isempty(token)
         continue
      end
      key = token{1};
      filename = string(fullfile(files(k).folder, files(k).name));
      if isKey(inventory, key)
         inventory(key) = [inventory(key); filename];
      else
         inventory(key) = filename;
      end
   end
end

function tf = preservesUnrelatedPayload(before, after)
   %PRESERVESUNRELATEDPAYLOAD Verify only approved source channels changed.
   tf = height(before) == height(after);
   if tf && istimetable(before)
      % Timetables may preserve a source-specific row-dimension name, so compare
      % the canonical row-time property instead of assuming it is named Time.
      tf = isequaln(before.Properties.RowTimes, after.Properties.RowTimes);
   end
   if ~tf
      return
   end

   allowed = ["precip", "ppt", "shf", "lhf", "subl", "modis", "runoff", ...
      "smb", "albedo", "snowd", "swe"];
   before_names = string(before.Properties.VariableNames);
   after_names = string(after.Properties.VariableNames);
   changed_names = setxor(before_names, after_names);
   if ~all(ismember(changed_names, allowed))
      tf = false;
      return
   end
   common = setdiff(intersect(before_names, after_names, "stable"), allowed);
   for name = reshape(common, 1, [])
      if ~isequaln(before.(name), after.(name))
         tf = false;
         return
      end
   end
end

function [T, metadata, actions, record, ok, reason] = ...
      repairMarDailyDiagnostics(T, metadata, record, state)
   %REPAIRMARDAILYDIAGNOSTICS Replace RUH/SMBH products from daily RU/SMB.
   actions = strings(1, 0);
   ok = true;
   reason = "";

   % A complete native-daily marker plus constant within-day values is enough
   % for an idempotent second dry run; no raw source variable is read.
   if marDailyMetadataCurrent(T, metadata)
      record.mar_qc_status = "current";
      record.mar_qc_sector = double(metadata.mar_qc_sector);
      record.mar_replaced_runoff_count = ...
         double(metadata.mar_qc_replaced_runoff_count);
      record.mar_replaced_smb_count = ...
         double(metadata.mar_qc_replaced_smb_count);
      return
   end

   % Native daily repair is an explicit source-backed operation. Metadata-only
   % repair callers remain source-light unless mar_dir (or ICEMODEL_MAR_DIR)
   % is supplied; the bounded production repair passes that source explicitly.
   if state.mar_dir == ""
      record.mar_qc_status = "not_requested";
      return
   end

   file_key = char(absolutePath(record.filename));
   if ~isfolder(state.mar_dir)
      ok = false;
      reason = "native MAR directory is unavailable for daily RU/SMB repair";
      return
   elseif ~isKey(state.mar_plan.by_file, file_key)
      ok = false;
      reason = "MAR artifact has no exact nearest-point daily-source plan";
      return
   end
   entry = state.mar_plan.by_file(file_key);
   if ~entry.valid
      ok = false;
      reason = entry.reason;
      return
   end

   % The former native_daily policy already discarded RUH/SMBH structure.
   % Daily files alone cannot reconstruct it, so never relabel those artifacts
   % as hybrid-preserved; a full source build must reread the hourly products.
   if isfield(metadata, 'mar_qc_method') ...
         && string(metadata.mar_qc_method) == "native_daily"
      ok = false;
      reason = "legacy constant-daily MAR artifact requires full restage " ...
         + "to recover RUH/SMBH";
      return
   end

   % Prime each source year once, then project the daily accumulation / 24
   % onto the artifact's exact hourly rows. Partial first/last days receive the
   % same known daily rate on their available rows.
   replacements = struct('runoff', nan(height(T), 1), ...
      'smb', nan(height(T), 1));
   years = unique(year(T.Time))';
   for yyyy = years
      [primed, source_reason] = primeMarYear(yyyy, state);
      if ~primed
         ok = false;
         reason = source_reason;
         return
      end
      key = char(marSeriesKey(yyyy, entry));
      if ~isKey(state.mar_series, key)
         ok = false;
         reason = "MAR daily source cache does not contain the artifact point";
         return
      end
      series = state.mar_series(key);
      inyear = year(T.Time) == yyyy;
      day = dateshift(T.Time(inyear), 'start', 'day');
      [found, index] = ismember(day, series.time);
      if ~all(found)
         ok = false;
         reason = "MAR daily source does not cover every artifact UTC day";
         return
      end
      replacements.runoff(inyear) = series.runoff(index) / 24;
      replacements.smb(inyear) = series.smb(index) / 24;
   end

   [T, metadata] = ...
      icemodel.forcing.helpers.applyMarDailyQualityControl( ...
      T, replacements, metadata, sector=entry.sector);
   actions(end + 1) = "apply_mar_daily_constrained_qc";
   record.mar_qc_status = "applied";
   record.mar_qc_sector = entry.sector;
   record.mar_replaced_runoff_count = ...
      double(metadata.mar_qc_replaced_runoff_count);
   record.mar_replaced_smb_count = ...
      double(metadata.mar_qc_replaced_smb_count);
end

function tf = marDailyMetadataCurrent(T, metadata)
   %MARDAILYMETADATACURRENT Validate a source-light idempotence marker.
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
   tf = all(isfield(metadata, required));
   if ~tf
      return
   end

   names = string(T.Properties.VariableNames);
   channels = string(metadata.mar_qc_channels);
   sector = double(metadata.mar_qc_sector);
   sector_name = string(metadata.mar_qc_sector_name);
   expected_name = "";
   if isscalar(sector) && sector == 1
      expected_name = "permanent_ice";
   elseif isscalar(sector) && sector == 2
      expected_name = "tundra";
   end
   expected_basis = ...
      "MAR hourly RUH/SMBH preserved where complete UTC-day sums agree " + ...
      "with native daily RU/SMB; missing, partial, or inconsistent days " + ...
      "use the native daily rate";
   tf = string(metadata.mar_qc_method) == "daily_constrained_hourly" ...
      && string(metadata.mar_qc_status) == "applied" ...
      && string(metadata.mar_qc_fallback) == "none" ...
      && all(ismember(["runoff", "smb"], names)) ...
      && isequal(channels(:), ["runoff"; "smb"]) ...
      && expected_name ~= "" && sector_name == expected_name ...
      && string(metadata.mar_qc_runoff_source) == "RU" ...
      && string(metadata.mar_qc_smb_source) == "SMB" ...
      && string(metadata.mar_qc_hourly_distribution) == ...
      "preserve_matching_hourly_else_daily_divided_by_24" ...
      && string(metadata.mar_qc_partial_day_policy) == ...
      "native_daily_rate_applied_to_available_rows_marked_replaced" ...
      && string(metadata.mar_qc_day_status_codes) == ...
      "1=preserved_hourly;2=replaced_daily;3=unverified" ...
      && string(metadata.mar_qc_daily_reference_units) == "mWE/day" ...
      && string(metadata.mar_qc_basis) == expected_basis;
   if ~tf
      return
   end

   % The per-day source reference makes the hybrid contract independently
   % checkable without reopening MAR. Only replaced days must be constant;
   % preserved days retain their validated native hourly structure.
   [groups, days] = findgroups(dateshift(T.Time, 'start', 'day'));
   complete = false(numel(days), 1);
   for day = 1:numel(days)
      complete(day) = isequal(T.Time(groups == day), ...
         days(day) + hours((0:23)'));
   end
   tf = double(metadata.mar_qc_complete_utc_day_count) ...
      == nnz(complete) ...
      && double(metadata.mar_qc_partial_utc_day_count) ...
      == nnz(~complete) ...
      && isscalar(metadata.mar_qc_replaced_runoff_count) ...
      && isscalar(metadata.mar_qc_replaced_smb_count);
   if ~tf
      return
   end
   tf = marLedgerCurrent(T.runoff, groups, complete, metadata, "runoff") ...
      && marLedgerCurrent(T.smb, groups, complete, metadata, "smb");
end

function tf = marLedgerCurrent(values, groups, complete, metadata, channel)
   %MARLEDGERCURRENT Verify status counts, daily totals, and replacement shape.
   status = uint8(metadata.("mar_qc_" + channel + "_day_status"));
   status = status(:);
   reference = double(metadata.( ...
      "mar_qc_" + channel + "_daily_reference_mwe"));
   reference = reference(:);
   complete = complete(:);
   ndays = max([groups; 0]);
   tf = numel(status) == ndays && numel(reference) == ndays ...
      && all(ismember(status, uint8([1 2 3]))) ...
      && double(metadata.("mar_qc_preserved_" + channel ...
      + "_day_count")) == nnz(status == 1) ...
      && double(metadata.("mar_qc_replaced_" + channel ...
      + "_day_count")) == nnz(status == 2) ...
      && double(metadata.("mar_qc_unverified_" + channel ...
      + "_day_count")) == nnz(status == 3);
   if tf
      % Unverified means the daily reference is unavailable, while preserved
      % hourly structure is possible only on an exact 24-row UTC day.
      tf = all(isnan(reference(status == 3))) ...
         && all(isfinite(reference(status ~= 3))) ...
         && ~any(status == 1 & ~complete);
   end
   if ~tf
      return
   end
   abs_tolerance = double(metadata.mar_qc_abs_tolerance_mwe_day);
   rel_tolerance = double(metadata.mar_qc_rel_tolerance);
   for day = 1:ndays
      if status(day) == 3
         continue
      end
      hourly = values(groups == day);
      if any(~isfinite(hourly)) || ~isfinite(reference(day))
         tf = false;
         return
      end
      if status(day) == 2
         scale = max(1, max(abs(hourly)));
         tf = max(hourly) - min(hourly) <= 16 * eps(scale);
      end
      if numel(hourly) == 24
         tolerance = abs_tolerance + rel_tolerance ...
            * max(abs([sum(hourly), reference(day)]));
         numeric_slack = 16 * eps(max(1, ...
            max(abs([sum(hourly), reference(day)]))));
         tf = tf && abs(sum(hourly) - reference(day)) ...
            <= tolerance + numeric_slack;
      end
      if ~tf
         return
      end
   end
end

function plan = marPointPlan(files, locations, mar_dir)
   %MARPOINTPLAN Resolve every exact MAR artifact to one source-grid sector.
   plan = struct();
   plan.by_file = containers.Map('KeyType', 'char', 'ValueType', 'any');
   prototype = struct('valid', false, 'reason', "", 'pixel_i', 0, ...
      'pixel_j', 0, 'pixel_linear', 0, 'sector', 0, 'cache_key', "");
   plan.points = repmat(prototype, 0, 1);
   if ~isfolder(mar_dir)
      return
   end
   source_files = dir(fullfile(mar_dir, '*.nc'));
   if isempty(source_files)
      return
   end

   % Reproduce buildMarData's requested-geographic -> native-grid nearest
   % mapping once. Every current PROMICE artifact declares nearest sampling.
   grid_file = string(fullfile(source_files(1).folder, source_files(1).name));
   grid = icemodel.forcing.marGridInfo(grid_file);
   to_native_x = scatteredInterpolant(grid.LON(:), grid.LAT(:), ...
      grid.Xnat(:), 'natural', 'nearest');
   to_native_y = scatteredInterpolant(grid.LON(:), grid.LAT(:), ...
      grid.Ynat(:), 'natural', 'nearest');
   seen = strings(0, 1);
   for k = 1:numel(files)
      record = parseArtifactFilename(files(k));
      if record.source_id ~= "mar3.11"
         continue
      end
      entry = prototype;
      [found, location, sample_method, ambiguous] = artifactLocation( ...
         files(k), record.alias, locations);
      if ambiguous || ~found
         entry.reason = "MAR artifact location is missing or ambiguous";
      elseif string(sample_method) ~= "nearest"
         entry.reason = "non-nearest MAR artifacts require full restage";
      else
         xq = to_native_x(double(location.lon_wgs84), ...
            double(location.lat_wgs84));
         yq = to_native_y(double(location.lon_wgs84), ...
            double(location.lat_wgs84));
         [~, index] = min(hypot(grid.Xnat(:) - xq, grid.Ynat(:) - yq));
         [entry.pixel_i, entry.pixel_j] = ind2sub(size(grid.LAT), index);
         entry.pixel_linear = index;
         entry.sector = 1 + double(grid.srf(index) ~= 4);
         entry.cache_key = compose('%03d_%03d_%d', ...
            entry.pixel_i, entry.pixel_j, entry.sector);
         entry.valid = true;
      end
      plan.by_file(char(absolutePath(files(k)))) = entry;
      if entry.valid && ~ismember(entry.cache_key, seen)
         plan.points(end + 1, 1) = entry;
         seen(end + 1, 1) = entry.cache_key; %#ok<AGROW>
      end
   end
end

function [ok, reason] = primeMarYear(yyyy, state)
   %PRIMEMARYEAR Read annual daily RU/SMB once and cache every planned point.
   ok = true;
   reason = "";
   points = state.mar_plan.points;
   if isempty(points)
      ok = false;
      reason = "MAR daily-source point plan is empty";
      return
   end
   first_key = char(marSeriesKey(yyyy, points(1)));
   if isKey(state.mar_series, first_key)
      return
   end
   matches = dir(fullfile(state.mar_dir, sprintf('*-%d.nc', yyyy)));
   if numel(matches) ~= 1
      ok = false;
      reason = sprintf('expected one MAR daily source for %d, found %d', ...
         yyyy, numel(matches));
      return
   end
   filename = string(fullfile(matches.folder, matches.name));
   info = ncinfo(filename, 'RU');
   smb_info = ncinfo(filename, 'SMB');
   if numel(info.Size) ~= 4 || info.Size(3) ~= 2 ...
         || numel(smb_info.Size) ~= 4 || smb_info.Size(3) ~= 2
      ok = false;
      reason = "MAR daily RU/SMB source lacks the two-sector layout";
      return
   end
   if string(ncreadatt(filename, 'RU', 'units')) ~= "mmWE/day" ...
         || string(ncreadatt(filename, 'SMB', 'units')) ~= "mmWE/day"
      ok = false;
      reason = "MAR daily RU/SMB source units are not mmWE/day";
      return
   end
   ndays = info.Size(end);
   t0 = datetime(yyyy, 1, 1, 'TimeZone', 'UTC');
   daily_time = (t0:days(1):(t0 + days(ndays - 1)))';

   % Read each full daily product once because its NetCDF chunks span the full
   % spatial grid. Extracting all planned pixels now avoids one decompression
   % of every day chunk per artifact.
   runoff = selectedMarDaily(filename, 'RU', points, info.Size);
   state.mar_reads(char(filename + '|RU')) = true;
   smb = selectedMarDaily(filename, 'SMB', points, smb_info.Size);
   state.mar_reads(char(filename + '|SMB')) = true;
   for k = 1:numel(points)
      series = struct('time', daily_time, 'runoff', runoff(:, k), ...
         'smb', smb(:, k));
      state.mar_series(char(marSeriesKey(yyyy, points(k)))) = series;
   end
end

function selected = selectedMarDaily(filename, variable, points, dims)
   %SELECTEDMARDAILY Retain planned point/sectors from one annual daily field.
   raw = ncread(filename, variable);
   raw = reshape(raw, prod(dims(1:2)), dims(3), dims(4));
   selected = nan(dims(4), numel(points));
   for k = 1:numel(points)
      values = double(reshape(raw(points(k).pixel_linear, ...
         points(k).sector, :), [], 1));
      values(abs(values) >= 1e30) = NaN;
      selected(:, k) = values / 1000;
   end
end

function key = marSeriesKey(yyyy, entry)
   %MARSERIESKEY Identify one native daily point/sector series in one year.
   key = string(yyyy) + '|' + entry.cache_key;
end

function mar_dir = resolveMarDir(mar_dir)
   %RESOLVEMARDIR Resolve an explicitly enabled native MAR repair source.
   if mar_dir ~= ""
      return
   end
   mar_dir = string(getenv("ICEMODEL_MAR_DIR"));
end

function [T, metadata, actions, orientation, ok, reason] = ...
      repairMerraFluxSigns(T, metadata, sample_method, state)
   %REPAIRMERRAFLUXSIGNS Normalize MERRA turbulent fluxes exactly once.
   actions = strings(1, 0);
   orientation = "not_applicable";
   ok = true;
   reason = "";
   channels = intersect(["shf", "lhf"], ...
      string(T.Properties.VariableNames), "stable");
   if isempty(channels)
      return
   end

   canonical = "positive_toward_surface";
   marker = "";
   if isfield(metadata, "merra_flux_sign_convention")
      marker = string(metadata.merra_flux_sign_convention);
   end
   if numel(marker) > 1
      ok = false;
      reason = "MERRA flux-sign marker is not scalar";
      return
   end

   if ~ismember(marker, ["", canonical, "positive_upward"])
      ok = false;
      reason = "unknown MERRA flux-sign marker: " + marker;
      return
   elseif marker ~= "" && ~state.verify_merra_sign
      orientation = marker;
   else
      % Verification mode deliberately ignores an existing marker and compares
      % numeric values to native HFLUX/EFLUX, catching a previously double-
      % flipped artifact whose metadata still claims the canonical convention.
      [orientation, reason] = inferMerraFluxOrientation( ...
         T, metadata, sample_method, state);
      if orientation == "unknown"
         ok = false;
         return
      end
   end

   % A positive-upward artifact matches native HFLUX/EFLUX and needs exactly
   % one sign inversion. A source-proven canonical artifact is only marked.
   if orientation == "positive_upward"
      for channel = reshape(channels, 1, [])
         T.(channel) = -T.(channel);
      end
      actions(end + 1) = "flip_merra_flux_sign";
   end
   if marker ~= canonical
      metadata.merra_flux_sign_convention = char(canonical);
      actions(end + 1) = "mark_merra_flux_sign";
   end
end

function [orientation, reason] = inferMerraFluxOrientation( ...
      T, metadata, sample_method, state)
   %INFERMERRAFLUXORIENTATION Compare unmarked data to native daily fluxes.
   orientation = "unknown";
   reason = "MERRA fluxes could not be compared with native source";
   if string(sample_method) ~= "nearest"
      reason = "non-nearest MERRA fluxes require full restage for sign audit";
      return
   end
   if state.merra_dir == "" || ~isfolder(state.merra_dir)
      reason = "native MERRA source is unavailable for sign audit";
      return
   end

   [found, source_lat, source_lon] = sourceLocation(T, metadata);
   if ~found
      reason = "MERRA source-grid coordinates are missing";
      return
   end

   % Try a bounded set of artifact days; zero-flux days can be inconclusive.
   days_present = unique(dateshift(T.Time, "start", "day"), "stable");
   checked = 0;
   for day = reshape(days_present, 1, [])
      matches = dir(fullfile(state.merra_dir, "flx", ...
         "*_Nx." + string(day, "yyyyMMdd") + ".nc4*"));
      if isempty(matches)
         continue
      end
      if numel(matches) > 1
         reason = "ambiguous native MERRA daily flux file";
         return
      end
      filename = string(fullfile(matches.folder, matches.name));
      raw = cachedMerraFlux(filename, source_lat, source_lon, state);
      [same_error, flipped_error, informative] = ...
         merraOrientationErrors(T, day, raw);
      checked = checked + 1;
      if informative && flipped_error <= 1e-8 ...
            && flipped_error * 10 < same_error
         orientation = "positive_toward_surface";
         reason = "";
         return
      elseif informative && same_error <= 1e-8 ...
            && same_error * 10 < flipped_error
         orientation = "positive_upward";
         reason = "";
         return
      end
      if checked >= 5
         reason = sprintf( ...
            'native MERRA sign comparison was inconclusive after %d days', ...
            checked);
         return
      end
   end
end

function raw = cachedMerraFlux(filename, source_lat, source_lon, state)
   %CACHEDMERRAFLUX Read one source-grid cell from one daily MERRA file once.
   lat = double(ncread(filename, "lat"));
   lon = double(ncread(filename, "lon"));
   [~, ilat] = min(abs(lat - source_lat));
   [~, ilon] = min(abs(lon - source_lon));
   key = char(filename + "|" + string(ilon) + "|" + string(ilat));
   if isKey(state.merra_flux, key)
      raw = state.merra_flux(key);
      return
   end

   raw = struct();
   raw.shf = reshape(icemodel.forcing.readMerra2( ...
      filename, "HFLUX", start=[ilon, ilat], count=[1, 1]), [], 1);
   raw.lhf = reshape(icemodel.forcing.readMerra2( ...
      filename, "EFLUX", start=[ilon, ilat], count=[1, 1]), [], 1);
   state.merra_flux(key) = raw;
end

function [same_error, flipped_error, informative] = ...
      merraOrientationErrors(T, day, raw)
   %MERRAORIENTATIONERRORS Compare hourly artifact values with native signs.
   on_day = dateshift(T.Time, "start", "day") == day ...
      & minute(T.Time) == 0 & second(T.Time) == 0;
   staged_values = zeros(0, 1);
   raw_values = zeros(0, 1);
   for channel = ["shf", "lhf"]
      if ~ismember(channel, string(T.Properties.VariableNames))
         continue
      end
      staged = double(T.(channel)(on_day));
      source = double(raw.(channel));
      n = min(numel(staged), numel(source));
      keep = isfinite(staged(1:n)) & isfinite(source(1:n));
      staged_values = [staged_values; staged(keep)]; %#ok<AGROW>
      raw_values = [raw_values; source(keep)]; %#ok<AGROW>
   end
   scale = max(1, norm(raw_values));
   informative = ~isempty(raw_values) && scale > 1;
   same_error = norm(staged_values - raw_values) / scale;
   flipped_error = norm(staged_values + raw_values) / scale;
end

function [found, source_lat, source_lon] = sourceLocation(T, metadata)
   %SOURCELOCATION Resolve sampled source-grid coordinates from an artifact.
   found = false;
   source_lat = NaN;
   source_lon = NaN;
   custom = T.Properties.CustomProperties;
   if isprop(custom, "Lat") && isprop(custom, "Lon")
      source_lat = double(custom.Lat);
      source_lon = double(custom.Lon);
      found = isscalar(source_lat) && isscalar(source_lon) ...
         && all(isfinite([source_lat, source_lon]));
      if found
         return
      end
   end
   if all(isfield(metadata, ["source_lat_wgs84", "source_lon_wgs84"]))
      source_lat = double(metadata.source_lat_wgs84);
      source_lon = double(metadata.source_lon_wgs84);
      found = isscalar(source_lat) && isscalar(source_lon) ...
         && all(isfinite([source_lat, source_lon]));
   end
end

function [T, metadata, actions, record, ok, reason] = repairModis( ...
      T, metadata, record, location, sample_method, state)
   %REPAIRMODIS Add or remove MODIS according to exact source-year coverage.
   actions = strings(1, 0);
   ok = true;
   reason = "";
   artifact_years = unique(year(T.Time))';
   [coverage, missing, ambiguous] = modisCoverage( ...
      artifact_years, state.modis_inventory);
   record.modis_coverage_years = coverage;
   record.modis_missing_years = missing;
   if ~isempty(ambiguous)
      record.modis_status = "repair_failed";
      ok = false;
      reason = "ambiguous GEUS MODIS files for years " ...
         + strjoin(string(ambiguous), ", ");
      return
   end

   modis_metadata = ...
      icemodel.forcing.helpers.geusModisCoverageMetadata( ...
      artifact_years, coverage);
   expected_status = string(modis_metadata.modis_status);
   record.modis_status = expected_status;
   if modisMetadataCurrent(T, metadata, modis_metadata)
      if ismember("modis", string(T.Properties.VariableNames))
         record.modis_finite_count = nnz(validModisValues(T.modis));
      end
      return
   end

   names = string(T.Properties.VariableNames);
   if isempty(coverage)
      % All-NaN temporary channels with no source coverage are absent data, not
      % successful augmentation. Remove them while retaining an explicit marker.
      if ismember("modis", names) && any(validModisValues(T.modis))
         record.modis_status = "repair_failed";
         ok = false;
         reason = "valid MODIS values exist without source-year coverage";
         return
      elseif ismember("modis", names)
         T = removevars(T, "modis");
         actions(end + 1) = "remove_uncovered_modis";
      end
      metadata = stampModisMetadata(metadata, modis_metadata);
      actions(end + 1) = "mark_modis_no_source_coverage";
      return
   end

   % Build one channel from cached daily point series. Every covered year must
   % contribute physical finite albedo; otherwise source coverage and repair
   % failure are kept distinct in the report and the artifact is not written.
   modis = nan(height(T), 1);
   try
      for yyyy = reshape(coverage, 1, [])
         [daily, daily_time] = cachedModisSeries( ...
            yyyy, location, sample_method, state);
         inyear = year(T.Time) == yyyy;
         modis(inyear) = icemodel.forcing.helpers.dailyToHourly( ...
            daily, daily_time, T.Time(inyear));
      end
   catch err
      record.modis_status = "repair_failed";
      ok = false;
      reason = "MODIS repair failed: " + string(err.message);
      return
   end
   valid_by_year = arrayfun(@(yyyy) ...
      any(validModisValues(modis(year(T.Time) == yyyy))), coverage);
   if ~all(valid_by_year)
      record.modis_status = "repair_failed";
      ok = false;
      reason = "GEUS MODIS source coverage produced no physical artifact values";
      return
   end

   if ~ismember("modis", names) || ~isequaln(T.modis, modis)
      T.modis = modis;
      actions(end + 1) = "write_modis";
   end
   metadata = stampModisMetadata(metadata, modis_metadata);
   actions(end + 1) = "mark_modis_coverage";
   record.modis_finite_count = nnz(validModisValues(modis));
end

function inventory = modisInventory(modis_dir)
   %MODISINVENTORY Index exact GEUS source files once per repair run.
   inventory = containers.Map('KeyType', 'char', 'ValueType', 'any');
   if modis_dir == ""
      return
   end
   if ~isfolder(modis_dir)
      error('icemodel:verification:repairRcmArtifactMetadata:modisNotFound', ...
         'GEUS MODIS source directory not found: %s', modis_dir)
   end
   files = dir(fullfile(modis_dir, "*.nc*"));
   for k = 1:numel(files)
      token = regexp(files(k).name, '_(\d{4})_', 'tokens', 'once');
      if isempty(token)
         continue
      end
      key = token{1};
      filename = string(fullfile(files(k).folder, files(k).name));
      if isKey(inventory, key)
         inventory(key) = [string(inventory(key)); filename];
      else
         inventory(key) = filename;
      end
   end
end

function [coverage, missing, ambiguous] = modisCoverage(years, inventory)
   %MODISCOVERAGE Classify available, missing, and ambiguous source years.
   coverage = zeros(1, 0);
   missing = zeros(1, 0);
   ambiguous = zeros(1, 0);
   for yyyy = reshape(years, 1, [])
      key = char(string(yyyy));
      if ~isKey(inventory, key)
         missing(end + 1) = yyyy; %#ok<AGROW>
      elseif isscalar(inventory(key))
         coverage(end + 1) = yyyy; %#ok<AGROW>
      else
         ambiguous(end + 1) = yyyy; %#ok<AGROW>
      end
   end
end

function tf = modisMetadataCurrent(T, metadata, expected)
   %MODISMETADATACURRENT True when a prior repair is complete and reusable.
   tf = all(isfield(metadata, ...
      ["modis_product", "modis_status", "modis_coverage_years"]));
   if ~tf
      return
   end
   coverage = double(expected.modis_coverage_years(:))';
   tf = string(metadata.modis_product) == string(expected.modis_product) ...
      && string(metadata.modis_status) == string(expected.modis_status) ...
      && isequal(sort(double(metadata.modis_coverage_years(:)))', ...
      sort(double(coverage(:)))');
   if ~tf
      return
   end
   names = string(T.Properties.VariableNames);
   if isempty(coverage)
      tf = ~ismember("modis", names);
      return
   end
   if ~ismember("modis", names)
      tf = false;
      return
   end
   valid = validModisValues(T.modis);
   invalid = ~isnan(T.modis) & ~valid;
   tf = ~any(invalid) && all(arrayfun(@(yyyy) ...
      any(valid(year(T.Time) == yyyy)), coverage));
end

function valid = validModisValues(values)
   %VALIDMODISVALUES Identify finite physical albedo used for coverage proof.
   valid = isfinite(values) & imag(values) == 0 ...
      & real(values) >= 0 & real(values) <= 1;
end

function metadata = stampModisMetadata(metadata, expected)
   %STAMPMODISMETADATA Overwrite stale fields with the shared exact contract.
   fields = fieldnames(expected);
   for k = 1:numel(fields)
      metadata.(fields{k}) = expected.(fields{k});
   end
end

function plan = modisPointPlan(files, locations, inventory)
   %MODISPOINTPLAN Group exact referenced point requests by GEUS source year.
   plan = containers.Map('KeyType', 'char', 'ValueType', 'any');
   if inventory.Count == 0
      return
   end
   prototype = struct('location', struct(), 'sample_method', "", ...
      'cache_key', "");
   for k = 1:numel(files)
      record = parseArtifactFilename(files(k));
      [found, location, sample_method, ambiguous] = artifactLocation( ...
         files(k), record.alias, locations);
      if ~found || ambiguous || record.status ~= "parsed"
         continue
      end
      first_year = str2double(extractBetween(record.window_start, 1, 4));
      last_year = str2double(extractBetween(record.window_end, 1, 4));
      for yyyy = first_year:last_year
         year_key = char(string(yyyy));
         if ~isKey(inventory, year_key) || ~isscalar(inventory(year_key))
            continue
         end
         filename = string(inventory(year_key));
         entry = prototype;
         entry.location = location;
         entry.sample_method = string(sample_method);
         entry.cache_key = modisCacheKey( ...
            filename, location, sample_method);
         if isKey(plan, year_key)
            entries = plan(year_key);
         else
            entries = repmat(prototype, 0, 1);
         end
         if isempty(entries) ...
               || ~any(string({entries.cache_key}) == entry.cache_key)
            entries(end + 1, 1) = entry; %#ok<AGROW>
            plan(year_key) = entries;
         end
      end
   end
end

function [daily, daily_time] = cachedModisSeries( ...
      yyyy, location, sample_method, state)
   %CACHEDMODISSERIES Serve one point/year from a lazily primed year cache.
   files = string(state.modis_inventory(char(string(yyyy))));
   filename = files(1);
   key = char(modisCacheKey(filename, location, sample_method));
   if isKey(state.modis_series, key)
      cached = state.modis_series(key);
      daily = cached.daily;
      daily_time = cached.time;
      return
   end

   % A nearest-point cache miss primes every referenced point for this source
   % year from one bounding albedo read. Non-nearest requests retain the shared
   % public reader as a narrow fallback.
   if string(sample_method) == "nearest"
      primeModisYear(yyyy, state);
      if isKey(state.modis_series, key)
         cached = state.modis_series(key);
         daily = cached.daily;
         daily_time = cached.time;
         return
      end
   end
   [daily, daily_time] = icemodel.forcing.readGeusModis( ...
      filename, [double(location.lat_wgs84), double(location.lon_wgs84)], ...
      string(sample_method));
   state.modis_series(key) = struct('daily', daily, 'time', daily_time);
   state.modis_reads(key) = true;
end

function primeModisYear(yyyy, state)
   %PRIMEMODISYEAR Read one source-year block for every planned nearest point.
   year_key = char(string(yyyy));
   if ~isKey(state.modis_plan, year_key)
      return
   end
   entries = state.modis_plan(year_key);
   entries = entries(string({entries.sample_method}) == "nearest");
   if isempty(entries)
      return
   end
   filename = string(state.modis_inventory(year_key));

   % Reconstruct the exact regular native grid used by readGeusModis, then
   % locate all nearest point cells before issuing one bounding NetCDF read.
   LON = double(ncread(filename, 'lon'));
   LAT = double(ncread(filename, 'lat'));
   geus = icemodel.forcing.helpers.geusModisProjection();
   [Xp, Yp] = projfwd(geus, LAT, LON);
   xax = linspace(mean(Xp(1, :)), mean(Xp(end, :)), size(Xp, 1)).';
   yax = linspace(mean(Yp(:, 1)), mean(Yp(:, end)), size(Yp, 2));
   ii = zeros(numel(entries), 1);
   jj = zeros(numel(entries), 1);
   for k = 1:numel(entries)
      location = entries(k).location;
      [xq, yq] = projfwd(geus, double(location.lat_wgs84), ...
         double(location.lon_wgs84));
      [~, ii(k)] = min(abs(xax - xq));
      [~, jj(k)] = min(abs(yax - yq));
   end

   i0 = min(ii);
   j0 = min(jj);
   ni = max(ii) - i0 + 1;
   nj = max(jj) - j0 + 1;
   info = ncinfo(filename, 'albedo');
   ndays = info.Size(end);
   block = ncread(filename, 'albedo', [i0, j0, 1], [ni, nj, ndays]);
   % Keep the batched nearest-point repair path identical to the public reader.
   block = icemodel.forcing.helpers.normalizeGeusModisAlbedo(block);
   % The inventory already resolved this exact source year. Using that value
   % avoids mistaking a four-digit token in a parent work-directory name for
   % the product year.
   t0 = datetime(yyyy, 1, 1, 'TimeZone', 'UTC');
   daily_time = (t0:days(1):(t0 + days(ndays - 1)))';
   for k = 1:numel(entries)
      daily = reshape(block(ii(k) - i0 + 1, ...
         jj(k) - j0 + 1, :), [], 1);
      state.modis_series(char(entries(k).cache_key)) = ...
         struct('daily', daily, 'time', daily_time);
   end
   state.modis_reads(char(filename)) = true;
end

function key = modisCacheKey(filename, location, sample_method)
   %MODISCACHEKEY Identify one requested point/method series in one year file.
   key = filename + "|" + compose("%.10f", location.lat_wgs84) ...
      + "|" + compose("%.10f", location.lon_wgs84) ...
      + "|" + string(sample_method);
end

function merra_dir = resolveMerraDir(merra_dir)
   %RESOLVEMERRADIR Use the same native-source fallback as buildMerraData.
   if merra_dir ~= ""
      return
   end
   merra_dir = string(getenv("ICEMODEL_MERRA_DIR"));
   if merra_dir == ""
      merra_dir = "/Volumes/S03/DATA/merra2/1hrly/ncfiles";
   end
end

function [found, location, sample_method, ambiguous] = artifactLocation( ...
      filename, alias, locations)
   %ARTIFACTLOCATION Resolve requested-site metadata for one artifact.
   found = false;
   ambiguous = false;
   location = struct();
   sample_method = "nearest";
   file_key = char(absolutePath(filename));
   if isKey(locations.by_file, file_key)
      ambiguous = locations.ambiguous_file(file_key);
      if ambiguous
         return
      end
      entry = locations.by_file(file_key);
      location = entry.location;
      sample_method = entry.sample_method;
      found = true;
      return
   end

   alias_key = char(matlab.lang.makeValidName(alias));
   if ~isKey(locations.by_alias, alias_key)
      return
   end
   ambiguous = locations.ambiguous_alias(alias_key);
   if ambiguous
      return
   end
   entry = locations.by_alias(alias_key);
   location = entry.location;
   sample_method = entry.sample_method;
   found = true;
end

function metadata = artifactMetadata(T, location, sample_method, metadata)
   %ARTIFACTMETADATA Build the writer-compatible artifact metadata struct.
   if isempty(metadata) || ~isstruct(metadata)
      metadata = struct();
   end
   if ~isfield(metadata, "sample_method") || strlength(string(metadata.sample_method)) == 0
      metadata.sample_method = char(sample_method);
   end
   metadata.lat_wgs84 = double(location.lat_wgs84);
   metadata.lon_wgs84 = double(location.lon_wgs84);

   % Existing CustomProperties are source-grid metadata from the original RCM
   % extraction. Keep them separate from the requested site location.
   cp = T.Properties.CustomProperties;
   if isprop(cp, "Lat")
      metadata.source_lat_wgs84 = cp.Lat;
   end
   if isprop(cp, "Lon")
      metadata.source_lon_wgs84 = cp.Lon;
   end
   metadata = icemodel.forcing.helpers.columnizeMetadata(metadata);
end

function record = parseArtifactFilename(filename)
   %PARSEARTIFACTFILENAME Decode alias/product/window from a current RCM path.
   record = emptyRecord();
   record.filename = filename;
   [~, stem, ~] = fileparts(filename);
   expr = "^(?:met_)?(.+?)_(mar3\.11|merra2|racmo2\.3p3)_(\d{8})_(\d{8})(?:_(?:1hr|15m))?$";
   tokens = regexp(stem, expr, "tokens", "once");
   if isempty(tokens)
      expr = "^(?:met_)?(.+?)_(mar3\.11|merra2|racmo2\.3p3)_(\d{4})(?:_1hr)?$";
      tokens = regexp(stem, expr, "tokens", "once");
      if ~isempty(tokens)
         tokens{4} = tokens{3};
      end
   end
   if isempty(tokens)
      record.status = "skipped";
      record.reason = "filename does not match current RCM artifact pattern";
      return
   end
   record.status = "parsed";
   record.alias = string(tokens{1});
   record.source_id = string(tokens{2});
   record.window_start = string(tokens{3});
   record.window_end = string(tokens{4});
end

function pathname = absolutePath(pathname)
   %ABSOLUTEPATH Return one canonical absolute key for path-map lookups.
   pathname = string(pathname);
   if ~startsWith(pathname, filesep)
      pathname = string(fullfile(pwd, pathname));
   end

   % Java canonical paths resolve existing symlink components and lexically
   % normalize dot segments in missing suffixes without creating the target.
   path_object = java.io.File(char(pathname));
   pathname = string(char(path_object.getCanonicalPath()));
end

function record = emptyRecord()
   %EMPTYRECORD Prototype report row.
   record = struct( ...
      "filename", "", ...
      "alias", "", ...
      "source_id", "", ...
      "window_start", "", ...
      "window_end", "", ...
      "variable", "", ...
      "has_artifact_metadata", false, ...
      "actions", {{}}, ...
      "changed_variables", {{}}, ...
      "changed_metadata_fields", {{}}, ...
      "hash_before", "", ...
      "hash_after", "", ...
      "unrelated_payload_preserved", false, ...
      "mar_qc_status", "not_applicable", ...
      "mar_qc_sector", 0, ...
      "mar_replaced_runoff_count", 0, ...
      "mar_replaced_smb_count", 0, ...
      "merra_flux_orientation", "not_applicable", ...
      "modis_status", "not_requested", ...
      "modis_coverage_years", zeros(1, 0), ...
      "modis_missing_years", zeros(1, 0), ...
      "modis_finite_count", 0, ...
      "status", "skipped", ...
      "reason", "");
end

function summary = summarizeRecords(records)
   %SUMMARIZERECORDS Count statuses for compact reporting.
   statuses = string({records.status});
   kinds = unique(statuses, "stable");
   summary = struct();
   summary.total = numel(records);
   for k = 1:numel(kinds)
      key = matlab.lang.makeValidName(kinds(k));
      summary.(key) = sum(statuses == kinds(k));
   end
end

function summary = summarizeActions(records)
   %SUMMARIZEACTIONS Count each planned or completed repair action.
   summary = struct();
   for k = 1:numel(records)
      actions = string(records(k).actions);
      for n = 1:numel(actions)
         if actions(n) == ""
            continue
         end
         key = matlab.lang.makeValidName(actions(n));
         if ~isfield(summary, key)
            summary.(key) = 0;
         end
         summary.(key) = summary.(key) + 1;
      end
   end
end
