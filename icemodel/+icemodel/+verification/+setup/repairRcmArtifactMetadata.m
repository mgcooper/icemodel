function report = repairRcmArtifactMetadata(input_root, kwargs)
   %REPAIRRCMARTIFACTMETADATA Classify and repair current-token RCM artifacts.
   %
   %  report = icemodel.verification.setup.repairRcmArtifactMetadata()
   %  report = ... repairRcmArtifactMetadata(input_root, eval_root=..., ...
   %     dataset_family="promice", source_id="mar3.11", dry_run=false)
   %
   % Unscoped calls inventory the current-token RCM cache; family-scoped calls
   % inspect exact manifest references. Mapped artifacts can be repaired without
   % reopening every source archive. With no callback, this function synchronizes
   % writer and requested-site metadata. A repair_function can add or replace
   % explicitly declared timetable variables, UserData fields, CustomProperties,
   % or the table Description. This function supplies manifest discovery,
   % dry-run classification, preservation checks, atomic replacement, hashes,
   % and idempotence reporting.
   %
   % The callback signature is
   %
   %    [artifact, actions] = repair_function(artifact, context)
   %
   % where context identifies the file, payload variable, source product, case
   % alias, requested location, sample method, and merged artifact metadata.
   % actions is a string array used in the report. The callback must preserve the
   % timetable time axis and all undeclared variables and metadata; new variables
   % must be appended after retained columns. Canonical
   % builders/importers must also implement every lasting contract change; the
   % callback is the bounded migration path for artifacts already accepted.
   %
   % Input
   %  input_root  current RCM input tree (default: <repo>/data/input).
   %
   % Name-value options
   %  eval_root                        manifest tree (default: <repo>/data/eval)
   %  dataset_family                   exact family-manifest selection
   %  source_id                        exact source-product inventory filter
   %  repair_function                  callback described above
   %  allowed_variable_changes         timetable variables the callback may change
   %  allowed_metadata_changes         UserData fields the callback may change
   %  allowed_custom_property_changes  CustomProperties the callback may change
   %  allowed_table_property_changes   table properties the callback may change;
   %                                   currently only Description is supported
   %  dry_run                          classify without writing (default: true)
   %
   % Return
   %  report  transaction settings, per-file records, status/action summaries,
   %          hashes, and preservation/idempotence evidence.
   %
   % dry_run defaults true. Unsupported filenames, unmapped or ambiguous
   % manifests, malformed payloads, and callback contract violations are
   % classified without mutation. This function does not change manifests,
   % observations.mat, source-grid CustomProperties, or artifact cadence.
   %
   % See also: icemodel.verification.setup.refreshManifestSourceLists,
   %  icemodel.verification.setup.repairMetTimeSupport,
   %  icemodel.verification.setup.stageRcmForcing

   arguments
      input_root (1, 1) string = string(fullfile( ...
         icemodel.internal.fullpath("data"), "input"))
      kwargs.eval_root (1, 1) string = string(fullfile( ...
         icemodel.internal.fullpath("data"), "eval"))
      kwargs.dataset_family (1, :) string = strings(1, 0)
      kwargs.source_id (1, :) string = strings(1, 0)
      kwargs.repair_function = []
      kwargs.allowed_variable_changes (1, :) string = strings(1, 0)
      kwargs.allowed_metadata_changes (1, :) string = strings(1, 0)
      kwargs.allowed_custom_property_changes (1, :) string = strings(1, 0)
      kwargs.allowed_table_property_changes (1, :) string ...
         {mustBeMember(kwargs.allowed_table_property_changes, ...
         "Description")} = strings(1, 0)
      kwargs.dry_run (1, 1) logical = true
   end

   % A callback is optional, but when supplied its mutation boundary must be
   % explicit so an ad-hoc repair cannot silently broaden its own scope.
   if ~isempty(kwargs.repair_function) ...
         && ~isa(kwargs.repair_function, 'function_handle')
      error('icemodel:verification:repairRcmArtifactMetadata:badRepairFunction', ...
         'repair_function must be empty or a function handle')
   end
   if ~isempty(kwargs.repair_function) ...
         && isempty(kwargs.allowed_variable_changes) ...
         && isempty(kwargs.allowed_metadata_changes) ...
         && isempty(kwargs.allowed_custom_property_changes) ...
         && isempty(kwargs.allowed_table_property_changes)
      error('icemodel:verification:repairRcmArtifactMetadata:missingBoundary', ...
         ['a repair_function requires at least one declared variable or ' ...
         'metadata/property change'])
   end
   forbidden = intersect(kwargs.allowed_custom_property_changes, ...
      protectedCustomProperties(), "stable");
   if ~isempty(forbidden)
      error( ...
         'icemodel:verification:repairRcmArtifactMetadata:protectedProperty', ...
         'canonical/provenance CustomProperties are not callback-owned: %s', ...
         strjoin(forbidden, ', '))
   end

   % Current manifests are the authority for both artifact selection and the
   % requested-site metadata attached to reused source-grid payloads.
   families = unique(reshape(kwargs.dataset_family, 1, []), "stable");
   locations = manifestLocations(kwargs.eval_root, input_root, families);
   files = rcmArtifactFiles(input_root, locations, kwargs.source_id);
   records = repmat(emptyRecord(), numel(files), 1);
   for k = 1:numel(files)
      records(k) = repairOne(files(k), locations, kwargs);
   end

   % Keep transaction settings and compact status/action summaries with the
   % per-file evidence so dry run, write, and pass-two results can be archived.
   report = struct();
   report.input_root = input_root;
   report.eval_root = kwargs.eval_root;
   report.dataset_family = families;
   report.source_id = kwargs.source_id;
   report.repair_function = repairFunctionName(kwargs.repair_function);
   report.allowed_variable_changes = kwargs.allowed_variable_changes;
   report.allowed_metadata_changes = kwargs.allowed_metadata_changes;
   report.allowed_custom_property_changes = ...
      kwargs.allowed_custom_property_changes;
   report.allowed_table_property_changes = ...
      kwargs.allowed_table_property_changes;
   report.dry_run = kwargs.dry_run;
   report.records = records;
   report.summary = summarizeRecords(records);
   report.actions = summarizeActions(records);
end

function locations = manifestLocations(eval_root, input_root, families)
   %MANIFESTLOCATIONS Build current artifact and alias location lookups.
   manifests = dir(fullfile(eval_root, "*", "manifest.json"));
   locations = struct();
   locations.by_file = containers.Map('KeyType', 'char', 'ValueType', 'any');
   locations.by_alias = containers.Map('KeyType', 'char', 'ValueType', 'any');
   locations.by_alias_source = containers.Map( ...
      'KeyType', 'char', 'ValueType', 'any');
   locations.ambiguous_file = containers.Map('KeyType', 'char', ...
      'ValueType', 'logical');
   locations.ambiguous_alias = containers.Map('KeyType', 'char', ...
      'ValueType', 'logical');
   locations.ambiguous_alias_source = containers.Map( ...
      'KeyType', 'char', 'ValueType', 'logical');
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
   key = aliasKey(alias);
   entry = manifestLocationEntry(loc, "nearest", "");
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

function key = aliasKey(alias)
   %ALIASKEY Preserve the exact alias text in a collision-free map key.
   key = char(jsonencode(string(alias)));
end

function locations = addAliasSourceMethod(locations, alias, source_id, method)
   %ADDALIASSOURCEMETHOD Record source-specific alias sampling agreement.
   key = aliasSourceKey(alias, source_id);
   if ~isKey(locations.by_alias_source, key)
      locations.by_alias_source(key) = string(method);
      locations.ambiguous_alias_source(key) = false;
      return
   end
   if locations.by_alias_source(key) ~= string(method)
      locations.ambiguous_alias_source(key) = true;
   end
end

function key = aliasSourceKey(alias, source_id)
   %ALIASSOURCEKEY Build a collision-free map key for one alias and product.
   key = char(jsonencode([string(alias), string(source_id)]));
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
      method = "nearest";
      if isfield(leg, "sample_method") ...
            && strlength(string(leg.sample_method)) > 0
         method = string(leg.sample_method);
      end
      if isfield(leg, "source_id") ...
            && strlength(string(leg.source_id)) > 0
         locations = addAliasSourceMethod(locations, ...
            string(c.case_id), string(leg.source_id), method);
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
   source_id = "";
   if isfield(leg, "source_id")
      source_id = string(leg.source_id);
   end
   entry = manifestLocationEntry(loc, method, source_id);
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

function entry = manifestLocationEntry(loc, method, source_id)
   %MANIFESTLOCATIONENTRY Bundle requested point and sampling method metadata.
   entry = struct('location', loc, 'sample_method', string(method), ...
      'source_id', string(source_id));
end

function tf = sameManifestLocationEntry(a, b)
   %SAMEMANIFESTLOCATIONENTRY True when exact-file metadata agrees.
   tf = sameLocation(a.location, b.location) ...
      && string(a.sample_method) == string(b.sample_method) ...
      && string(a.source_id) == string(b.source_id);
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


function record = repairOne(filename, locations, kwargs)
   %REPAIRONE Classify and optionally repair one staged MAT file.
   record = parseArtifactFilename(filename);
   if record.status ~= "parsed"
      return
   end
   record.hash_before = ...
      icemodel.verification.setup.fileSha256(filename);

   % Exact manifest references win; alias fallback is allowed only when every
   % current family using the alias agrees on the requested coordinates.
    [found, location, sample_method, ambiguous] = artifactLocation( ...
       filename, record.alias, record.source_id, locations);
   if ambiguous
      record.status = "ambiguous";
      record.reason = "current manifests disagree for this artifact";
      return
   end
    if ~found
       record.status = "unmapped";
       record.reason = "artifact alias/source is not present in current manifests";
       return
   end

   try
      inventory = whos('-file', filename);
      names = string({inventory.name});
      payload = intersect(["met", "Data"], names, "stable");
      if numel(payload) ~= 1
         record.status = "error";
         record.reason = "artifact must contain exactly one timetable payload";
         return
      end

      % Load only the payload and optional metadata sidecar. Atomic writes copy
      % the original MAT before replacing these variables, so unrelated top-level
      % data never enter memory or get reserialized.
      selected = payload;
      if ismember("artifact_metadata", names)
         selected(end + 1) = "artifact_metadata";
      end
      contents = load(filename, selected{:}, '-mat');
      if ~istimetable(contents.(payload))
         record.status = "error";
         record.reason = "artifact must contain exactly one timetable payload";
         return
      end

      % Merge the two established metadata locations before invoking an optional
      % repair so the callback sees the complete accepted-artifact context.
      variable = payload(1);
      record.variable = variable;
      before = contents.(variable);
      record.has_artifact_metadata = ismember("artifact_metadata", names);
      existing_metadata = struct();
      if record.has_artifact_metadata && ~isstruct(contents.artifact_metadata)
         record.status = "error";
         record.reason = "top-level artifact_metadata must be a struct";
         return
      elseif record.has_artifact_metadata
         existing_metadata = contents.artifact_metadata;
      end
      userdata_metadata = metadataStruct(before.Properties.UserData);
      userdata_method = "";
      artifact_method = "";
      if isfield(userdata_metadata, "sample_method")
         userdata_method = string(userdata_metadata.sample_method);
      end
      if isfield(existing_metadata, "sample_method")
         artifact_method = string(existing_metadata.sample_method);
      end
      if strlength(userdata_method) > 0 && strlength(artifact_method) > 0 ...
            && userdata_method ~= artifact_method
         record.status = "restage_required";
         record.reason = ...
            "artifact metadata copies disagree on sample_method";
         return
      end
      metadata_before = mergeArtifactMetadata( ...
         userdata_metadata, existing_metadata);
      if isfield(metadata_before, "sample_method") ...
            && strlength(string(metadata_before.sample_method)) > 0 ...
            && string(metadata_before.sample_method) ~= sample_method
         record.status = "restage_required";
         record.reason = ...
            "artifact sample_method conflicts with the current manifest";
         return
      end
      if variable == "met" && isfield(metadata_before, ...
            "met_resample_policy") ...
            && string(metadata_before.met_resample_policy) ...
            == "linear_adjacent_finite_only"
         record.status = "repair_required";
         record.reason = ...
            "legacy linear met support requires the explicit " + ...
            "repairMetTimeSupport path before metadata repair";
         return
      end
      canonical_before = icemodel.forcing.helpers.stampMetadata( ...
         before, strict=false);
      context = repairContext(filename, record, variable, location, ...
         sample_method, metadata_before);

      % The callback owns only its declared variables/UserData fields. The
      % coordinator verifies that time, source-grid properties, and every
      % undeclared payload remain unchanged.
      repaired = before;
      actions = strings(1, 0);
      if ~isempty(kwargs.repair_function)
         [repaired, callback_actions] = ...
            kwargs.repair_function(repaired, context);
         actions = normalizeActions(callback_actions);
         validateRepairResult(before, repaired, ...
            kwargs.allowed_variable_changes, ...
            kwargs.allowed_metadata_changes, ...
            kwargs.allowed_custom_property_changes, ...
            kwargs.allowed_table_property_changes);
         if isempty(actions) && ~isequaln(before, repaired)
            actions = "custom_repair";
         end
      end

      % Shared metadata stamping is the durable default repair. It can add newly
      % canonical variable metadata without duplicating source-specific builders.
      transformed = repaired;
      repaired = icemodel.forcing.helpers.stampMetadata( ...
         repaired, strict=false);
      if ~isempty(kwargs.repair_function) ...
            && ~preservesCommonVariableMetadata( ...
            canonical_before, repaired)
         record.status = "error";
         record.reason = ...
            "callback changed canonical metadata for a retained variable";
         return
      end
      metadata = mergeArtifactMetadata( ...
         repaired.Properties.UserData, existing_metadata);
      metadata = artifactMetadata( ...
         repaired, location, sample_method, metadata);
      repaired.Properties.UserData = metadata;
      if ~isequaln(transformed, repaired)
         actions(end + 1) = "restamp_metadata";
      end

      top_level_changed = ~record.has_artifact_metadata ...
         || ~isequaln(existing_metadata, metadata);
      if top_level_changed
         actions(end + 1) = "sync_artifact_metadata";
      end
      record.actions = cellstr(unique(actions, "stable"));
      record.changed_variables = cellstr(changedVariables(before, repaired));
      record.changed_metadata_fields = cellstr( ...
         changedMetadataFields(metadata_before, metadata));
      record.unrelated_payload_preserved = preservesUnrelatedPayload( ...
         before, repaired, kwargs.allowed_variable_changes);
      if ~record.unrelated_payload_preserved
         record.status = "error";
         record.reason = "undeclared payload changed during artifact repair";
         return
      end

      changed = ~isequaln(before, repaired) || top_level_changed;
      if ~changed
         record.status = "unchanged";
         record.reason = "";
         record.hash_after = record.hash_before;
      elseif kwargs.dry_run
         record.status = "would_repair";
         record.reason = "";
      else
         saveRepairedArtifact(filename, variable, repaired, metadata);
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

function context = repairContext(filename, record, variable, location, ...
      sample_method, metadata)
   %REPAIRCONTEXT Describe one bounded repair target to a custom callback.
   context = struct();
   context.filename = filename;
   context.payload_variable = variable;
   context.alias = record.alias;
   context.source_id = record.source_id;
   context.window_start = record.window_start;
   context.window_end = record.window_end;
   context.requested_location = location;
   context.sample_method = sample_method;
   context.artifact_metadata = metadata;
end

function name = repairFunctionName(repair_function)
   %REPAIRFUNCTIONNAME Return a stable readable callback name for the report.
   name = "";
   if ~isempty(repair_function)
      name = string(func2str(repair_function));
   end
end

function actions = normalizeActions(actions)
   %NORMALIZEACTIONS Make callback action labels compact and deterministic.
   actions = reshape(string(actions), 1, []);
   actions = actions(strlength(actions) > 0);
   actions = unique(actions, "stable");
end

function validateRepairResult(before, after, allowed_variables, ...
      allowed_metadata, allowed_custom_properties, allowed_table_properties)
   %VALIDATEREPAIRRESULT Enforce the callback's declared mutation boundary.
   if ~istimetable(after)
      error('icemodel:verification:repairRcmArtifactMetadata:badRepairResult', ...
         'repair_function must return a timetable')
   end
   if ~isequaln(before.Time, after.Time)
      error('icemodel:verification:repairRcmArtifactMetadata:timeChanged', ...
         'repair_function must preserve the timetable time axis')
   end

   changed_variables = changedVariables(before, after);
   if ~withinBoundary(changed_variables, allowed_variables)
      error('icemodel:verification:repairRcmArtifactMetadata:variableBoundary', ...
         'repair_function changed undeclared variables: %s', ...
         strjoin(setdiff(changed_variables, allowed_variables), ', '))
   end
   before_metadata = metadataStruct(before.Properties.UserData);
   after_metadata = metadataStruct(after.Properties.UserData);
   changed_metadata = changedMetadataFields(before_metadata, after_metadata);
   if ~withinBoundary(changed_metadata, allowed_metadata)
      error('icemodel:verification:repairRcmArtifactMetadata:metadataBoundary', ...
         'repair_function changed undeclared UserData fields: %s', ...
         strjoin(setdiff(changed_metadata, allowed_metadata), ', '))
   end
   if ~preservesUnrelatedTableProperties( ...
         before, after, allowed_custom_properties, ...
         allowed_table_properties)
      error('icemodel:verification:repairRcmArtifactMetadata:propertyBoundary', ...
         ['repair_function changed source-grid or undeclared variable ' ...
         'properties'])
   end
end

function tf = withinBoundary(changed, allowed)
   %WITHINBOUNDARY Test whether all observed changes were explicitly declared.
   tf = ismember("*", allowed) || all(ismember(changed, allowed));
end

function metadata = metadataStruct(value)
   %METADATASTRUCT Normalize empty UserData while rejecting opaque metadata.
   if isempty(value)
      metadata = struct();
   elseif isstruct(value)
      metadata = value;
   else
      error('icemodel:verification:repairRcmArtifactMetadata:badMetadata', ...
         'artifact timetable UserData must be empty or a struct')
   end
end

function tf = preservesUnrelatedPayload(before, after, allowed)
   %PRESERVESUNRELATEDPAYLOAD Prove undeclared values and time are unchanged.
   if ~isequaln(before.Time, after.Time)
      tf = false;
      return
   end
   before_names = string(before.Properties.VariableNames);
   after_names = string(after.Properties.VariableNames);
   unrelated = setdiff(union(before_names, after_names, "stable"), allowed, ...
      "stable");
   if ismember("*", allowed)
      unrelated = strings(1, 0);
   end
   tf = all(ismember(unrelated, before_names)) ...
      && all(ismember(unrelated, after_names));
   for name = reshape(unrelated, 1, [])
      tf = tf && isequaln(before.(name), after.(name));
   end
end

function tf = preservesUnrelatedTableProperties(before, after, ...
      allowed_custom_properties, allowed_table_properties)
   %PRESERVESUNRELATEDTABLEPROPERTIES Protect source-grid and variable metadata.
   if ~isequaln(before.Properties.DimensionNames, ...
         after.Properties.DimensionNames)
      tf = false;
      return
   end
   if ~any(allowed_table_properties == "Description") ...
         && ~isequaln(before.Properties.Description, ...
         after.Properties.Description)
      tf = false;
      return
   end
   changed_custom_properties = changedMetadataFields( ...
      customProperties(before), customProperties(after));
   if ~isempty(intersect(changed_custom_properties, ...
         protectedCustomProperties()))
      tf = false;
      return
   end
   if ~withinBoundary( ...
         changed_custom_properties, allowed_custom_properties)
      tf = false;
      return
   end
   % Variable-associated properties for undeclared columns must remain exact.
   before_names = string(before.Properties.VariableNames);
   after_names = string(after.Properties.VariableNames);
   before_retained = before_names(ismember(before_names, after_names));
   after_existing = after_names(ismember(after_names, before_names));
   if ~isequal(before_retained, after_existing) ...
         || ~isequal(after_names(1:numel(before_retained)), before_retained)
      tf = false;
      return
   end
   names = intersect(before_names, after_names, "stable");
   tf = true;
   properties_to_check = [ ...
      "VariableUnits", "VariableDescriptions", "VariableContinuity"];
   for property = properties_to_check
      if ~isprop(before.Properties, property) ...
            || ~isprop(after.Properties, property)
         continue
      end
      before_values = before.Properties.(property);
      after_values = after.Properties.(property);
      for name = reshape(names, 1, [])
         before_index = before_names == name;
         after_index = after_names == name;
         before_value = variablePropertyValue( ...
            before_values, before_index, numel(before_names));
         after_value = variablePropertyValue( ...
            after_values, after_index, numel(after_names));
         tf = tf && isequaln(before_value, after_value);
      end
   end
end

function value = variablePropertyValue(values, index, variable_count)
   %VARIABLEPROPERTYVALUE Read optional per-variable metadata without over-indexing.
   value = [];
   if isempty(values)
      return
   end
   if numel(values) ~= variable_count
      error('icemodel:verification:repairRcmArtifactMetadata:badProperties', ...
         'timetable variable property width does not match its variables')
   end
   value = values(index);
end

function values = customProperties(T)
   %CUSTOMPROPERTIES Convert dynamic source-grid properties to comparable data.
   values = struct();
   names = sort(string(properties(T.Properties.CustomProperties)));
   for name = reshape(names, 1, [])
      values.(char(name)) = T.Properties.CustomProperties.(name);
   end
end

function names = protectedCustomProperties()
   %PROTECTEDCUSTOMPROPERTIES Metadata owned by extraction or canonical stamping.
   names = [ ...
      "X", "Y", "Lat", "Lon", "Elev", "Slope", "ScalarUnits", ...
      "StandardNames"];
end

function tf = preservesCommonVariableMetadata(expected, actual)
   %PRESERVESCOMMONVARIABLEMETADATA Keep retained fields on canonical metadata.
   expected_names = string(expected.Properties.VariableNames);
   actual_names = string(actual.Properties.VariableNames);
   names = intersect(expected_names, actual_names, "stable");
   tf = true;
   properties_to_check = [ ...
      "VariableUnits", "VariableDescriptions", "VariableContinuity"];
   for property = properties_to_check
      if ~isprop(expected.Properties, property) ...
            || ~isprop(actual.Properties, property)
         continue
      end
      expected_values = expected.Properties.(property);
      actual_values = actual.Properties.(property);
      for name = reshape(names, 1, [])
         expected_value = variablePropertyValue(expected_values, ...
            expected_names == name, numel(expected_names));
         actual_value = variablePropertyValue(actual_values, ...
            actual_names == name, numel(actual_names));
         tf = tf && isequaln(expected_value, actual_value);
      end
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

function saveRepairedArtifact(filename, variable, repaired, metadata)
   %SAVEREPAIREDARTIFACT Atomically replace one MAT while preserving top-level data.
   [folder, stem, suffix] = fileparts(filename);
   temp_file = string(fullfile(folder, ...
      "." + string(stem) + ".repair-" + string(char(java.util.UUID.randomUUID)) ...
      + string(suffix)));
   cleanup = onCleanup(@() deleteIfPresent(temp_file));
   [ok, message] = copyfile(filename, temp_file, 'f');
   if ~ok
      error('icemodel:verification:repairRcmArtifactMetadata:copyFailed', ...
         'failed to copy %s for atomic repair: %s', filename, message)
   end
   updates = struct();
   updates.(variable) = repaired;
   updates.artifact_metadata = metadata;
   save(temp_file, '-struct', 'updates', '-append')
   [ok, message] = movefile(temp_file, filename, 'f');
   if ~ok
      error('icemodel:verification:repairRcmArtifactMetadata:replaceFailed', ...
         'failed to replace %s: %s', filename, message)
   end
   clear cleanup
end

function deleteIfPresent(filename)
   %DELETEIFPRESENT Remove an abandoned sibling temporary file.
   if isfile(filename)
      delete(filename)
   end
end

function metadata = mergeArtifactMetadata(userdata, artifact_metadata)
   %MERGEARTIFACTMETADATA Preserve top-level fields absent from UserData.
   metadata = metadataStruct(userdata);
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

function [found, location, sample_method, ambiguous] = artifactLocation( ...
      filename, alias, source_id, locations)
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
       if entry.source_id == "" || entry.source_id ~= string(source_id)
          ambiguous = true;
          return
       end
       location = entry.location;
      sample_method = entry.sample_method;
      found = true;
      return
   end

    alias_key = aliasKey(alias);
   if ~isKey(locations.by_alias, alias_key)
      return
   end
   ambiguous = locations.ambiguous_alias(alias_key);
   if ambiguous
      return
   end
    entry = locations.by_alias(alias_key);
    location = entry.location;
    method_key = aliasSourceKey(alias, source_id);
    if ~isKey(locations.by_alias_source, method_key)
       return
    end
    ambiguous = locations.ambiguous_alias_source(method_key);
    if ambiguous
       return
    end
    sample_method = locations.by_alias_source(method_key);
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
      expr = "^(?:met_)?(.+?)_(mar3\.11|merra2|racmo2\.3p3)_(\d{4})(?:_(?:1hr|15m))?$";
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
   %EMPTYRECORD Define deterministic per-file repair evidence.
   record = struct( ...
      'filename', "", ...
      'status', "", ...
      'reason', "", ...
      'alias', "", ...
      'source_id', "", ...
      'window_start', "", ...
      'window_end', "", ...
      'variable', "", ...
      'has_artifact_metadata', false, ...
      'hash_before', "", ...
      'hash_after', "", ...
      'actions', {cell(0, 1)}, ...
      'changed_variables', {cell(0, 1)}, ...
      'changed_metadata_fields', {cell(0, 1)}, ...
      'unrelated_payload_preserved', false);
end

function summary = summarizeRecords(records)
   %SUMMARIZERECORDS Count statuses for compact reporting.
   statuses = string({records.status});
   kinds = unique(statuses, "stable");
   summary = struct( ...
      'total', numel(records), ...
      'unchanged', 0, ...
       'would_repair', 0, ...
       'repaired', 0, ...
       'repair_required', 0, ...
       'restage_required', 0, ...
       'unmapped', 0, ...
      'ambiguous', 0, ...
      'error', 0, ...
      'skipped', 0);
   for k = 1:numel(kinds)
      if strlength(kinds(k)) == 0
         continue
      end
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
