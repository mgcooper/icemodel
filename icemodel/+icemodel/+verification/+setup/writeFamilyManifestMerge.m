function manifest = writeFamilyManifestMerge(manifest_file, manifest, kwargs)
   %WRITEFAMILYMANIFESTMERGE Merge new case entries into a family manifest.
   %
   %  manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
   %     manifest_file, manifest)
   %  manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
   %     manifest_file, manifest, overwrite_family=true)
   %  manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
   %     manifest_file, manifest, requested_ids=["kanl","kanm"])
   %
   %  Incremental staging primitive shared by the firn importers
   %  (importPromiceSites, importSumup). Staging one site must NOT churn or drop
   %  the other sites' committed case entries. This helper ADDS or UPDATES only
   %  the requested cases and PRESERVES every other existing case entry byte for
   %  byte, then writes the merged manifest.
   %
   %  Merge semantics (the DEFAULT)
   %    * Existing manifest is read from manifest_file when present (raw decode,
   %      so untouched cases re-encode identically - no field reordering).
   %    * Each NEW case in manifest.cases PATCHES the existing case with the
   %      same case_id, or is APPENDED when new. Omitted colocation legs remain
   %      unchanged because source selectors describe only the current call.
   %      A narrower ordinary patch cannot contract the durable case period or a
   %      retained artifact leg window; wider patches expand those intervals.
   %      Existing cases whose case_id is not in the new set are PRESERVED
   %      unchanged, in their original order.
   %      Stale duplicate existing case_ids are collapsed to one entry because
   %      a manifest case_id is the stable identity of one staged case folder.
   %      Newly added cases are appended after preserved existing cases.
   %    * Family-level fields the new manifest does not carry (e.g. a hand-added
   %      "schema" descriptor) are PRESERVED from the existing manifest, so a
   %      re-stage never silently drops them.
   %    * skipped[]: skip records for the REQUESTED ids are recomputed from the
   %      new manifest; skip records for OTHER ids are preserved. Re-staging a
   %      site that now succeeds clears its stale skip entry.
   %
   %  overwrite_family=true forces a full rewrite from manifest.cases alone
   %  (legacy whole-family behavior), discarding any prior cases and family
   %  fields. When that rewrite actually removes prior cases, coverage, sources,
   %  artifact references, skipped records, or extension fields, this helper
   %  emits an overwriteFamily warning. Use only to deliberately rebuild a
   %  family root from scratch.
   %
   %  Inputs
   %    manifest_file : string  destination eval/<family>/manifest.json path.
   %    manifest      : struct   family manifest for the REQUESTED cases only
   %                    (from makeFamilyManifest); manifest.skipped optional.
   %
   %  Name-value
   %    requested_ids   : string vector  case_ids this stage requested. Defaults
   %                      to the case_ids present in manifest.cases plus the
   %                      sites in manifest.skipped, so an idempotent re-stage of
   %                      the same site updates exactly that entry.
   %    overwrite_family: logical (default false)  force a full rewrite.
   %
   %  Returns
   %    manifest : struct  the merged manifest, also written to manifest_file.
   %
   % See also: icemodel.verification.setup.makeFamilyManifest,
   %  icemodel.verification.setup.writeManifest,
   %  icemodel.verification.setup.importPromiceSites,
   %  icemodel.verification.setup.importSumup

   arguments
      manifest_file (1, 1) string
      manifest (1, 1) struct
      kwargs.requested_ids (1, :) string = strings(1, 0)
      kwargs.overwrite_family (1, 1) logical = false
   end

   new_cases = caseCell(manifest.cases);
   new_ids = caseIds(new_cases);
   validateUniqueNewCaseIds(new_ids);

   new_skipped = skipCell(manifestSkipped(manifest));

   % The requested-id set scopes which prior cases/skips this stage may touch.
   % It defaults to exactly what the new manifest carries, so an unspecified
   % re-stage is still scoped to its own cases.
   requested = kwargs.requested_ids;
   if isempty(requested)
      requested = unique([new_ids, skipIds(new_skipped)], 'stable');
   end

   % Full-rewrite escape hatch: ignore any prior manifest entirely. Warn only
   % when replacement actually removes prior cases, sources, skipped records,
   % or family fields; repeated kill-safe persistence while building a family
   % is not destructive.
   if kwargs.overwrite_family || ~isfile(manifest_file)
      if kwargs.overwrite_family && isfile(manifest_file)
         existing = jsondecode(fileread(manifest_file));
         old_cases = caseCell(existingField(existing, 'cases'));
         old_ids = caseIds(old_cases);
         removes_state = ~isempty(setdiff(old_ids, new_ids));

         % A skipped record is durable family state too. Compare normalized
         % station identities so legacy underscore spelling does not create a
         % false removal warning.
         old_skip_ids = compactSkipKey(skipIds( ...
            skipCell(existingField(existing, 'skipped'))));
         new_skip_ids = compactSkipKey(skipIds(new_skipped));
         removes_state = removes_state ...
            || ~isempty(setdiff(old_skip_ids, new_skip_ids));

         % A retained case can still lose source state during replacement.
         % Compare both public source lists and concrete colocation legs so the
         % warning remains useful even for older manifests with incomplete lists.
         for k = 1:numel(old_cases)
            hit = find(new_ids == old_ids(k), 1);
            if isempty(hit)
               continue
            end
            old_case = old_cases{k};
            new_case = new_cases{hit};
            for source_field = ["forcing_sources", "eval_sources"]
               name = char(source_field);
               old_sources = strings(1, 0);
               new_sources = strings(1, 0);
               if isfield(old_case, name)
                  old_sources = reshape(string(old_case.(name)), 1, []);
               end
               if isfield(new_case, name)
                  new_sources = reshape(string(new_case.(name)), 1, []);
               end
               removes_state = removes_state ...
                  || ~isempty(setdiff(old_sources, new_sources));
            end
            old_legs = strings(1, 0);
            new_legs = strings(1, 0);
            if isfield(old_case, 'colocation') && isstruct(old_case.colocation)
               old_legs = reshape(string(fieldnames(old_case.colocation)), 1, []);
            end
            if isfield(new_case, 'colocation') && isstruct(new_case.colocation)
               new_legs = reshape(string(fieldnames(new_case.colocation)), 1, []);
            end
            removes_state = removes_state ...
               || ~isempty(setdiff(old_legs, new_legs));

            % Whole-case replacement can be destructive without changing a
            % source label or leg name. Detect narrowed periods and lost staged
            % artifact state on the case and on every retained colocation leg.
            removes_state = removes_state ...
               || replacementRemovesArtifactState(old_case, new_case);
            shared_legs = intersect(old_legs, new_legs, 'stable');
            for leg_name = shared_legs
               old_leg = old_case.colocation.(char(leg_name));
               new_leg = new_case.colocation.(char(leg_name));
               if isstruct(old_leg) && isstruct(new_leg)
                  removes_state = removes_state ...
                     || replacementRemovesArtifactState(old_leg, new_leg);
               end
            end
         end

         % overwrite_family also discards family-level extension fields.
         removes_state = removes_state || ~isempty(setdiff( ...
            fieldnames(existing), fieldnames(manifest)));
         if removes_state
            warning( ...
               'icemodel:verification:writeFamilyManifestMerge:overwriteFamily', ...
               ['Replacing the existing family manifest at %s; prior cases, ' ...
               'coverage, sources, artifacts, skipped records, or family ' ...
               'fields will be removed.'], ...
               manifest_file);
         end
      end
      merged = manifest;
      merged.cases = cell2caseArray(new_cases);
      merged.skipped = cell2skipArray(new_skipped);
      icemodel.verification.setup.writeManifest(manifest_file, merged);
      manifest = merged;
      return
   end

   % Merge into the existing manifest. Raw decode preserves every untouched
   % case struct exactly as written so its re-encoding is byte-stable.
   existing = jsondecode(fileread(manifest_file));
   old_cases = caseCell(existingField(existing, 'cases'));
   old_ids = caseIds(old_cases);
   old_skipped = skipCell(existingField(existing, 'skipped'));

   % Cases: keep prior cases not superseded by this stage (preserve order),
   % update in place where a requested id already exists, append the rest. The
   % merged length is the prior count plus the new ids not already present, so
   % preallocate to that size (no growing in the loop).
   appended = false(1, numel(new_cases));
   n_appended_new = sum(~ismember(new_ids, old_ids));
   merged_cases = cell(1, numel(old_cases) + n_appended_new);
   seen_old_ids = strings(1, numel(old_cases));
   n_seen_old = 0;
   pos = 0;
   for k = 1:numel(old_cases)
      id = old_ids(k);
      if ismember(id, seen_old_ids(1:n_seen_old))
         continue
      end
      n_seen_old = n_seen_old + 1;
      seen_old_ids(n_seen_old) = id;
      hit = find(new_ids == id, 1);
      pos = pos + 1;
      if ~isempty(hit)
         merged_cases{pos} = mergeCasePatch(old_cases{k}, new_cases{hit});
         appended(hit) = true;
      else
         merged_cases{pos} = old_cases{k};
      end
   end
   for k = 1:numel(new_cases)
      if ~appended(k)
         pos = pos + 1;
         merged_cases{pos} = new_cases{k};
      end
   end
   merged_cases = merged_cases(1:pos);

   % skipped: drop prior skips for the requested ids (they were re-evaluated
   % this stage), keep all other prior skips, add the new skips. Mask the prior
   % skips to keep, then concatenate (no growing in the loop).
   keep_old = false(1, numel(old_skipped));
   for k = 1:numel(old_skipped)
      keep_old(k) = ~isRequestedSkip(skipIds(old_skipped(k)), requested);
   end
   merged_skipped = [old_skipped(keep_old), new_skipped];

   % Family-level fields: take the new manifest's provenance values, but keep
   % any extra fields (e.g. "schema") the existing manifest carried and the new
   % one does not, so a re-stage never drops them.
   merged = manifest;
   extra = setdiff(fieldnames(existing), fieldnames(manifest), 'stable');
   for f = reshape(string(extra), 1, [])
      if f == "cases" || f == "skipped"
         continue
      end
      merged.(char(f)) = existing.(char(f));
   end

   merged.cases = cell2caseArray(merged_cases);
   merged.skipped = cell2skipArray(merged_skipped);

   icemodel.verification.setup.writeManifest(manifest_file, merged);
   manifest = merged;
end

%% Local helpers
function patched = mergeCasePatch(existing, incoming)
   %MERGECASEPATCH Preserve omitted legs and enclosing additive coverage.
   patched = incoming;
   case_identity_matches = concreteIdentityMatches(existing, incoming);
   if case_identity_matches
      patched = mergeWindowField(existing, patched, 'period');
   end
   if ~isfield(existing, 'colocation') || ~isstruct(existing.colocation) ...
         || ~isfield(incoming, 'colocation') || ~isstruct(incoming.colocation)
      return
   end

   % Incoming fields replace explicitly refreshed legs; every omitted leg stays
   % byte-equivalent even when the case-level point or period changed. Concrete
   % leg conflicts are destructive only for names explicitly present below.
   incoming_names = fieldnames(incoming.colocation);
   patched.colocation = icemodel.verification.setup.mergeColocation( ...
      existing.colocation, incoming.colocation);

   % A surgical refresh may update one retained source leg, but it does not
   % declare that leg's complete historical coverage. Preserve the enclosing
   % interval and the artifacts that support it. A partial extension keeps the
   % stable union of old and new references; an incoming enclosing/equal
   % rebuild replaces the old references. Concrete identity conflicts are a
   % replacement, never an invented union across products, points, or methods.
   shared_names = intersect(fieldnames(existing.colocation), ...
      incoming_names, 'stable');
   for k = 1:numel(shared_names)
      name = shared_names{k};
      old_leg = existing.colocation.(name);
      new_leg = patched.colocation.(name);
      if ~isstruct(old_leg) || ~isstruct(new_leg)
         continue
      end
      replace_prior_artifacts = false;
      if isfield(new_leg, 'replace_prior_artifacts')
         % This transient signal is set only after a requested refresh proves
         % prior files missing or concretely incompatible. Consume it here so
         % additive coverage merging cannot resurrect invalid references and
         % the implementation detail never enters the durable manifest.
         replace_prior_artifacts = ...
            scalarTrue(new_leg.replace_prior_artifacts);
         new_leg = rmfield(new_leg, 'replace_prior_artifacts');
         patched.colocation.(name) = new_leg;
      end
      if replace_prior_artifacts
         continue
      end
      if ~case_identity_matches || ~concreteIdentityMatches(old_leg, new_leg)
         continue
      end
      [new_leg, period_relation] = mergeWindowField( ...
         old_leg, new_leg, 'period');
      [new_leg, window_relation] = mergeWindowField( ...
         old_leg, new_leg, 'window');
      coverage_relation = preferredCoverageRelation( ...
         old_leg, new_leg, period_relation, window_relation);
      patched.colocation.(name) = mergeArtifactState( ...
         old_leg, new_leg, coverage_relation);
   end

   % Derive source labels from the final graph so retained/unioned legs cannot
   % drift from forcing_sources or eval_sources. Keep incoming extension labels
   % that are not represented by a standard colocation leg.
   if isfield(existing, 'forcing_sources') ...
         || isfield(incoming, 'forcing_sources') ...
         || isfield(existing, 'eval_sources') ...
         || isfield(incoming, 'eval_sources')
      [derived_forcing, derived_eval] = ...
         icemodel.verification.setup.colocationSourceLists( ...
         patched.colocation);
      incoming_forcing = sourceList(incoming, 'forcing_sources');
      incoming_eval = sourceList(incoming, 'eval_sources');
      patched.forcing_sources = cellstr(union( ...
         incoming_forcing, derived_forcing, 'stable'));
      patched.eval_sources = cellstr(union( ...
         incoming_eval, derived_eval, 'stable'));
   end
end

function [patched, relation] = mergeWindowField(existing, patched, fieldname)
   %MERGEWINDOWFIELD Keep the union and report which artifact state supports it.
   relation = "incoming";
   if ~isfield(existing, fieldname)
      return
   end
   if ~isfield(patched, fieldname) || ~isstruct(patched.(fieldname))
      patched.(fieldname) = existing.(fieldname);
      relation = "existing";
      return
   end

   old_window = existing.(fieldname);
   if ~isstruct(old_window) || ~all(isfield(old_window, ["start", "end"]))
      return
   end
   try
      [old_start, old_end, old_enabled] = icemodel.internal.pairedWindow( ...
         old_window.start, old_window.end);
   catch
      % Legacy non-datetime sentinels cannot define enclosing coverage.
      return
   end
   if ~old_enabled
      % Blank bounds mean all available data. A bounded surgical call cannot
      % narrow that durable declaration.
      patched.(fieldname) = old_window;
      relation = "existing";
      return
   end

   new_window = patched.(fieldname);
   try
      [new_start, new_end, new_enabled] = icemodel.internal.pairedWindow( ...
         new_window.start, new_window.end);
   catch
      % Do not let an ordinary malformed patch destroy valid durable coverage.
      patched.(fieldname) = old_window;
      relation = "existing";
      return
   end
   if ~new_enabled
      % An incoming all-available build encloses any prior bounded interval.
      relation = "incoming";
      return
   end

   % One scalar window cannot represent separated support. Preserve the durable
   % interval instead of inventing continuous coverage across an unproven gap;
   % callers can rebuild an enclosing artifact or use overwrite_family explicitly.
   if old_end < new_start || new_end < old_start
      patched.(fieldname) = old_window;
      relation = "existing";
      warning( ...
         'icemodel:verification:writeFamilyManifestMerge:disjointCoverage', ...
         ['Ordinary manifest patch has disjoint %s coverage; preserving the ' ...
         'existing interval and artifact state instead of spanning the gap.'], ...
         fieldname)
      return
   end

   % Select the original serialized endpoint that contributes the union so an
   % idempotent contained refresh retains the existing bytes and text format.
   old_encloses = old_start <= new_start && old_end >= new_end;
   incoming_encloses = new_start <= old_start && new_end >= old_end;
   if old_encloses && ~incoming_encloses
      relation = "existing";
   elseif incoming_encloses
      relation = "incoming";
   else
      relation = "union";
   end
   if old_start < new_start
      patched.(fieldname).start = old_window.start;
   end
   if old_end > new_end
      patched.(fieldname).end = old_window.end;
   end
end

function relation = preferredCoverageRelation(existing, incoming, ...
      period_relation, window_relation)
   %PREFERREDCOVERAGERELATION Use the most specific declared artifact window.
   if isfield(existing, 'window') || isfield(incoming, 'window')
      relation = window_relation;
   else
      relation = period_relation;
   end
end

function patched = mergeArtifactState(existing, patched, relation)
   %MERGEARTIFACTSTATE Keep references consistent with the merged coverage.
   incoming = patched;
   fields = artifactFieldNames();
   if relation == "existing"
      % The incoming refresh does not span the durable envelope. Remove its
      % narrower artifact claims, then restore the exact prior references.
      for field = fields
         name = char(field);
         if isfield(patched, name)
            patched = rmfield(patched, name);
         end
         if isfield(existing, name)
            patched.(name) = existing.(name);
         end
      end
      if isfield(patched, 'staged')
         patched = rmfield(patched, 'staged');
      end
      if isfield(existing, 'staged')
         patched.staged = existing.staged;
      end
      patched = mergeModelOutputState(existing, patched, incoming);
      return
   end
   if relation ~= "union"
      patched = mergeModelOutputState(existing, patched, incoming);
      return
   end

   % Partial/disjoint extensions require both artifact sets. Preserve the
   % original representation when only one value exists and use a cell list
   % only when the field genuinely needs multiple references.
   for field = fields
      name = char(field);
      old_refs = artifactReferences(existing, field);
      new_refs = artifactReferences(patched, field);
      refs = union(old_refs, new_refs, 'stable');
      if isempty(refs)
         continue
      elseif isempty(new_refs) && isfield(existing, name)
         patched.(name) = existing.(name);
      elseif numel(refs) > 1
         patched.(name) = cellstr(reshape(refs, [], 1));
      end
   end
   if (isfield(existing, 'staged') && scalarTrue(existing.staged)) ...
         || (isfield(patched, 'staged') && scalarTrue(patched.staged))
      patched.staged = true;
   end
   patched = mergeModelOutputState(existing, patched, incoming);
end

function patched = mergeModelOutputState(existing, patched, incoming)
   %MERGEMODELOUTPUTSTATE Merge a sidecar and its metadata as one artifact.
   if hasModelOutputReference(incoming)
      % A newly generated sidecar supersedes every field describing the old
      % artifact; mixing their metadata would make ownership ambiguous.
      source = incoming;
   elseif hasModelOutputReference(existing)
      % A forcing-only refresh cannot rebuild the optional sidecar, so keep
      % its prior file, status, variables, format, and provenance together.
      source = existing;
   else
      % With no durable sidecar to protect, retain the incoming diagnostic
      % state (for example, profile_not_available) exactly as reported.
      source = incoming;
   end

   % Replace the complete model-output field group so stale metadata cannot
   % survive when the chosen artifact changes.
   fields = string(fieldnames(patched));
   fields = fields(startsWith(fields, "model_output_"));
   if ~isempty(fields)
      patched = rmfield(patched, cellstr(fields));
   end
   fields = string(fieldnames(source));
   fields = fields(startsWith(fields, "model_output_"));
   for field = reshape(fields, 1, [])
      name = char(field);
      patched.(name) = source.(name);
   end
end

function tf = hasModelOutputReference(entry)
   %HASMODELOUTPUTREFERENCE True when a leg owns a concrete model sidecar.
   tf = ~isempty(artifactReferences(entry, "model_output_file")) ...
      || ~isempty(artifactReferences(entry, "model_output_files"));
end

function tf = concreteIdentityMatches(existing, incoming)
   %CONCRETEIDENTITYMATCHES Reject unions across conflicting known identity.
   tf = icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
      existing, incoming);
   if ~tf
      return
   end

   % Point identity appears either directly, as a numeric point, or inside a
   % location struct depending on the family/import generation.
   point_fields = ["point", "site_location", "location"];
   for field = point_fields
      if concretePointConflicts(existing, incoming, field)
         tf = false;
         return
      end
   end

   % Met and Data cadence are distinct artifact identities: production legs may
   % legitimately pair 15-minute model met with hourly userdata. Compare each
   % class independently before any coverage/reference union.
   tf = concreteArtifactCadenceMatches(existing, incoming);
end

function tf = concretePointConflicts(existing, incoming, field)
   %CONCRETEPOINTCONFLICTS Compare coordinates only when both points are finite.
   name = char(field);
   tf = false;
   if ~isfield(existing, name) || ~isfield(incoming, name)
      return
   end
   old = pointValue(existing.(name));
   new = pointValue(incoming.(name));
   if numel(old) == 2 && numel(new) == 2 ...
         && all(isfinite(old)) && all(isfinite(new))
      tf = any(abs(old - new) > 1e-8);
   end
end

function point = pointValue(value)
   %POINTVALUE Normalize a numeric point or location struct to [lat lon].
   point = [];
   if isnumeric(value) && numel(value) == 2
      point = reshape(double(value), 1, 2);
   elseif isstruct(value) ...
         && all(isfield(value, ["lat_wgs84", "lon_wgs84"]))
      point = [double(value.lat_wgs84), double(value.lon_wgs84)];
   end
end

function tf = concreteArtifactCadenceMatches(existing, incoming)
   %CONCRETEARTIFACTCADENCEMATCHES Compare known met/Data cadence independently.
   tf = true;
   classes = { ...
      ["met_file", "met_files"], ...
      ["data_file", "data_files"]};
   labels = ["met", "data"];
   for k = 1:numel(classes)
      [old_cadence, old_known, old_valid] = ...
         manifestArtifactCadence(existing, classes{k}, labels(k));
      [new_cadence, new_known, new_valid] = ...
         manifestArtifactCadence(incoming, classes{k}, labels(k));
      if ~old_valid || ~new_valid ...
            || (old_known && new_known && old_cadence ~= new_cadence)
         tf = false;
         return
      end
   end
end

function [cadence, known, valid] = manifestArtifactCadence(entry, fields, kind)
   %MANIFESTARTIFACTCADENCE Resolve one class from explicit or durable evidence.
   refs = artifactClassReferences(entry, fields);
   cadence = NaN;
   known = false;
   valid = true;
   if isempty(refs)
      return
   end

   % A generic explicit cadence is authoritative only when the leg contains one
   % artifact class. On mixed met/Data legs it cannot identify which class it
   % describes, so the durable per-class filenames remain authoritative.
   if kind == "met"
      other = artifactClassReferences(entry, ["data_file", "data_files"]);
   else
      other = artifactClassReferences(entry, ["met_file", "met_files"]);
   end
   if isempty(other)
      [cadence, known, valid] = explicitArtifactCadence(entry);
      if ~valid || known
         return
      end
   end

   % Filename inference is accepted only when every recorded reference follows
   % the durable class grammar and every recognized cadence agrees.
   values = NaN(numel(refs), 1);
   for n = 1:numel(refs)
      [values(n), parsed] = artifactFilenameCadence(refs(n), kind);
      if ~parsed
         return
      end
   end
   values = unique(values);
   valid = isscalar(values);
   if valid
      cadence = values(1);
      known = true;
   end
end

function refs = artifactClassReferences(entry, fields)
   %ARTIFACTCLASSREFERENCES Collect stable nonblank refs for one artifact class.
   refs = strings(0, 1);
   for field = fields
      refs = union(refs, artifactReferences(entry, field), 'stable');
   end
end

function [cadence, known, valid] = explicitArtifactCadence(entry)
   %EXPLICITARTIFACTCADENCE Read a generic cadence from a single-class leg.
   cadence = NaN;
   known = false;
   valid = true;
   values = NaN(2, 1);
   n_values = 0;
   if isfield(entry, 'artifact_cadence_seconds')
      [value, ok] = cadenceScalar(entry.artifact_cadence_seconds);
      if ~ok
         valid = false;
         return
      end
      n_values = n_values + 1;
      values(n_values) = value;
   end
   if isfield(entry, 'artifact_metadata') ...
         && isstruct(entry.artifact_metadata) ...
         && isfield(entry.artifact_metadata, 'artifact_cadence_seconds')
      [value, ok] = cadenceScalar( ...
         entry.artifact_metadata.artifact_cadence_seconds);
      if ~ok
         valid = false;
         return
      end
      n_values = n_values + 1;
      values(n_values) = value;
   end
   if n_values == 0
      return
   end
   values = unique(values(1:n_values));
   valid = isscalar(values);
   if valid
      cadence = values(1);
      known = true;
   end
end

function [cadence, valid] = cadenceScalar(raw)
   %CADENCESCALAR Normalize one positive explicit output cadence to seconds.
   cadence = NaN;
   if isduration(raw) && isscalar(raw)
      cadence = seconds(raw);
   elseif (isnumeric(raw) || islogical(raw)) && isscalar(raw)
      cadence = double(raw);
   elseif (ischar(raw) || isstring(raw)) && isscalar(string(raw))
      token = lower(strtrim(string(raw)));
      switch token
         case {"15m", "15 min", "15-minute"}
            cadence = 900;
         case {"30m", "30 min", "30-minute"}
            cadence = 1800;
         case {"1h", "1hr", "hourly"}
            cadence = 3600;
         otherwise
            cadence = str2double(token);
      end
   end
   valid = isscalar(cadence) && isfinite(cadence) && cadence > 0;
end

function [cadence, known] = artifactFilenameCadence(reference, kind)
   %ARTIFACTFILENAMECADENCE Parse one durable met or userdata filename grammar.
   cadence = NaN;
   known = false;
   [~, stem, ext] = fileparts(reference);
   name = string(stem) + string(ext);
   if kind == "met"
      token = regexp(name, ...
         "^met_.+_(?:\d{4}|\d{8}_\d{8})_(15m|30m|1hr)\.mat$", ...
         "tokens", "once");
      if isempty(token)
         return
      end
      switch string(token{1})
         case "15m"
            cadence = 900;
         case "30m"
            cadence = 1800;
         case "1hr"
            cadence = 3600;
      end
      known = true;
      return
   end

   % Hourly userdata intentionally retain the suffix-free legacy grammar;
   % native variants carry minute or raw-second suffixes.
   if ~isempty(regexp(name, ...
         "^.+_(?:\d{4}|\d{8}_\d{8})\.mat$", "once"))
      cadence = 3600;
      known = true;
      return
   end
   token = regexp(name, ...
      "^.+_(?:\d{4}|\d{8}_\d{8})_(15m|30m|(\d+)s)\.mat$", ...
      "tokens", "once");
   if isempty(token)
      return
   end
   switch string(token{1})
      case "15m"
         cadence = 900;
      case "30m"
         cadence = 1800;
      otherwise
         cadence = str2double(string(token{2}));
   end
   known = isfinite(cadence) && cadence > 0;
end

function values = sourceList(entry, field)
   %SOURCELIST Normalize one optional manifest source list to a column string.
   values = strings(0, 1);
   if isfield(entry, field)
      values = reshape(string(entry.(field)), [], 1);
   end
end

function tf = replacementRemovesArtifactState(existing, incoming)
   %REPLACEMENTREMOVESARTIFACTSTATE Detect narrowed coverage or lost file state.
   tf = windowNarrows(existingField(existing, 'period'), ...
      existingField(incoming, 'period')) ...
      || windowNarrows(existingField(existing, 'window'), ...
      existingField(incoming, 'window')) ...
      || stagedStateDrops(existing, incoming);
   if tf
      return
   end

   % Compare only staged-artifact references, not raw provenance source_file
   % fields whose replacement is expected during a source refresh.
   artifact_fields = artifactFieldNames();
   for field = artifact_fields
      old_refs = artifactReferences(existing, field);
      if isempty(old_refs)
         continue
      end
      new_refs = artifactReferences(incoming, field);
      if ~isempty(setdiff(old_refs, new_refs))
         tf = true;
         return
      end
   end
end

function fields = artifactFieldNames()
   %ARTIFACTFIELDNAMES Manifest fields that identify staged output artifacts.
   fields = ["evaluation_file", "reference_file", "obs_file", ...
      "observations_file", "forcing_file", "forcing_files", ...
      "met_file", "met_files", "data_file", "data_files", ...
      "profile_file", "profile_files", "surface_file", ...
      "model_output_file", "model_output_files"];
end

function tf = windowNarrows(existing_window, incoming_window)
   %WINDOWNARROWS True when a valid prior interval is no longer fully covered.
   tf = false;
   if ~isstruct(existing_window) ...
         || ~all(isfield(existing_window, ["start", "end"]))
      return
   end
   try
      [old_start, old_end, old_enabled] = icemodel.internal.pairedWindow( ...
         existing_window.start, existing_window.end);
   catch
      % Legacy non-datetime sentinel periods are not comparable coverage.
      return
   end
   if ~old_enabled
      if ~isstruct(incoming_window) ...
            || ~all(isfield(incoming_window, ["start", "end"]))
         tf = true;
         return
      end
      try
         [~, ~, new_enabled] = icemodel.internal.pairedWindow( ...
            incoming_window.start, incoming_window.end);
         tf = new_enabled;
      catch
         tf = true;
      end
      return
   end
   if ~isstruct(incoming_window) ...
         || ~all(isfield(incoming_window, ["start", "end"]))
      tf = true;
      return
   end
   try
      [new_start, new_end, new_enabled] = icemodel.internal.pairedWindow( ...
         incoming_window.start, incoming_window.end);
      tf = ~new_enabled || new_start > old_start || new_end < old_end;
   catch
      tf = true;
   end
end

function tf = stagedStateDrops(existing, incoming)
   %STAGEDSTATEDROPS True only for an explicit prior staged=true loss.
   tf = isfield(existing, 'staged') && scalarTrue(existing.staged) ...
      && (~isfield(incoming, 'staged') || ~scalarTrue(incoming.staged));
end

function tf = scalarTrue(value)
   %SCALARTRUE Normalize JSON/logical scalar truth without accepting vectors.
   tf = false;
   if islogical(value) || isnumeric(value)
      tf = isscalar(value) && isfinite(double(value)) && logical(value);
   elseif ischar(value) || isstring(value)
      text = lower(strtrim(string(value)));
      tf = isscalar(text) && text == "true";
   end
end

function refs = artifactReferences(entry, field)
   %ARTIFACTREFERENCES Normalize one artifact file field to unique nonblank refs.
   refs = strings(1, 0);
   name = char(field);
   if ~isfield(entry, name)
      return
   end
   try
      refs = reshape(strtrim(string(entry.(name))), 1, []);
   catch
      return
   end
   refs = unique(refs(strlength(refs) > 0), 'stable');
end

function c = caseCell(cases)
   %CASECELL Normalize a cases value (struct array / struct([]) / cell) to cell.
   if isempty(cases)
      c = {};
   elseif iscell(cases)
      c = reshape(cases, 1, []);
   else
      c = arrayfun(@(s) s, cases(:).', 'UniformOutput', false);
   end
end

function ids = caseIds(case_cell)
   %CASEIDS Extract case_id strings from a cell of case structs.
   ids = strings(1, numel(case_cell));
   for k = 1:numel(case_cell)
      ids(k) = string(case_cell{k}.case_id);
   end
end

function validateUniqueNewCaseIds(ids)
   %VALIDATEUNIQUENEWCASEIDS Require one new entry per staged case identity.
   duplicate = strings(1, numel(ids));
   seen = strings(1, numel(ids));
   n_duplicate = 0;
   n_seen = 0;
   for k = 1:numel(ids)
      if ismember(ids(k), seen(1:n_seen))
         n_duplicate = n_duplicate + 1;
         duplicate(n_duplicate) = ids(k);
      else
         n_seen = n_seen + 1;
         seen(n_seen) = ids(k);
      end
   end
   duplicate = unique(duplicate(1:n_duplicate), 'stable');
   if isempty(duplicate)
      return
   end
   error('icemodel:verification:writeFamilyManifestMerge:duplicateCaseId', ...
      'new manifest cases contain duplicate case_id values: %s', ...
      strjoin(string(duplicate), ', '))
end

function s = skipCell(skipped)
   %SKIPCELL Normalize a skipped value to a cell of scalar structs.
   if isempty(skipped) || (isstruct(skipped) && isempty(fieldnames(skipped)))
      s = {};
   elseif iscell(skipped)
      s = reshape(skipped, 1, []);
   else
      s = arrayfun(@(x) x, skipped(:).', 'UniformOutput', false);
   end
end

function ids = skipIds(skip_cell)
   %SKIPIDS Extract the site id from skip records (handles cell or scalar).
   if iscell(skip_cell)
      ids = strings(1, numel(skip_cell));
      for k = 1:numel(skip_cell)
         ids(k) = string(skip_cell{k}.site);
      end
   else
      ids = string(skip_cell.site);
   end
end

function tf = isRequestedSkip(skip_id, requested)
   %ISREQUESTEDSKIP Match skip ids using exact and compact station-style keys.
   tf = ismember(skip_id, requested) ...
      || ismember(compactSkipKey(skip_id), compactSkipKey(requested));
end

function key = compactSkipKey(ids)
   %COMPACTSKIPKEY Normalize legacy station labels used in skipped records.
   key = lower(erase(string(ids), "_"));
end

function skipped = manifestSkipped(manifest)
   %MANIFESTSKIPPED Return manifest.skipped or an empty skip struct.
   if isfield(manifest, 'skipped')
      skipped = manifest.skipped;
   else
      skipped = struct('site', {}, 'reason', {});
   end
end

function v = existingField(s, name)
   %EXISTINGFIELD Read a field from a decoded manifest, [] when absent.
   if isfield(s, name)
      v = s.(name);
   else
      v = [];
   end
end

function arr = cell2caseArray(case_cell)
   %CELL2CASEARRAY Concatenate scalar case structs into one struct array.
   %
   % JSON-decoded and freshly built case structs can carry their fields in a
   % different ORDER (decode preserves the on-disk order; builders use the
   % canonical schema order). Reorder every entry's fields to match the first
   % entry before concatenating so the array assembles without error and the
   % canonical order wins for the touched entries.
   if isempty(case_cell)
      arr = struct([]);
      return
   end
   ref = fieldnames(case_cell{1});
   for k = 2:numel(case_cell)
      case_cell{k} = orderfields(alignFields(case_cell{k}, ref), ref);
   end
   arr = [case_cell{:}];
end

function arr = cell2skipArray(skip_cell)
   %CELL2SKIPARRAY Concatenate scalar skip structs into one struct array.
   if isempty(skip_cell)
      arr = struct('site', {}, 'reason', {});
      return
   end
   for k = 1:numel(skip_cell)
      skip_cell{k} = struct('site', string(skip_cell{k}.site), ...
         'reason', string(skip_cell{k}.reason));
   end
   arr = [skip_cell{:}];
end

function s = alignFields(s, ref)
   %ALIGNFIELDS Ensure struct s carries exactly the ref field set.
   %
   % A mismatch in the field SET (not just order) between a preserved case and a
   % touched case would otherwise break concatenation; surface it instead of
   % silently fabricating fields.
   have = fieldnames(s);
   if ~isempty(setxor(have, ref))
      error('icemodel:verification:writeFamilyManifestMerge:fieldMismatch', ...
         'case entries carry different field sets; cannot merge')
   end
end
