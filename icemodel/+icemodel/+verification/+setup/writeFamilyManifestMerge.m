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
   %    * Each NEW case in manifest.cases REPLACES the existing case with the
   %      same case_id, or is APPENDED when new. Existing cases whose case_id is
   %      not in the new set are PRESERVED unchanged, in their original order;
   %      newly added cases are appended after them.
   %    * Family-level fields the new manifest does not carry (e.g. a hand-added
   %      "schema" descriptor) are PRESERVED from the existing manifest, so a
   %      re-stage never silently drops them.
   %    * skipped[]: skip records for the REQUESTED ids are recomputed from the
   %      new manifest; skip records for OTHER ids are preserved. Re-staging a
   %      site that now succeeds clears its stale skip entry.
   %
   %  overwrite_family=true forces a full rewrite from manifest.cases alone
   %  (legacy whole-family behavior), discarding any prior cases and family
   %  fields. Use only to deliberately rebuild a family root from scratch.
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

   new_skipped = skipCell(manifestSkipped(manifest));

   % The requested-id set scopes which prior cases/skips this stage may touch.
   % It defaults to exactly what the new manifest carries, so an unspecified
   % re-stage is still scoped to its own cases.
   requested = kwargs.requested_ids;
   if isempty(requested)
      requested = unique([new_ids, skipIds(new_skipped)], 'stable');
   end

   % Full-rewrite escape hatch: ignore any prior manifest entirely.
   if kwargs.overwrite_family || ~isfile(manifest_file)
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
   pos = 0;
   for k = 1:numel(old_cases)
      id = old_ids(k);
      hit = find(new_ids == id, 1);
      pos = pos + 1;
      if ~isempty(hit)
         merged_cases{pos} = new_cases{hit};
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

   % skipped: drop prior skips for the requested ids (they were re-evaluated
   % this stage), keep all other prior skips, add the new skips. Mask the prior
   % skips to keep, then concatenate (no growing in the loop).
   keep_old = false(1, numel(old_skipped));
   for k = 1:numel(old_skipped)
      keep_old(k) = ~ismember(skipIds(old_skipped(k)), requested);
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
      error(['icemodel:verification:writeFamilyManifestMerge:fieldMismatch ' ...
         'case entries carry different field sets; cannot merge'])
   end
end
