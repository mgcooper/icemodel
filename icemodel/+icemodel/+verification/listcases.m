function cases = listcases(kwargs)
   %LISTCASES Enumerate staged verification cases from family manifests.
   %
   %  cases = icemodel.verification.listcases()
   %  cases = icemodel.verification.listcases(dataset_family="esm_snowmip")
   %
   % Inputs
   %  evaluation_data_root       Base evaluation-data root. When blank, the
   %                             path is resolved from icemodel.config.
   %  icemodel_config_casename   Config casename used to resolve the default
   %                             evaluation-data root without mutating config.
   %  dataset_family             Optional family filter, for example
   %                             "esm_snowmip" or "laugh_tests".
   %
   % Outputs
   %  cases   Struct array of resolved case manifest entries. File paths are
   %          absolute paths under the staged verification data tree.
   %
   % Role
   %  This is a normal verification workflow entry point. It reads committed
   %  manifests and does not create or refresh staged data.

   arguments
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.dataset_family (1, 1) string ...
         {icemodel.verification.validators.mustBeDatasetFamilyFilter} = ""
   end

   % Discover family manifests first; all later filters operate on manifest
   % contents so the function stays independent of any hard-coded case list.
   files = icemodel.verification.helpers.familyManifestFiles( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);

   % Read each family manifest, apply optional filters, and keep only families
   % that contribute at least one case.
   selected_family_cases = cell(numel(files), 1);
   n_families = 0;
   for i = 1:numel(files)
      family = icemodel.verification.helpers.readFamilyManifest(files(i));
      if skipFamily(family, kwargs.dataset_family)
         continue
      end

      group = resolveFamilyCases(family);
      if isempty(group)
         continue
      end

      n_families = n_families + 1;
      selected_family_cases{n_families} = group;
   end

   % Return an empty struct when no manifests match. This makes caller-side
   % missing-data handling explicit without requiring try/catch.
   if n_families == 0
      cases = struct([]);
      return
   end

   % Sort by case id so runner and test output are stable across filesystems.
   % All families share one canonical case-entry schema, so vertcat works
   % without harmonization.
   cases = vertcat(selected_family_cases{1:n_families});
   [~, order] = sort([cases.case_id]);
   cases = cases(order);
end

function tf = skipFamily(family, dataset_family)
   %SKIPFAMILY True when a family manifest does not match the requested filter.

   tf = ~isblanktext(dataset_family) && ...
      family.dataset_family ~= dataset_family;
end

function cases = resolveFamilyCases(family)
   %RESOLVEFAMILYCASES Resolve one family's case entries.

   entries = reshape(family.cases, [], 1);
   resolved_rows = repmat(resolveCase(entries(1), family), numel(entries), 1);
   for i = 1:numel(entries)
      resolved_rows(i) = resolveCase(entries(i), family);
   end
   cases = resolved_rows;
end

function resolved = resolveCase(entry, family)
   %RESOLVECASE Combine family metadata with one case entry.

   % Copy case-local metadata first, then add family-level provenance.
   resolved = entry;
   resolved.dataset_family = family.dataset_family;
   resolved.source_doi = family.source_doi;
   resolved.source_url = family.source_url;
   resolved.source_version = family.source_version;
   resolved.retrieval_date = family.retrieval_date;
   resolved.manifest_path = family.manifest_path;
   resolved.family_root = family.family_root;

   % Resolve relative artifact paths at read time so manifests stay portable
   % inside demo/data while workflow functions receive absolute paths.
   resolved.forcing_path = resolveCasePath(family.family_root, entry, ...
      'forcing_file');
   resolved.evaluation_path = resolveCasePath(family.family_root, entry, ...
      'evaluation_file');
   resolved.reference_path = resolveCasePath(family.family_root, entry, ...
      'reference_file');
end

function pathname = resolveCasePath(family_root, entry, fieldname)
   %RESOLVECASEPATH Resolve one optional relative manifest path.

   if isfield(entry, fieldname)
      pathname = fullfile(family_root, entry.(fieldname));
   else
      pathname = "";
   end
end
