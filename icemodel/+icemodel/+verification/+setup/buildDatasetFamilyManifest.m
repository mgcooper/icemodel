function manifest = buildDatasetFamilyManifest(state, alive, kwargs)
   %BUILDDATASETFAMILYMANIFEST Build and merge-write a family manifest.
   %
   %  manifest = icemodel.verification.setup.buildDatasetFamilyManifest( ...
   %     state, alive, dataset_family=..., entry_callback=@makeEntry)
   %
   % The callback receives one family-specific state record and returns one
   % schema-valid case entry, usually from makeFirnCaseManifestEntry. This helper
   % keeps manifest creation and writeFamilyManifestMerge semantics neutral so
   % importers do not duplicate manifest assembly. Case-entry callbacks should
   % derive forcing_sources/eval_sources with colocationSourceLists.
   %
   % See also: icemodel.verification.setup.runDatasetFamilyImport,
   %  icemodel.verification.setup.writeFamilyManifestMerge

   arguments
      state (1, :) struct
      alive (1, :) logical
      kwargs.dataset_family (1, 1) string
      kwargs.manifest_file (1, 1) string
      kwargs.requested_ids (1, :) string = strings(1, 0)
      kwargs.skipped (1, :) struct = struct('site', {}, 'reason', {})
      kwargs.source_doi (1, 1) string = ""
      kwargs.source_url (1, 1) string = ""
      kwargs.source_version (1, 1) string = ""
      kwargs.retrieval_date (1, 1) string = string(datetime('today'))
      kwargs.citation (1, 1) string = ""
      kwargs.license (1, 1) string = ""
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.entry_callback = []
   end

   if isempty(kwargs.entry_callback)
      error('icemodel:verification:buildDatasetFamilyManifest:missingCallback', ...
         'entry_callback is required')
   end

   case_entries = cell(1, nnz(alive));
   n_cases = 0;
   for n = 1:numel(state)
      if ~alive(n)
         continue
      end
      n_cases = n_cases + 1;
      case_entries{n_cases} = kwargs.entry_callback(state(n));
   end
   case_entries = case_entries(1:n_cases);

   if isempty(case_entries)
      cases = struct([]);
   else
      cases = vertcat(case_entries{:});
   end

   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      kwargs.dataset_family, kwargs.source_doi, kwargs.source_url, ...
      kwargs.source_version, kwargs.retrieval_date, cases, kwargs.skipped, ...
      citation=kwargs.citation, license=kwargs.license);

   manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
      kwargs.manifest_file, manifest, requested_ids=kwargs.requested_ids, ...
      overwrite_family=kwargs.overwrite_family);
end
