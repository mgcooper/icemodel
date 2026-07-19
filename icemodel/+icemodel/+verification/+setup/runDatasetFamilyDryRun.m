function manifest = runDatasetFamilyDryRun(state, alive, kwargs)
   %RUNDATASETFAMILYDRYRUN Build a dry-run manifest through shared orchestration.
   %
   %  manifest = icemodel.verification.setup.runDatasetFamilyDryRun( ...
   %     state, alive, dataset_family=..., entry_callback=@makeEntry)
   %
   % Dry-run importers should return the same manifest shape as real imports
   % without writing to the staged eval/input tree. This helper routes them through
   % runDatasetFamilyImport using a temporary manifest file that is removed before
   % returning to the caller.

   arguments
      state (1, :) struct
      alive (1, :) logical
      kwargs.dataset_family (1, 1) string
      kwargs.requested_ids (1, :) string = strings(1, 0)
      kwargs.skipped (1, :) struct = struct('site', {}, 'reason', {})
      kwargs.source_doi (1, 1) string = ""
      kwargs.source_url (1, 1) string = ""
      kwargs.source_version (1, 1) string = ""
      kwargs.retrieval_date (1, 1) string = string(datetime('today'))
      kwargs.entry_callback = []
   end

   manifest_file = string(tempname) + ".json";
   cleanup = onCleanup(@() deleteIfFile(manifest_file));

   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family=kwargs.dataset_family, ...
      manifest_file=manifest_file, requested_ids=kwargs.requested_ids, ...
      skipped=kwargs.skipped, source_doi=kwargs.source_doi, ...
      source_url=kwargs.source_url, source_version=kwargs.source_version, ...
      retrieval_date=kwargs.retrieval_date, overwrite_family=true, ...
      entry_callback=kwargs.entry_callback);
end

function deleteIfFile(filename)
   %DELETEIFFILE Remove the temporary dry-run manifest file when present.
   if isfile(filename)
      delete(filename)
   end
end
