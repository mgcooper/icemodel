function manifest = makeFamilyManifest(dataset_family, source_doi, source_url, ...
      source_version, retrieval_date, cases, skipped)
   %MAKEFAMILYMANIFEST Build one verification family manifest struct.
   %
   %  manifest = icemodel.verification.setup.makeFamilyManifest( ...
   %     dataset_family, source_doi, source_url, source_version, ...
   %     retrieval_date, cases, skipped)
   %
   % Inputs
   %  dataset_family   Dataset family id, for example "esm_snowmip".
   %  source_doi       Source DOI string, or blank if none exists.
   %  source_url       Source URL used for provenance.
   %  source_version   Source version or bundle identifier.
   %  retrieval_date   Date the staged subset was generated or refreshed.
   %  cases            Struct array of case manifest entries.
   %  skipped          Struct array of requested cases that were not staged.
   %
   % Outputs
   %  manifest   Family-level manifest struct written to manifest.json.
   %
   % Role
   %  Setup helper used by dataset importers while creating or refreshing
   %  staged manifest files.

   % Keep family manifest fields in one canonical order for stable JSON output.
   names = icemodel.verification.setup.familyManifestFieldNames();
   values = {dataset_family, source_doi, source_url, source_version, ...
      retrieval_date, cases, skipped};

   % Match makeCaseManifestEntry: fail early if the schema and values diverge
   % instead of writing a malformed or partially shifted JSON manifest.
   if numel(values) ~= numel(names)
      error('family manifest expects %d values', numel(names))
   end

   manifest = cell2struct(values, names, 2);
end
