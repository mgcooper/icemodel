function manifest = makeFamilyManifest(dataset_family, source_doi, source_url, ...
      source_version, retrieval_date, cases, skipped, kwargs)
   %MAKEFAMILYMANIFEST Build one verification family manifest struct.
   %
   %  manifest = icemodel.verification.setup.makeFamilyManifest( ...
   %     dataset_family, source_doi, source_url, source_version, ...
   %     retrieval_date, cases, skipped)
   %  manifest = icemodel.verification.setup.makeFamilyManifest( ...
   %     ..., citation="Smeets et al. (2022) ...", license="CC-BY-4.0")
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
   % Name-value
   %  citation   Optional required-attribution citation for sources whose
   %             license mandates it (first used by ktransect / PANGAEA).
   %  license    Optional source license identifier, e.g. "CC-BY-4.0".
   %             Blank optional fields are omitted from the manifest so
   %             existing families keep their exact schema.
   %
   % Outputs
   %  manifest   Family-level manifest struct written to manifest.json.
   %
   % Role
   %  Setup helper used by dataset importers while creating or refreshing
   %  staged manifest files.

   arguments
      dataset_family
      source_doi
      source_url
      source_version
      retrieval_date
      cases
      skipped
      kwargs.citation (1, 1) string = ""
      kwargs.license (1, 1) string = ""
   end

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

   % Attach the optional attribution fields after the provenance block and
   % before cases/skipped so the JSON stays readable; blank values are omitted
   % entirely to preserve existing families' schemas.
   optional = ["citation", "license"];
   present = optional([kwargs.citation ~= "", kwargs.license ~= ""]);
   for name = present
      manifest.(char(name)) = kwargs.(char(name));
   end
   if ~isempty(present)
      manifest = orderfields(manifest, ...
         [names(1:5); cellstr(present(:)); names(6:7)]);
   end
end
