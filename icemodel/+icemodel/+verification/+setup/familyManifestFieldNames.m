function names = familyManifestFieldNames()
   %FAMILYMANIFESTFIELDNAMES Return canonical family-manifest fields.
   %
   %  names = icemodel.verification.setup.familyManifestFieldNames()
   %
   % Outputs
   %  names   Cell array in the exact field order written to manifest.json.
   %
   % Role
   %  Setup helper used while constructing family manifests. These names remain
   %  centralized because both supported importers write the same manifest
   %  schema and should fail together if that schema changes.

   names = { ...
      'dataset_family'
      'source_doi'
      'source_url'
      'source_version'
      'retrieval_date'
      'cases'};
end
