function names = firnCaseManifestFieldNames()
   %FIRNCASEMANIFESTFIELDNAMES Return canonical firn case-manifest fields.
   %
   %  names = icemodel.verification.setup.firnCaseManifestFieldNames()
   %
   % Outputs
   %  names   String array in the exact field order written for each
   %          co-located firn evaluation case.
   %
   % Role
   %  Setup helper used while constructing firn case entries. This is the
   %  firn-side counterpart of caseManifestFieldNames; it is a separate
   %  schema because a firn case records a co-located multi-model forcing
   %  bundle, projected site coordinates, and per-model comparison windows
   %  that the snow case schema does not carry. Keeping them distinct means
   %  the snow manifests never have to absorb firn-only fields and the two
   %  schemas can evolve independently.
   %
   % See also: icemodel.verification.setup.caseManifestFieldNames,
   %  icemodel.verification.setup.makeFirnCaseManifestEntry,
   %  icemodel.verification.setup.importPromiceSites

   names = [ ...
      "case_id"
      "case_type"
      "site_id"
      "site_name"
      "site_location"
      "evaluation_file"
      "reference_file"
      "native_timestep"
      "comparison_window"
      "comparison_variables"
      "observation_variables"
      "colocated_forcing"
      "notes"];
end
