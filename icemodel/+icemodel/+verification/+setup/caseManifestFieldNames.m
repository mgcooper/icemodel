function names = caseManifestFieldNames()
   %CASEMANIFESTFIELDNAMES Return canonical case-manifest fields.
   %
   %  names = icemodel.verification.setup.caseManifestFieldNames()
   %
   % Outputs
   %  names   String array in the exact field order written for each case.
   %
   % Role
   %  Setup helper for building case entries. The operational manifest reader
   %  can inspect the JSON fields directly. Every importer reads the field
   %  order from here, so ESM-SnowMIP and Laugh-Tests entries stay consistent.

   names = [ ...
      "case_id"
      "case_type"
      "site_id"
      "site_name"
      "surface_zone"
      "eval_target"
      "permafrost_zone"
      "evaluation_file"
      "reference_file"
      "native_timestep"
      "period"
      "comparison_variables"
      "observation_variables"
      "notes"];
end
