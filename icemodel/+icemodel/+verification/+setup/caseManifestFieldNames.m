function names = caseManifestFieldNames()
   %CASEMANIFESTFIELDNAMES Return canonical case-manifest fields.
   %
   %  names = icemodel.verification.setup.caseManifestFieldNames()
   %
   % Outputs
   %  names   String array in the exact field order written for each case.
   %
   % Role
   %  Setup helper used while constructing case entries. The operational
   %  manifest reader can inspect JSON fields directly, but importers need this
   %  single source of truth so ESM-SnowMIP and Laugh-Tests entries do not drift.

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
