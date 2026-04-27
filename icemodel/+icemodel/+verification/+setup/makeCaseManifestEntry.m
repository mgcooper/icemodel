function entry = makeCaseManifestEntry(values)
   %MAKECASEMANIFESTENTRY Build one case manifest entry from canonical fields.
   %
   %  entry = icemodel.verification.setup.makeCaseManifestEntry(values)
   %
   % Inputs
   %  values   Cell array matching the canonical case manifest field order.
   %
   % Outputs
   %  entry    Struct with the canonical case manifest schema.
   %
   % Role
   %  Setup helper used by dataset importers to make manifest schemas explicit
   %  and shared across verification families.

   % Get the canonical case manifest field names from the local schema helper.
   names = caseManifestFieldNames();

   % Fail early if an importer adds or removes a field without updating the
   % shared schema. This prevents silent per-family manifest drift.
   if numel(values) ~= numel(names)
      error('case manifest entry expects %d values', numel(names))
   end

   % Build the manifest entry in the canonical field order.
   entry = cell2struct(values(:), names, 1);
end

function names = caseManifestFieldNames()
   %CASEMANIFESTFIELDNAMES Return the canonical manifest-case field names.

   names = [ ...
      "case_id"
      "case_type"
      "tier"
      "site_id"
      "site_name"
      "forcing_file"
      "evaluation_file"
      "reference_file"
      "native_timestep"
      "comparison_window"
      "comparison_variables"
      "observation_variables"
      "notes"];
end
