function manifest = readFamilyManifest(pathname)
   %READFAMILYMANIFEST Read one verification family manifest JSON file.
   %
   %  manifest = icemodel.verification.helpers.readFamilyManifest(pathname)
   %
   % Inputs
   %  pathname   Path to `eval/<family>/manifest.json`.
   %
   % Outputs
   %  manifest   Decoded family manifest with `manifest_path` and `family_root`
   %             fields added for downstream path resolution.
   %
   % Role
   %  Operational helper used while listing cases. It does not create or update
   %  staged setup artifacts.

   arguments
      pathname (1, 1) string
   end

   % Decode JSON first, then normalize text fields once at the manifest
   % boundary. Downstream workflow code should not need repeated string/char
   % conversions just because JSON decodes text as char vectors.
   manifest = jsondecode(fileread(pathname));
   manifest.dataset_family = string(manifest.dataset_family);
   manifest.source_doi = string(manifest.source_doi);
   manifest.source_url = string(manifest.source_url);
   manifest.source_version = string(manifest.source_version);
   manifest.retrieval_date = string(manifest.retrieval_date);
   manifest.cases = normalizeCaseEntries(manifest.cases);

   % Add local path context that should not be stored in the portable committed
   % manifest.
   manifest.manifest_path = pathname;
   manifest.family_root = string(fileparts(pathname));
end

function cases = normalizeCaseEntries(cases)
   %NORMALIZECASEENTRIES Convert JSON-decoded case text fields to strings.

   text_fields = ["case_id"; "case_type"; "site_id"; "site_name"; ...
      "surface_zone"; "evaluation_file"; "reference_file"; ...
      "native_timestep"; "notes"];
   for icase = 1:numel(cases)
      for jfield = 1:numel(text_fields)
         fieldname = text_fields(jfield);
         if isfield(cases(icase), fieldname)
            cases(icase).(fieldname) = string(cases(icase).(fieldname));
         end
      end

      % String-array fields. comparison_variables is shared by both schemas;
      % eval_target / forcing_sources / eval_sources are the metadata-only firn
      % descriptors. jsondecode renders a single-element array as a scalar char
      % and a multi-element array as a cellstr, so coerce both to a string array.
      array_fields = ["comparison_variables", "eval_target", ...
         "forcing_sources", "eval_sources"];
      for jfield = 1:numel(array_fields)
         fieldname = array_fields(jfield);
         if isfield(cases(icase), fieldname)
            cases(icase).(fieldname) = string(cases(icase).(fieldname));
         end
      end

      % Window fields. The snow schema names it comparison_window; the
      % metadata-only firn schema names it period. Both carry {start, end}.
      for window_field = ["comparison_window", "period"]
         if isfield(cases(icase), window_field)
            cases(icase).(window_field).start = ...
               string(cases(icase).(window_field).start);
            cases(icase).(window_field).end = ...
               string(cases(icase).(window_field).end);
         end
      end
   end
end
