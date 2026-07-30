function dataset_families = datasetfamily()
   %DATASETFAMILY Return the supported snow-verification dataset families.
   %
   %  dataset_families = icemodel.verification.namelists.datasetfamily()
   %
   % Outputs
   %  dataset_families   Supported manifest dataset_family values. These are
   %                     the staged source-family folder names directly under
   %                     the selected evaluation root (family-flat taxonomy).
   %
   % Role
   %  Canonical source-family list shared by setup importers, validators, and
   %  normal verification workflow filters. This is not a case type.

   % Keep family ids aligned with manifests under the selected evaluation root.
   % Families are flat: source families get their own folder and manifest, while
   % the physical regime (snow/firn) lives in per-case surface_zone metadata.
   dataset_families = [ ...
      "esm_snowmip"
      "laugh_tests"
      "promice"
      "sumup"
      "retmip"
      "imau"
      "research_site"
      "ktransect"];
end
