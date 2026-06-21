function dataset_families = datasetfamily()
   %DATASETFAMILY Return the supported snow-verification dataset families.
   %
   %  dataset_families = icemodel.verification.namelists.datasetfamily()
   %
   % Outputs
   %  dataset_families   Supported manifest dataset_family values. These are
   %                     the staged source-family folder names directly under
   %                     demo/data/eval (family-flat taxonomy).
   %
   % Role
   %  Canonical source-family list shared by setup importers, validators, and
   %  normal verification workflow filters. This is not a case type.

   % Keep family ids aligned with committed folders under demo/data/eval.
   % Families are flat: esm_snowmip, laugh_tests, promice, sumup (the SUMup
   % family is data-gated; its importer stages into sumup/ when the cache is
   % present, so the filter must accept it ahead of staged data). The physical
   % regime (snow/firn) now lives in per-case surface_zone metadata, not the
   % directory layout.
   dataset_families = ["esm_snowmip"; "laugh_tests"; "promice"; "sumup"];
end
