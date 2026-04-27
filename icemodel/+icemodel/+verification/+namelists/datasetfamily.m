function dataset_families = datasetfamily()
   %DATASETFAMILY Return the supported snow-verification dataset families.
   %
   %  dataset_families = icemodel.verification.namelists.datasetfamily()
   %
   % Outputs
   %  dataset_families   Supported manifest dataset_family values. These are
   %                     the staged source-family folder names under
   %                     demo/data/eval/snow.
   %
   % Role
   %  Canonical source-family list shared by setup importers, validators, and
   %  normal verification workflow filters. This is not a case type.

   % Keep family ids aligned with committed folders under demo/data/eval/snow.
   dataset_families = ["esm_snowmip"; "laugh_tests"];
end
