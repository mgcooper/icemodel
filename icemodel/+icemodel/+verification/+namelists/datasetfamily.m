function dataset_families = datasetfamily()
   %DATASETFAMILY Return the supported snow-verification dataset families.
   %
   %  dataset_families = icemodel.verification.namelists.datasetfamily()
   %
   % Outputs
   %  dataset_families   Supported manifest dataset_family values. These are
   %                     the staged source-family folder names under
   %                     demo/data/eval/snow and demo/data/eval/firn.
   %
   % Role
   %  Canonical source-family list shared by setup importers, validators, and
   %  normal verification workflow filters. This is not a case type.

   % Keep family ids aligned with committed folders under demo/data/eval.
   % snow/: esm_snowmip, laugh_tests. firn/: promice, sumup (the SUMup family
   % is data-gated; its importer stages into firn/sumup/ when the cache is
   % present, so the filter must accept it ahead of staged data).
   dataset_families = ["esm_snowmip"; "laugh_tests"; "promice"; "sumup"];
end
