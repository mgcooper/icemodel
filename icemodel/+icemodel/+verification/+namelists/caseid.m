function case_ids = caseid(dataset_family)
   %CASEID Return supported runnable snow-verification case ids.
   %
   %  case_ids = icemodel.verification.namelists.caseid()
   %  case_ids = icemodel.verification.namelists.caseid("esm_snowmip")
   %  case_ids = icemodel.verification.namelists.caseid("laugh_tests")
   %
   % Inputs
   %  dataset_family   Source dataset family selector or "all".
   %
   % Outputs
   %  case_ids         Runnable case ids for the selected source family.
   %
   % Role
   %  Canonical runnable-case list shared by setup importers, validators, and
   %  normal verification workflow filters. This selector is keyed by
   %  dataset_family, not by case_type.

   arguments
      dataset_family (1, 1) string = "all"
   end

   validateDatasetFamily(dataset_family)

   % Keep each dataset family explicit so adding cases remains
   % reviewable. The ESM-SnowMIP roster is sourced from the canonical
   % site catalog (icemodel.verification.namelists.snowmipsite) so
   % adding a new site only requires updating the catalog.
   switch dataset_family
      case "all"
         case_ids = [esmSnowmipCaseIds(); "colbeck1976"];

      case "esm_snowmip"
         case_ids = esmSnowmipCaseIds();

      case "laugh_tests"
         case_ids = "colbeck1976";
   end
end

function ids = esmSnowmipCaseIds()
   %ESMSNOWMIPCASEIDS Lift the ESM-SnowMIP site list from the catalog.
   info = icemodel.verification.namelists.snowmipsite();
   ids = reshape([info.sitename], [], 1);
end

function validateDatasetFamily(dataset_family)
   %VALIDATEDATASETFAMILY Validate the requested source-family selector.

   % The only supported selectors are "all" plus the dataset-family namelist.
   valid = ["all"; icemodel.verification.namelists.datasetfamily()];
   if ~ismember(dataset_family, valid)
      error('unsupported snow-verification dataset family: %s', dataset_family)
   end
end
