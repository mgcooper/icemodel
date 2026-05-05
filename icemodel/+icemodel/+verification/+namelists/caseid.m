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
   %  dataset_family, not by case_type. Family-specific case ids are sourced
   %  from each family's dedicated namelist (snowmipsite for "esm_snowmip",
   %  laughtests for "laugh_tests") so adding cases requires only updating
   %  the family namelist.

   arguments
      dataset_family (1, 1) string = "all"
   end

   validateDatasetFamily(dataset_family)

   % Dispatch by dataset family. Each family's case-id list is owned by its
   % namelist function; this dispatcher just concatenates them.
   switch dataset_family
      case "all"
         case_ids = [ ...
            icemodel.verification.namelists.snowmipsite(); ...
            icemodel.verification.namelists.laughtests()];

      case "esm_snowmip"
         case_ids = icemodel.verification.namelists.snowmipsite();

      case "laugh_tests"
         case_ids = icemodel.verification.namelists.laughtests();
   end
end

function validateDatasetFamily(dataset_family)
   %VALIDATEDATASETFAMILY Validate the requested source-family selector.

   % The only supported selectors are "all" plus the dataset-family namelist.
   valid = ["all"; icemodel.verification.namelists.datasetfamily()];
   if ~ismember(dataset_family, valid)
      error('unsupported snow-verification dataset family: %s', dataset_family)
   end
end
