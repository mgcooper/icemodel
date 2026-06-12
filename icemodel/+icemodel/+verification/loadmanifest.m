function manifest = loadmanifest(case_id, kwargs)
   %LOADMANIFEST Return one resolved verification case manifest.
   %
   %  manifest = icemodel.verification.loadmanifest("cdp")
   %  manifest = icemodel.verification.loadmanifest("colbeck1976")
   %
   % Inputs
   %  case_id                    Case id to resolve from the staged manifests.
   %  evaluation_data_root       Base evaluation-data root. When blank, the
   %                             path is resolved from icemodel.config.
   %  icemodel_config_casename   Config casename used to resolve the default
   %                             evaluation-data root without mutating config.
   %
   % Outputs
   %  manifest   One resolved case-entry struct with family provenance and
   %             absolute artifact paths.
   %
   % Role
   %  This is a normal verification workflow entry point. It reads committed
   %  setup artifacts and does not mutate the staged data tree.

   arguments
      case_id (1, :) string
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
   end

   % Reuse listcases so manifest path resolution and filtering live in one
   % operational path.
   cases = icemodel.verification.listcases( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);

   % Give a path-aware error when no staged manifests are available.
   if isempty(cases)
      snow_data_root = icemodel.verification.helpers.snowDataRoot( ...
         "evaluation_data_root", kwargs.evaluation_data_root, ...
         "icemodel_config_casename", kwargs.icemodel_config_casename);
      error('no snow-verification cases found under %s', snow_data_root)
   end

   % Match case ids case-insensitively for interactive convenience while still
   % returning the canonical manifest case id.
   ids = [cases.case_id];
   idx = find(ids == case_id, 1);
   if isempty(idx)
      error('snow-verification case not found: %s', case_id)
   end
   manifest = cases(idx);
end
