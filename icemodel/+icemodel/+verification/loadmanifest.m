function manifest = loadmanifest(case_id, kwargs)
   %LOADMANIFEST Return one resolved verification case manifest.
   %
   %  manifest = icemodel.verification.loadmanifest("cdp")
   %  manifest = icemodel.verification.loadmanifest("colbeck1976")
   %  manifest = icemodel.verification.loadmanifest("kanl", ...
   %     dataset_family="sumup")
   %
   % Inputs
   %  case_id                    Case id to resolve from the staged manifests.
%  evaluation_data_root       Base evaluation-data root. When blank, the
%                             repo-local data/eval tree is used.
%  input_data_root            Optional paired input-data root. When blank,
%                             resolved cases infer input/ beside eval/.
%  icemodel_config_casename   Config casename used to resolve the default
%                             evaluation-data root without mutating config.
   %  dataset_family             Optional family filter to disambiguate case
   %                             ids shared across families. The firn families
   %                             promice and sumup both publish kanl/kanm/kanu
   %                             (distinguished only by family folder), so a
   %                             bare loadmanifest("kanl") returns the first
   %                             match; pass dataset_family to select one.
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
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
      kwargs.dataset_family (1, 1) string = ""
   end

   % Reuse listcases so manifest path resolution and filtering live in one
   % operational path. An optional dataset_family narrows the search to one
   % family so the shared firn case ids (kanl/kanm/kanu) resolve unambiguously.
   cases = icemodel.verification.listcases( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "input_data_root", kwargs.input_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename, ...
      "dataset_family", kwargs.dataset_family);

   % Give a path-aware error when no staged manifests are available.
   if isempty(cases)
      evaluation_data_root = icemodel.verification.helpers.evaluationDataRoot( ...
         "evaluation_data_root", kwargs.evaluation_data_root, ...
         "icemodel_config_casename", kwargs.icemodel_config_casename);
      error('no verification cases found under %s', evaluation_data_root)
   end

   % Match the requested case id against the canonical manifest case ids and
   % return the first match (ids are the compact lowercase aliases the
   % importers write, e.g. kanl/kanm).
   ids = [cases.case_id];
   idx = find(ids == case_id, 1);
   if isempty(idx)
      error('snow-verification case not found: %s', case_id)
   end
   manifest = cases(idx);
end
