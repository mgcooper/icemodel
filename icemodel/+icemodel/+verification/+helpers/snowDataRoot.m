function root = snowDataRoot(kwargs)
   %SNOWDATAROOT Return the staged snow-verification data root.
   %
   %  root = icemodel.verification.helpers.snowDataRoot()
   %  root = icemodel.verification.helpers.snowDataRoot(icemodel_config_casename="test")
   % 
   % Inputs
   %  evaluation_data_root       Base evaluation-data root. When blank, the
   %                             path is resolved from icemodel.config.
   %  icemodel_config_casename   Config casename used to resolve the default
   %                             evaluation-data root without mutating config.
   %
   % Outputs
   %  root   Full path to the snow-verification data root, normally
   %         `<repo>/demo/data/eval/snow` for test/demo workflows.
   %
   % Role
   %  Operational path helper shared by setup and normal workflow functions. It
   %  is the only helper that appends the verification-specific `snow` folder.

   arguments
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
   end

   % Resolve the base eval root first so the appended snow folder remains an
   % explicit verification-boundary decision.
   evaluation_data_root = icemodel.verification.helpers.evaluationDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);

   % Append the snow family root used by both ESM-SnowMIP and Laugh-Tests.
   root = fullfile(evaluation_data_root, "snow");
end
