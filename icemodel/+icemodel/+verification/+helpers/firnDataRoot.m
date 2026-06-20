function root = firnDataRoot(kwargs)
   %FIRNDATAROOT Return the staged firn-verification data root.
   %
   %  root = icemodel.verification.helpers.firnDataRoot()
   %  root = icemodel.verification.helpers.firnDataRoot(icemodel_config_casename="test")
   %
   % Inputs
   %  evaluation_data_root       Base evaluation-data root. When blank, the
   %                             path is resolved from icemodel.config.
   %  icemodel_config_casename   Config casename used to resolve the default
   %                             evaluation-data root without mutating config.
   %
   % Outputs
   %  root   Full path to the firn-verification data root, normally
   %         `<repo>/demo/data/eval/firn` for test/demo workflows. Per-source
   %         families (promice, sumup) are appended by callers as the next
   %         path segment, mirroring how the snow root holds esm_snowmip and
   %         laugh_tests.
   %
   % Role
   %  Operational path helper shared by setup and normal workflow functions. It
   %  is the firn-side counterpart of snowDataRoot: the only helper that
   %  appends the verification-specific `firn` folder.
   %
   % See also: icemodel.verification.helpers.snowDataRoot,
   %  icemodel.verification.helpers.evaluationDataRoot,
   %  icemodel.verification.setup.importPromiceSites

   arguments
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
   end

   % Resolve the base eval root first so the appended firn folder remains an
   % explicit verification-boundary decision.
   evaluation_data_root = icemodel.verification.helpers.evaluationDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);

   % Append the firn family root (parallel to snow/). Source families such as
   % promice/ and sumup/ are appended by callers.
   root = fullfile(evaluation_data_root, "firn");
end
