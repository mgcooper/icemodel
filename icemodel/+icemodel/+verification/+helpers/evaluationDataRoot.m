function root = evaluationDataRoot(kwargs)
   %EVALUATIONDATAROOT Resolve the base evaluation-data root.
   %
   %  root = icemodel.verification.helpers.evaluationDataRoot()
   %  root = ...
   %     icemodel.verification.helpers.evaluationDataRoot(evaluation_data_root=path)
   %  root = ...
   %     icemodel.verification.helpers.evaluationDataRoot(icemodel_config_casename="test")
   %
   % Inputs
   %  evaluation_data_root       Explicit base evaluation-data root. When this
   %                             is provided, it is returned unchanged.
   %  icemodel_config_casename   Config casename used to resolve the default
   %                             evaluation-data root without mutating config.
   %
   % Outputs
   %  root   Base evaluation-data root. Dataset-family subfolders such as
   %         `esm_snowmip/` or `promice/` are appended by callers as the next
   %         path segment. The taxonomy is dataset-family-flat: there is no
   %         intermediate `snow/` or `firn/` process-split level.
   %
   % Role
   %  Operational path helper shared by setup and normal workflow functions.
   %  It owns base eval-root resolution; callers append the family folder.

   arguments
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
   end

   % Return the caller-supplied path as a string. This allows any root folder
   % location to be passed in from upstream callers.
   if ~isblanktext(kwargs.evaluation_data_root)
      root = kwargs.evaluation_data_root;
      return
   end

   % Parse either the active-config ICEMODEL_EVAL_PATH or the default one.
   if isblanktext(kwargs.icemodel_config_casename)

      % Return the current value of ICEMODEL_EVAL_PATH. Unless the caller has
      % set ICEMODEL_EVAL_PATH to a custom location, this returns the default
      % demo/data/eval path.
      root = string(icemodel.getpath('eval'));
   else

      % Return the casename-specific value of ICEMODEL_EVAL_PATH without setting
      % the ICEMODEL_EVAL_PATH variable. This enables retrieval of the
      % default ICEMODEL_EVAL_PATH corresponding to a specific icemodel.config
      % casename without mutating the currently active config.
      cfg = icemodel.config( ...
         "casename", kwargs.icemodel_config_casename, ...
         'setenv', false);
      root = string(cfg.ICEMODEL_EVAL_PATH);
   end
end
