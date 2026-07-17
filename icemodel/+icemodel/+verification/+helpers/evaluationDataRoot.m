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
   %  icemodel_config_casename   Optional config casename used to resolve the
   %                             evaluation-data root without mutating config.
   %                             Leave blank to use the repo top-level data
   %                             tree, independent of ICEMODEL_DATA_PATH.
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

   % Use a repo-root default for staging workflows so importer output does not
   % depend on the user's active ICEMODEL_DATA_PATH. A nonblank casename keeps
   % the older committed-fixture path available for tests and demo reads.
   if isblanktext(kwargs.icemodel_config_casename)
      root = string(fullfile(icemodel.internal.fullpath('data'), 'eval'));
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
