function root = inputDataRoot(kwargs)
   %INPUTDATAROOT Resolve the base icemodel input-data root.
   %
   %  root = icemodel.verification.helpers.inputDataRoot()
   %  root = icemodel.verification.helpers.inputDataRoot( ...
   %     input_data_root=path)
   %  root = icemodel.verification.helpers.inputDataRoot( ...
   %     icemodel_config_casename="test")
   %
   % Inputs
   %  input_data_root            Explicit base input-data root. When this
   %                             is provided, it is returned unchanged.
   %  icemodel_config_casename   Optional config casename used to resolve the
   %                             input-data root without mutating config.
   %                             Leave blank to use the repo top-level data
   %                             tree, independent of ICEMODEL_DATA_PATH.
   %
   % Outputs
   %  root   Base input-data root (the parent of `met/`). The verification
   %         setup writes staged ESM-SnowMIP met files into `<root>/met/`
   %         following the standard icemodel naming convention.
   %
   % Role
   %  Mirrors evaluationDataRoot for the input side. Setup tooling uses
   %  this to stage met files in the standard icemodel input layout so
   %  configureRun + createMetFileNames + loadmet resolve them without
   %  verification-only branches.

   arguments
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
   end

   if ~isblanktext(kwargs.input_data_root)
      root = kwargs.input_data_root;
      return
   end

   % Use a repo-root default for staging workflows so importer output does not
   % depend on the user's active ICEMODEL_DATA_PATH. A nonblank casename keeps
   % the older committed-fixture path available for tests and demo reads.
   if isblanktext(kwargs.icemodel_config_casename)
      root = string(fullfile(icemodel.internal.fullpath('data'), 'input'));
   else
      cfg = icemodel.config( ...
         "casename", kwargs.icemodel_config_casename, ...
         'setenv', false);
      root = string(cfg.ICEMODEL_INPUT_PATH);
   end
end
