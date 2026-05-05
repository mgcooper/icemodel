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
   %  icemodel_config_casename   Config casename used to resolve the default
   %                             input-data root without mutating config.
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

   if isblanktext(kwargs.icemodel_config_casename)
      root = string(icemodel.getpath('input'));
   else
      cfg = icemodel.config( ...
         "casename", kwargs.icemodel_config_casename, ...
         'setenv', false);
      root = string(cfg.ICEMODEL_INPUT_PATH);
   end
end
