function [evaluation_data_root, input_root] = resolveStagingRoots(kwargs)
   %RESOLVESTAGINGROOTS Resolve paired eval/input staging roots for importers.
   %
   %  [evaluation_data_root, input_root] = ...
   %     icemodel.verification.setup.resolveStagingRoots( ...
   %     output_root=..., evaluation_data_root=..., input_data_root=...)
   %
   % Importers share one root convention: output_root, when set, owns both
   % <output_root>/eval and <output_root>/input. Otherwise explicit lower-level
   % roots are honored, falling back through the verification root helpers.

   arguments
      kwargs.output_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
   end

   if kwargs.output_root ~= ""
      evaluation_data_root = fullfile(kwargs.output_root, 'eval');
      input_root = fullfile(kwargs.output_root, 'input');
      return
   end

   evaluation_data_root = icemodel.verification.helpers.evaluationDataRoot( ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
   if kwargs.input_data_root == "" && kwargs.evaluation_data_root ~= ""
      input_root = fullfile(fileparts(evaluation_data_root), 'input');
   else
      input_root = icemodel.verification.helpers.inputDataRoot( ...
         input_data_root=kwargs.input_data_root, ...
         icemodel_config_casename=kwargs.icemodel_config_casename);
   end
end
