function [evaluation_data_root, input_root] = resolveStagingRoots(kwargs)
   %RESOLVESTAGINGROOTS Resolve paired eval/input staging roots for importers.
   %
   %  [evaluation_data_root, input_root] = ...
   %     icemodel.verification.setup.resolveStagingRoots( ...
   %     data_root=..., evaluation_data_root=..., input_data_root=...)
   %
   % DATA_ROOT owns both <data_root>/eval and <data_root>/input. OUTPUT_ROOT is
   % retained as a legacy alias and cannot be combined with DATA_ROOT. Otherwise
   % explicit leaf roots are honored with two-way sibling inference before a
   % nonmutating config-case or repo-local verification default is resolved.

   arguments
      kwargs.output_root (1, 1) string = ""
      kwargs.data_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
   end

   % DATA_ROOT supersedes the legacy alias, but accepting both would make the
   % selected tree ambiguous and is therefore rejected at the boundary.
   if ~isblanktext(kwargs.data_root) && ~isblanktext(kwargs.output_root)
      error('icemodel:verification:resolveStagingRoots:conflictingRoots', ...
         'data_root and output_root cannot both be supplied')
   end

   % Resolve a whole-tree selector before any lower-level leaf roots.
   data_root = kwargs.data_root;
   if isblanktext(data_root)
      data_root = kwargs.output_root;
   end
   if ~isblanktext(data_root)
      evaluation_data_root = fullfile(data_root, 'eval');
      input_root = fullfile(data_root, 'input');
      return
   end

   % Explicit leaf roots are authoritative. When only one is supplied, derive
   % the other from the same parent so a scoped call cannot cross-compose trees.
   if ~isblanktext(kwargs.evaluation_data_root)
      evaluation_data_root = kwargs.evaluation_data_root;
      if ~isblanktext(kwargs.input_data_root)
         input_root = kwargs.input_data_root;
      else
         input_root = fullfile(fileparts(evaluation_data_root), 'input');
      end
      return
   elseif ~isblanktext(kwargs.input_data_root)
      input_root = kwargs.input_data_root;
      evaluation_data_root = fullfile(fileparts(input_root), 'eval');
      return
   end

   % With no explicit roots, resolve both leaves through the same nonmutating
   % case selector. A blank selector preserves the top-level verification tree.
   evaluation_data_root = icemodel.verification.helpers.evaluationDataRoot( ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
   input_root = icemodel.verification.helpers.inputDataRoot( ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
end
