function cleanupSyntheticWorkspace(workspace)
%CLEANUPSYNTHETICWORKSPACE Restore env vars and remove a temp workspace.
%
%  icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace)

   envnames = {'ICEMODEL_DATA_PATH', 'ICEMODEL_INPUT_PATH', ...
      'ICEMODEL_OUTPUT_PATH', 'ICEMODEL_EVAL_PATH', 'ICEMODEL_USERDATA_PATH'};
   for n = 1:numel(envnames)
      if isfield(workspace, 'oldenv') && isfield(workspace.oldenv, envnames{n})
         setenv(envnames{n}, workspace.oldenv.(envnames{n}));
      end
   end

   if isfield(workspace, 'rootdir') && exist(workspace.rootdir, 'dir') == 7
      rmdir(workspace.rootdir, 's');
   end
end
