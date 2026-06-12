function cleanupSyntheticWorkspace(workspace)
   %CLEANUPSYNTHETICWORKSPACE Restore env vars and remove a temp workspace.
   %
   %  icemodel.test.fixtures.cleanupSyntheticWorkspace(workspace)

   % Restore whichever icemodel config fields the fixture snapshot captured.
   if isfield(workspace, 'oldenv')
      envnames = fieldnames(workspace.oldenv);
      for n = 1:numel(envnames)
         setenv(envnames{n}, workspace.oldenv.(envnames{n}));
      end
   end

   % Remove the temporary workspace tree after the caller has finished.
   if isfield(workspace, 'rootdir') && isfolder(workspace.rootdir)
      rmdir(workspace.rootdir, 's');
   end
end
