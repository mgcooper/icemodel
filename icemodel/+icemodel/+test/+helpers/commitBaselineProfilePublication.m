function commitBaselineProfilePublication(transaction)
   %COMMITBASELINEPROFILEPUBLICATION Remove a retained sidecar backup.
   %
   %  icemodel.test.helpers.commitBaselineProfilePublication(transaction)

   arguments
      transaction (1, 1) struct
   end

   if isfield(transaction, 'backup_dir')
      icemodel.test.helpers.removeBaselineProfileStage( ...
         string(transaction.backup_dir));
   end
end
