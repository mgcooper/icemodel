function rollbackBaselineProfilePublication(transaction)
   %ROLLBACKBASELINEPROFILEPUBLICATION Restore the prior profiler sidecar.
   %
   %  icemodel.test.helpers.rollbackBaselineProfilePublication(transaction)

   arguments
      transaction (1, 1) struct
   end

   if ~isfield(transaction, 'final_dir') ...
         || isblanktext(string(transaction.final_dir))
      return
   end

   final_dir = string(transaction.final_dir);
   backup_dir = string(transaction.backup_dir);
   icemodel.test.helpers.removeBaselineProfileStage(final_dir);
   if ~isblanktext(backup_dir)
      [ok, message] = movefile(backup_dir, final_dir);
      if ~ok
         error('icemodel:test:perf:profileRollbackFailed', '%s', message)
      end
   end
end
