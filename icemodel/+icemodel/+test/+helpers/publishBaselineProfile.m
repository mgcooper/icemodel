function [profile_artifacts, transaction] = publishBaselineProfile( ...
      stage_dir, final_dir, profile_artifacts)
   %PUBLISHBASELINEPROFILE Replace one managed profiler sidecar.
   %
   %  [profile_artifacts, transaction] = ...
   %     icemodel.test.helpers.publishBaselineProfile( ...
   %     stage_dir, final_dir, profile_artifacts)

   arguments
      stage_dir string {mustBeTextScalarOrEmpty}
      final_dir string {mustBeTextScalarOrEmpty}
      profile_artifacts struct
   end

   transaction = struct('final_dir', "", 'backup_dir', "");
   if isblanktext(final_dir)
      return
   end

   final_parent = string(fileparts(final_dir));
   if ~isfolder(final_parent)
      mkdir(final_parent);
   end

   backup_dir = string.empty();
   if isfolder(final_dir)
      backup_dir = string(tempname(final_parent));
      [ok, message] = movefile(final_dir, backup_dir);
      if ~ok
         error('icemodel:test:perf:profilePublishFailed', '%s', message)
      end
   end

   % A blank stage means the accepted baseline has no profiler sidecar.
   if isblanktext(stage_dir)
      if ~isblanktext(backup_dir)
         transaction.final_dir = final_dir;
         transaction.backup_dir = backup_dir;
      end
      return
   end

   [ok, message] = movefile(stage_dir, final_dir);
   if ~ok
      publish_error = MException( ...
         'icemodel:test:perf:profilePublishFailed', '%s', message);
      icemodel.test.helpers.removeBaselineProfileStage(final_dir);
      if ~isblanktext(backup_dir)
         [restored, restore_message] = movefile(backup_dir, final_dir);
         if ~restored
            rollback_error = MException( ...
               'icemodel:test:perf:profileRollbackFailed', ...
               '%s', restore_message);
            publish_error = addCause(publish_error, rollback_error);
         end
      end
      throw(publish_error)
   end
   transaction.final_dir = final_dir;
   transaction.backup_dir = backup_dir;

   % The recorded paths still name the staging directory the capture wrote
   % to. Point them at the published location. These stay absolute; the
   % tracked snapshot converts them to repo-relative POSIX form in
   % snapshotBaseline.
   fields = ["dir", "index_file", "info_file"];
   for field = fields
      if isfield(profile_artifacts, field)
         profile_artifacts.(field) = replace( ...
            string(profile_artifacts.(field)), stage_dir, final_dir);
      end
   end
end
