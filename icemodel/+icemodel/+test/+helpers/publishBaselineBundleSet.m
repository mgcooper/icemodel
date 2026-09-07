function bundles = publishBaselineBundleSet(kind, bundles)
   %PUBLISHBASELINEBUNDLESET Publish a complete model baseline set.
   %
   %  bundles = icemodel.test.helpers.publishBaselineBundleSet(kind, bundles)
   %
   % See also: build_regression_baseline, build_perf_baseline

   arguments
      kind (1, 1) string {mustBeMember(kind, ["regression", "perf"])}
      bundles cell
   end

   % Archive the prior rolling set before changing any managed file.
   for k = 1:numel(bundles)
      bundle = bundles{k};
      if bundle.baseline_type == "rolling"
         icemodel.test.helpers.archiveManagedBaseline( ...
            bundle.output_file, kind);
      end
   end

   profile_transactions = cell(numel(bundles), 1);
   file_transactions = cell(numel(bundles), 1);
   n_profiles = 0;
   n_files = 0;
   try
      for k = 1:numel(bundles)
         bundle = bundles{k};
         final_profile_dir = string( ...
            icemodel.test.helpers.baselineProfilerDir(bundle.output_file));
         [bundle.profile_artifacts, profile_transactions{k}] = ...
            icemodel.test.helpers.publishBaselineProfile( ...
            bundle.profile_stage_dir, final_profile_dir, ...
            bundle.profile_artifacts);
         n_profiles = k;

         stage_file = saveBundle(kind, bundle);
         try
            file_transactions{k} = publishFile( ...
               stage_file, bundle.output_file);
         catch err
            deleteIfExists(stage_file);
            rethrow(err)
         end
         n_files = k;
         bundles{k} = bundle;
      end
   catch err
      for k = n_files:-1:1
         try
            rollbackFile(file_transactions{k});
         catch cleanup_err
            err = addCause(err, cleanup_err);
         end
      end
      for k = n_profiles:-1:1
         try
            icemodel.test.helpers.rollbackBaselineProfilePublication( ...
               profile_transactions{k});
         catch cleanup_err
            err = addCause(err, cleanup_err);
         end
      end
      rethrow(err)
   end

   % Remove retained backups only after the complete set is public.
   for k = 1:numel(bundles)
      try
         icemodel.test.helpers.commitBaselineProfilePublication( ...
            profile_transactions{k});
      catch cleanup_err
         warning('icemodel:test:baselineBackupCleanupFailed', ...
            'Could not remove a committed profile backup: %s', ...
            cleanup_err.message)
      end
      try
         commitFile(file_transactions{k});
      catch cleanup_err
         warning('icemodel:test:baselineBackupCleanupFailed', ...
            'Could not remove a committed baseline backup: %s', ...
            cleanup_err.message)
      end
   end
end

function stage_file = saveBundle(kind, bundle)
   %SAVEBUNDLE Write one candidate beside its managed destination.

   % Stage in the destination folder so the later move is a rename on the
   % same volume, which no partial copy can interrupt.
   outdir = destinationParent(bundle.output_file);
   if ~isfolder(outdir)
      mkdir(outdir);
   end
   stage_file = string(tempname(outdir)) + ".mat";
   try
      switch kind
         case "regression"
            save(char(stage_file), '-struct', 'bundle', ...
               'RegressionBaseline', 'case_opts', 'meta', ...
               'profile_summary', 'profile_meta', 'profile_artifacts');
         case "perf"
            save(char(stage_file), '-struct', 'bundle', ...
               'PerfBaseline', 'case_opts', 'meta', 'BenchmarkBaseline', ...
               'benchmark_meta', 'profile_summary', 'profile_meta', ...
               'profile_artifacts');
      end
   catch err
      % A failed save can still leave a partial stage file behind.
      deleteIfExists(stage_file);
      rethrow(err)
   end
end

function transaction = publishFile(stage_file, final_file)
   %PUBLISHFILE Replace one baseline file and retain its prior bytes.

   % A directory at the destination would make every move below ambiguous.
   if isfolder(final_file)
      error('icemodel:test:baselinePublishFailed', ...
         'Baseline destination is a directory: %s', final_file)
   end

   % Move the prior file aside instead of overwriting it, so a later bundle's
   % failure can restore the whole set. The backup sits beside the
   % destination, not in tempdir, because a cross-volume move can fail
   % halfway and leave neither copy in place.
   backup_file = "";
   if isfile(final_file)
      backup_file = string(tempname(destinationParent(final_file))) + ".mat";
      [ok, message] = movefile(final_file, backup_file);
      if ~ok
         error('icemodel:test:baselinePublishFailed', '%s', message)
      end
   end

   % When the candidate cannot take its place, put the prior file back before
   % reporting, and report both failures when the restore also fails.
   [ok, message] = movefile(stage_file, final_file);
   if ~ok
      publish_error = MException( ...
         'icemodel:test:baselinePublishFailed', '%s', message);
      if ~isblanktext(backup_file)
         [restored, restore_message] = movefile(backup_file, final_file);
         if ~restored
            publish_error = addCause(publish_error, MException( ...
               'icemodel:test:baselineRollbackFailed', ...
               '%s', restore_message));
         end
      end
      throw(publish_error)
   end
   transaction = struct( ...
      'final_file', string(final_file), 'backup_file', backup_file);
end

function parent = destinationParent(filename)
   %DESTINATIONPARENT Return the absolute parent for an output file.

   % A bare filename has no parent, and MATLAB resolves it against the
   % current folder, so stage and back up there too.
   parent = string(fileparts(char(filename)));
   if isblanktext(parent)
      parent = string(pwd);
   end
end

function rollbackFile(transaction)
   %ROLLBACKFILE Restore one prior managed baseline file.

   % Remove the published candidate first; movefile will not replace it.
   % A blank backup means this publish created the file, so removal is the
   % whole rollback.
   deleteIfExists(transaction.final_file);
   if ~isblanktext(transaction.backup_file)
      [ok, message] = movefile( ...
         transaction.backup_file, transaction.final_file);
      if ~ok
         error('icemodel:test:baselineRollbackFailed', '%s', message)
      end
   end
end

function commitFile(transaction)
   %COMMITFILE Remove one retained baseline backup.

   % Every bundle published, so no rollback can need this copy again.
   deleteIfExists(transaction.backup_file);
end

function deleteIfExists(pathname)
   %DELETEIFEXISTS Delete one file when its path is nonblank.

   % Callers pass a blank path for a transaction that created no backup.
   if ~isblanktext(pathname) && isfile(pathname)
      delete(pathname)
   end
end
