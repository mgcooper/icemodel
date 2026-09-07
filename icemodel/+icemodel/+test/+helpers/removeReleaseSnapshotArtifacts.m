function removeReleaseSnapshotArtifacts(pathname)
   %REMOVERELEASESNAPSHOTARTIFACTS Remove a snapshot MAT file and sidecar.
   %
   %  icemodel.test.helpers.removeReleaseSnapshotArtifacts(pathname)
   %
   % A failed snapshot leaves a baseline file and its profiler sidecar behind.
   % Remove both, and report every removal that failed.
   %
   % See also: icemodel.test.helpers.baselineProfilerDir,
   %  icemodel.test.helpers.snapshotBaseline

   % Hold the first failure instead of rethrowing it, so a failed baseline
   % removal cannot leave the sidecar in place.
   cleanup_error = [];
   try
      if isfile(pathname)
         delete(pathname)
      end
   catch err
      cleanup_error = err;
   end

   % The sidecar lives in the managed profiler tree keyed by baseline stem,
   % not beside the baseline file.
   profile_dir = icemodel.test.helpers.baselineProfilerDir(pathname);
   try
      if isfolder(profile_dir)
         rmdir(profile_dir, 's')
      end
   catch err
      % Chain the second failure onto the first so the caller sees both.
      if isempty(cleanup_error)
         cleanup_error = err;
      else
         cleanup_error = addCause(cleanup_error, err);
      end
   end

   % Both removals ran. Report now if either failed.
   if ~isempty(cleanup_error)
      throw(cleanup_error)
   end
end
