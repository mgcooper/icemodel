function revision = sourceRevisionGuard(expected_revision, revision_reader)
   %SOURCEREVISIONGUARD Capture or verify the current source revision.
   %
   %  revision = icemodel.test.helpers.sourceRevisionGuard()
   %  revision = icemodel.test.helpers.sourceRevisionGuard(expected_revision)
   %
   % With no expected revision, this function returns the current worktree
   % identity. With an expected revision, it requires the current identity to
   % match. Tests can supply REVISION_READER.
   %
   % Input
   %  expected_revision - worktree identity captured before measurement
   %  revision_reader   - function that returns the current worktree identity
   %
   % Output
   %  revision - current nonblank worktree identity
   %
   % See also: icemodel.test.helpers.worktreeRevision, build_perf_baseline

   arguments
      expected_revision string = string.empty()
      revision_reader (1, 1) function_handle = ...
         @icemodel.test.helpers.worktreeRevision
   end

   revision = string(revision_reader());
   if isblanktext(revision)
      error('icemodel:test:baseline:sourceRevisionMissing', ...
         'The baseline source revision is unavailable.')
   end
   if ~isempty(expected_revision) && revision ~= expected_revision
      error('icemodel:test:baseline:sourceRevisionChanged', ...
         'The source tree changed during baseline generation.')
   end
end
