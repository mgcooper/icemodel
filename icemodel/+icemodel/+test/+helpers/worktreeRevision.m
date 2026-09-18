function revision = worktreeRevision(repo_root, command_runner)
   %WORKTREEREVISION Return the Git description of the source tree.
   %
   %  revision = icemodel.test.helpers.worktreeRevision()
   %  revision = icemodel.test.helpers.worktreeRevision(repo_root)
   %  revision = icemodel.test.helpers.worktreeRevision(repo_root, ...
   %     command_runner)
   %
   % The value is `git describe --always --dirty`. A clean tree returns its
   % description, such as v1.0.0-458-g0f869e8a. A tree with a tracked change
   % adds the -dirty suffix. Untracked files do not change the value. A Git
   % failure returns "" so a caller that requires a nonblank identity rejects
   % the run.
   %
   % A baseline built on a dirty tree records the -dirty suffix, and
   % snapshotBaseline refuses to freeze a rolling file whose recorded
   % revision carries it. Git is the rollback for tracked baselines.
   %
   % Tests can pass COMMAND_RUNNER to exercise the failure path. Production
   % callers use system.
   %
   % See also: build_perf_baseline, build_regression_baseline,
   %  run_perf_suite, icemodel.test.helpers.snapshotBaseline

   arguments
      repo_root (1, 1) string = icemodel.internal.fullpath()
      command_runner (1, 1) function_handle = @system
   end

   % Quote the repository once because system() passes each command through
   % the host shell.
   git = "git --no-pager -C " + icemodel.shellQuote(repo_root);
   [describe_status, describe_text] = ...
      command_runner(git + " describe --always --dirty");
   if describe_status ~= 0
      revision = "";
      return
   end
   revision = string(strtrim(describe_text));
end
