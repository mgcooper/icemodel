function revision = worktreeRevision(repo_root, file_hasher, ...
      command_runner, kwargs)
   %WORKTREEREVISION Return a content-sensitive Git source identity.
   %
   %  revision = icemodel.test.helpers.worktreeRevision()
   %  revision = icemodel.test.helpers.worktreeRevision(repo_root)
   %  revision = icemodel.test.helpers.worktreeRevision(repo_root, file_hasher)
   %  revision = icemodel.test.helpers.worktreeRevision( ...
   %     repo_root, file_hasher, command_runner)
   %  revision = icemodel.test.helpers.worktreeRevision( ...
   %     repo_root, ...
   %     ignored_paths="test/baselines/perf_baseline_2016_rolling_icemodel.mat")
   %
   % A clean tree returns its Git description. A dirty tree adds a digest of
   % tracked changes, status entries, and untracked file bytes. Listing
   % untracked files uses git --exclude-standard, so a file excluded by
   % .gitignore, .git/info/exclude, or the global excludes does not affect the
   % result. A Git or file-read failure returns "" so the performance A/A gate
   % rejects the run.
   %
   % Name-value
   %  ignored_paths  Paths this call must also leave out of the digest, so a
   %                 builder's own outputs cannot change the identity of the
   %                 source that produced them. Relative paths resolve against
   %                 REPO_ROOT; paths outside REPO_ROOT are dropped.
   %
   % Tests can pass FILE_HASHER or COMMAND_RUNNER to exercise failures.
   % Production callers use fileSha256 and system.

   arguments
      repo_root (1, 1) string = icemodel.internal.fullpath()
      file_hasher (1, 1) function_handle = ...
         @icemodel.verification.setup.fileSha256
      command_runner (1, 1) function_handle = @system
      kwargs.ignored_paths string = strings(0, 1)
   end

   % Quote the repository once because system() passes each command through
   % the host shell.
   git = "git --no-pager -C " + icemodel.shellQuote(repo_root);
   pathspec = gitPathspec(repo_root, kwargs.ignored_paths);
   [describe_status, describe_text] = ...
      command_runner(git + " describe --always");
   [status_status, status_text] = ...
      command_runner(git + " status --porcelain" + pathspec);
   if describe_status ~= 0 || status_status ~= 0
      revision = "";
      return
   end

   revision = string(strtrim(describe_text));
   if isblanktext(status_text)
      return
   end

   % Git renders a binary-safe tracked patch as text. Hash each untracked file
   % separately so its bytes do not pass through system().
   [diff_status, tracked_diff] = command_runner( ...
      git + " diff --binary HEAD" + pathspec);
   [untracked_status, untracked_text] = command_runner( ...
      git + " -c core.quotepath=false ls-files --others --exclude-standard" ...
      + pathspec);
   if diff_status ~= 0 || untracked_status ~= 0
      revision = "";
      return
   end

   untracked_files = sort(splitlines(string(untracked_text)));
   untracked_files = untracked_files(strlength(untracked_files) > 0);
   records = strings(2 + numel(untracked_files), 2);
   records(1, :) = ["tracked-diff", ...
      icemodel.verification.setup.textSha256(string(tracked_diff))];
   records(2, :) = ["status", ...
      icemodel.verification.setup.textSha256(string(status_text))];
   try
      for n = 1:numel(untracked_files)
         pathname = fullfile(repo_root, untracked_files(n));
         records(n + 2, :) = ["untracked:" + untracked_files(n), ...
            file_hasher(pathname)];
      end
   catch
      revision = "";
      return
   end

   worktree_hash = icemodel.verification.setup.textSha256( ...
      strjoin(records(:, 1) + ":" + records(:, 2), newline));
   revision = revision + "-dirty-" + extractBetween(worktree_hash, 1, 12);
end

function pathspec = gitPathspec(repo_root, ignored_paths)
   %GITPATHSPEC Exclude exact builder-owned outputs from source identity.
   pathspec = "";
   repo_root = icemodel.helpers.canonicalPath(repo_root);
   ignored_paths = reshape(string(ignored_paths), [], 1);
   for n = 1:numel(ignored_paths)
      % java.io.File classifies absolute paths on every supported platform.
      % icemodel.helpers.absolutePath uses the same call, so the two agree.
      if ~java.io.File(char(ignored_paths(n))).isAbsolute()
         ignored_paths(n) = fullfile(repo_root, ignored_paths(n));
      end
      ignored_paths(n) = icemodel.helpers.canonicalPath(ignored_paths(n));
   end
   ignored_paths = ignored_paths(arrayfun(@(pathname) ...
      icemodel.isPathInside(pathname, repo_root), ignored_paths));
   if isempty(ignored_paths)
      return
   end

   pathspec = " -- .";
   for pathname = ignored_paths.'
      relative = icemodel.verification.setup.fixtureRelativePosix( ...
         repo_root, pathname);
      pathspec = pathspec + " " ...
         + icemodel.shellQuote(":(exclude,literal)" + relative);
   end
end
