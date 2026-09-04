function revision = worktreeRevision(repo_root, file_hasher, command_runner)
   %WORKTREEREVISION Return a content-sensitive Git source identity.
   %
   %  revision = icemodel.test.helpers.worktreeRevision()
   %  revision = icemodel.test.helpers.worktreeRevision(repo_root)
   %  revision = icemodel.test.helpers.worktreeRevision(repo_root, file_hasher)
   %  revision = icemodel.test.helpers.worktreeRevision( ...
   %     repo_root, file_hasher, command_runner)
   %
   % A clean tree returns its Git description. A dirty tree adds a digest of
   % tracked changes, status entries, and untracked file bytes. Ignored files
   % do not affect the result. A Git or file-read failure returns "" so the
   % performance A/A gate rejects the run.
   %
   % Tests can pass FILE_HASHER or COMMAND_RUNNER to exercise failures.
   % Production callers use fileSha256 and system.

   arguments
      repo_root (1, 1) string = icemodel.internal.fullpath()
      file_hasher (1, 1) function_handle = ...
         @icemodel.verification.setup.fileSha256
      command_runner (1, 1) function_handle = @system
   end

   % Quote the repository once because system() passes each command through
   % the host shell.
   git = "git --no-pager -C " + icemodel.shellQuote(repo_root);
   [describe_status, describe_text] = ...
      command_runner(git + " describe --always");
   [status_status, status_text] = ...
      command_runner(git + " status --porcelain");
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
   [diff_status, tracked_diff] = command_runner(git + " diff --binary HEAD");
   [untracked_status, untracked_text] = command_runner( ...
      git + " -c core.quotepath=false ls-files --others --exclude-standard");
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
