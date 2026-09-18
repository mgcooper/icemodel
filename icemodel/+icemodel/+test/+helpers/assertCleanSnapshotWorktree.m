function assertCleanSnapshotWorktree(baseline_tag, kwargs)
   %ASSERTCLEANSNAPSHOTWORKTREE Require a clean worktree for a release snapshot.
   %
   %  icemodel.test.helpers.assertCleanSnapshotWorktree("v1.3")
   %
   % A versioned snapshot must record a revision that one commit names, so
   % it requires no tracked change and no untracked file. The regression and
   % perf snapshot tools run as separate calls, and the first call writes
   % versioned baseline files that the second call would otherwise see as
   % untracked. This check therefore ignores an untracked file whose name is
   % a release baseline file of BASELINE_TAG for a formal model, of either
   % kind and of any benchmark year. A tracked release file that is
   % modified or deleted is not an output of the sequence: it is a change to
   % an immutable file and refuses the snapshot.
   %
   % Name-value
   %  repo_root       Repository to inspect. Defaults to the icemodel root.
   %                  Another root must be an icemodel checkout, because the
   %                  own-output names are relative to the checkout layout.
   %  command_runner  Runs the git command. Tests pass a stub that returns
   %                  porcelain text; production uses system.
   %
   % See also: snapshot_regression_baseline, snapshot_perf_baseline,
   %  icemodel.test.helpers.managedBaselineSiblings

   arguments
      baseline_tag (1, 1) string
      kwargs.repo_root (1, 1) string = icemodel.internal.fullpath()
      kwargs.command_runner (1, 1) function_handle = @system
   end

   % Quote the repository once because system() passes the command through
   % the host shell. --untracked-files=all lists every file, not only the
   % folder that contains it, so a path comparison below is exact. -z ends
   % each record with NUL and never quotes a path, so a file name that
   % contains a newline, a quote, or " -> " cannot be misread.
   git = "git --no-pager -C " + icemodel.shellQuote(kwargs.repo_root);
   [status, text] = kwargs.command_runner( ...
      git + " status --porcelain -z --untracked-files=all");
   if status ~= 0
      error('icemodel:test:snapshot:worktreeStatusUnavailable', ...
         'git status failed for %s, so the worktree state is unknown.', ...
         char(kwargs.repo_root))
   end
   [entries, statuses, paths] = parsePorcelainRecords(string(text));
   if isempty(entries)
      return
   end

   % The snapshot sequence itself writes these release files; a run that
   % already froze one kind must not refuse the other kind. The perf file
   % name carries the benchmark year, and the regression tool has no year,
   % so the perf names are matched for any four-digit year. Only an
   % untracked file counts as an output: a tracked release file that is
   % modified or deleted is a change to an immutable file.
   own_patterns = ownOutputPatterns(baseline_tag);
   is_own = false(numel(entries), 1);
   for k = 1:numel(entries)
      is_own(k) = statuses(k) == "??" ...
         && any(~cellfun(@isempty, regexp(paths(k), own_patterns, 'once')));
   end
   foreign = entries(~is_own);
   if ~isempty(foreign)
      error('icemodel:test:snapshot:dirtyWorktree', ...
         ['A versioned %s snapshot requires a clean worktree. Commit or ', ...
         'remove:\n%s'], char(baseline_tag), ...
         char(strjoin(foreign, newline)))
   end
end

function patterns = ownOutputPatterns(baseline_tag)
   %OWNOUTPUTPATTERNS Anchored regexps of the release files the sequence writes.
   %
   % The names come from managedBaselineSiblings for the default benchmark
   % year, made relative to this checkout's root because that is where
   % baselineFilePath resolves them; git prints the same relative names for
   % any checkout with the icemodel layout. The year digits in each perf
   % name are then generalized so a sequence run with another simyear is
   % recognized too.
   checkout_root = icemodel.internal.fullpath();
   regression_outputs = icemodel.test.helpers.managedBaselineSiblings( ...
      "regression", baseline_tag, "");
   perf_outputs = icemodel.test.helpers.managedBaselineSiblings( ...
      "perf", baseline_tag, "");
   own_outputs = [regression_outputs(:); perf_outputs(:)];
   patterns = strings(numel(own_outputs), 1);
   for k = 1:numel(own_outputs)
      relative = icemodel.verification.setup.fixtureRelativePosix( ...
         checkout_root, own_outputs(k));
      escaped = regexptranslate('escape', relative);
      patterns(k) = "^" + regexprep(escaped, '\d{4}', '\\d{4}') + "$";
   end
end

function [entries, statuses, paths] = parsePorcelainRecords(text)
   %PARSEPORCELAINRECORDS Split NUL-delimited porcelain into entries and paths.
   %
   % Each record is two status characters, one space, and the path. A rename
   % or copy record (first status character R or C) is followed by one extra
   % NUL-terminated record that holds the original path; the first path is
   % the one on disk.
   records = split(text, char(0));
   records = records(strlength(records) > 0);
   % At most one entry per record; a rename consumes two records.
   entries = strings(numel(records), 1);
   statuses = strings(numel(records), 1);
   paths = strings(numel(records), 1);
   n = 0;
   k = 1;
   while k <= numel(records)
      record = records(k);
      status_code = extractBefore(record, 3);
      pathname = extractAfter(record, 3);
      n = n + 1;
      entries(n) = status_code + " " + pathname;
      statuses(n) = status_code;
      paths(n) = pathname;
      if startsWith(status_code, "R") || startsWith(status_code, "C")
         k = k + 1;
      end
      k = k + 1;
   end
   entries = entries(1:n);
   statuses = statuses(1:n);
   paths = paths(1:n);
end
