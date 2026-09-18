function baseline = snapshotBaseline(kind, baseline_tag, smbmodel, ...
      overwrite, output_file, simyear, source_file)
   %SNAPSHOTBASELINE Save a release snapshot from the rolling test baseline.
   %
   %  baseline = icemodel.test.helpers.snapshotBaseline( ...
   %     "perf", "v1.1", "skinmodel", false, string.empty(), 2016)
   %  baseline = icemodel.test.helpers.snapshotBaseline( ...
   %     "regression", "v1.1", "icemodel", false)
   %
   % Existing release files are immutable. New snapshots also require the
   % rolling source rows to match the forcing identity registered for the
   % requested release tag, and a rolling source whose recorded revision
   % carries the -dirty suffix is refused. A managed release file, whether
   % OUTPUT_FILE is blank or names the managed path, also requires a clean
   % worktree and, for perf, the release-only measurement conditions of
   % assertReleasePerfBaselineSource. A custom OUTPUT_FILE outside the
   % managed tree is a diagnostic copy and skips those two checks.
   % SOURCE_FILE is a test seam that replaces the managed rolling file;
   % production callers omit it.

   arguments
      kind (1, :) string {mustBeMember(kind, ["perf", "regression"])}
      baseline_tag (1, :) string
      smbmodel (1, :) string ...
         {icemodel.validators.mustBeFormalSmbmodelName(smbmodel)}
      overwrite (1, 1) logical = false
      output_file string = string.empty()
      simyear double = NaN
      source_file string = string.empty()
   end

   % Use the default release baseline path for this target. A caller that
   % names the managed release path explicitly gets the managed checks too,
   % so the release gates cannot be skipped by spelling the path out.
   managed_file = string(icemodel.test.helpers.baselineFilePath(kind, ...
      smbmodel=smbmodel, baseline_type="release", ...
      baseline_tag=baseline_tag, simyear=simyear));
   if isblanktext(output_file)
      output_file = managed_file;
   end
   managed_output = icemodel.helpers.canonicalPath(string(output_file)) ...
      == icemodel.helpers.canonicalPath(managed_file);

   % Copy the current rolling baseline bundle without rerunning the model.
   if isblanktext(source_file)
      source_file = icemodel.test.helpers.baselineFilePath(kind, ...
         smbmodel=smbmodel, simyear=simyear);
   end
   if ~isfile(char(source_file))
      error('rolling %s baseline is missing: %s', kind, char(source_file))
   end

   % Select and validate the saved baseline before checking the target state.
   % An out-of-date rolling source is then the first error the caller sees.
   S = load(char(source_file));
   switch kind
      case "perf"
         if ~isfield(S, 'PerfBaseline')
            error('rolling perf baseline file is malformed: %s', char(source_file))
         end
         baseline = S.PerfBaseline;

      case "regression"
         if ~isfield(S, 'RegressionBaseline')
            error('rolling regression baseline file is malformed: %s', ...
               char(source_file))
         end
         baseline = S.RegressionBaseline;
   end

   if ~isfield(S, 'meta') || ~isfield(S.meta, 'git_revision') ...
         || isblanktext(string(S.meta.git_revision))
      error('icemodel:test:releaseBaselineProvenanceMissing', ...
         'The rolling %s baseline has no source revision.', kind)
   end
   % A rolling file built on a dirty tree records a revision no commit
   % names, so a release file frozen from it could not be reproduced.
   if endsWith(string(S.meta.git_revision), "-dirty")
      error('icemodel:test:releaseBaselineSourceDirty', ...
         ['The rolling %s baseline was built on a dirty worktree ', ...
         '(%s). Rebuild it on a clean tree before the snapshot.'], ...
         kind, char(string(S.meta.git_revision)))
   end

   icemodel.test.helpers.assertFormalBaselineForcing(baseline, "rolling");
   icemodel.test.helpers.assertFormalBaselineForcing(baseline, baseline_tag);
   % A managed perf release file freezes a measurement, so the rolling
   % source must pass the release-only measurement conditions read from its
   % metadata. A custom output file is a diagnostic copy outside the managed
   % tree and skips the release conditions, as it skips the clean-tree check.
   if kind == "perf" && managed_output
      icemodel.test.helpers.assertReleasePerfBaselineSource(baseline, S.meta);
   end
   % A managed release file must record a revision one commit names, so the
   % worktree must be clean apart from the release files this snapshot
   % sequence writes. The check lives here so every caller of this helper,
   % not only the snapshot tools, meets it.
   if managed_output
      icemodel.test.helpers.assertCleanSnapshotWorktree(baseline_tag);
   end
   icemodel.test.helpers.assertNewReleaseBaselineTarget( ...
      output_file, overwrite);

   % Rewrite the validated in-memory baseline before saving the new release.
   baseline = rewriteBaselineTag(baseline, baseline_tag);
   switch kind
      case "perf"
         S.PerfBaseline = baseline;
      case "regression"
         S.RegressionBaseline = baseline;
   end

   % Keep the managed benchmark timing bundle aligned with the snapshot.
   if isfield(S, 'BenchmarkBaseline')
      S.BenchmarkBaseline = rewriteBaselineTag(S.BenchmarkBaseline, baseline_tag);
   end
   if isfield(S, 'meta')
      S.meta.baseline_type = "release";
      S.meta.baseline_tag = baseline_tag;
      S.meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   end
   if isfield(S, 'benchmark_meta')
      S.benchmark_meta.baseline_type = "release";
      S.benchmark_meta.baseline_tag = baseline_tag;
      S.benchmark_meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   end
   if isfield(S, 'profile_meta')
      S.profile_meta.baseline_type = "release";
      S.profile_meta.baseline_tag = baseline_tag;
      S.profile_meta.snapshot_utc = datetime('now', 'TimeZone', 'UTC');
   end

   % Copy the profiler sidecar and save both snapshot parts as one operation.
   snapshot_profdir = icemodel.test.helpers.baselineProfilerDir(output_file);
   try
      [S, snapshot_profdir] = copyProfilerArtifacts( ...
         S, source_file, output_file);
      save(char(output_file), '-struct', 'S');
   catch err
      try
         removeSnapshotProfiler(snapshot_profdir);
      catch cleanup_err
         err = addCause(err, cleanup_err);
      end
      try
         if isfile(output_file)
            delete(output_file)
         end
      catch cleanup_err
         err = addCause(err, cleanup_err);
      end
      rethrow(err)
   end
end

function baseline = rewriteBaselineTag(baseline, baseline_tag)
   %REWRITEBASELINETAG Rewrite baseline_type/baseline_tag for a snapshot.

   if isempty(baseline)
      return
   end

   n_rows = height(baseline);
   if ismember('baseline_type', baseline.Properties.VariableNames)
      baseline.baseline_type = repmat("release", n_rows, 1);
   end
   if ismember('baseline_tag', baseline.Properties.VariableNames)
      baseline.baseline_tag = repmat(baseline_tag, n_rows, 1);
   end
end

function [S, dst_profdir] = copyProfilerArtifacts( ...
      S, source_file, output_file)
   %COPYPROFILERARTIFACTS Copy and relink a snapshot's profiler sidecar.

   src_profdir = icemodel.test.helpers.baselineProfilerDir(source_file);
   dst_profdir = "";
   if ~isfolder(src_profdir)
      if isfield(S, 'profile_artifacts')
         S.profile_artifacts = struct();
      end
      return
   end

   dst_profdir = icemodel.test.helpers.baselineProfilerDir(output_file);
   if isfolder(dst_profdir)
      rmdir(dst_profdir, 's');
   end
   copyfile(src_profdir, dst_profdir);

   % A snapshot is tracked, so record the profiler paths relative to the
   % repository root in POSIX form. An absolute path would pin the snapshot to
   % the machine that created it. A destination outside the repository keeps
   % its absolute path, because no relative form exists.
   if isfield(S, 'profile_artifacts')
      repo_root = string(icemodel.internal.fullpath());
      if icemodel.isPathInside(dst_profdir, repo_root)
         recorded_dir = ...
            icemodel.verification.setup.fixtureRelativePosix( ...
            repo_root, dst_profdir);
      else
         recorded_dir = replace(string(dst_profdir), filesep, "/");
      end
      fields = ["dir", "index_file", "info_file"];
      for field = fields
         if isfield(S.profile_artifacts, field)
            if field == "dir"
               S.profile_artifacts.(field) = recorded_dir;
            else
               [~, name, extension] = fileparts(char( ...
                  S.profile_artifacts.(field)));
               S.profile_artifacts.(field) = recorded_dir + "/" ...
                  + string(name) + string(extension);
            end
         end
      end
   end
end

function removeSnapshotProfiler(pathname)
   %REMOVESNAPSHOTPROFILER Remove an unresolved snapshot sidecar.

   if ~isblanktext(pathname) && isfolder(pathname)
      rmdir(pathname, 's');
   end
end
