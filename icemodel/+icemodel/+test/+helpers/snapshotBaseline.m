function baseline = snapshotBaseline(kind, baseline_tag, smbmodel, overwrite, output_file, simyear)
   %SNAPSHOTBASELINE Save a release snapshot from the rolling test baseline.
   %
   %  baseline = icemodel.test.helpers.snapshotBaseline( ...
   %     "perf", "v1.1", "skinmodel", false, string.empty(), 2016)
   %  baseline = icemodel.test.helpers.snapshotBaseline( ...
   %     "regression", "v1.1", "icemodel", false)
   %
   % Existing release files are immutable. New snapshots also require the
   % rolling source rows to match the forcing identity registered for the
   % requested release tag.

   arguments
      kind (1, :) string {mustBeMember(kind, ["perf", "regression"])}
      baseline_tag (1, :) string
      smbmodel (1, :) string ...
         {icemodel.validators.mustBeFormalSmbmodelName(smbmodel)}
      overwrite (1, 1) logical = false
      output_file string = string.empty()
      simyear double = NaN
   end

   % Use the default release baseline path for this target.
   if isblanktext(output_file)
      output_file = icemodel.test.helpers.baselineFilePath(kind, ...
         smbmodel=smbmodel, baseline_type="release", ...
         baseline_tag=baseline_tag, simyear=simyear);
   end

   % Copy the current rolling baseline bundle without rerunning the model.
   source_file = icemodel.test.helpers.baselineFilePath(kind, ...
      smbmodel=smbmodel, simyear=simyear);
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

   icemodel.test.helpers.assertFormalBaselineForcing(baseline, "rolling");
   icemodel.test.helpers.assertFormalBaselineForcing(baseline, baseline_tag);
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
