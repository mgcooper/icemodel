function baseline = snapshotBaseline(kind, baseline_tag, smbmodel, overwrite, output_file, simyear)
   %SNAPSHOTBASELINE Save a release snapshot from the rolling test baseline.
   %
   %  baseline = icemodel.test.helpers.snapshotBaseline("perf", "v1.1", "skinmodel", true, string.empty(), 2016)
   %  baseline = icemodel.test.helpers.snapshotBaseline("regression", "v1.1", "icemodel", true)

   arguments
      kind (1, :) string {mustBeMember(kind, ["perf", "regression"])}
      baseline_tag (1, :) string
      smbmodel (1, :) string ...
         {icemodel.validators.mustBeFormalSmbmodelName(smbmodel)}
      overwrite (1, 1) logical = false
      output_file string = string.empty()
      simyear double = NaN
   end

   % Default to the canonical release-baseline path for this target.
   if isblanktext(output_file)
      output_file = icemodel.test.helpers.baselineFilePath(kind, ...
         smbmodel=smbmodel, baseline_type="release", ...
         baseline_tag=baseline_tag, simyear=simyear);
   end

   if isfile(char(output_file)) && ~overwrite
      error('release %s baseline already exists: %s', kind, char(output_file))
   end

   % Copy from the current rolling baseline bundle, not from a rerun.
   source_file = icemodel.test.helpers.baselineFilePath(kind, ...
      smbmodel=smbmodel, simyear=simyear);
   if ~isfile(char(source_file))
      error('rolling %s baseline is missing: %s', kind, char(source_file))
   end

   % Rewrite the saved baseline metadata in-memory before saving the new
   % release file.
   S = load(char(source_file));
   switch kind
      case "perf"
         if ~isfield(S, 'PerfBaseline')
            error('rolling perf baseline file is malformed: %s', char(source_file))
         end
         S.PerfBaseline = rewriteBaselineTag(S.PerfBaseline, baseline_tag);
         baseline = S.PerfBaseline;

      case "regression"
         if ~isfield(S, 'RegressionBaseline')
            error('rolling regression baseline file is malformed: %s', ...
               char(source_file))
         end
         S.RegressionBaseline = rewriteBaselineTag( ...
            S.RegressionBaseline, baseline_tag);
         baseline = S.RegressionBaseline;
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

   % Save the rewritten bundle and copy the matching profiler sidecar.
   save(char(output_file), '-struct', 'S');
   copyProfilerArtifacts(source_file, output_file);
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

function copyProfilerArtifacts(source_file, output_file)
   %COPYPROFILERARTIFACTS Copy the profiler sidecar folder for a snapshot.

   src_profdir = icemodel.test.helpers.baselineProfilerDir(source_file);
   if ~isfolder(src_profdir)
      return
   end

   dst_profdir = icemodel.test.helpers.baselineProfilerDir(output_file);
   if isfolder(dst_profdir)
      rmdir(dst_profdir, 's');
   end
   copyfile(src_profdir, dst_profdir);
end
