function pathname = baselineFilePath(kind, kwargs)
   %BASELINEFILEPATH Return the canonical baseline file path.
   %
   %  pathname = icemodel.test.helpers.baselineFilePath("perf")
   %  pathname = icemodel.test.helpers.baselineFilePath("perf", smbmodel="skinmodel")
   %  pathname = icemodel.test.helpers.baselineFilePath("perf", baseline_tag="v1.1")
   %  pathname = icemodel.test.helpers.baselineFilePath("perf", baseline_type="release")
   %  pathname = icemodel.test.helpers.baselineFilePath("regression", smbmodel="icemodel")
   %
   % If BASELINE_TAG is provided without BASELINE_TYPE, the type is inferred
   % as "release". If BASELINE_TYPE is "release" without a tag, the most
   % recent release baseline is resolved by scanning the baselines directory.

   arguments
      kind (1, :) string {mustBeMember( ...
         kind, ["perf", "regression"])}

      kwargs.smbmodel (1, :) string ...
         {icemodel.validators.mustBeSmbmodelName} ...
         = "icemodel"

      kwargs.baseline_type (1, :) string {mustBeMember( ...
         kwargs.baseline_type, ["rolling", "release"])} ...
         = "rolling"

      kwargs.baseline_tag string {mustBeTextScalarOrEmpty} ...
         = string.empty()

      kwargs.simyear double ...
         = 2016
   end

   [smbmodel, baseline_type, baseline_tag, simyear] = deal( ...
      kwargs.smbmodel, kwargs.baseline_type, kwargs.baseline_tag, ...
      kwargs.simyear);

   % A tag implies a release baseline.
   if ~isblanktext(baseline_tag) && baseline_type == "rolling"
      baseline_type = "release";
   end

   % Resolve the shared test-root and model tag once for both baseline types.
   testdir = icemodel.getpath('test');
   model_tag = icemodel.test.helpers.smbmodelTag(smbmodel);
   baselines_dir = fullfile(testdir, 'baselines');

   % If release is requested but no tag is specified, find the latest one.
   if baseline_type == "release" && isblanktext(baseline_tag)
      baseline_tag = findLatestReleaseTag( ...
         baselines_dir, kind, model_tag, simyear);
   end

   % Perf baselines are keyed by year; regression baselines are not.
   switch kind
      case "perf"
         if ~(isscalar(simyear) && isfinite(simyear))
            error('simyear is required for perf baseline paths')
         end
         if baseline_type == "rolling"
            pathname = fullfile(baselines_dir, ...
               sprintf('perf_baseline_%d_rolling_%s.mat', ...
               simyear, model_tag));
         else
            pathname = fullfile(baselines_dir, ...
               sprintf('perf_baseline_%d_%s_%s.mat', ...
               simyear, icemodel.test.helpers.sanitizeTag( ...
               baseline_tag), model_tag));
         end

      case "regression"
         if baseline_type == "rolling"
            pathname = fullfile(baselines_dir, ...
               "regression_baseline_rolling_" + model_tag + ".mat");
         else
            pathname = fullfile(baselines_dir, ...
               "regression_baseline_" + ...
               icemodel.test.helpers.sanitizeTag(baseline_tag) + ...
               "_" + model_tag + ".mat");
         end
   end
end

function tag = findLatestReleaseTag(baselines_dir, kind, model_tag, simyear)
   %FINDLATESTRELEASETAG Scan the baselines directory for the latest release.

   switch kind
      case "perf"
         pattern = sprintf('perf_baseline_%d_*_%s.mat', simyear, model_tag);
         prefix = sprintf('perf_baseline_%d_', simyear);
      case "regression"
         pattern = sprintf('regression_baseline_*_%s.mat', model_tag);
         prefix = 'regression_baseline_';
   end

   files = dir(fullfile(baselines_dir, pattern));
   if isempty(files)
      error('icemodel:test:noReleaseBaseline', ...
         'No release baseline files found matching: %s', pattern)
   end

   % Filter out rolling baselines.
   names = string({files.name});
   names = names(~contains(names, "_rolling_"));
   if isempty(names)
      error('icemodel:test:noReleaseBaseline', ...
         'No release baseline files found matching: %s', pattern)
   end

   % Extract the sanitized tag from each filename.
   suffix = "_" + model_tag + ".mat";
   tags = strings(size(names));
   for i = 1:numel(names)
      stem = extractBefore(names(i), suffix);
      tags(i) = extractAfter(stem, prefix);
   end

   % Sort version tags and return the latest one. Convert sanitized tags
   % back to dotted form (v1_1 -> v1.1) for natural version comparison,
   % then return the dotted form so the caller sanitizes it normally.
   dotted = strrep(tags, "_", ".");
   sorted = sort(dotted, 'descend');
   tag = sorted(1);
end
