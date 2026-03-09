function pathname = defaultBaselinePath(kind, baseline_type, baseline_tag, smbmodel, simyear)
%DEFAULTBASELINEPATH Return the canonical baseline path for test tooling.
%
%  pathname = test.helpers.defaultBaselinePath("perf", "rolling", "", "skinmodel", 2016)
%  pathname = test.helpers.defaultBaselinePath("regression", "release", "v1.1", "icemodel")

   arguments
      kind (1, :) string {mustBeMember(kind, ["perf", "regression"])}
      baseline_type (1, :) string {mustBeMember(baseline_type, ["rolling", "release"])}
      baseline_tag string = string.empty()
      smbmodel string = "all"
      simyear double = NaN
   end

   rootdir = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
   model_tag = test.helpers.smbmodelTag(smbmodel);

   switch kind
      case "perf"
         if ~(isscalar(simyear) && isfinite(simyear))
            error('simyear is required for perf baseline paths')
         end
         if baseline_type == "rolling"
            pathname = fullfile(rootdir, 'test', 'baselines', ...
               sprintf('perf_baseline_%d_rolling_%s.mat', simyear, model_tag));
         else
            pathname = fullfile(rootdir, 'test', 'baselines', ...
               sprintf('perf_baseline_%d_%s_%s.mat', ...
               simyear, test.helpers.sanitizeTag(baseline_tag), model_tag));
         end

      case "regression"
         if baseline_type == "rolling"
            pathname = fullfile(rootdir, 'test', 'baselines', ...
               "regression_baseline_rolling_" + model_tag + ".mat");
         else
            pathname = fullfile(rootdir, 'test', 'baselines', ...
               "regression_baseline_" + test.helpers.sanitizeTag(baseline_tag) ...
               + "_" + model_tag + ".mat");
         end
   end
end
