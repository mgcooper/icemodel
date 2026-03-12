function [baseline_type, baseline_tag, output_file, rootdir, input_path, ...
      output_path, cases] = prepareBaselineBuild(kind, baseline, ...
      baseline_tag, tier, smbmodel, output_file, simyear, solver)
%PREPAREBASELINEBUILD Resolve shared setup for perf/regression baseline builds.
%
%  [baseline_type, baseline_tag, output_file, rootdir, input_path, ...
%     output_path, cases] = test.helpers.prepareBaselineBuild(kind, ...
%     baseline, baseline_tag, tier, smbmodel, output_file, simyear)
%
% This helper centralizes the common build-time setup shared by
% `build_perf_baseline` and `build_regression_baseline`:
%  1. resolve rolling vs release baseline naming
%  2. locate the repo root and configure the public test paths
%  3. load the canonical formal case matrix for the requested suite
%
% Inputs:
%  kind - "perf" or "regression"
%  baseline, baseline_tag, smbmodel, output_file, simyear, solver
%    forwarded to the baseline-path/case selection logic
%
% Output:
%  baseline_type, baseline_tag, output_file - resolved baseline metadata
%  rootdir - repo root
%  input_path, output_path - configured public test paths
%  cases - canonical formal case table for the requested suite

   arguments
      kind (1, :) string {mustBeMember(kind, ["perf", "regression"])}
      baseline string = string.empty()
      baseline_tag string = string.empty()
      tier (1, :) string {mustBeMember(tier, ["smoke", "full", "all"])} = "full"
      smbmodel (1, :) string {mustBeMember(smbmodel, ...
         ["all", "icemodel", "skinmodel"])} = "all"
      output_file string = string.empty()
      simyear double = NaN
      solver = []
   end

   [baseline_type, baseline_tag, output_file] = ...
      test.helpers.resolveBaselineBuild(kind, baseline, baseline_tag, ...
      smbmodel, output_file, simyear);

   rootdir = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
   addpath(rootdir);
   addpath(fullfile(rootdir, 'test'));
   [input_path, output_path] = test.helpers.configureModelPaths(rootdir);

   switch kind
      case "perf"
         cases = test.helpers.getCaseMatrix(tier, smbmodel, solver);
         if isempty(cases)
            error('no performance cases matched tier=%s smbmodel=%s', ...
               tier, smbmodel)
         end
         cases = cases(cases.simyear == simyear, :);
         if isempty(cases)
            error('no performance cases matched simyear=%d', simyear)
         end

      case "regression"
         cases = test.helpers.getRegressionCaseMatrix(tier, smbmodel, solver);
         if isempty(cases)
            error('no regression cases matched tier=%s smbmodel=%s', ...
               tier, smbmodel)
         end
   end
end
