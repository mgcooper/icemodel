function [baseline_type, baseline_tag, output_file, rootdir, input_path, ...
      output_path, cases] = prepareBaselineBuild(kind, baseline, ...
      baseline_tag, tier, smbmodel, output_file, simyear, solver, ...
      smoke_sites, full_sites)
%PREPAREBASELINEBUILD Resolve shared setup for perf/regression baseline builds.
%
%  [baseline_type, baseline_tag, output_file, rootdir, input_path, ...
%     output_path, cases] = icemodel.test.helpers.prepareBaselineBuild(kind, ...
%     baseline, baseline_tag, tier, smbmodel, output_file, simyear, ...
%     solver, smoke_sites, full_sites)
%
% This helper centralizes the common build-time setup shared by
% `build_perf_baseline` and `build_regression_baseline`:
%  1. resolve rolling vs release baseline naming
%  2. locate the repo root and configure the public test paths
%  3. load the canonical formal case matrix for the requested suite
%
% Inputs:
%  kind - "perf" or "regression"
%  baseline, baseline_tag, smbmodel, output_file, simyear, solver,
%  smoke_sites, full_sites
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
      smoke_sites string = "kanm"
      full_sites string = ["kanm"; "kanl"]
   end

   [baseline_type, baseline_tag, output_file] = ...
      icemodel.test.helpers.resolveBaselineBuild(kind, baseline, baseline_tag, ...
      smbmodel, output_file, simyear);

   rootdir = icemodel.internal.fullpath();
   addpath(rootdir);
   addpath(fullfile(rootdir, 'test'));
   [input_path, output_path] = icemodel.test.helpers.configureModelPaths(rootdir);

   switch kind
      case "perf"
         cases = icemodel.test.helpers.getPerfCaseMatrix( ...
            tier=tier, smbmodel=smbmodel, solver=solver, simyear=simyear, ...
            smoke_sites=smoke_sites, full_sites=full_sites);
         if isempty(cases)
            error('no performance cases matched tier=%s smbmodel=%s', ...
               tier, smbmodel)
         end

      case "regression"
         cases = icemodel.test.helpers.getRegressionCaseMatrix( ...
            tier=tier, smbmodel=smbmodel, solver=solver, simyear=simyear, ...
            smoke_sites=smoke_sites, full_sites=full_sites);
         if isempty(cases)
            error('no regression cases matched tier=%s smbmodel=%s', ...
               tier, smbmodel)
         end
   end
end
