function cases = getFormalTestSuiteCases()
%GETFORMALTESTSUITECASES Return the canonical formal test-suite cases.
%
%  cases = test.helpers.getFormalTestSuiteCases()
%
% Use this helper when you need the ordered list of top-level formal suite
% cases. These are suite-level build/snapshot/run cases such as
% `perf_full_release` or `regression_smoke_rolling`.
%
% Output:
%  cases - table with one row per formal suite case and columns:
%    suite         : regression | perf
%    action        : build | snapshot | run
%    tier          : smoke | full for compare runs; full/blank for build/snapshot
%    baseline_mode : rolling | release
%    result_field  : generated field name used in run_test_bootstrap results
%
% Notes:
%  - This is the single source of truth for the possible formal suite cases
%    used by `run_test_bootstrap`.
%  - These are not the underlying model run cases. Model/site/year/solver
%    combinations are defined separately by `test.helpers.getCaseMatrix` and
%    `test.helpers.getRegressionCaseMatrix`.
%  - The helper returns regression cases first, then perf cases. Within each
%    suite, the ordering is:
%       1. rolling baseline build
%       2. release baseline snapshot
%       3. smoke/full compare against rolling baseline
%       4. smoke/full compare against release baseline
%
% Why this exists:
%  - `run_test_bootstrap` needs one canonical definition of the formal suite
%    lifecycle so adding or removing a suite-level case later is a one-file
%    edit.
%  - The generated `result_field` values keep bootstrap result names
%    consistent with the case metadata instead of maintaining a second
%    hand-written list.

   % Concatenate the regression and perf lifecycle cases into one ordered
   % table. Each suite contributes the same build/snapshot/run pattern.
   cases = [composeCases("regression"); composeCases("perf")];
end

function cases = composeCases(suite)
   % Each suite first defines two non-compare lifecycle cases:
   %  - build rolling baseline
   %  - snapshot rolling baseline into a release baseline
   build_actions = ["build"; "snapshot"];
   build_tier = ["full"; ""];
   build_baseline = ["rolling"; "release"];
   build_cases = table( ...
      repmat(suite, numel(build_actions), 1), ...
      build_actions, ...
      build_tier, ...
      build_baseline, ...
      'VariableNames', {'suite', 'action', 'tier', 'baseline_mode'});

   % Then define the compare runs. `ndgrid` gives the Cartesian product of:
   %  - tier = smoke/full
   %  - baseline_mode = rolling/release
   %
   % This yields four compare cases per suite:
   %  - smoke vs rolling
   %  - full vs rolling
   %  - smoke vs release
   %  - full vs release
   [tiers, baselines] = ndgrid(["smoke"; "full"], ["rolling"; "release"]);
   run_cases = table( ...
      repmat(suite, numel(tiers), 1), ...
      repmat("run", numel(tiers), 1), ...
      tiers(:), ...
      baselines(:), ...
      'VariableNames', {'suite', 'action', 'tier', 'baseline_mode'});

   % Combine lifecycle and compare cases for this suite, then generate the
   % canonical `results.<field>` names used by run_test_bootstrap.
   cases = [build_cases; run_cases];
   cases.result_field = composeResultFields(cases);
end

function fields = composeResultFields(cases)
   % Build the bootstrap results field names directly from the case metadata.
   %
   % Naming rules:
   %  - build/snapshot cases omit tier because they are not compare runs:
   %      regression_rolling
   %      perf_release
   %  - compare cases include tier and baseline mode:
   %      regression_smoke_rolling
   %      perf_full_release
   fields = strings(height(cases), 1);
   for i = 1:height(cases)
      if cases.action(i) == "run"
         fields(i) = cases.suite(i) + "_" + cases.tier(i) + "_" ...
            + cases.baseline_mode(i);
      else
         fields(i) = cases.suite(i) + "_" + cases.baseline_mode(i);
      end
   end
end
