% Run the focused MAR semantics, grid, and density-profile validation gate.
%
% Run from the repository root with:
%   matlab -batch "run('test/tools/run_mar_profile_validation.m')"
%
% The gate runs the complete profile-reader/stager tests, the four shared-path
% regressions touched by profile integration, the source/grid audit, and the
% bounded two-case real-source artifact proof. It does not broadly restage RCM
% products.

% RUN temporarily enters this script's directory, so bootstrap the package path
% from the tracked file location before using the canonical root helper.
bootstrap_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(bootstrap_root, 'icemodel'))
repo_root = string(icemodel.internal.fullpath());

% The two dedicated files are narrow enough to run completely.
complete_files = [ ...
   "test/unit/test_profile_groups.m"; ...
   "test/unit/test_mar_density_profiles.m"];
suite = reshape(testsuite(fullfile(repo_root, complete_files(1))), 1, []);
for file = complete_files(2:end)'
   candidates = testsuite(fullfile(repo_root, file));
   suite = [suite, reshape(candidates, 1, [])]; %#ok<AGROW>
end

% Larger shared test files contribute only the procedures proving the changed
% callback, discovery, comparison/plot, and audit contracts.
selected = { ...
   "test/unit/test_dataset_family_import_helpers.m", ...
   "test_stageDatasetRcmForcing_runs_after_source_callback_before_persist"; ...
   "test/unit/test_comparison_compatibility.m", ...
   "test_colocated_loader_combines_profile_sidecar_with_userdata"; ...
   "test/unit/test_retmip_imau_sources.m", ...
   "test_dated_profiles_pair_exact_dates_without_depth_collapse"; ...
   "test/unit/test_verification_artifact_audit.m", ...
   "test_audits_mar_profile_model_output_contract"};
for n = 1:size(selected, 1)
   candidates = testsuite(fullfile(repo_root, selected{n, 1}));
   names = string({candidates.Name});
   procedure = string(selected{n, 2});
   hit = names == procedure | endsWith(names, "/" + procedure);
   assert(nnz(hit) == 1, ...
      'Focused test procedure was not uniquely resolved: %s', selected{n, 2})
   suite(end + 1) = candidates(hit); %#ok<SAGROW>
end

% Fail before touching durable QA outputs when any focused regression fails.
results = run(suite);
assert(all([results.Passed]), ...
   'Focused MAR profile validation has one or more failing tests.')

% Both replay scripts validate fully before promoting compact evidence under
% data/preview/qa; the profile proof deletes its isolated stage on completion.
run(fullfile(repo_root, 'test', 'tools', 'audit_mar_semantics_and_grid.m'))
run(fullfile(repo_root, 'test', 'tools', 'run_mar_profile_bundle_proof.m'))

fprintf('MAR validation passed: %d focused tests plus both source-backed replays.\n', ...
   numel(results))
