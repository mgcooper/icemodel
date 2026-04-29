function tests = test_verification_data_rebuild
   %TEST_VERIFICATION_DATA_REBUILD End-to-end validation that the staged
   %  verification artifacts can be rebuilt from local native source data.
   %
   %  Both ESM-SnowMIP and Laugh-Tests source caches are optional. Each
   %  test case probes the corresponding fetch helper in non-strict mode
   %  and skips with a clear assumption when the cache is absent. With
   %  the cache present, the test rebuilds the staged artifacts to a
   %  temporary directory and validates their schema (manifest fields,
   %  evaluation.mat top-level keys, expected experiments).
   %
   %  This is the operational shape of icemodel-79n.5: a documented
   %  command that proves the rebuild workflow works without forcing
   %  every fresh clone to fetch external data.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   testCase.TestData.tmp = tempname;
   icemodel.helpers.ensureDirExists(testCase.TestData.tmp);
end

function teardown(testCase)
   if exist(testCase.TestData.tmp, 'dir') == 7
      rmdir(testCase.TestData.tmp, 's');
   end
end

function test_rebuild_laugh_tests_colbeck1976(testCase)
   % Skip when the Laugh-Tests checkout is not available (sibling
   % fallback covers the common Matt-machine case).
   src = icemodel.verification.setup.fetchLaughTests( ...
      strict=false, silent=true);
   has_source = laughTestsCheckoutComplete(src);
   testCase.assumeTrue(has_source, ...
      sprintf('Laugh-Tests source cache not present at %s', src));

   % Rebuild into an isolated evaluation_data_root so the live
   % demo/data tree is not perturbed by the test.
   eval_root = fullfile(testCase.TestData.tmp, 'eval');
   manifest = icemodel.verification.setup.importLaughTests(src, ...
      evaluation_data_root=string(eval_root), overwrite=true);

   % Manifest schema check.
   verifyEqual(testCase, manifest.dataset_family, "laugh_tests");
   verifyTrue(testCase, isscalar(manifest.cases) || ...
      ~isempty(manifest.cases));

   % Evaluation.mat carries both target sources keyed at the top level.
   eval_path = fullfile(eval_root, 'snow', 'laugh_tests', ...
      'colbeck1976', 'evaluation.mat');
   verifyTrue(testCase, exist(eval_path, 'file') == 2);
   loaded = load(eval_path, 'targets');
   verifyTrue(testCase, isfield(loaded.targets, 'numerical_summa'));
   verifyTrue(testCase, isfield(loaded.targets, 'analytical_clark2017'));

   % Both bundles expose the three Colbeck experiments.
   for src_name = ["numerical_summa", "analytical_clark2017"]
      bundle = loaded.targets.(char(src_name));
      verifyEqual(testCase, ...
         sort(string(fieldnames(bundle.experiments))), ...
         ["exp1"; "exp2"; "exp3"], ...
         sprintf('%s missing expected experiments', src_name));
   end
end

function test_rebuild_esm_snowmip_smoke_sites(testCase)
   src = icemodel.verification.setup.fetchEsmSnowmip( ...
      strict=false, silent=true);
   has_source = esmSnowmipCacheComplete(src);
   testCase.assumeTrue(has_source, ...
      sprintf('ESM-SnowMIP source cache not present at %s', src));

   eval_root = fullfile(testCase.TestData.tmp, 'eval');
   manifest = icemodel.verification.setup.importEsmSnowmip(src, ...
      evaluation_data_root=string(eval_root), overwrite=true);

   verifyEqual(testCase, manifest.dataset_family, "esm_snowmip");

   for case_id = ["cdp", "wfj"]
      eval_path = fullfile(eval_root, 'snow', 'esm_snowmip', ...
         char(case_id), 'evaluation.mat');
      verifyTrue(testCase, exist(eval_path, 'file') == 2, ...
         sprintf('%s evaluation.mat missing after rebuild', case_id));
      loaded = load(eval_path, 'targets');
      verifyEqual(testCase, loaded.targets.format, 'timeseries');
      verifyTrue(testCase, ...
         istimetable(loaded.targets.data) || ...
         isstruct(loaded.targets.data));
   end
end

% ----- helpers ---------------------------------------------------------

function tf = laughTestsCheckoutComplete(src)
   required = ["test_cases/input_data/colbeck1976/colbeck1976_forcing.nc";
               "validation_data/m2_mac_Sept23/colbeck1976/colbeck1976-exp1_G1-1_timestep.nc";
               "validation_data/m2_mac_Sept23/colbeck1976/colbeck1976-exp2_G1-1_timestep.nc";
               "validation_data/m2_mac_Sept23/colbeck1976/colbeck1976-exp3_G1-1_timestep.nc"];
   tf = exist(char(src), 'dir') == 7;
   for i = 1:numel(required)
      tf = tf && exist(char(fullfile(src, required(i))), 'file') == 2;
   end
end

function tf = esmSnowmipCacheComplete(src)
   required = ["met_insitu_cdp_1994_2014.nc";
               "obs_insitu_cdp_1994_2014.nc";
               "met_insitu_wfj_1996_2016.nc";
               "obs_insitu_wfj_1996_2016.nc"];
   tf = exist(char(src), 'dir') == 7;
   for i = 1:numel(required)
      tf = tf && exist(char(fullfile(src, required(i))), 'file') == 2;
   end
end
