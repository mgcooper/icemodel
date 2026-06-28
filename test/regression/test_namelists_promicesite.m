function tests = test_namelists_promicesite
   %TEST_NAMELISTS_PROMICESITE Auto-discovered PROMICE station namelist (vq5).
   %
   % icemodel.verification.namelists.promicesite is the single source of truth
   % for the full PROMICE station list, discovered from the on-disk hourly L3
   % product. Verifies it discovers stations from the committed verification
   % cache and errors cleanly when the product is absent.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
   testCase.TestData.promice_dir = string(fullfile(repo_root, ...
      'data', 'verification', 'promice'));
end

function teardown(testCase)
   testCase.TestData.cleanup = [];
end

function test_discovers_from_cache(testCase)
   % Discovers a non-empty station row from the committed product, including the
   % KAN_L anchor, with no duplicates.
   testCase.assumeTrue(isfolder(fullfile(testCase.TestData.promice_dir, 'hour')), ...
      'PROMICE hourly product absent; skipping discovery test.');

   sites = icemodel.verification.namelists.promicesite( ...
      testCase.TestData.promice_dir);

   testCase.verifyClass(sites, 'string');
   testCase.verifyEqual(size(sites, 1), 1, 'sites must be a row');
   testCase.verifyNotEmpty(sites);
   testCase.verifyTrue(ismember("KAN_L", sites));
   testCase.verifyEqual(numel(sites), numel(unique(sites)), 'no duplicates');
end

function test_errors_when_absent(testCase)
   % A directory with no <STATION>_hour.nc files errors with the namelist id.
   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));

   testCase.verifyError( ...
      @() icemodel.verification.namelists.promicesite(string(tmp)), ...
      'icemodel:verification:namelists:promicesite:noStations');
end
