function tests = test_import_sumup
   %TEST_IMPORT_SUMUP Smoke-test the SUMup 2025 observation parser.
   %
   % Exercises icemodel.verification.setup.buildSumupObservations against the
   % real SUMup 2025 Greenland NetCDF files when they are present in the
   % committed verification cache (data/verification/sumup). The test parses
   % the grouped /DATA + /METADATA NetCDF layout near a PROMICE KAN anchor and
   % asserts that the returned firn-observation bundle is well-formed: a
   % density profile table with depth/value/error columns, a subsurface-
   % temperature table, and an SMB table, each linked to a resolvable
   % name_key. No eval artifacts are committed; the test reads only the
   % verification source cache.
   %
   % When the cache is absent (e.g. a checkout without the verification data),
   % the test is assumed-failed (skipped) rather than erroring, so it stays a
   % light smoke check that never demands committing the large source files.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   % The committed SUMup verification source cache lives at the repo-root
   % data/verification/sumup (whitelisted in .gitignore), independent of the
   % test bootstrap's demo/data redirect. Resolve it from this test file's
   % location: test/unit/<file>.m -> repo root is three levels up.
   repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
   testCase.TestData.cache = fullfile(repo_root, ...
      'data', 'verification', 'sumup');
   % KAN_U: densest SUMup co-location (FirnCover / GC-Net KANU cores & pits).
   testCase.TestData.point = [67.00045840765831, -47.02720909947265];
end

function teardown(testCase)
   testCase.TestData.cleanup = [];
end

function assumeCachePresent(testCase)
   nc = dir(fullfile(testCase.TestData.cache, ...
      'SUMup_2025_density_greenland.nc'));
   testCase.assumeNotEmpty(nc, ...
      'SUMup 2025 Greenland cache absent; skipping parser smoke test.');
end

function test_density_profile_parsed(testCase)
   % buildSumupObservations parses the grouped density NetCDF into a profile
   % table carrying the depth / density / error channels near KAN_U.

   assumeCachePresent(testCase);

   [obs, meta] = icemodel.verification.setup.buildSumupObservations( ...
      testCase.TestData.point, ...
      source_dir=string(testCase.TestData.cache), radius_km=7.5);

   testCase.verifyEqual(obs.format, 'subsurface_profile_bundle');
   testCase.verifyTrue(istable(obs.density));
   testCase.verifyGreaterThan(height(obs.density), 0);

   % Density is a depth profile TABLE carrying a real profile datetime column.
   needed = ["density", "start_depth", "stop_depth", "midpoint", ...
      "error", "latitude", "longitude", "name_key", "name", "datetime"];
   present = string(obs.density.Properties.VariableNames);
   testCase.verifyTrue(isdatetime(obs.density.datetime));
   for v = needed
      testCase.verifyTrue(ismember(v, present), ...
         sprintf('density profile missing column %s', v));
   end

   % name_key must resolve to a non-empty core/site label for at least one row.
   testCase.verifyTrue(any(strlength(obs.density.name) > 0));
   testCase.verifyTrue(contains(meta.density_note, "within 7.5 km"));
end

function test_import_driver_obs_staged_when_forcing_throws(testCase)
   % importSumup driver guard: the SUMup observations are the PRIMARY target, so
   % a throwing buildSumupForcing (here forced by a non-existent MAR/RACMO source
   % dir) must degrade to SKIPPED forcing legs, NOT a skipped SUMup point. The
   % obs case still stages (observations.mat + a manifest entry) with mar/racmo
   % recorded as skipped legs. Exercises the per-leg try/catch guard added to
   % importSumup so the obs are never blocked by a forcing failure.

   assumeCachePresent(testCase);

   % Private eval / input roots so the committed demo tree is untouched.
   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');
   mkdir(eval_root);
   mkdir(fullfile(input_root, 'met'));
   mkdir(fullfile(input_root, 'userdata'));

   % Bogus RCM source dirs so buildSumupForcing throws (source not found); the
   % obs cache is real so the observation leg succeeds.
   bogus = fullfile(tmp, 'no_such_rcm_source');

   manifest = icemodel.verification.setup.importSumup( ...
      string(testCase.TestData.cache), ...
      points=testCase.TestData.point, case_ids="kanu", ...
      mar_dir=bogus, racmo_dir=bogus, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      overwrite=true);

   % The point is NOT skipped: it stages as a case despite the forcing failure.
   testCase.verifyEqual(numel(manifest.cases), 1);
   c = manifest.cases(1);
   testCase.verifyEqual(string(c.case_id), "kanu");
   testCase.verifyEqual(string(c.case_type), "firn_observational");

   % The observations bundle is staged on disk.
   obs_file = fullfile(eval_root, 'sumup', 'kanu', 'observations.mat');
   testCase.verifyTrue(isfile(obs_file));

   % sumup_obs is recorded as an eval source; the failed forcing legs are NOT.
   testCase.verifyTrue(ismember("sumup_obs", string(c.eval_sources)));
   testCase.verifyFalse(ismember("mar", string(c.forcing_sources)));
   testCase.verifyFalse(ismember("racmo", string(c.eval_sources)));

   % The colocation record marks mar/racmo as skipped (staged=false) legs.
   testCase.verifyTrue(isfield(c.colocation, 'mar'));
   testCase.verifyFalse(logical(c.colocation.mar.staged));
   testCase.verifyTrue(isfield(c.colocation, 'racmo'));
   testCase.verifyFalse(logical(c.colocation.racmo.staged));

   % The SUMup observation leg is still recorded.
   testCase.verifyTrue(isfield(c.colocation, 'sumup'));
   testCase.verifyEqual(string(c.colocation.sumup.kind), "firn_profile_obs");
end

function test_temperature_and_smb_parsed(testCase)
   % Subsurface-temperature and SMB groups parse into tables with their value
   % and depth/time channels.

   assumeCachePresent(testCase);

   obs = icemodel.verification.setup.buildSumupObservations( ...
      testCase.TestData.point, ...
      source_dir=string(testCase.TestData.cache), radius_km=7.5);

   % Subsurface temperature is a datetime-indexed TIMETABLE (a time series).
   testCase.verifyTrue(istimetable(obs.subsurface_temperature));
   testCase.verifyTrue(isdatetime(obs.subsurface_temperature.Time));
   tvars = string(obs.subsurface_temperature.Properties.VariableNames);
   testCase.verifyTrue(ismember("temperature", tvars));
   testCase.verifyTrue(ismember("depth", tvars));
   % The raw days-since-1900 timestamp column is consumed into Time.
   testCase.verifyFalse(ismember("timestamp", tvars));

   % Accumulation is a period TABLE with real start_date/end_date datetimes.
   testCase.verifyTrue(istable(obs.accumulation));
   svars = string(obs.accumulation.Properties.VariableNames);
   testCase.verifyTrue(ismember("smb", svars));
   testCase.verifyTrue(ismember("start_date", svars));
   testCase.verifyTrue(isdatetime(obs.accumulation.start_date));
   testCase.verifyTrue(isdatetime(obs.accumulation.end_date));
end
