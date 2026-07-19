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
   testCase.TestData.cache = fullfile(icemodel.internal.fullpath(), ...
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
      "error", "latitude", "longitude", "name_key", "name", ...
      "measurement_id", "datetime"];
   present = string(obs.density.Properties.VariableNames);
   testCase.verifyTrue(isdatetime(obs.density.datetime));
   for v = needed
      testCase.verifyTrue(ismember(v, present), ...
         sprintf('density profile missing column %s', v));
   end
   density_units = string(obs.density.Properties.VariableUnits);
   testCase.verifyEqual(density_units(present == "error"), "kg m-3");

   % name_key must resolve to a non-empty core/site label for at least one row.
   testCase.verifyTrue(any(strlength(obs.density.name) > 0));
   testCase.verifyTrue(contains(meta.density_note, "within 7.5 km"));
   verifySumupIdentityBundle(testCase, obs, meta, "density", "density");
   testCase.verifyTrue(isdatetime(meta.observation_period_start));
   testCase.verifyTrue(isdatetime(meta.observation_period_end));
   testCase.verifyFalse(isnat(meta.observation_period_start));
   testCase.verifyFalse(isnat(meta.observation_period_end));
   testCase.verifyLessThanOrEqual(meta.observation_period_start, ...
      meta.observation_period_end);
end

function test_import_driver_obs_staged_when_forcing_throws(testCase)
   % importSumup driver guard: the SUMup observations are the PRIMARY target, so
   % absent RCM sources (here non-existent MAR/MERRA/RACMO source dirs) must
   % degrade to SKIPPED forcing legs, NOT a skipped SUMup point. The obs case
   % still stages (observations.mat + a manifest entry) with mar/merra/racmo
   % recorded as skipped legs. Exercises the decoupled obs/forcing flow so the
   % observations are never blocked by a forcing failure.

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

   % Bogus RCM source dirs so the forcing legs are skipped (source not found);
   % the obs cache is real so the observation leg succeeds.
   bogus = fullfile(tmp, 'no_such_rcm_source');

   manifest = icemodel.verification.setup.importSumup( ...
      string(testCase.TestData.cache), ...
      points=testCase.TestData.point, case_ids="kanu", build_forcing=true, ...
      forcing_sources="mar", ...
      mar_dir=bogus, merra_dir=bogus, racmo_dir=bogus, ...
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
   testCase.verifyFalse(ismember("mar3.11", string(c.forcing_sources)));
   testCase.verifyFalse(ismember("racmo2.3p3", string(c.eval_sources)));

   % Unbounded SUMup imports report the actual observation coverage, not blanks.
   testCase.verifyGreaterThan(strlength(string(c.period.start)), 0);
   testCase.verifyGreaterThan(strlength(string(c.period.end)), 0);
   testCase.verifyTrue(endsWith(string(c.period.start), "00:00"));
   testCase.verifyTrue(contains(string(c.period.end), ":"));

   % The colocation record marks only the requested MAR leg as skipped.
   testCase.verifyTrue(isfield(c.colocation, 'mar'));
   testCase.verifyFalse(logical(c.colocation.mar.staged));
   testCase.verifyFalse(isfield(c.colocation, 'merra'));
   testCase.verifyFalse(isfield(c.colocation, 'racmo'));

   % The SUMup observation leg is still recorded.
   testCase.verifyTrue(isfield(c.colocation, 'sumup'));
   testCase.verifyEqual(string(c.colocation.sumup.kind), "firn_profile_obs");
end

function test_import_driver_stages_forcing_legs_when_present(testCase)
   % With staged MAR + MERRA + RACMO fixtures present, the co-located forcing
   % legs stage. MAR + MERRA each write both a met file and a userdata (Data)
   % file; RACMO writes Data only.

   assumeCachePresent(testCase);
   mar_dir = string(fullfile(icemodel.internal.fullpath('data'), 'test', 'forcing', 'mar'));
   merra_dir = string(fullfile(icemodel.internal.fullpath('data'), 'test', 'forcing', 'merra2'));
   racmo_dir = string(fullfile(icemodel.internal.fullpath('data'), 'test', 'forcing', 'racmo'));
   testCase.assumeTrue( ...
      isfolder(mar_dir) && isfolder(merra_dir) && isfolder(racmo_dir), ...
      'MAR/MERRA/RACMO fixtures absent; skipping green-path test.');

   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');
   mkdir(eval_root);
   mkdir(fullfile(input_root, 'met'));
   mkdir(fullfile(input_root, 'userdata'));

   manifest = icemodel.verification.setup.importSumup( ...
      string(testCase.TestData.cache), ...
      points=testCase.TestData.point, case_ids="kanu", years=2012, ...
      startdate="2012-01-01", enddate="2012-01-31 23:00:00", ...
      build_forcing=true, ...
      mar_dir=mar_dir, merra_dir=merra_dir, racmo_dir=racmo_dir, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      overwrite=true);

   c = manifest.cases(1);
   % MAR + MERRA each write a met file and a userdata (Data) file.
   for src = ["mar", "merra"]
      product_id = icemodel.verification.namelists.rcmProductIds(src);
      leg = c.colocation.(char(src));
      testCase.verifyTrue(logical(leg.staged), ...
         sprintf('%s leg should be staged', src));
      testCase.verifyTrue(isfield(leg, 'met_files') && ~isempty(leg.met_files), ...
         sprintf('%s must write a met file', src));
      testCase.verifyTrue(isfield(leg, 'data_files') && ~isempty(leg.data_files), ...
         sprintf('%s must ALSO write a userdata (Data) file', src));
      testCase.verifyNotEmpty( ...
         dir(fullfile(input_root, 'met', char(product_id), '*.mat')));
      testCase.verifyNotEmpty( ...
         dir(fullfile(input_root, 'userdata', char(product_id), '*.mat')));
      met_file = fullfile(input_root, 'met', string(leg.met_files(1)));
      data_file = fullfile(input_root, 'userdata', string(leg.data_files(1)));
      saved_met = load(met_file, 'met');
      saved_data = load(data_file, 'Data');
      testCase.verifyTrue(endsWith(string(met_file), "_15m.mat"));
      testCase.verifyEqual(seconds(median(diff(saved_met.met.Time))), 900);
      testCase.verifyEqual(seconds(median(diff(saved_data.Data.Time))), 3600);
   end

   testCase.verifyTrue(ismember("mar3.11", string(c.forcing_sources)));
   testCase.verifyTrue(ismember("merra2", string(c.forcing_sources)));
   testCase.verifyTrue(all(ismember(["sumup_obs", "mar3.11", "merra2", ...
      "racmo2.3p3"], string(c.eval_sources))));
   % RACMO: Data only, no met.
   testCase.verifyTrue(isfield(c.colocation.racmo, 'data_files') ...
      && ~isempty(c.colocation.racmo.data_files));
   testCase.verifyFalse(isfield(c.colocation.racmo, 'met_files'));
   testCase.verifyNotEmpty( ...
      dir(fullfile(input_root, 'userdata', 'racmo2.3p3', '*racmo*.mat')));
   prior_merra = jsonencode(c.colocation.merra);
   prior_racmo = jsonencode(c.colocation.racmo);

   % A MAR-only refresh is additive and keeps omitted existing RCM legs.
   absent_cache = fullfile(tmp, 'absent-sumup-cache');
   manifest = icemodel.verification.setup.importSumup( ...
      string(absent_cache), ...
      case_ids="kanu", years=2012, ...
      startdate="2012-01-01", enddate="2012-01-31 23:00:00", ...
      forcing_sources="mar", build_forcing=true, ...
      build_observations=false, ...
      mar_dir=mar_dir, merra_dir=merra_dir, racmo_dir=racmo_dir, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      overwrite=true);
   c = manifest.cases(1);
   testCase.verifyTrue(logical(c.colocation.merra.staged));
   testCase.verifyTrue(logical(c.colocation.racmo.staged));
   testCase.verifyEqual(jsonencode(c.colocation.merra), prior_merra);
   testCase.verifyEqual(jsonencode(c.colocation.racmo), prior_racmo);
   testCase.verifyTrue(ismember("sumup_obs", string(c.eval_sources)));
   testCase.verifyTrue(ismember("merra2", string(c.forcing_sources)));
   testCase.verifyTrue(ismember("racmo2.3p3", string(c.eval_sources)));
   testCase.verifyFalse(isfolder(absent_cache));
end

function test_import_skips_pre_rcm_observations_before_writing(testCase)
   % WEG-B observations predate every 2012 fixture. The importer must reject
   % all three legs before delegated staging so no unreferenced RCM file can be
   % left in the isolated input tree.
   assumeCachePresent(testCase);
   mar_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
      'test', 'forcing', 'mar'));
   merra_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
      'test', 'forcing', 'merra2'));
   racmo_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
      'test', 'forcing', 'racmo'));
   testCase.assumeTrue( ...
      isfolder(mar_dir) && isfolder(merra_dir) && isfolder(racmo_dir), ...
      'MAR/MERRA/RACMO fixtures absent; skipping no-overlap test.');

   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');
   mkdir(eval_root);
   mkdir(fullfile(input_root, 'met'));
   mkdir(fullfile(input_root, 'userdata'));

   manifest = icemodel.verification.setup.importSumup( ...
      string(testCase.TestData.cache), ...
      points=[71.14145799276248, -51.222045560858525], ...
      case_ids="wegb", years=2012, ...
      forcing_sources=["mar", "merra", "racmo"], build_forcing=true, ...
      mar_dir=mar_dir, merra_dir=merra_dir, racmo_dir=racmo_dir, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      overwrite=true);

   c = manifest.cases(1);
   testCase.verifyLessThan( ...
      icemodel.verification.setup.ensureUtc(c.period.end), ...
      icemodel.verification.setup.ensureUtc('2012-01-01'));
   for src = ["mar", "merra", "racmo"]
      leg = c.colocation.(char(src));
      testCase.verifyFalse(logical(leg.staged));
      testCase.verifySubstring(string(leg.reason), ...
         "does not overlap SUMup observations");
      testCase.verifyFalse(isfield(leg, 'met_files'));
      testCase.verifyFalse(isfield(leg, 'data_files'));
   end
   testCase.verifyEmpty(dir(fullfile(input_root, '**', '*.mat')));
end

function test_forcing_only_fast_path_requires_explicit_case_ids(testCase)
   % Point discovery is observation work and is unavailable on the fast path.
   root = fullfile(tempname, 'sumup-fast-path');
   testCase.verifyError(@() ...
      icemodel.verification.setup.importSumup( ...
      "", output_root=root, build_observations=false), ...
      'icemodel:verification:importSumup:fastPathRequiresCaseIds');
end

function test_forcing_only_default_uses_each_staged_case_period(testCase)
   % An omitted forcing window belongs to each staged observation case. A
   % pre-RCM WEG-B case remains byte-logically unchanged while an overlapping
   % KAN-U case attaches MAR without reading the SUMup source cache.
   mar_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
      'test', 'forcing', 'mar'));
   testCase.assumeTrue(isfolder(mar_dir), ...
      'MAR source-light fixture absent; skipping forcing-only period test.');

   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));
    eval_root = fullfile(tmp, 'eval');
    input_root = fullfile(tmp, 'input');
    writeSumupForcingOnlyPeriodManifest(eval_root);
    prior = jsondecode(fileread( ...
       fullfile(eval_root, 'sumup', 'manifest.json')));
    prior_ids = string({prior.cases.case_id});
    prior_wegb = prior.cases(prior_ids == "wegb");

   manifest = icemodel.verification.setup.importSumup( ...
      "", case_ids=["wegb", "kanu"], build_observations=false, ...
      build_forcing=true, forcing_sources="mar", mar_dir=mar_dir, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      overwrite=false);

   ids = string({manifest.cases.case_id});
    wegb = manifest.cases(ids == "wegb");
    testCase.verifyEqual(jsonencode(wegb), jsonencode(prior_wegb));
    testCase.verifyEqual(string(wegb.period.start), ...
       "1929-07-31 00:00:00");
   testCase.verifyEqual(string(wegb.period.end), ...
      "1931-10-04 00:00:00");
   testCase.verifyFalse(isfield(wegb.colocation, 'mar'), ...
      jsonencode(wegb.colocation));
   testCase.verifyEmpty(dir(fullfile(input_root, '**', '*wegb*.mat')));

   kanu = manifest.cases(ids == "kanu");
   testCase.verifyEqual(string(kanu.period.start), ...
      "2012-01-01 00:00:00");
   testCase.verifyEqual(string(kanu.period.end), ...
      "2012-01-31 23:00:00");
   testCase.verifyTrue(logical(kanu.colocation.mar.staged));
    testCase.verifyTrue(ismember("mar3.11", string(kanu.forcing_sources)));
    testCase.verifyNotEmpty(kanu.colocation.mar.met_files);
    testCase.verifyNotEmpty(kanu.colocation.mar.data_files);
    testCase.verifyEqual(string(kanu.colocation.mar.window.start), ...
       "2012-01-01 00:00:00");
    testCase.verifyEqual(string(kanu.colocation.mar.window.end), ...
       "2012-12-31 23:00:00");
 end

function test_reuse_path_preserves_existing_staged_and_diagnostic_legs(testCase)
   % A metadata-only repeat call preserves every existing case leg verbatim,
   % including staged disjoint artifacts and a diagnostic skipped leg.
   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
    input_root = fullfile(tmp, 'input');
    writeSumupReuseOverlapManifest(eval_root);
    prior = jsondecode(fileread( ...
       fullfile(eval_root, 'sumup', 'manifest.json')));

   manifest = icemodel.verification.setup.importSumup( ...
      "", case_ids="kanu", build_observations=false, ...
      build_forcing=false, evaluation_data_root=eval_root, ...
      input_data_root=input_root);
    c = manifest.cases(1);

   testCase.verifyEqual(jsonencode(c), jsonencode(prior.cases));
   testCase.verifyTrue(logical(c.colocation.merra.staged));
   testCase.verifyEqual(string(c.colocation.merra.data_files), ...
      "merra2/kanu_merra2_disjoint.mat");
   testCase.verifyFalse(logical(c.colocation.racmo.staged));
   testCase.verifyEqual(string(c.colocation.racmo.reason), ...
      "diagnostic leg retained");
   testCase.verifyFalse(isfolder(input_root));
end

function test_reuse_path_restores_requested_prior_zero_overlap_leg(testCase)
   % A zero-overlap refresh cannot replace an existing staged leg with a skip.
   mar_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
      'test', 'forcing', 'mar'));
   testCase.assumeTrue(isfolder(mar_dir), ...
      'MAR source-light fixture absent; skipping prior-leg restore test.');

   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');
   writeSumupReuseOverlapManifest(eval_root);
   manifest_file = fullfile(eval_root, 'sumup', 'manifest.json');
   prior = jsondecode(fileread(manifest_file));
   prior.cases.period = struct('start', '1929-07-31 00:00:00', ...
      'end', '1931-10-04 00:00:00');
   icemodel.verification.setup.writeManifest(manifest_file, prior);

   manifest = icemodel.verification.setup.importSumup( ...
      "", case_ids="kanu", build_observations=false, ...
      build_forcing=true, forcing_sources="mar", mar_dir=mar_dir, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      overwrite=false);

   testCase.verifyEqual(jsonencode(manifest.cases), ...
      jsonencode(prior.cases));
   testCase.verifyEmpty(dir(fullfile(input_root, '**', '*.mat')));
end

function test_import_rejects_half_comparison_window(testCase)
   % Invalid windows win before a non-dry import can create its staging roots.
   root = fullfile(tempname, 'sumup-half-window');
   testCase.verifyError(@() ...
      icemodel.verification.setup.importSumup( ...
      "", points=testCase.TestData.point, case_ids="kanu", ...
      output_root=root, startdate="2012-01-01"), ...
      'icemodel:internal:pairedWindow:invalidWindow');
   testCase.verifyFalse(isfolder(root));

   end_only_root = fullfile(tempname, 'sumup-end-only-window');
   testCase.verifyError(@() ...
      icemodel.verification.setup.importSumup( ...
      "", points=testCase.TestData.point, case_ids="kanu", ...
      output_root=end_only_root, enddate="2012-01-02"), ...
      'icemodel:internal:pairedWindow:invalidWindow');
   testCase.verifyFalse(isfolder(end_only_root));

   reversed_root = fullfile(tempname, 'sumup-reversed-window');
   testCase.verifyError(@() ...
      icemodel.verification.setup.importSumup( ...
      "", points=testCase.TestData.point, case_ids="kanu", ...
      output_root=reversed_root, startdate="2012-01-02", ...
      enddate="2012-01-01"), ...
      'icemodel:internal:pairedWindow:invalidWindow');
   testCase.verifyFalse(isfolder(reversed_root));
end

function test_reuse_and_normal_paths_share_final_provenance(testCase)
   % Reuse-only and freshly staged state converge on the same family finalizer.
   normal = icemodel.verification.setup.importSumup( ...
      "", points=testCase.TestData.point, case_ids="normal", ...
      dry_run=true, build_forcing=false);

   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');
   writeSumupReuseOverlapManifest(eval_root);
   reused = icemodel.verification.setup.importSumup( ...
      "", case_ids="kanu", build_observations=false, ...
      build_forcing=false, evaluation_data_root=eval_root, ...
      input_data_root=input_root);

   for field = ["dataset_family", "source_doi", "source_url", "source_version"]
      testCase.verifyEqual(string(reused.(char(field))), ...
         string(normal.(char(field))));
   end
   testCase.verifyTrue(isfield(reused.cases.colocation, 'mar'));
   testCase.verifyEqual(string(reused.cases.evaluation_file), ...
      "kanu/observations.mat");
   testCase.verifyFalse(isfolder(input_root));
end

function test_unbounded_sumup_advertises_overlapping_short_rcm_legs(testCase)
   % Unbounded SUMup imports report actual all-available observations. A shorter
   % convenience RCM build is still comparable over its overlap and should stay
   % advertised with a clipped window.
   assumeCachePresent(testCase);
   mar_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
      'test', 'forcing', 'mar'));
   merra_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
      'test', 'forcing', 'merra2'));
   testCase.assumeTrue(isfolder(mar_dir) && isfolder(merra_dir), ...
      'MAR/MERRA fixtures absent; skipping source-list guard.');

   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');
   mkdir(eval_root);
   mkdir(fullfile(input_root, 'met'));
   mkdir(fullfile(input_root, 'userdata'));

   manifest = icemodel.verification.setup.importSumup( ...
      string(testCase.TestData.cache), ...
      points=testCase.TestData.point, case_ids="kanu", years=2012, ...
      forcing_sources=["mar", "merra"], build_forcing=true, ...
      mar_dir=mar_dir, merra_dir=merra_dir, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      overwrite=true);

   c = manifest.cases(1);
   testCase.verifyTrue(ismember("sumup_obs", string(c.eval_sources)));
   testCase.verifyTrue(ismember("mar3.11", string(c.forcing_sources)));
   testCase.verifyTrue(ismember("merra2", string(c.forcing_sources)));
   testCase.verifyTrue(ismember("mar3.11", string(c.eval_sources)));
   testCase.verifyTrue(ismember("merra2", string(c.eval_sources)));
   testCase.verifyTrue(logical(c.colocation.mar.staged));
   testCase.verifyTrue(logical(c.colocation.merra.staged));
   testCase.verifyTrue(isfield(c.colocation.mar, 'met_files'));
   testCase.verifyTrue(isfield(c.colocation.merra, 'data_files'));
   testCase.verifyEqual(string(c.colocation.mar.window.start), ...
      "2012-01-01 00:00:00");
   testCase.verifyEqual(string(c.colocation.mar.window.end), ...
      "2012-12-31 23:00:00");
end

function test_temperature_and_smb_parsed(testCase)
   % Subsurface-temperature and SMB groups parse into tables with their value
   % and depth/time channels.

   assumeCachePresent(testCase);

   [obs, meta] = icemodel.verification.setup.buildSumupObservations( ...
      testCase.TestData.point, ...
      source_dir=string(testCase.TestData.cache), radius_km=7.5);

   % Subsurface temperature is a datetime-indexed TIMETABLE (a time series).
   testCase.verifyTrue(istimetable(obs.subsurface_temperature));
   testCase.verifyTrue(isdatetime(obs.subsurface_temperature.Time));
   tvars = string(obs.subsurface_temperature.Properties.VariableNames);
   testCase.verifyTrue(ismember("subsurface_temperature", tvars));
   testCase.verifyTrue(ismember("depth", tvars));
   testCase.verifyTrue(ismember("measurement_id", tvars));
   % The raw days-since-1900 timestamp column is consumed into Time.
   testCase.verifyFalse(ismember("timestamp", tvars));
   tunits = string(obs.subsurface_temperature.Properties.VariableUnits);
   testCase.verifyEqual(tunits(tvars == "error"), "degC");
   verifySumupIdentityBundle(testCase, obs, meta, ...
      "subsurface_temperature", "temperature");

   % SMB is a period TABLE with real start_date/end_date datetimes.
   testCase.verifyTrue(istable(obs.smb));
   svars = string(obs.smb.Properties.VariableNames);
   testCase.verifyTrue(ismember("smb", svars));
   testCase.verifyTrue(ismember("start_date", svars));
   testCase.verifyTrue(ismember("measurement_id", svars));
   testCase.verifyTrue(isdatetime(obs.smb.start_date));
   testCase.verifyTrue(isdatetime(obs.smb.end_date));
   sunits = string(obs.smb.Properties.VariableUnits);
   testCase.verifyEqual(sunits(svars == "error"), "m w.e.");
   verifySumupIdentityBundle(testCase, obs, meta, "smb", "SMB");
end

function test_duplicate_bearing_artifact_reports_exact_row_provenance(testCase)
   % The documented Humphrey duplicate selection should stage one identity-unique
   % artifact whose flat metadata and human note agree on every row count.
   assumeCachePresent(testCase);
   site = icemodel.verification.setup.researchSiteCatalog("humphrey");
   point = [site.lat_wgs84, site.lon_wgs84];
   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');

   manifest = icemodel.verification.setup.importSumup( ...
      string(testCase.TestData.cache), points=point, ...
      case_ids="humphrey_identity", evaluation_data_root=eval_root, ...
      input_data_root=fullfile(tmp, 'input'), build_forcing=false, ...
      overwrite_family=true);
   artifact = fullfile(eval_root, 'sumup', manifest.cases.evaluation_file);
   loaded = load(artifact, 'targets');
   observations = loaded.targets.data;
   meta = loaded.targets.metadata;

   testCase.verifyEqual( ...
      meta.subsurface_temperature_duplicate_rows_removed, 119);
   for bundle = ["density", "subsurface_temperature", "smb"]
      raw_name = char(bundle + "_raw_rows");
      unique_name = char(bundle + "_unique_rows");
      removed_name = char(bundle + "_duplicate_rows_removed");
      testCase.verifyTrue(isfield(meta, raw_name));
      testCase.verifyEqual(meta.(raw_name), ...
         meta.(unique_name) + meta.(removed_name));
      if ~isempty(observations.(bundle))
         source_variable = bundle;
         if bundle == "subsurface_temperature"
            source_variable = "temperature";
         end
         verifySumupIdentityBundle(testCase, observations, meta, ...
            bundle, source_variable);
      end
   end
end

function test_skip_missing_suppresses_fetch_banner(testCase)
   % Import-level skip_missing follows sibling importers and suppresses banners.
   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));
   missing_cache = fullfile(tmp, 'missing-sumup');
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');

   output = evalc([ ...
      'manifest = icemodel.verification.setup.importSumup(' ...
      '"' missing_cache '", points=[' num2str(testCase.TestData.point) '], ' ...
      'case_ids="quiet", build_forcing=false, skip_missing=true, ' ...
      'evaluation_data_root="' eval_root '", input_data_root="' ...
      input_root '", overwrite_family=true);']);

   testCase.verifyFalse(contains(output, ...
      "SUMup firn source cache incomplete"));
   testCase.verifyEmpty(manifest.cases);
   testCase.verifyEqual(string(manifest.skipped.site), "quiet");
end

function test_strict_missing_sumup_prints_retrieval_guidance(testCase)
   % Fatal cache misses keep the actionable SUMup DOI banner visible to callers.
   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));
   missing_cache = fullfile(tmp, 'missing-sumup-strict');
   eval_root = fullfile(tmp, 'eval');
   f = @() icemodel.verification.setup.importSumup( ...
      missing_cache, points=testCase.TestData.point, case_ids="strict", ...
      build_forcing=false, skip_missing=false, ...
      evaluation_data_root=eval_root); %#ok<NASGU>

   output = evalc([ ...
      'testCase.verifyError(f, ' ...
      '''icemodel:verification:fetchSumup:missingSources'');']);

   testCase.verifyTrue(contains(output, "10.18739/A2M61BR5M"));
end

function verifySumupIdentityBundle(testCase, observations, meta, bundle, ...
      source_variable)
   %VERIFYSUMUPIDENTITYBUNDLE Recheck staged rows and count provenance together.
   record = observations.(bundle);
   raw_rows = meta.(char(bundle + "_raw_rows"));
   unique_rows = meta.(char(bundle + "_unique_rows"));
   removed_rows = meta.(char(bundle + "_duplicate_rows_removed"));
   note = string(meta.(char(bundle + "_note")));

   testCase.verifyEqual(raw_rows - unique_rows, removed_rows);
   testCase.verifyEqual(height(record), unique_rows);
   testCase.verifyTrue(contains(note, "rows raw="));
   testCase.verifyTrue(contains(note, "duplicates_removed="));

   % Reverse only datetime shaping, then prove the saved table has no remaining
   % duplicate under the same pre-datetime scientific identity.
   numeric_record = numericSumupRecord(record, bundle);
   [~, checked_raw, checked_unique, checked_removed] = ...
      icemodel.verification.setup.deduplicateSumupRecords( ...
      numeric_record, source_variable);
   testCase.verifyEqual(checked_raw, unique_rows);
   testCase.verifyEqual(checked_unique, unique_rows);
   testCase.verifyEqual(checked_removed, 0);
end

function record = numericSumupRecord(record, bundle)
   %NUMERICSUMUPRECORD Reverse presentation-only datetime shaping for identity QA.
   epoch = datetime(1900, 1, 1, 'TimeZone', 'UTC');
   switch bundle
      case "density"
         record.timestamp = days(record.datetime - epoch);
         record.datetime = [];
      case "subsurface_temperature"
         record = timetable2table(record, 'ConvertRowTimes', true);
         time_name = string(record.Properties.VariableNames{1});
         record.timestamp = days(record.(time_name) - epoch);
         record.(time_name) = [];
         record = renamevars(record, "subsurface_temperature", "temperature");
      case "smb"
         record.start_date = days(record.start_date - epoch);
         record.end_date = days(record.end_date - epoch);
      otherwise
         error('icemodel:test:sumup:badBundle', ...
            'unsupported SUMup bundle %s', bundle)
   end
end

function writeSumupReuseOverlapManifest(eval_root)
   %WRITESUMUPREUSEOVERLAPMANIFEST Seed one prior case with two RCM windows.
   family_root = fullfile(eval_root, 'sumup');
   mkdir(family_root)
   sumup = struct('kind', 'firn_profile_obs', 'staged', true, ...
      'obs_file', 'kanu/observations.mat');
   mar = struct('kind', 'point_met', 'staged', true, ...
      'source', 'mar', 'source_id', 'mar3.11', ...
      'met_files', {{'mar3.11/met_kanu_mar3.11_partial_15m.mat'}}, ...
      'data_files', {{'mar3.11/kanu_mar3.11_partial.mat'}}, ...
      'sample_method', 'nearest', ...
      'window', struct('start', '2012-01-01 00:00:00', ...
      'end', '2012-06-15 23:00:00'));
    merra = struct('kind', 'point_met', 'staged', true, ...
      'source', 'merra', 'source_id', 'merra2', ...
      'met_files', {{'merra2/met_kanu_merra2_disjoint_15m.mat'}}, ...
      'data_files', {{'merra2/kanu_merra2_disjoint.mat'}}, ...
      'sample_method', 'nearest', ...
       'window', struct('start', '2011-01-01 00:00:00', ...
       'end', '2011-12-31 23:00:00'));
   racmo = struct('kind', 'point_data', 'staged', false, ...
      'source', 'racmo', 'source_id', 'racmo2.3p3', ...
      'reason', 'diagnostic leg retained');
   colocation = struct( ...
      'sumup', sumup, 'mar', mar, 'merra', merra, 'racmo', racmo);
   values = { ...
      'kanu'
      'firn_observational'
      'KAN_U'
      'KAN_U'
      ''
      {'firn'}
      ''
      struct('lat_wgs84', 67, 'lon_wgs84', -47, ...
      'x_epsg3413', 0, 'y_epsg3413', 0, 'elev_m', NaN)
      struct('start', '2012-06-01 00:00:00', ...
      'end', '2012-06-30 23:00:00')
      'kanu/observations.mat'
      {'mar3.11', 'merra2'}
      {'sumup_obs', 'mar3.11', 'merra2'}
      {'density'}
      struct()
      colocation
      'irregular'
      'SUMup reuse-overlap fixture'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
   manifest = struct('dataset_family', 'sumup', ...
      'source_doi', '10.18739/A2M61BR5M', ...
      'source_url', 'https://nsidc.org/data/g02288', ...
      'source_version', 'sumup2025', 'retrieval_date', '', ...
      'cases', entry, 'skipped', []);
   fid = fopen(fullfile(family_root, 'manifest.json'), 'w');
   cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '%s\n', jsonencode(manifest));
   clear cleanup
end

function writeSumupForcingOnlyPeriodManifest(eval_root)
   %WRITESUMUPFORCINGONLYPERIODMANIFEST Seed disjoint and overlapping cases.
   family_root = fullfile(eval_root, 'sumup');
   mkdir(family_root)
   entries(1) = sumupForcingOnlyEntry("wegb", "WEG_B", ...
      [71.14145799276248, -51.222045560858525], ...
      struct('start', '1929-07-31 00:00:00', ...
      'end', '1931-10-04 00:00:00'));
   entries(2) = sumupForcingOnlyEntry("kanu", "KAN_U", ...
      [67.00045840765831, -47.02720909947265], ...
      struct('start', '2012-01-01 00:00:00', ...
      'end', '2012-01-31 23:00:00'));
   manifest = struct('dataset_family', 'sumup', ...
      'source_doi', '10.18739/A2M61BR5M', ...
      'source_url', 'https://nsidc.org/data/g02288', ...
      'source_version', 'sumup2025', 'retrieval_date', '', ...
      'cases', entries, 'skipped', []);
   fid = fopen(fullfile(family_root, 'manifest.json'), 'w');
   cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '%s\n', jsonencode(manifest));
   clear cleanup
end

function entry = sumupForcingOnlyEntry(case_id, site_id, point, period)
   %SUMUPFORCINGONLYENTRY Build one observation-only manifest fixture entry.
   alias = char(case_id);
   sumup = struct('kind', 'firn_profile_obs', 'staged', true, ...
      'obs_file', char(fullfile(alias, 'observations.mat')));
   values = { ...
      alias
      'firn_observational'
      char(site_id)
      char(site_id)
      ''
      {'firn'}
      ''
      struct('lat_wgs84', point(1), 'lon_wgs84', point(2), ...
      'x_epsg3413', 0, 'y_epsg3413', 0, 'elev_m', NaN)
      period
      char(fullfile(alias, 'observations.mat'))
      cell(1, 0)
      {'sumup_obs'}
      {'density'}
      struct()
      struct('sumup', sumup)
      'irregular'
      'SUMup forcing-only period fixture'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
end
