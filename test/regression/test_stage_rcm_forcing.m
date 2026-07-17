function tests = test_stage_rcm_forcing
   %TEST_STAGE_RCM_FORCING Decoupled RCM forcing/Data staging.
   %
   % Covers the observation-import / RCM-forcing decoupling and the contract
   % that every supported RCM writes its complete Data ("userdata") timetable:
   %   * MAR / MERRA write BOTH a met file AND a userdata (Data) file.
   %   * RACMO writes a userdata (Data) file ONLY (no met - it lacks the
   %     near-surface state channels).
   %
   % Two lanes run WITHOUT the external RCM archives (they only need the
   % committed PROMICE verification cache):
   %   - importPromiceSites(build_forcing=false) imports observations only and
   %     writes the manifest, staging NO runtime forcing files.
   %   - stageRcmForcing in manifest-convenience mode resolves legs from a staged
   %     manifest and degrades every source to skip-with-reason when the sources
   %     are absent, never throwing (the cheap fail-early gate).
   % The Data-write contract + completed-source-survives-failure guard use the
   % staged fast RCM fixtures under data/test/forcing.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;

   testCase.TestData.promice_dir = string(fullfile(icemodel.internal.fullpath(), ...
      'data', 'verification', 'promice'));
   testCase.TestData.site = "KAN_L";
   % A single January window keeps the staged fixture-backed green-path builds quick.
   testCase.TestData.startdate = "2012-01-01";
   testCase.TestData.enddate = "2012-01-31";
   testCase.TestData.mar = firstWithData( ...
      string(fullfile(icemodel.internal.fullpath('data'), 'test', 'forcing', 'mar')), ...
      @(p) ~isempty(dir(fullfile(p, "MARv3.11*.nc"))));
   testCase.TestData.merra = firstWithData( ...
      string(fullfile(icemodel.internal.fullpath('data'), 'test', 'forcing', 'merra2')), ...
      @(p) ~isempty(dir(fullfile(p, "slv", "*_Nx.*.nc4*"))));
   testCase.TestData.racmo = firstWithData( ...
      string(fullfile(icemodel.internal.fullpath('data'), 'test', 'forcing', 'racmo')), ...
      @(p) ~isempty(dir(fullfile(p, "*.RACMO*.nc"))));
end

function setup(testCase)
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   tmp = tempname;
   mkdir(tmp);
   testCase.addTeardown(@() rmdir(tmp, 's'));
   testCase.TestData.root = tmp;
end

function teardown(testCase)
   testCase.TestData.cleanup = [];
end

function assumePromicePresent(testCase)
   testCase.assumeTrue(isfile(fullfile(testCase.TestData.promice_dir, ...
      'hour', char(testCase.TestData.site) + "_hour.nc")), ...
      'PROMICE verification cache absent; skipping.');
end

function test_obs_only_when_build_forcing_false(testCase)
   % build_forcing=false imports the PROMICE observations and writes the
   % manifest, but stages NO PROMICE or gridded-RCM runtime files.
   assumePromicePresent(testCase);
   root = testCase.TestData.root;

   manifest = icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, ...
      output_root=string(root), build_forcing=false, overwrite=true);

   testCase.verifyEqual(numel(manifest.cases), 1);
   testCase.verifyEqual(string(manifest.source_version), ...
      "pypromice-L3-hour");
   c = manifest.cases(1);

   % The eval target is staged; the manifest references it.
   testCase.verifyTrue(isfile(fullfile(root, 'eval', 'promice', ...
      'kanl', 'observations.mat')));
   testCase.verifyEqual(string(c.evaluation_file), "kanl/observations.mat");

   % The PROMICE eval leg is present; NO gridded-RCM legs were staged.
   testCase.verifyTrue(isfield(c.colocation, 'promice'));
   testCase.verifyFalse(isfield(c.colocation, 'mar'));
   testCase.verifyFalse(isfield(c.colocation, 'merra'));
   testCase.verifyFalse(isfield(c.colocation, 'racmo'));

   % No runtime forcing/userdata files are written when build_forcing=false.
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'met', ...
      'promice', '*.mat')));
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'userdata', ...
      'promice', '*.mat')));
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'met', ...
      'mar3.11', '*.mat')));
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'met', ...
      'merra2', '*.mat')));
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'userdata', ...
      'mar3.11', '*.mat')));
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'userdata', ...
      'racmo2.3p3', '*.mat')));
end

function test_promice_blank_dt_preserves_native_model_met(testCase)
   % The public PROMICE importer forwards the explicit native-cadence escape
   % hatch while leaving hourly userdata untouched.
   assumePromicePresent(testCase);
   root = testCase.TestData.root;

   manifest = icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, output_root=string(root), ...
      forcing_sources="promice", dt_out="", overwrite=true);

   c = manifest.cases(1);
   met_file = fullfile(root, 'input', 'met', ...
      string(c.colocation.promice.met_files(1)));
   data_file = fullfile(root, 'input', 'userdata', ...
      string(c.colocation.promice.data_files(1)));
   testCase.verifyTrue(endsWith(met_file, "_1hr.mat"));
   saved_met = load(met_file, 'met');
   saved_data = load(data_file, 'Data');
   testCase.verifyEqual(seconds(median(diff(saved_met.met.Time))), 3600);
   testCase.verifyEqual(seconds(median(diff(saved_data.Data.Time))), 3600);
end

function test_disabled_forcing_ignores_inactive_source_selection(testCase)
   % build_forcing=false writes observations only and ignores source selectors.
   assumePromicePresent(testCase);
   root = testCase.TestData.root;

   manifest = icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, ...
      forcing_sources=["mar", "merra", "racmo"], ...
      output_root=string(root), build_forcing=false, overwrite=true);

   c = manifest.cases(1);
   testCase.verifyFalse(isfield(c.colocation, 'mar'));
   testCase.verifyFalse(isfield(c.colocation, 'merra'));
   testCase.verifyFalse(isfield(c.colocation, 'racmo'));
   testCase.verifyEmpty(c.forcing_sources);
   testCase.verifyEqual(string(c.eval_sources), "promice_obs");
   testCase.verifyFalse(c.colocation.promice.staged);
   testCase.verifyTrue(c.colocation.promice.eval_staged);
   testCase.verifyEmpty(c.colocation.promice.met_files);
   testCase.verifyEmpty(c.colocation.promice.data_files);
   testCase.verifyFalse(isfolder(fullfile(root, 'input')));
end

function test_invalid_direct_method_fails_without_output_mutation(testCase)
   % Direct mode rejects unknown sampling methods before discovery or writes.
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'direct-invalid', 'met');
   userdata_outdir = fullfile(root, 'direct-invalid', 'userdata');

   testCase.verifyError( ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      [], forcing_sources="mar", method="linear", ...
      met_outdir=met_outdir, userdata_outdir=userdata_outdir), ...
      'MATLAB:validators:mustBeMember');

   testCase.verifyFalse(isfolder(met_outdir));
   testCase.verifyFalse(isfolder(userdata_outdir));
end

function test_invalid_manifest_method_preserves_manifest_bytes(testCase)
   % Manifest mode applies the same fail-before-read/write method contract.
   root = testCase.TestData.root;
   manifest_file = fullfile(root, 'invalid-method-manifest.json');
   writeCorruptFile(manifest_file);
   before = fileBytes(manifest_file);
   met_outdir = fullfile(root, 'manifest-invalid', 'met');
   userdata_outdir = fullfile(root, 'manifest-invalid', 'userdata');

   testCase.verifyError( ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=string(manifest_file), forcing_sources="mar", ...
      method="linear", met_outdir=met_outdir, ...
      userdata_outdir=userdata_outdir), ...
      'MATLAB:validators:mustBeMember');

   testCase.verifyEqual(fileBytes(manifest_file), before);
   testCase.verifyFalse(isfolder(met_outdir));
   testCase.verifyFalse(isfolder(userdata_outdir));
end

function test_empty_forcing_sources_require_observation_only_mode(testCase)
   % Empty forcing_sources is valid only when forcing builds are disabled.
   assumePromicePresent(testCase);
   root = testCase.TestData.root;

   testCase.verifyError(@() icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, ...
      forcing_sources="", output_root=string(root), overwrite=true), ...
      'icemodel:verification:importPromiceSites:emptyForcingSources');
end

function test_manifest_convenience_skips_absent_sources(testCase)
   % stageRcmForcing manifest-convenience mode: after an observations-only
   % import, it resolves the legs from the staged manifest and, when the RCM
   % sources are absent (bogus dirs), degrades EVERY source to a skip-with-reason
   % WITHOUT throwing - validating the manifest-mode plumbing (read cases ->
   % resolve points/windows -> stage -> merge -> persist) off the fail-early gate.
   % An unrelated pre-existing skipped record is preserved exactly once.
   assumePromicePresent(testCase);
   root = testCase.TestData.root;

   icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, ...
      output_root=string(root), build_forcing=false, overwrite=true);

   manifest_file = fullfile(root, 'eval', 'promice', 'manifest.json');
   prior = jsondecode(fileread(manifest_file));
   prior.skipped = struct('site', "missing_promice", 'reason', "missing");
   icemodel.verification.setup.writeManifest(manifest_file, prior);
   bogus = fullfile(root, 'no_such_rcm_source');

   manifest = icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=string(manifest_file), manifest_file=string(manifest_file), ...
      forcing_sources=["mar", "merra", "racmo"], ...
      met_outdir=fullfile(root, 'input', 'met'), ...
      userdata_outdir=fullfile(root, 'input', 'userdata'), ...
      mar_dir=string(bogus), merra_dir=string(bogus), racmo_dir=string(bogus));

   c = manifest.cases(1);
   testCase.verifyTrue(isfield(c.colocation, 'mar'));
   testCase.verifyFalse(logical(c.colocation.mar.staged));
   testCase.verifyTrue(isfield(c.colocation, 'merra'));
   testCase.verifyFalse(logical(c.colocation.merra.staged));
   testCase.verifyTrue(isfield(c.colocation, 'racmo'));
   testCase.verifyFalse(logical(c.colocation.racmo.staged));
   testCase.verifyEqual(numel(manifest.skipped), 1);
   testCase.verifyEqual(string(manifest.skipped.site), "missing_promice");

   % Nothing was written for the absent sources.
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'met', 'mar3.11', '*.mat')));
   testCase.verifyEmpty( ...
      dir(fullfile(root, 'input', 'userdata', 'racmo2.3p3', '*.mat')));
end

function test_existing_enclosing_files_skip_rebuild(testCase)
   % A previously staged full-period file should satisfy a narrower requested
   % leg. This guards the write-side equivalent of createMetFileNames/loadmet
   % enclosing-file resolution without needing the external RCM archives.
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');
   alias = "kanl";

   point = [67.0, -50.0];
   writeReusableRcmMat(fullfile(met_outdir, 'mar3.11', ...
      'met_kanl_mar3.11_20100101_20200101_1hr.mat'), point);
   writeReusableRcmMat(fullfile(met_outdir, 'merra2', ...
      'met_kanl_merra2_20100101_20200101_1hr.mat'), point);
   writeReusableRcmMat(fullfile(userdata_outdir, 'mar3.11', ...
      'kanl_mar3.11_20100101_20200101.mat'), point);
   writeReusableRcmMat(fullfile(userdata_outdir, 'merra2', ...
      'kanl_merra2_20100101_20200101.mat'), point);
   writeReusableRcmMat(fullfile(userdata_outdir, 'racmo2.3p3', ...
      'kanl_racmo2.3p3_20100101_20200101.mat'), point);

   L = struct('staged', true, 'years', 2013, ...
      'start', icemodel.verification.setup.ensureUtc('2013-01-01'), ...
      'end', icemodel.verification.setup.ensureUtc('2013-12-31 23:00:00'), ...
      'reason', "");
   legspec = struct('alias', alias, 'mar', L, 'merra', L, 'racmo', L);

   colocation = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      point, legspec=legspec, ...
      forcing_sources=["mar", "merra", "racmo"], ...
      met_outdir=string(met_outdir), ...
      userdata_outdir=string(userdata_outdir), ...
      mar_dir=string(fullfile(root, 'absent_mar')), ...
      merra_dir=string(fullfile(root, 'absent_merra')), ...
      racmo_dir=string(fullfile(root, 'absent_racmo')), dt_out=""), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   c = colocation{1};
   testCase.verifyTrue(logical(c.mar.staged));
   testCase.verifyEqual(string(c.mar.met_files), ...
      "mar3.11/met_kanl_mar3.11_20100101_20200101_1hr.mat");
   testCase.verifyEqual(string(c.mar.data_files), ...
      "mar3.11/kanl_mar3.11_20100101_20200101.mat");
   verifyLegReadiness(testCase, c.mar, fullfile(met_outdir, ...
      string(c.mar.met_files)));
   testCase.verifyTrue(logical(c.merra.staged));
   testCase.verifyEqual(string(c.merra.met_files), ...
      "merra2/met_kanl_merra2_20100101_20200101_1hr.mat");
   verifyLegReadiness(testCase, c.merra, fullfile(met_outdir, ...
      string(c.merra.met_files)));
   testCase.verifyTrue(logical(c.racmo.staged));
   testCase.verifyEqual(string(c.racmo.data_files), ...
      "racmo2.3p3/kanl_racmo2.3p3_20100101_20200101.mat");

   % No requested-window duplicates were created because every required output
   % was already covered by an existing staged file.
   testCase.verifyFalse(isfile(fullfile(met_outdir, 'mar3.11', ...
      'met_kanl_mar3.11_20130101_20131231_1hr.mat')));
   testCase.verifyFalse(isfile(fullfile(userdata_outdir, 'racmo2.3p3', ...
      'kanl_racmo2.3p3_20130101_20131231.mat')));
end

function test_partial_cache_fallback_is_visible_when_rebuild_unavailable(testCase)
   % A requested wider window may retain a clipped cache only after rebuild fails.
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');
   alias = "kanl";

   point = [67.0, -50.0];
   writeReusableRcmMat(fullfile(met_outdir, 'mar3.11', ...
      'met_kanl_mar3.11_20120101_20121231_1hr.mat'), point);
   writeReusableRcmMat(fullfile(userdata_outdir, 'mar3.11', ...
      'kanl_mar3.11_20120101_20121231.mat'), point);

   L = struct('staged', true, 'years', 2011:2014, ...
      'start', icemodel.verification.setup.ensureUtc('2011-01-01'), ...
      'end', icemodel.verification.setup.ensureUtc('2014-12-31 23:00:00'), ...
      'reason', "");
   legspec = struct('alias', alias, 'mar', L);

   colocation = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      point, legspec=legspec, forcing_sources="mar", ...
      met_outdir=string(met_outdir), ...
      userdata_outdir=string(userdata_outdir), ...
      mar_dir=string(fullfile(root, 'absent_mar')), dt_out=""), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   c = colocation{1};
   testCase.verifyTrue(logical(c.mar.staged));
   testCase.verifyEqual(string(c.mar.met_files), ...
      "mar3.11/met_kanl_mar3.11_20120101_20121231_1hr.mat");
   testCase.verifyEqual(string(c.mar.data_files), ...
      "mar3.11/kanl_mar3.11_20120101_20121231.mat");
   testCase.verifyEqual(string(c.mar.window.start), ...
      "2012-01-01 00:00:00");
   testCase.verifyEqual(string(c.mar.window.end), ...
      "2012-12-31 23:59:59");
   testCase.verifyTrue(contains(string(c.mar.note), ...
      "not fully cover requested window/output requirements"));
   testCase.verifyTrue(contains(string(c.mar.note), "Raw refresh failed"));

   testCase.verifyFalse(isfile(fullfile(met_outdir, 'mar3.11', ...
      'met_kanl_mar3.11_20110101_20141231_1hr.mat')));
end

function test_corrupt_source_preserves_partial_met_cache(testCase)
   % A corrupt raw NetCDF is a source-local build failure, not a reason to
   % abort the import or discard a compatible partial forcing cache.
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');
   mar_dir = fullfile(root, 'corrupt_mar');
   alias = "kanl";
   point = [67.0, -50.0];
   cached_met = "met_kanl_mar3.11_20120601_20120601_1hr.mat";

   writeReusableRcmMat(fullfile(met_outdir, 'mar3.11', cached_met), point);
   writeCorruptFile(fullfile(mar_dir, 'MARv3.11-2012.nc'));

   L = struct('staged', true, 'years', 2012, ...
      'start', icemodel.verification.setup.ensureUtc('2012-06-01'), ...
      'end', icemodel.verification.setup.ensureUtc('2012-06-02 23:00:00'), ...
      'reason', "");
   legspec = struct('alias', alias, 'mar', L);

   colocation = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      point, legspec=legspec, forcing_sources="mar", ...
      met_outdir=string(met_outdir), ...
      userdata_outdir=string(userdata_outdir), ...
      mar_dir=string(mar_dir), dt_out=""), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   leg = colocation{1}.mar;
   testCase.verifyTrue(logical(leg.staged));
   testCase.verifyEqual(string(leg.met_files), ...
      "mar3.11/" + cached_met);
   testCase.verifyFalse(isfield(leg, 'data_files'));
   testCase.verifyTrue(contains(string(leg.note), "Raw refresh failed"));
end

function test_partial_cache_rebuilds_full_requested_window(testCase)
   % A partial cache cannot suppress a wider requested-source rebuild.
   testCase.assumeTrue(strlength(testCase.TestData.mar) > 0, ...
      'MAR fixture absent; skipping requested-window widening guard.');
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');
   alias = "kanl";
   point = [67.0, -50.0];
   partial_met = "met_kanl_mar3.11_20120110_20120120_1hr.mat";
   partial_data = "kanl_mar3.11_20120110_20120120.mat";
   writeReusableRcmMat(fullfile(met_outdir, 'mar3.11', partial_met), point);
   writeReusableRcmMat(fullfile(userdata_outdir, 'mar3.11', partial_data), point);

   L = struct('staged', true, 'years', 2012, ...
      'start', icemodel.verification.setup.ensureUtc('2012-01-01'), ...
      'end', icemodel.verification.setup.ensureUtc('2012-01-31 23:00:00'), ...
      'discovery_start', ...
      icemodel.verification.setup.ensureUtc('2011-01-01'), ...
      'discovery_end', ...
      icemodel.verification.setup.ensureUtc('2013-12-31 23:00:00'), ...
      'reason', "");
   legspec = struct('alias', alias, 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      point, legspec=legspec, forcing_sources="mar", ...
      met_outdir=string(met_outdir), ...
      userdata_outdir=string(userdata_outdir), ...
      mar_dir=testCase.TestData.mar, dt_out="");

   c = colocation{1}.mar;
   testCase.verifyTrue(logical(c.staged));
   testCase.verifyEqual(string(c.window.start), "2012-01-01 00:00:00");
   testCase.verifyEqual(string(c.window.end), "2012-01-31 23:00:00");
   testCase.verifyFalse(any(endsWith(string(c.met_files), partial_met)));
   testCase.verifyFalse(any(endsWith(string(c.data_files), partial_data)));
   testCase.verifyTrue(all(endsWith(string(c.met_files), "_1hr.mat")));
   native_met = load(fullfile(met_outdir, string(c.met_files(1))), 'met');
   testCase.verifyEqual(seconds(median(diff(native_met.met.Time))), 3600);
end

function test_mar_midday_window_aligns_userdata_and_15m_met_ledgers(testCase)
   % Whole-year MAR metadata must be clipped to the exact three retained UTC
   % days before both hourly userdata and the default 15-minute met are written.
   testCase.assumeTrue(strlength(testCase.TestData.mar) > 0, ...
      'MAR fixture absent; skipping staged daily-ledger alignment guard.');
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');
   alias = "kanl";
   point = [67.0, -50.0];
   L = struct('staged', true, 'years', 2012, ...
      'start', icemodel.verification.setup.ensureUtc('2012-01-10 12:00:00'), ...
      'end', icemodel.verification.setup.ensureUtc('2012-01-12 11:00:00'), ...
      'reason', "");
   legspec = struct('alias', alias, 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      point, legspec=legspec, forcing_sources="mar", ...
      met_outdir=string(met_outdir), ...
      userdata_outdir=string(userdata_outdir), ...
      mar_dir=testCase.TestData.mar);

   % Writers must synchronize Properties.UserData with top-level metadata and
   % preserve the same aligned day ledger across the derived-met boundary.
   leg = colocation{1}.mar;
   data_file = fullfile(userdata_outdir, string(leg.data_files(1)));
   met_file = fullfile(met_outdir, string(leg.met_files(1)));
   data = load(data_file, 'Data', 'artifact_metadata');
   derived = load(met_file, 'met', 'artifact_metadata');
   testCase.verifyEqual(data.Data.Properties.UserData, data.artifact_metadata);
   testCase.verifyEqual(derived.met.Properties.UserData, ...
      derived.artifact_metadata);
   testCase.verifyEqual(seconds(median(diff(data.Data.Time))), 3600);
   testCase.verifyEqual(seconds(median(diff(derived.met.Time))), 900);

   daily_fields = ["mar_qc_runoff_day_status", ...
      "mar_qc_smb_day_status", "mar_qc_runoff_daily_reference_mwe", ...
      "mar_qc_smb_daily_reference_mwe", ...
      "mar_qc_complete_utc_day_count", "mar_qc_partial_utc_day_count", ...
      "mar_diagnostic_melt_day_status", ...
      "mar_diagnostic_melt_daily_reference_mwe", ...
      "mar_diagnostic_melt_residual_mwe_day"];
   for field = daily_fields
      testCase.verifyEqual(data.artifact_metadata.(field), ...
         derived.artifact_metadata.(field));
   end
   testCase.verifyEqual(numel( ...
      data.artifact_metadata.mar_qc_runoff_day_status), 3);
   testCase.verifyEqual( ...
      data.artifact_metadata.mar_qc_complete_utc_day_count, 1);
   testCase.verifyEqual( ...
      data.artifact_metadata.mar_qc_partial_utc_day_count, 2);
   testCase.verifyEqual( ...
      data.artifact_metadata.mar_qc_runoff_day_status([1 3]), uint8([3; 3]));
   testCase.verifyTrue(all(isnan( ...
      data.artifact_metadata.mar_qc_runoff_daily_reference_mwe([1 3]))));
end

function test_overwrite_build_failure_preserves_existing_fallback(testCase)
   % overwrite=true forces a fresh source build attempt, but a raw-source
   % failure must preserve the compatible cache rather than drop the leg.
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');
   alias = "kanl";
   point = [67.0, -50.0];

   writeReusableRcmMat(fullfile(met_outdir, 'mar3.11', ...
      'met_kanl_mar3.11_20100101_20200101_1hr.mat'), point);
   writeReusableRcmMat(fullfile(userdata_outdir, 'mar3.11', ...
      'kanl_mar3.11_20100101_20200101.mat'), point);

   L = struct('staged', true, 'years', 2013, ...
      'start', icemodel.verification.setup.ensureUtc('2013-01-01'), ...
      'end', icemodel.verification.setup.ensureUtc('2013-12-31 23:00:00'), ...
      'reason', "");
   legspec = struct('alias', alias, 'mar', L);

   colocation = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.stageRcmForcing(point, ...
      legspec=legspec, forcing_sources="mar", ...
      met_outdir=string(met_outdir), ...
      userdata_outdir=string(userdata_outdir), ...
      mar_dir=string(fullfile(root, 'absent_mar')), dt_out="", ...
      overwrite=true), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   c = colocation{1};
   testCase.verifyTrue(logical(c.mar.staged));
   testCase.verifyEqual(string(c.mar.met_files), ...
      "mar3.11/met_kanl_mar3.11_20100101_20200101_1hr.mat");
   testCase.verifyEqual(string(c.mar.data_files), ...
      "mar3.11/kanl_mar3.11_20100101_20200101.mat");
   testCase.verifyTrue(contains(string(c.mar.note), "Raw refresh failed"));
end

function test_overwrite_success_replaces_compatible_cache(testCase)
   % A compatible cache is fallback only: successful overwrite writes fresh
   % source-derived met and Data payloads at the same requested filenames.
   testCase.assumeTrue(strlength(testCase.TestData.mar) > 0, ...
      'MAR fixture absent; skipping overwrite-success guard.');
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');
   alias = "kanl";
   point = [67.0, -50.0];
   met_file = fullfile(met_outdir, 'mar3.11', ...
      'met_kanl_mar3.11_20120101_20120131_1hr.mat');
   data_file = fullfile(userdata_outdir, 'mar3.11', ...
      'kanl_mar3.11_20120101_20120131.mat');
   writeReusableRcmMat(met_file, point);
   writeReusableRcmMat(data_file, point);
   old_met = fileBytes(met_file);
   old_data = fileBytes(data_file);

   L = struct('staged', true, 'years', 2012, ...
      'start', icemodel.verification.setup.ensureUtc('2012-01-01'), ...
      'end', icemodel.verification.setup.ensureUtc('2012-01-31 23:00:00'), ...
      'reason', "");
   legspec = struct('alias', alias, 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing(point, ...
      legspec=legspec, forcing_sources="mar", ...
      met_outdir=string(met_outdir), ...
      userdata_outdir=string(userdata_outdir), ...
      mar_dir=testCase.TestData.mar, dt_out="", overwrite=true);

   testCase.verifyTrue(logical(colocation{1}.mar.staged));
   testCase.verifyFalse(isequal(fileBytes(met_file), old_met));
   testCase.verifyFalse(isequal(fileBytes(data_file), old_data));
   testCase.verifyTrue(any(string({whos('-file', met_file).name}) == "met"));
   testCase.verifyTrue(any(string({whos('-file', data_file).name}) == "Data"));
end

function test_concrete_cache_conflict_rebuilds_when_source_available(testCase)
   % A concrete point/method conflict is stale cache, not a reason to skip a
   % requested leg when the raw source can rebuild the correct artifacts.
   testCase.assumeTrue(strlength(testCase.TestData.mar) > 0, ...
      'MAR fixture absent; skipping conflict-rebuild guard.');
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');
   alias = "kanl";
   point = [67.0, -50.0];
   met_file = fullfile(met_outdir, 'mar3.11', ...
      'met_kanl_mar3.11_20120101_20120131_1hr.mat');
   data_file = fullfile(userdata_outdir, 'mar3.11', ...
      'kanl_mar3.11_20120101_20120131.mat');
   writeReusableRcmMat(met_file, [68.0, -49.0], "natural");
   writeReusableRcmMat(data_file, [68.0, -49.0], "natural");

   L = struct('staged', true, 'years', 2012, ...
      'start', icemodel.verification.setup.ensureUtc('2012-01-01'), ...
      'end', icemodel.verification.setup.ensureUtc('2012-01-31 23:00:00'), ...
      'reason', "");
   legspec = struct('alias', alias, 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing(point, ...
      legspec=legspec, forcing_sources="mar", ...
      met_outdir=string(met_outdir), ...
      userdata_outdir=string(userdata_outdir), ...
      mar_dir=testCase.TestData.mar, dt_out="", method="nearest");
   saved_met = load(met_file, 'artifact_metadata');
   saved_data = load(data_file, 'artifact_metadata');

   testCase.verifyTrue(logical(colocation{1}.mar.staged));
   testCase.verifyEqual(string(colocation{1}.mar.sample_method), "nearest");
   testCase.verifyEqual(string(saved_met.artifact_metadata.sample_method), ...
      "nearest");
   testCase.verifyEqual(string(saved_data.artifact_metadata.sample_method), ...
      "nearest");
   testCase.verifyEqual([saved_met.artifact_metadata.lat_wgs84, ...
      saved_met.artifact_metadata.lon_wgs84], point);
   testCase.verifyEqual([saved_data.artifact_metadata.lat_wgs84, ...
      saved_data.artifact_metadata.lon_wgs84], point);
end

function test_manifest_mode_preserves_case_discovery_for_partial_raw_coverage(testCase)
   % Manifest mode pairs eval/input roots and searches the original case period,
   % even when raw coverage clips the required/staged leg to one middle year.
   assumePromicePresent(testCase);
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');

   icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, ...
      output_root=string(root), build_forcing=false, overwrite=true);
   manifest_file = fullfile(root, 'eval', 'promice', 'manifest.json');
   obs_manifest = jsondecode(fileread(manifest_file));
   obs_manifest.cases.period = struct( ...
      'start', '2011-01-01 00:00:00', 'end', '2013-12-31 23:00:00');
   icemodel.verification.setup.writeManifest(manifest_file, obs_manifest);
   point = [obs_manifest.cases.site_location.lat_wgs84, ...
      obs_manifest.cases.site_location.lon_wgs84];

   writeReusableRcmMat(fullfile(met_outdir, 'mar3.11', ...
      'met_kanl_mar3.11_20110101_20131231_1hr.mat'), point);
   writeReusableRcmMat(fullfile(userdata_outdir, 'mar3.11', ...
      'kanl_mar3.11_20110101_20131231.mat'), point);
   % A source-clipped exact candidate coexists, but full-case discovery/ranking
   % must continue to prefer the broader compatible artifact.
   writeReusableRcmMat(fullfile(met_outdir, 'mar3.11', ...
      'met_kanl_mar3.11_20120101_20121231_1hr.mat'), point);
   writeReusableRcmMat(fullfile(userdata_outdir, 'mar3.11', ...
      'kanl_mar3.11_20120101_20121231.mat'), point);
   partial_mar = fullfile(root, 'partial_mar');
   writeCorruptFile(fullfile(partial_mar, 'MARv3.11-test-2012.nc'));

   manifest = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=string(manifest_file), forcing_sources="mar", ...
      mar_dir=string(partial_mar), dt_out=""), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   c = manifest.cases(1);
   testCase.verifyEqual(string(c.colocation.mar.met_files), ...
      "mar3.11/met_kanl_mar3.11_20110101_20131231_1hr.mat");
   testCase.verifyEqual(string(c.colocation.mar.data_files), ...
      "mar3.11/kanl_mar3.11_20110101_20131231.mat");
   testCase.verifyEqual(string(c.colocation.mar.window.start), ...
      "2012-01-01 00:00:00");
   testCase.verifyEqual(string(c.colocation.mar.window.end), ...
      "2012-12-31 23:00:00");
end

function test_manifest_mode_exact_clipped_cache_satisfies_required_window(testCase)
   % Full-case discovery must not widen a source-clipped cache requirement.
   assumePromicePresent(testCase);
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');

   icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, ...
      output_root=string(root), build_forcing=false, overwrite=true);
   manifest_file = fullfile(root, 'eval', 'promice', 'manifest.json');
   obs_manifest = jsondecode(fileread(manifest_file));
   obs_manifest.cases.period = struct( ...
      'start', '2011-01-01 00:00:00', 'end', '2013-12-31 23:00:00');
   icemodel.verification.setup.writeManifest(manifest_file, obs_manifest);
   point = [obs_manifest.cases.site_location.lat_wgs84, ...
      obs_manifest.cases.site_location.lon_wgs84];

   writeReusableRcmMat(fullfile(met_outdir, 'mar3.11', ...
      'met_kanl_mar3.11_20120101_20121231_1hr.mat'), point);
   writeReusableRcmMat(fullfile(userdata_outdir, 'mar3.11', ...
      'kanl_mar3.11_20120101_20121231.mat'), point);
   partial_mar = fullfile(root, 'partial_mar');
   writeCorruptFile(fullfile(partial_mar, 'MARv3.11-test-2012.nc'));

   manifest = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=string(manifest_file), forcing_sources="mar", ...
      mar_dir=string(partial_mar), dt_out=""), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   c = manifest.cases(1).colocation.mar;
   testCase.verifyEqual(string(c.met_files), ...
      "mar3.11/met_kanl_mar3.11_20120101_20121231_1hr.mat");
   testCase.verifyEqual(string(c.data_files), ...
      "mar3.11/kanl_mar3.11_20120101_20121231.mat");
   testCase.verifyEqual(string(c.window.start), "2012-01-01 00:00:00");
   testCase.verifyEqual(string(c.window.end), "2012-12-31 23:00:00");
   testCase.verifyTrue(contains(string(c.note), "fully cover"));
   testCase.verifyFalse(contains(string(c.note), "Raw refresh failed"));
   verifyLegReadiness(testCase, c, fullfile(met_outdir, string(c.met_files)));
end

function test_existing_met_only_does_not_block_userdata_build(testCase)
   % A MAR/MERRA leg is complete only when both met and userdata bracket the
   % requested window. A met-only artifact must not block creation of the
   % missing userdata file.
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');
   alias = "kanl";

   writeEmptyMat(fullfile(met_outdir, 'mar3.11', ...
      'met_kanl_mar3.11_20100101_20200101_1hr.mat'));

   L = struct('staged', true, 'years', 2013, ...
      'start', icemodel.verification.setup.ensureUtc('2013-01-01'), ...
      'end', icemodel.verification.setup.ensureUtc('2013-12-31 23:00:00'), ...
      'reason', "");
   legspec = struct('alias', alias, 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      [67.0, -50.0], legspec=legspec, forcing_sources="mar", ...
      met_outdir=string(met_outdir), ...
      userdata_outdir=string(userdata_outdir), ...
      mar_dir=string(fullfile(root, 'absent_mar')), dt_out="");

   testCase.verifyFalse(logical(colocation{1}.mar.staged));
end

function test_existing_met_only_rcm_is_forcing_only(testCase)
   % A valid MAR/MERRA met-only cache is reusable forcing even without Data.
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');
   alias = "kanl";
   point = [67.0, -50.0];

   writeReusableRcmMat(fullfile(met_outdir, 'mar3.11', ...
      'met_kanl_mar3.11_20100101_20200101_1hr.mat'), point);

   L = struct('staged', true, 'years', 2013, ...
      'start', icemodel.verification.setup.ensureUtc('2013-01-01'), ...
      'end', icemodel.verification.setup.ensureUtc('2013-12-31 23:00:00'), ...
      'reason', "");
   legspec = struct('alias', alias, 'mar', L);

   colocation = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      point, legspec=legspec, forcing_sources="mar", ...
      met_outdir=string(met_outdir), ...
      userdata_outdir=string(userdata_outdir), ...
      mar_dir=string(fullfile(root, 'absent_mar')), dt_out=""), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   testCase.verifyTrue(logical(colocation{1}.mar.staged));
   testCase.verifyTrue(isfield(colocation{1}.mar, 'met_files'));
   testCase.verifyFalse(isfield(colocation{1}.mar, 'data_files'));
   testCase.verifyEqual(string(colocation{1}.mar.window.start), ...
      "2013-01-01 00:00:00");
   testCase.verifyEqual(string(colocation{1}.mar.window.end), ...
      "2013-12-31 23:00:00");
end

function test_mar_merra_write_met_and_userdata(testCase)
   % MAR and MERRA each write both a met file and a
   % userdata (Data) file; RACMO writes a userdata file ONLY (no met).
   assumePromicePresent(testCase);
   testCase.assumeTrue(strlength(testCase.TestData.mar) > 0 ...
      && strlength(testCase.TestData.merra) > 0 ...
      && isfolder(testCase.TestData.racmo), ...
      'MAR/MERRA/RACMO fixtures not available; skipping Data-write guard.');
   root = testCase.TestData.root;

   manifest = icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, output_root=string(root), ...
      forcing_sources=["promice", "mar", "merra", "racmo"], ...
      mar_dir=testCase.TestData.mar, merra_dir=testCase.TestData.merra, ...
      racmo_dir=string(testCase.TestData.racmo), overwrite=true);

   c = manifest.cases(1);

   native_met = load(fullfile(root, 'input', 'met', ...
      string(c.colocation.promice.met_files(1))), 'met');
   native_data = load(fullfile(root, 'input', 'userdata', ...
      string(c.colocation.promice.data_files(1))), 'Data');
   testCase.verifyTrue(endsWith( ...
      string(c.colocation.promice.met_files(1)), "_15m.mat"));
   testCase.verifyEqual(seconds(median(diff(native_met.met.Time))), 900);
   testCase.verifyEqual(seconds(median(diff(native_data.Data.Time))), 3600);

   % MAR + MERRA: a met file AND a Data (userdata) file (the fix).
   for src = ["mar", "merra"]
      product = icemodel.verification.namelists.rcmProductIds(src);
      leg = c.colocation.(char(src));
      testCase.verifyTrue(logical(leg.staged), ...
         sprintf('%s leg should be staged', src));
      testCase.verifyTrue(isfield(leg, 'met_files') && ~isempty(leg.met_files), ...
         sprintf('%s must write a met file', src));
      testCase.verifyTrue(isfield(leg, 'data_files') && ~isempty(leg.data_files), ...
         sprintf('%s must ALSO write a userdata (Data) file', src));
      verifyLegReadiness(testCase, leg, fullfile(root, 'input', 'met', ...
         string(leg.met_files)));
      testCase.verifyNotEmpty( ...
         dir(fullfile(root, 'input', 'met', char(product), '*.mat')));
      testCase.verifyNotEmpty( ...
       dir(fullfile(root, 'input', 'userdata', char(product), '*.mat')));
      met_match = dir(fullfile(root, 'input', 'met', char(product), '*.mat'));
      loaded = load(fullfile(met_match(1).folder, met_match(1).name), ...
         'met', 'artifact_metadata');
      testCase.verifyTrue(endsWith(string(met_match(1).name), "_15m.mat"));
      testCase.verifyEqual(seconds(median(diff(loaded.met.Time))), 900);
      testCase.verifyTrue(ismember("modis", ...
         string(loaded.met.Properties.VariableNames)), ...
         sprintf('%s met should include the default MODIS albedo channel', src));
      testCase.verifyGreaterThan(nnz(isfinite(loaded.met.modis)), 0, ...
         sprintf('%s met MODIS channel should contain finite samples', src));
      verifyModisProvenance(testCase, loaded, "met");

      data_match = dir(fullfile(root, 'input', 'userdata', ...
         char(product), '*.mat'));
      native_data = load(fullfile(data_match(1).folder, ...
         data_match(1).name), 'Data', 'artifact_metadata');
      testCase.verifyEqual(seconds(median(diff(native_data.Data.Time))), 3600);
      verifyModisProvenance(testCase, native_data, "Data");
      if src == "merra"
         % Case-window staging must narrow the builder's full-year tavg3 ledger
         % and carry that exact proof into both saved artifact cadences.
         testCase.verifyTrue( ...
            icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid( ...
            native_data.Data, native_data.artifact_metadata));
         testCase.verifyTrue( ...
            icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid( ...
            loaded.met, loaded.artifact_metadata));
         metadata = native_data.artifact_metadata;
         expected = native_data.Data.Time( ...
            minute(native_data.Data.Time) == 0 ...
            & second(native_data.Data.Time) == 0 ...
            & mod(hour(native_data.Data.Time), 3) == 0);
         testCase.verifyEqual( ...
            metadata.merra_tavg3_expected_source_row_count, numel(expected));
      end
   end

   % RACMO: Data only, no met.
   testCase.verifyTrue(logical(c.colocation.racmo.staged));
   testCase.verifyTrue(isfield(c.colocation.racmo, 'data_files') ...
      && ~isempty(c.colocation.racmo.data_files));
   testCase.verifyFalse(isfield(c.colocation.racmo, 'met_files'));
   testCase.verifyNotEmpty( ...
      dir(fullfile(root, 'input', 'userdata', 'racmo2.3p3', '*.mat')));
   testCase.verifyEmpty( ...
      dir(fullfile(root, 'input', 'met', 'racmo2.3p3', '*.mat')));
   racmo_match = dir(fullfile(root, 'input', 'userdata', ...
      'racmo2.3p3', '*.mat'));
   racmo_data = load(fullfile(racmo_match(1).folder, ...
      racmo_match(1).name), 'Data', 'artifact_metadata');
   verifyModisProvenance(testCase, racmo_data, "Data");

   % Repeating the same ordinary import reuses every current observation,
   % native, and indirect RCM artifact without changing its bytes.
   artifacts = stagedCaseArtifacts(root, c);
   before = cellfun(@fileBytes, cellstr(artifacts), 'UniformOutput', false);
   icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, output_root=string(root), ...
      forcing_sources=["promice", "mar", "merra", "racmo"], ...
      mar_dir=testCase.TestData.mar, merra_dir=testCase.TestData.merra, ...
      racmo_dir=string(testCase.TestData.racmo), overwrite=false);
   after = cellfun(@fileBytes, cellstr(artifacts), 'UniformOutput', false);
   testCase.verifyEqual(after, before);
end

function test_mar_fresh_stage_marks_no_modis_source_coverage(testCase)
   % A valid but nonmatching optional cache writes explicit absence to Data/met,
   % omits an all-NaN pseudo-channel, and passes first-class MODIS artifact QA.
   assumePromicePresent(testCase);
   testCase.assumeTrue(strlength(testCase.TestData.mar) > 0, ...
      'MAR fixture not available; skipping no-MODIS-coverage staging guard.');
   root = testCase.TestData.root;
   modis_dir = fullfile(root, 'empty-modis-cache');
   mkdir(modis_dir)

   manifest = icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, output_root=string(root), ...
      forcing_sources=["promice", "mar"], ...
      mar_dir=testCase.TestData.mar, modis_dir=string(modis_dir), ...
      overwrite=true);
   leg = manifest.cases(1).colocation.mar;

   met_loaded = load(fullfile(root, 'input', 'met', ...
      string(leg.met_files(1))), 'met', 'artifact_metadata');
   data_loaded = load(fullfile(root, 'input', 'userdata', ...
      string(leg.data_files(1))), 'Data', 'artifact_metadata');
   verifyNoModisProvenance(testCase, met_loaded, "met");
   verifyNoModisProvenance(testCase, data_loaded, "Data");

   % The reusable artifact auditor must accept the fresh no-source encoding
   % without a post-hoc metadata repair.
   report = icemodel.verification.auditArtifacts( ...
      evaluation_data_root=fullfile(root, 'eval'), ...
      input_data_root=fullfile(root, 'input'), families="promice", ...
      report_dir="");
   codes = string({report.findings.code});
   modis_codes = codes(startsWith(codes, "modis_"));
   testCase.verifyEmpty(modis_codes);
end

function test_completed_sources_survive_racmo_failure(testCase)
   % Per-source isolation: with MAR/MERRA present but RACMO absent, the
   % MAR/MERRA met+userdata files AND their manifest legs are staged; only the
   % RACMO leg degrades to skip-with-reason. A failing/absent source never rolls
   % back a completed one.
   assumePromicePresent(testCase);
   testCase.assumeTrue(strlength(testCase.TestData.mar) > 0 ...
      && strlength(testCase.TestData.merra) > 0, ...
      'MAR/MERRA fixtures not available; skipping isolation guard.');
   root = testCase.TestData.root;
   bogus = fullfile(root, 'no_such_racmo');

   manifest = icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, output_root=string(root), ...
      forcing_sources=["promice", "mar", "merra", "racmo"], ...
      mar_dir=testCase.TestData.mar, merra_dir=testCase.TestData.merra, ...
      racmo_dir=string(bogus), overwrite=true);

   c = manifest.cases(1);
   testCase.verifyTrue(logical(c.colocation.mar.staged));
   testCase.verifyTrue(logical(c.colocation.merra.staged));
   testCase.verifyFalse(logical(c.colocation.racmo.staged));
   testCase.verifyNotEmpty( ...
      dir(fullfile(root, 'input', 'met', 'mar3.11', '*.mat')));
   testCase.verifyNotEmpty( ...
      dir(fullfile(root, 'input', 'userdata', 'merra2', '*.mat')));
end

function test_manifest_checkpoint_survives_later_source_write_failure(testCase)
   % A non-skippable MERRA write failure must leave the completed MAR artifacts
   % linked from disk while preserving the manifest state that predated staging.
   assumePromicePresent(testCase);
   testCase.assumeTrue(strlength(testCase.TestData.mar) > 0 ...
      && strlength(testCase.TestData.merra) > 0, ...
      'MAR/MERRA fixtures not available; skipping manifest checkpoint guard.');
   root = testCase.TestData.root;

   % Seed an observation-only manifest and an unrelated skip record that the
   % first source checkpoint must retain byte-logically through ordinary merge.
   icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, output_root=string(root), ...
      build_forcing=false, overwrite=true);
   manifest_file = fullfile(root, 'eval', 'promice', 'manifest.json');
   prior = jsondecode(fileread(manifest_file));
   prior.skipped = struct('site', "prior_missing_site", ...
      'reason', "preserve this state");
   icemodel.verification.setup.writeManifest(manifest_file, prior);

   % A regular file at userdata/merra2 deterministically injects a later-source
   % mkdir failure after MAR has completed and checkpointed its own artifacts.
   userdata_outdir = fullfile(root, 'input', 'userdata');
   mkdir(userdata_outdir);
   writeCorruptFile(fullfile(userdata_outdir, 'merra2'));
   caught = MException.empty;
   try
      icemodel.verification.setup.stageRcmForcing( ...
         obs_manifest=string(manifest_file), ...
         forcing_sources=["mar", "merra"], ...
         mar_dir=testCase.TestData.mar, merra_dir=testCase.TestData.merra, ...
         dt_out="", overwrite=true);
   catch err
      caught = err;
   end

   testCase.verifyNotEmpty(caught);
   testCase.verifyTrue(contains(string(caught.message), "merra2"));
   persisted = jsondecode(fileread(manifest_file));
   c = persisted.cases(1);
   testCase.verifyTrue(isfield(c.colocation, 'promice'));
   testCase.verifyTrue(logical(c.colocation.mar.staged));
   testCase.verifyFalse(isfield(c.colocation, 'merra'));
   testCase.verifyEqual(numel(persisted.skipped), 1);
   testCase.verifyEqual(string(persisted.skipped.site), "prior_missing_site");

   % Both first-source artifact classes must exist at the paths now recorded by
   % the durable checkpoint, not merely as unreferenced files on disk.
   testCase.verifyTrue(all(isfile(fullfile(root, 'input', 'met', ...
      string(c.colocation.mar.met_files)))));
   testCase.verifyTrue(all(isfile(fullfile(root, 'input', 'userdata', ...
      string(c.colocation.mar.data_files)))));
end

%% Local helpers
function writeEmptyMat(filename)
   %WRITEEMPTYMAT Create a placeholder .mat file for filename-resolution tests.
   folder = fileparts(filename);
   if ~isfolder(folder)
      mkdir(folder)
   end
   placeholder = true;
   save(filename, 'placeholder')
end

function writeReusableRcmMat(filename, point, method)
   %WRITEREUSABLERCMMAT Create a metadata-bearing placeholder RCM artifact.
   if nargin < 3
      % Most cache-reuse tests exercise the public nearest-sampling default.
      method = "nearest";
   end
   folder = fileparts(filename);
   if ~isfolder(folder)
      mkdir(folder)
   end
   artifact_metadata = struct('sample_method', char(method), ...
      'lat_wgs84', point(1), 'lon_wgs84', point(2), ...
      'racmo_ice_mask_applied', true, ...
      'racmo_point_max_distance_m', 1);
   [~, name] = fileparts(filename);
   if startsWith(string(name), "met_")
      % Reuse tests select the exact saved met bytes, so the fixture must carry
      % a real structurally valid payload rather than metadata alone.
      stamps = regexp(name, '\d{8}', 'match');
      assert(numel(stamps) >= 2, 'reusable met fixture needs a window name')
      t1 = datetime(stamps{1}, InputFormat='yyyyMMdd', TimeZone='UTC');
      t2 = datetime(stamps{2}, InputFormat='yyyyMMdd', TimeZone='UTC');
      if t2 == t1
         % A one-day filename still needs two postings for validatemet.
         t2 = t1 + hours(1);
      end
      required = icemodel.forcing.helpers.metvariables();
      met = array2timetable(ones(2, numel(required)), ...
         RowTimes=[t1; t2], VariableNames=cellstr(required));
      met = icemodel.forcing.helpers.stampMetadata(met);
      met.Properties.UserData = artifact_metadata;
      save(filename, 'met', 'artifact_metadata')
   else
      % Data-only cache tests exercise discovery metadata, not payload loading.
      save(filename, 'artifact_metadata')
   end
end

function verifyLegReadiness(testCase, leg, met_file)
   %VERIFYLEGREADINESS Compare manifest diagnostics with selected saved bytes.
   [ready, reason, windows] = ...
      icemodel.verification.setup.metArtifactReadiness(string(met_file));
   testCase.verifyEqual(logical(leg.forcing_ready), logical(ready));
   testCase.verifyEqual(string(leg.forcing_ready_reason), string(reason));
   testCase.verifyEqual(leg.forcing_complete_windows, windows);
end

function writeCorruptFile(filename)
   %WRITECORRUPTFILE Write bytes that cannot be decoded as a NetCDF source.
   folder = fileparts(filename);
   if ~isfolder(folder)
      mkdir(folder)
   end
   fid = fopen(filename, 'w');
   cleanup = onCleanup(@() fclose(fid));
   fwrite(fid, uint8('not a NetCDF file'), 'uint8');
end

function verifyModisProvenance(testCase, loaded, variable)
   %VERIFYMODISPROVENANCE Check fresh physical data and synchronized metadata.
   T = loaded.(char(variable));
   metadata = loaded.artifact_metadata;
   testCase.verifyTrue(ismember("modis", string(T.Properties.VariableNames)));
   testCase.verifyGreaterThan(nnz(isfinite(T.modis)), 0);
   testCase.verifyEqual(string(metadata.modis_product), ...
      "GEUS Greenland Reflectivity 5km C6");
   testCase.verifyEqual(string(metadata.modis_status), "source_coverage");
   testCase.verifyEqual(double(metadata.modis_coverage_years(:))', 2012);
   testCase.verifyEqual(T.Properties.UserData, metadata);
end

function verifyNoModisProvenance(testCase, loaded, variable)
   %VERIFYNOMODISPROVENANCE Check the canonical explicit no-source encoding.
   T = loaded.(char(variable));
   metadata = loaded.artifact_metadata;
   testCase.verifyFalse(ismember("modis", string(T.Properties.VariableNames)));
   testCase.verifyEqual(string(metadata.modis_product), ...
      "GEUS Greenland Reflectivity 5km C6");
   testCase.verifyEqual(string(metadata.modis_status), "no_source_coverage");
   testCase.verifyEmpty(metadata.modis_coverage_years);
   testCase.verifyEqual(T.Properties.UserData, metadata);
end

function filenames = stagedCaseArtifacts(root, c)
   %STAGEDCASEARTIFACTS Resolve every artifact referenced by one PROMICE case.
   filenames = string(fullfile(root, 'eval', 'promice', c.evaluation_file));
   for src = ["promice", "mar", "merra", "racmo"]
      leg = c.colocation.(char(src));
      if isfield(leg, 'met_files')
         filenames = [filenames; fullfile(root, 'input', 'met', ...
            string(leg.met_files(:)))]; %#ok<AGROW>
      end
      if isfield(leg, 'data_files')
         filenames = [filenames; fullfile(root, 'input', 'userdata', ...
            string(leg.data_files(:)))]; %#ok<AGROW>
      end
   end
end

function bytes = fileBytes(filename)
   %FILEBYTES Read one staged binary artifact for no-churn assertions.
   fid = fopen(filename, 'r');
   cleanup = onCleanup(@() fclose(fid));
   bytes = fread(fid, Inf, '*uint8');
end

function p = firstWithData(candidates, hasData)
   %FIRSTWITHDATA First candidate dir that exists and satisfies hasData, else "".
   p = "";
   for c = candidates
      if isfolder(c) && hasData(c)
         p = c;
         return
      end
   end
end
