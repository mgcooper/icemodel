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
   %   - importPromiceSites(build_forcing=false) imports observations + station
   %     met and writes the manifest, staging NO gridded-RCM files.
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
   % build_forcing=false imports the PROMICE observations + station met and
   % writes the manifest, but stages NO gridded-RCM (mar/merra/racmo) files or
   % colocation legs - the observation import is decoupled from RCM forcing.
   assumePromicePresent(testCase);
   root = testCase.TestData.root;

   manifest = icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, ...
      output_root=string(root), build_forcing=false, overwrite=true);

   testCase.verifyEqual(numel(manifest.cases), 1);
   c = manifest.cases(1);

   % The eval target is staged; the manifest references it.
   testCase.verifyTrue(isfile(fullfile(root, 'eval', 'promice', ...
      'kanl', 'observations.mat')));
   testCase.verifyEqual(string(c.evaluation_file), "kanl/observations.mat");

   % The PROMICE leg is present; NO gridded-RCM legs were staged.
   testCase.verifyTrue(isfield(c.colocation, 'promice'));
   testCase.verifyFalse(isfield(c.colocation, 'mar'));
   testCase.verifyFalse(isfield(c.colocation, 'merra'));
   testCase.verifyFalse(isfield(c.colocation, 'racmo'));

   % NO gridded-RCM files on disk (the PROMICE met/userdata may exist).
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'met', 'mar', '*.mat')));
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'met', 'merra', '*.mat')));
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'userdata', 'mar', '*.mat')));
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'userdata', 'racmo', '*.mat')));
end

function test_manifest_convenience_skips_absent_sources(testCase)
   % stageRcmForcing manifest-convenience mode: after an observations-only
   % import, it resolves the legs from the staged manifest and, when the RCM
   % sources are absent (bogus dirs), degrades EVERY source to a skip-with-reason
   % WITHOUT throwing - validating the manifest-mode plumbing (read cases ->
   % resolve points/windows -> stage -> merge -> persist) off the fail-early gate.
   assumePromicePresent(testCase);
   root = testCase.TestData.root;

   icemodel.verification.setup.importPromiceSites( ...
      sites=testCase.TestData.site, promice_dir=testCase.TestData.promice_dir, ...
      startdate=testCase.TestData.startdate, ...
      enddate=testCase.TestData.enddate, ...
      output_root=string(root), build_forcing=false, overwrite=true);

   manifest_file = fullfile(root, 'eval', 'promice', 'manifest.json');
   bogus = fullfile(root, 'no_such_rcm_source');

   manifest = icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=string(manifest_file), manifest_file=string(manifest_file), ...
      models=["mar", "merra", "racmo"], ...
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

   % Nothing was written for the absent sources.
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'met', 'mar', '*.mat')));
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'userdata', 'racmo', '*.mat')));
end

function test_existing_enclosing_files_skip_rebuild(testCase)
   % A previously staged full-period file should satisfy a narrower requested
   % leg. This guards the write-side equivalent of createMetFileNames/loadmet
   % enclosing-file resolution without needing the external RCM archives.
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');
   alias = "kanl";

   writeEmptyMat(fullfile(met_outdir, 'mar', ...
      'met_kanl_mar_20100101_20200101_1hr.mat'));
   writeEmptyMat(fullfile(met_outdir, 'merra', ...
      'met_kanl_merra_20100101_20200101_1hr.mat'));
   writeEmptyMat(fullfile(userdata_outdir, 'mar', ...
      'kanl_mar_20100101_20200101.mat'));
   writeEmptyMat(fullfile(userdata_outdir, 'merra', ...
      'kanl_merra_20100101_20200101.mat'));
   writeEmptyMat(fullfile(userdata_outdir, 'racmo', ...
      'kanl_racmo_20100101_20200101.mat'));

   L = struct('staged', true, 'years', 2013, ...
      'start', icemodel.verification.setup.ensureUtc('2013-01-01'), ...
      'end', icemodel.verification.setup.ensureUtc('2013-12-31 23:00:00'), ...
      'reason', "");
   legspec = struct('alias', alias, 'mar', L, 'merra', L, 'racmo', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      [67.0, -50.0], legspec=legspec, ...
      models=["mar", "merra", "racmo"], ...
      met_outdir=string(met_outdir), ...
      userdata_outdir=string(userdata_outdir), ...
      mar_dir=string(fullfile(root, 'absent_mar')), ...
      merra_dir=string(fullfile(root, 'absent_merra')), ...
      racmo_dir=string(fullfile(root, 'absent_racmo')));

   c = colocation{1};
   testCase.verifyTrue(logical(c.mar.staged));
   testCase.verifyEqual(string(c.mar.met_files), ...
      "mar/met_kanl_mar_20100101_20200101_1hr.mat");
   testCase.verifyEqual(string(c.mar.data_files), ...
      "mar/kanl_mar_20100101_20200101.mat");
   testCase.verifyTrue(logical(c.merra.staged));
   testCase.verifyEqual(string(c.merra.met_files), ...
      "merra/met_kanl_merra_20100101_20200101_1hr.mat");
   testCase.verifyTrue(logical(c.racmo.staged));
   testCase.verifyEqual(string(c.racmo.data_files), ...
      "racmo/kanl_racmo_20100101_20200101.mat");

   % No requested-window duplicates were created because every required output
   % was already bracketed by a wider staged file.
   testCase.verifyFalse(isfile(fullfile(met_outdir, 'mar', ...
      'met_kanl_mar_20130101_20131231_1hr.mat')));
   testCase.verifyFalse(isfile(fullfile(userdata_outdir, 'racmo', ...
      'kanl_racmo_20130101_20131231.mat')));
end

function test_existing_met_only_does_not_block_userdata_build(testCase)
   % A MAR/MERRA leg is complete only when both met and userdata bracket the
   % requested window. A pre-refactor met-only artifact must not block creation
   % of the missing userdata file.
   root = testCase.TestData.root;
   met_outdir = fullfile(root, 'input', 'met');
   userdata_outdir = fullfile(root, 'input', 'userdata');
   alias = "kanl";

   writeEmptyMat(fullfile(met_outdir, 'mar', ...
      'met_kanl_mar_20100101_20200101_1hr.mat'));

   L = struct('staged', true, 'years', 2013, ...
      'start', icemodel.verification.setup.ensureUtc('2013-01-01'), ...
      'end', icemodel.verification.setup.ensureUtc('2013-12-31 23:00:00'), ...
      'reason', "");
   legspec = struct('alias', alias, 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      [67.0, -50.0], legspec=legspec, models="mar", ...
      met_outdir=string(met_outdir), ...
      userdata_outdir=string(userdata_outdir), ...
      mar_dir=string(fullfile(root, 'absent_mar')));

   testCase.verifyFalse(logical(colocation{1}.mar.staged));
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
      build_forcing=true, ...
      mar_dir=testCase.TestData.mar, merra_dir=testCase.TestData.merra, ...
      racmo_dir=string(testCase.TestData.racmo), overwrite=true);

   c = manifest.cases(1);

   % MAR + MERRA: a met file AND a Data (userdata) file (the fix).
   for src = ["mar", "merra"]
      leg = c.colocation.(char(src));
      testCase.verifyTrue(logical(leg.staged), ...
         sprintf('%s leg should be staged', src));
      testCase.verifyTrue(isfield(leg, 'met_files') && ~isempty(leg.met_files), ...
         sprintf('%s must write a met file', src));
      testCase.verifyTrue(isfield(leg, 'data_files') && ~isempty(leg.data_files), ...
         sprintf('%s must ALSO write a userdata (Data) file', src));
      testCase.verifyNotEmpty(dir(fullfile(root, 'input', 'met', char(src), '*.mat')));
      testCase.verifyNotEmpty( ...
         dir(fullfile(root, 'input', 'userdata', char(src), '*.mat')));
   end

   % RACMO: Data only, no met.
   testCase.verifyTrue(logical(c.colocation.racmo.staged));
   testCase.verifyTrue(isfield(c.colocation.racmo, 'data_files') ...
      && ~isempty(c.colocation.racmo.data_files));
   testCase.verifyFalse(isfield(c.colocation.racmo, 'met_files'));
   testCase.verifyNotEmpty( ...
      dir(fullfile(root, 'input', 'userdata', 'racmo', '*.mat')));
   testCase.verifyEmpty(dir(fullfile(root, 'input', 'met', 'racmo', '*.mat')));
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
      build_forcing=true, ...
      mar_dir=testCase.TestData.mar, merra_dir=testCase.TestData.merra, ...
      racmo_dir=string(bogus), overwrite=true);

   c = manifest.cases(1);
   testCase.verifyTrue(logical(c.colocation.mar.staged));
   testCase.verifyTrue(logical(c.colocation.merra.staged));
   testCase.verifyFalse(logical(c.colocation.racmo.staged));
   testCase.verifyNotEmpty(dir(fullfile(root, 'input', 'met', 'mar', '*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(root, 'input', 'userdata', 'merra', '*.mat')));
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
