function tests = test_firn_staging_promice
   %TEST_FIRN_STAGING_PROMICE Verify the co-located firn staging driver.
   %
   % Exercises icemodel.verification.setup.importPromiceSites end to end:
   % the co-located PROMICE/MAR/MERRA/RACMO bundle is staged for a PROMICE
   % anchor site and the per-site manifest entry resolves. Reads the raw
   % multi-model sources from the S03 external-drive layout (or local
   % caches) and skips cleanly when any required source is unavailable.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Resolve every raw source the bundle needs; skip the whole file when
   % any is absent (the bundle is co-located, so a partial mount cannot
   % stage a faithful bundle).

   promice = firstWithData([ ...
      string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'promice')), ...
      "/Volumes/S03/DATA/greenland/geus/aws/l3"], ...
      @(p) ~isempty(dir(fullfile(p, "hour", "*_hour.nc"))));
   mar = firstWithData([ ...
      "/Volumes/S03/DATA/greenland/mar3p11/RUH2", ...
      string(fullfile(icemodel.getpath('data'), 'forcing', 'mar'))], ...
      @(p) ~isempty(dir(fullfile(p, "MARv3.11*.nc"))));
   merra = firstWithData([ ...
      "/Volumes/S03/DATA/merra2/1hrly/ncfiles", ...
      string(fullfile(icemodel.getpath('data'), 'forcing', 'merra2'))], ...
      @(p) ~isempty(dir(fullfile(p, "slv", "*_Nx.*.nc4*"))));
   racmo = firstWithData([ ...
      "/Volumes/S03/DATA/greenland/racmo2p3/surface", ...
      string(fullfile(icemodel.getpath('data'), 'forcing', 'racmo'))], ...
      @(p) ~isempty(dir(fullfile(p, "*.RACMO23p3_*.nc"))));

   testCase.assumeTrue(all([ ...
      strlength(promice) > 0, strlength(mar) > 0, ...
      strlength(merra) > 0, strlength(racmo) > 0]), ...
      'co-located firn sources not all available (S03 unmounted/spun down)');

   testCase.TestData.promice = promice;
   testCase.TestData.mar = mar;
   testCase.TestData.merra = merra;
   testCase.TestData.racmo = racmo;

   % Stage into a private eval / input root so the committed demo tree is
   % untouched by the test.
   tmp = tempname;
   mkdir(tmp)
   testCase.TestData.tmp = tmp;
   testCase.TestData.eval_root = fullfile(tmp, 'eval');
   testCase.TestData.input_root = fullfile(tmp, 'input');
   mkdir(testCase.TestData.eval_root)
   mkdir(fullfile(testCase.TestData.input_root, 'met'))
   mkdir(fullfile(testCase.TestData.input_root, 'userdata'))
end

function teardownOnce(testCase)
   if isfield(testCase.TestData, 'tmp') && isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
end

function test_colocated_bundle_and_manifest_resolve(testCase)
   % Stage a short single-site bundle and assert every model artifact and
   % the firn manifest entry resolve.

   % A short window inside the RACMO 2012-2015 archive keeps the test fast.
   manifest = stageKanm(testCase, "2013-06-01", "2013-06-30");

   % Family manifest shape.
   testCase.verifyEqual(string(manifest.dataset_family), "promice");
   testCase.verifyEqual(numel(manifest.cases), 1);
   testCase.verifyEmpty(manifest.skipped);

   c = manifest.cases(1);
   testCase.verifyEqual(string(c.case_id), "kanm");
   testCase.verifyEqual(string(c.case_type), "firn_observational");
   testCase.verifyEqual(string(c.site_id), "KAN_M");

   % surface_zone is single-sourced from promicesiteinfo (KAN_M is the upper
   % ablation zone) and must validate against the canonical namelist.
   testCase.verifyEqual(string(c.surface_zone), "upper_ablation");
   testCase.verifyTrue(ismember(string(c.surface_zone), ...
      icemodel.verification.namelists.surfacezone()));

   % Site location: WGS84 + EPSG:3413 recorded.
   testCase.verifyEqual(c.site_location.lat_wgs84, 67.067, 'AbsTol', 1e-2);
   testCase.verifyTrue(isfinite(c.site_location.x_epsg3413));
   testCase.verifyTrue(isfinite(c.site_location.y_epsg3413));

   % Co-located forcing: all four models present with the right kinds.
   cf = c.colocated_forcing;
   testCase.verifyTrue(all(isfield(cf, {'promice', 'mar', 'merra', 'racmo'})));
   testCase.verifyEqual(string(cf.racmo.kind), "point_data_smb_eval");
   testCase.verifyTrue(~isempty(cf.promice.met_files));
   testCase.verifyTrue(~isempty(cf.mar.met_files));
   testCase.verifyTrue(~isempty(cf.merra.met_files));
   testCase.verifyTrue(~isempty(cf.racmo.data_files));
end

function test_staged_files_exist_on_disk(testCase)
   % The bundle's met / userdata / eval files must be written to disk.

   stageKanm(testCase, "2013-06-01", "2013-06-30");

   met_dir = fullfile(testCase.TestData.input_root, 'met');
   ud_dir = fullfile(testCase.TestData.input_root, 'userdata');
   eval_dir = fullfile(testCase.TestData.eval_root, 'promice', 'kanm');

   testCase.verifyNotEmpty(dir(fullfile(met_dir, 'met_kanm_promice_*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(met_dir, 'met_kanm_mar_*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(met_dir, 'met_kanm_merra_*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(ud_dir, 'kanm_promice_*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(ud_dir, 'kanm_racmo_*.mat')));

   testCase.verifyTrue(isfile(fullfile(eval_dir, 'evaluation.mat')));
   testCase.verifyTrue(isfile(fullfile(eval_dir, 'reference.mat')));
   testCase.verifyTrue(isfile(fullfile(testCase.TestData.eval_root, ...
      'promice', 'manifest.json')));
end

function test_evaluation_and_reference_artifacts_load(testCase)
   % evaluation.mat (PROMICE targets) and reference.mat (RACMO reference)
   % round-trip and carry the expected timeseries payloads.

   stageKanm(testCase, "2013-06-01", "2013-06-30");
   eval_dir = fullfile(testCase.TestData.eval_root, 'promice', 'kanm');

   ev = load(fullfile(eval_dir, 'evaluation.mat'), 'targets');
   testCase.verifyEqual(string(ev.targets.format), "timeseries");
   testCase.verifyTrue(istimetable(ev.targets.data));

   rf = load(fullfile(eval_dir, 'reference.mat'), 'reference');
   testCase.verifyEqual(string(rf.reference.metadata.reference_kind), ...
      "colocated_rcm");
   testCase.verifyTrue(istimetable(rf.reference.data));
end

function test_manifest_json_is_valid(testCase)
   % The written manifest.json must parse back to a struct with the firn
   % case schema fields.

   stageKanm(testCase, "2013-06-01", "2013-06-30");
   manifest_file = fullfile(testCase.TestData.eval_root, ...
      'promice', 'manifest.json');

   decoded = jsondecode(fileread(manifest_file));
   testCase.verifyEqual(string(decoded.dataset_family), "promice");

   needed = icemodel.verification.setup.firnCaseManifestFieldNames();
   have = string(fieldnames(decoded.cases));
   testCase.verifyTrue(all(ismember(needed, have)));
end

function test_data_gated_site_is_skipped_not_fabricated(testCase)
   % A nonexistent station id must be recorded in manifest.skipped with a
   % reason and must NOT create a fabricated case folder.

   manifest = icemodel.verification.setup.importPromiceSites( ...
      sites="ZZZ_NOPE", ...
      promice_dir=testCase.TestData.promice, ...
      mar_dir=testCase.TestData.mar, ...
      merra_dir=testCase.TestData.merra, ...
      racmo_dir=testCase.TestData.racmo, ...
      startdate="2013-06-01", enddate="2013-06-30", ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      icemodel_config_casename="", overwrite=true);

   testCase.verifyEqual(numel(manifest.cases), 0);
   testCase.verifyEqual(numel(manifest.skipped), 1);
   testCase.verifyEqual(string(manifest.skipped(1).site), "ZZZ_NOPE");
   testCase.verifyTrue(strlength(manifest.skipped(1).reason) > 0);
   testCase.verifyFalse(isfolder(fullfile(testCase.TestData.eval_root, ...
      'promice', 'zzznope')));
end

%% Local helpers
function manifest = stageKanm(testCase, startdate, enddate)
   %STAGEKANM Stage the KAN_M bundle into the private test roots.
   manifest = icemodel.verification.setup.importPromiceSites( ...
      sites="KAN_M", ...
      promice_dir=testCase.TestData.promice, ...
      mar_dir=testCase.TestData.mar, ...
      merra_dir=testCase.TestData.merra, ...
      racmo_dir=testCase.TestData.racmo, ...
      startdate=startdate, enddate=enddate, ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      icemodel_config_casename="", overwrite=true);
end

function path = firstWithData(candidates, hasdata)
   %FIRSTWITHDATA First candidate folder that actually contains data files.
   path = "";
   for c = candidates
      if isfolder(c) && hasdata(c)
         path = c;
         return
      end
   end
end
