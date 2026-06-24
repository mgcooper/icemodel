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

   % surface_zone is single-sourced from promicesiteinfo (KAN_M is in the
   % ablation zone) and must validate against the canonical namelist.
   testCase.verifyEqual(string(c.surface_zone), "ablation");
   testCase.verifyTrue(ismember(string(c.surface_zone), ...
      icemodel.verification.namelists.surfacezone()));

   % permafrost_zone is single-sourced from promicesiteinfo: KAN_M is an
   % ice-sheet anchor (not permafrost ground) so it is "none" and validates
   % against the canonical permafrost namelist.
   testCase.verifyEqual(string(c.permafrost_zone), "none");
   testCase.verifyTrue(ismember(string(c.permafrost_zone), ...
      icemodel.verification.namelists.permafrostzone()));

   % eval_target: KAN_M exercises seasonal snow + bare ice.
   testCase.verifyEqual(sort(string(c.eval_target(:))), ...
      sort(["bare_ice"; "seasonal_snow"]));

   % Metadata-only source records (ids, not bundled data).
   testCase.verifyTrue(all(ismember(["promice", "mar", "merra"], ...
      string(c.forcing_sources))));
   testCase.verifyTrue(all(ismember(["promice_obs", "racmo"], ...
      string(c.eval_sources))));

   % Site location: WGS84 + EPSG:3413 recorded.
   testCase.verifyEqual(c.site_location.lat_wgs84, 67.067, 'AbsTol', 1e-2);
   testCase.verifyTrue(isfinite(c.site_location.x_epsg3413));
   testCase.verifyTrue(isfinite(c.site_location.y_epsg3413));

   % Colocation metadata: all four model legs present with the right kinds.
   cf = c.colocation;
   testCase.verifyTrue(all(isfield(cf, {'promice', 'mar', 'merra', 'racmo'})));
   testCase.verifyEqual(string(cf.racmo.kind), "point_data_smb_eval");
   testCase.verifyTrue(~isempty(cf.promice.met_files));
   testCase.verifyTrue(~isempty(cf.mar.met_files));
   testCase.verifyTrue(~isempty(cf.merra.met_files));
   testCase.verifyTrue(~isempty(cf.racmo.data_files));
end

function test_staged_files_exist_on_disk(testCase)
   % The individual met / userdata files must be written to disk. The eval is a
   % data-only observations.mat; NO evaluation.mat / reference.mat bundle (the
   % forcing/reference side is never bundled with the eval target).

   stageKanm(testCase, "2013-06-01", "2013-06-30");

   met_dir = fullfile(testCase.TestData.input_root, 'met');
   ud_dir = fullfile(testCase.TestData.input_root, 'userdata');
   eval_dir = fullfile(testCase.TestData.eval_root, 'promice', 'kanm');

   % Forcing/userdata stage into per-source subfolders (met/<source>/, userdata/<source>/).
   testCase.verifyNotEmpty(dir(fullfile(met_dir, 'promice', 'met_kanm_promice_*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(met_dir, 'mar', 'met_kanm_mar_*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(met_dir, 'merra', 'met_kanm_merra_*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(ud_dir, 'promice', 'kanm_promice_*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(ud_dir, 'racmo', 'kanm_racmo_*.mat')));

   % Eval is the data-only observations.mat (no bundled evaluation/reference copy).
   testCase.verifyTrue(isfile(fullfile(eval_dir, 'observations.mat')));
   testCase.verifyFalse(isfile(fullfile(eval_dir, 'evaluation.mat')));
   testCase.verifyFalse(isfile(fullfile(eval_dir, 'reference.mat')));

   testCase.verifyTrue(isfile(fullfile(testCase.TestData.eval_root, ...
      'promice', 'manifest.json')));
end

function test_colocated_data_reconstitutes_from_individual_files(testCase)
   % The eval target (PROMICE obs) and the RCM reference (RACMO Data) must
   % reconstitute on demand from the individual per-year userdata files the
   % forcing-agnostic manifest declares - no committed forcing bundle needed.

   manifest = stageKanm(testCase, "2013-06-01", "2013-06-30");
   c = manifest.cases(1);
   c.manifest_path = fullfile(testCase.TestData.eval_root, 'promice', ...
      'manifest.json');

   promice = icemodel.verification.helpers.loadColocatedData(c, "promice", ...
      input_data_root=testCase.TestData.input_root);
   testCase.verifyEqual(string(promice.format), "timeseries");
   testCase.verifyTrue(istimetable(promice.data));
   testCase.verifyGreaterThan(height(promice.data), 0);

   racmo = icemodel.verification.helpers.loadColocatedData(c, "racmo", ...
      input_data_root=testCase.TestData.input_root);
   testCase.verifyEqual(string(racmo.format), "timeseries");
   testCase.verifyTrue(istimetable(racmo.data));
   testCase.verifyGreaterThan(height(racmo.data), 0);
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
   % reason and must NOT create a fabricated case folder. Under merge, any
   % previously staged cases are preserved, so assert the gated site itself is
   % skipped and never fabricated rather than that the manifest is empty.

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

   % No fabricated case for the gated site.
   case_ids = string(arrayfun(@(c) string(c.case_id), manifest.cases));
   testCase.verifyFalse(ismember("zzznope", case_ids));

   % The gated site is recorded in skipped with a reason.
   skip_sites = string(arrayfun(@(s) string(s.site), manifest.skipped));
   testCase.verifyTrue(ismember("ZZZ_NOPE", skip_sites));
   idx = find(skip_sites == "ZZZ_NOPE", 1);
   testCase.verifyTrue(strlength(manifest.skipped(idx).reason) > 0);

   testCase.verifyFalse(isfolder(fullfile(testCase.TestData.eval_root, ...
      'promice', 'zzznope')));
end

function test_staging_second_site_does_not_churn_first(testCase)
   % Staging a SECOND site into a family root that already holds a FIRST site
   % must ADD the second case and leave the first site's case entry + its
   % staged files byte for byte unchanged (the KAN no-churn guarantee). This is
   % the file-level counterpart to test_firn_manifest_merge.

   % Stage KAN_L first (a short window keeps it fast).
   stageSite(testCase, "KAN_L", "2013-06-01", "2013-06-30");

   manifest_file = fullfile(testCase.TestData.eval_root, 'promice', ...
      'manifest.json');
   kanl_manifest_before = fileread(manifest_file);

   % Fingerprint the always-staged KAN_L artifacts: the data-only eval bundle
   % (eval/promice/kanl/observations.mat) and the promice forcing/Data files.
   % writemet/writeuserdata stage into a per-source subfolder (met/<source>/,
   % userdata/<source>/), so glob there. The PROMICE MET leg is independently
   % skippable (a station/window with no longwave sensor degrades to a skipped
   % met leg while the eval obs still stage - #15 / RR1), so it is folded in
   % only when present rather than required.
   eval_case_dir = fullfile(testCase.TestData.eval_root, 'promice', 'kanl');
   met_dir = fullfile(testCase.TestData.input_root, 'met', 'promice');
   ud_dir = fullfile(testCase.TestData.input_root, 'userdata', 'promice');
   kanl_eval = dir(fullfile(eval_case_dir, 'observations.mat'));
   kanl_met = dir(fullfile(met_dir, 'met_kanl_*'));
   kanl_ud = dir(fullfile(ud_dir, 'kanl_*'));
   testCase.assertNotEmpty(kanl_ud, ...
      'KAN_L staged no promice userdata (eval leg failed)');
   before = fileFingerprints([kanl_eval; kanl_met; kanl_ud], ...
      {eval_case_dir, met_dir, ud_dir});

   % Snapshot the KAN_L case JSON region from the merged manifest.
   m_before = jsondecode(kanl_manifest_before);
   kanl_case_before = encodeCaseById(m_before, "kanl");

   % Stage KAN_M into the SAME roots (merge, default).
   stageSite(testCase, "KAN_M", "2013-06-01", "2013-06-30");

   m_after = jsondecode(fileread(manifest_file));
   ids = string(arrayfun(@(c) string(c.case_id), m_after.cases));
   testCase.verifyTrue(all(ismember(["kanl", "kanm"], ids)));

   % KAN_L case entry is byte-identical after KAN_M was added.
   testCase.verifyEqual(encodeCaseById(m_after, "kanl"), kanl_case_before, ...
      'KAN_L case churned when KAN_M was staged');

   % KAN_L's staged files are byte-identical (size + content hash).
   after = fileFingerprints([dir(fullfile(eval_case_dir, 'observations.mat')); ...
      dir(fullfile(met_dir, 'met_kanl_*')); ...
      dir(fullfile(ud_dir, 'kanl_*'))], {eval_case_dir, met_dir, ud_dir});
   testCase.verifyEqual(after, before, ...
      'KAN_L staged files churned when KAN_M was staged');
end

function test_missing_rcm_legs_still_yield_a_case(testCase)
   % #15 / RR1 regression: PROMICE met+eval must NEVER be gated by RCM
   % coverage. A site whose MAR/MERRA/RACMO source dirs are absent (the
   % builders THROW) must still stage a PROMICE case - the RCM legs degrade
   % to skipped legs in colocation, and the site is NOT wholesale skipped.
   % This guards the staging bug where an erroring RCM builder, called inside
   % the per-site try/catch, skipped the entire site.

   % Resolve PROMICE on its own; this test does not need the RCM sources, so
   % it must run even when S03 is unmounted (unlike the colocated tests).
   promice = firstWithData([ ...
      string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'promice')), ...
      "/Volumes/S03/DATA/greenland/geus/aws/l3"], ...
      @(p) ~isempty(dir(fullfile(p, "hour", "*_hour.nc"))));
   testCase.assumeTrue(strlength(promice) > 0, ...
      'PROMICE source not available (no *_hour.nc on disk)');

   % Private roots so nothing committed is touched.
   tmp = tempname;
   mkdir(tmp)
   cleanup = onCleanup(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');
   mkdir(eval_root)
   mkdir(fullfile(input_root, 'met'))
   mkdir(fullfile(input_root, 'userdata'))

   % Point every RCM dir at a path that cannot exist, forcing each builder to
   % throw the "source directory not found" error - the exact failure mode
   % that previously skipped the whole site.
   missing = fullfile(tmp, 'no_such_rcm_source');

   manifest = icemodel.verification.setup.importPromiceSites( ...
      sites="KAN_M", ...
      promice_dir=promice, ...
      mar_dir=missing, merra_dir=missing, racmo_dir=missing, ...
      startdate="2013-06-01", enddate="2013-06-30", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      icemodel_config_casename="", overwrite=true);

   % The site STAGED as a case (not skipped) despite every RCM leg failing.
   testCase.verifyEqual(numel(manifest.cases), 1, ...
      'PROMICE site was wrongly skipped when its RCM legs were absent');
   c = manifest.cases(1);
   testCase.verifyEqual(string(c.case_id), "kanm");

   % KAN_M is NOT in skipped - only its RCM legs are.
   if ~isempty(manifest.skipped)
      skip_sites = string(arrayfun(@(s) string(s.site), manifest.skipped));
      testCase.verifyFalse(ismember("KAN_M", skip_sites));
   end

   % PROMICE leg fully staged: it remains the forcing + eval source.
   testCase.verifyTrue(ismember("promice", string(c.forcing_sources)));
   testCase.verifyTrue(ismember("promice_obs", string(c.eval_sources)));

   % Each RCM leg recorded as a skipped leg with a reason, never fabricated.
   cf = c.colocation;
   testCase.verifyFalse(cf.mar.staged);
   testCase.verifyFalse(cf.merra.staged);
   testCase.verifyFalse(cf.racmo.staged);
   testCase.verifyGreaterThan(strlength(string(cf.mar.reason)), 0);

   % RCM legs are absent from the forcing/eval source lists (no fabrication).
   testCase.verifyFalse(ismember("mar", string(c.forcing_sources)));
   testCase.verifyFalse(ismember("merra", string(c.forcing_sources)));
   testCase.verifyFalse(ismember("racmo", string(c.eval_sources)));
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

function manifest = stageSite(testCase, site, startdate, enddate)
   %STAGESITE Stage one site's bundle into the private test roots (merge).
   manifest = icemodel.verification.setup.importPromiceSites( ...
      sites=site, ...
      promice_dir=testCase.TestData.promice, ...
      mar_dir=testCase.TestData.mar, ...
      merra_dir=testCase.TestData.merra, ...
      racmo_dir=testCase.TestData.racmo, ...
      startdate=startdate, enddate=enddate, ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      icemodel_config_casename="", overwrite=true);
end

function fp = fileFingerprints(entries, ~)
   %FILEFINGERPRINTS Map each file to {bytes, checksum} for a churn comparison.
   fp = struct();
   for k = 1:numel(entries)
      e = entries(k);
      f = fullfile(e.folder, e.name);
      key = matlab.lang.makeValidName(e.name);
      fid = fopen(f, 'r');
      raw = fread(fid, Inf, '*uint8');
      fclose(fid);
      fp.(key) = struct('bytes', e.bytes, 'sum', sum(double(raw)));
   end
end

function s = encodeCaseById(manifest, case_id)
   %ENCODECASEBYID Stable JSON encoding of the case with the given id.
   ids = string(arrayfun(@(c) string(c.case_id), manifest.cases));
   s = jsonencode(manifest.cases(ids == case_id), PrettyPrint=true);
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
