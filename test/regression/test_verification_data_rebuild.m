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
   % Skip when the Laugh-Tests checkout is not available.
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
   eval_path = fullfile(eval_root, 'laugh_tests', ...
      'colbeck1976', 'evaluation.mat');
   reference_path = fullfile(eval_root, 'laugh_tests', ...
      'colbeck1976', 'reference.mat');
   verifyTrue(testCase, exist(eval_path, 'file') == 2);
   loaded = load(eval_path, 'targets');
   verifyTrue(testCase, isfield(loaded.targets, 'numerical_summa'));
   verifyTrue(testCase, isfield(loaded.targets, 'analytical_clark2017'));

   % The two-source Laugh bundle is atomic, but an identical ordinary import
   % still reuses both current artifacts byte-for-byte.
   original = {fileBytes(eval_path), fileBytes(reference_path)};
   icemodel.verification.setup.importLaughTests(src, ...
      evaluation_data_root=string(eval_root), overwrite=false);
   verifyEqual(testCase, ...
      {fileBytes(eval_path), fileBytes(reference_path)}, original);

   % Both bundles expose the three Colbeck experiments.
   for src_name = ["numerical_summa", "analytical_clark2017"]
      bundle = loaded.targets.(char(src_name));
      verifyEqual(testCase, ...
         sort(string(fieldnames(bundle.experiments))), ...
         ["exp1"; "exp2"; "exp3"], ...
         sprintf('%s missing expected experiments', src_name));
   end
end

function test_laugh_tests_output_root_stages_eval_tree(testCase)
   % output_root should isolate Laugh-Tests eval artifacts like other importers.
   src = icemodel.verification.setup.fetchLaughTests( ...
      strict=false, silent=true);
   has_source = laughTestsCheckoutComplete(src);
   testCase.assumeTrue(has_source, ...
      sprintf('Laugh-Tests source cache not present at %s', src));

   output_root = fullfile(testCase.TestData.tmp, 'laugh-output-root');
   manifest = icemodel.verification.setup.importLaughTests(src, ...
      output_root=string(output_root), overwrite=true);

   verifyEqual(testCase, manifest.dataset_family, "laugh_tests");
   verifyTrue(testCase, isfile(fullfile(output_root, 'eval', ...
      'laugh_tests', 'manifest.json')));
end

function test_laugh_tests_dry_run_does_not_write_staging_tree(testCase)
   % Laugh-Tests dry runs should use the shared non-mutating import path without
   % reading a source checkout.
   src = fullfile(testCase.TestData.tmp, 'missing-laugh-tests');

   output_root = fullfile(testCase.TestData.tmp, 'laugh-dry-run');
   manifest = icemodel.verification.setup.importLaughTests(src, ...
      output_root=string(output_root), dry_run=true);

   verifyEqual(testCase, string(manifest.dataset_family), "laugh_tests");
   verifyEqual(testCase, string(manifest.cases.case_id), "colbeck1976");
   verifyFalse(testCase, isfolder(src));
   verifyFalse(testCase, isfolder(fullfile(output_root, 'eval')));
end

function test_laugh_tests_skip_missing_reaches_case_staging(testCase)
   % skip_missing=true should use the fetch validator and then record a skipped
   % case through the shared importer path instead of throwing.
   src = fullfile(testCase.TestData.tmp, 'empty-laugh-cache');
   eval_root = fullfile(testCase.TestData.tmp, 'empty-laugh-eval');
   mkdir(src)

   manifest = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.importLaughTests(src, ...
      evaluation_data_root=string(eval_root), skip_missing=true), ...
      'icemodel:verification:importLaughTests:caseSkipped');

   verifyEmpty(testCase, manifest.cases);
   verifyEqual(testCase, string(manifest.skipped.site), "colbeck1976");
   verifyFalse(testCase, isfolder(fullfile(eval_root, ...
      'laugh_tests', 'colbeck1976')));
end

function test_laugh_tests_invalid_case_is_side_effect_free(testCase)
   % Case membership is an argument-boundary contract in dry and real modes.
   source_dir = fullfile(testCase.TestData.tmp, 'invalid-laugh-source');
   for dry_run = [true, false]
      output_root = fullfile(testCase.TestData.tmp, ...
         "invalid-laugh-output-" + string(dry_run));
      verifyError(testCase, @() ...
         icemodel.verification.setup.importLaughTests(source_dir, ...
         case_id="unknown_case", output_root=output_root, dry_run=dry_run), ...
         'icemodel:verification:validators:mustBeLaughTestCase');
      verifyFalse(testCase, isfolder(output_root));
      verifyFalse(testCase, isfolder(source_dir));
   end
end

function test_rebuild_esm_snowmip_smoke_sites(testCase)
   src = icemodel.verification.setup.fetchEsmSnowmip( ...
      strict=false, silent=true);
   has_source = esmSnowmipCacheComplete(src);
   testCase.assumeTrue(has_source, ...
      sprintf('ESM-SnowMIP source cache not present at %s', src));

   eval_root = fullfile(testCase.TestData.tmp, 'eval');
   input_root = fullfile(testCase.TestData.tmp, 'input');
   manifest = icemodel.verification.setup.importEsmSnowmip(src, ...
      evaluation_data_root=string(eval_root), ...
      input_data_root=string(input_root), overwrite=true);

   verifyEqual(testCase, manifest.dataset_family, "esm_snowmip");
   verifyTrue(testCase, all(string({manifest.cases.native_timestep}) == "15m"));

   for case_id = ["cdp", "wfj"]
      case_dir = fullfile(eval_root, 'esm_snowmip', char(case_id));
      obs_path = fullfile(case_dir, 'observations.mat');
      verifyTrue(testCase, exist(obs_path, 'file') == 2, ...
         sprintf('%s observations.mat missing after rebuild', case_id));
      loaded = load(obs_path, 'targets');
      verifyEqual(testCase, loaded.targets.format, 'timeseries');
      verifyTrue(testCase, ...
         istimetable(loaded.targets.data) || ...
         isstruct(loaded.targets.data));
      if isfield(loaded.targets.metadata, 'soil_depths_m')
         soil_depths = double(loaded.targets.metadata.soil_depths_m(:));
         verifyEqual(testCase, soil_depths, round(soil_depths, 4), ...
            'AbsTol', 1e-12);
      end

      % Metadata-only: the redundant smoke reference.mat is no longer written.
      verifyTrue(testCase, ...
         exist(fullfile(case_dir, 'reference.mat'), 'file') == 0, ...
         sprintf('%s reference.mat should not be written (metadata-only)', ...
         case_id));
   end

   met_files = dir(fullfile(input_root, 'met', 'esm_snowmip', ...
      'met_cdp_esm_snowmip_*_15m.mat'));
   verifyNotEmpty(testCase, met_files);
   met_bundle = load(fullfile(met_files(1).folder, met_files(1).name), ...
      'artifact_metadata', 'met');
   verifyEqual(testCase, string(met_bundle.artifact_metadata.sitename), ...
      "cdp");
   verifyEqual(testCase, size(met_bundle.artifact_metadata.met_file, 2), 1);
   verifyEqual(testCase, seconds(median(diff(met_bundle.met.Time))), 900);

   % The manifest cadence must drive the standard runtime resolver to the
   % staged per-source 15-minute file, not synthesize the legacy flat 1hr name.
   cdp = manifest.cases(string({manifest.cases.case_id}) == "cdp");
   cdp.dataset_family = manifest.dataset_family;
   cdp.input_data_root = string(input_root);
   opts = icemodel.test.helpers.setModelOptsForCase(cdp);
   verifyEqual(testCase, opts.dt, 900);
   verifyTrue(testCase, isfile(opts.metfname{1}));

   % ESM observation/forcing source pairs remain atomic, while repeated
   % ordinary imports preserve the already-current staged bytes.
   obs_path = fullfile(eval_root, 'esm_snowmip', 'cdp', 'observations.mat');
   met_path = fullfile(met_files(1).folder, met_files(1).name);
   original = {fileBytes(obs_path), fileBytes(met_path)};
   icemodel.verification.setup.importEsmSnowmip(src, ...
      case_ids="cdp", evaluation_data_root=string(eval_root), ...
      input_data_root=string(input_root), overwrite=false);
   verifyEqual(testCase, {fileBytes(obs_path), fileBytes(met_path)}, original);
end

function test_esm_builders_enforce_paired_utc_window(testCase)
   % Both direct builders reject malformed pairs before source discovery.
   missing_source = fullfile(testCase.TestData.tmp, 'missing-esm-window-source');
   error_id = 'icemodel:internal:pairedWindow:invalidWindow';
   missing_builders = { ...
      @(a, b) icemodel.verification.setup.buildEsmSnowmipForcing( ...
      "cdp", source_dir=missing_source, startdate=a, enddate=b)
      @(a, b) icemodel.verification.setup.buildEsmSnowmipObservations( ...
      "cdp", source_dir=missing_source, startdate=a, enddate=b)};
   for k = 1:numel(missing_builders)
      build = missing_builders{k};
      verifyError(testCase, @() build("2000-01-01", ""), error_id);
      verifyError(testCase, @() build("", "2000-01-02"), error_id);
      verifyError(testCase, @() build("2000-01-02", "2000-01-01"), error_id);
   end
   verifyFalse(testCase, isfolder(missing_source));

   % With the optional cache present, zoned bounds preserve their instants and
   % every returned source axis is represented in UTC.
   src = icemodel.verification.setup.fetchEsmSnowmip( ...
      strict=false, silent=true);
   testCase.assumeTrue(esmSnowmipCacheComplete(src), ...
      sprintf('ESM-SnowMIP source cache not present at %s', src));
   builders = { ...
      @(a, b) icemodel.verification.setup.buildEsmSnowmipForcing( ...
      "cdp", source_dir=src, startdate=a, enddate=b)
      @(a, b) icemodel.verification.setup.buildEsmSnowmipObservations( ...
      "cdp", source_dir=src, startdate=a, enddate=b)};
   for k = 1:numel(builders)
      full_product = builders{k}("", "");
      expected = full_product.Time(1:min(3, height(full_product)));
      window_start = expected(1);
      window_end = expected(end);
      window_start.TimeZone = 'America/New_York';
      window_end.TimeZone = 'America/New_York';
      subset = builders{k}(window_start, window_end);
      verifyEqual(testCase, subset.Time.TimeZone, 'UTC');
      verifyEqual(testCase, subset.Time, expected);
   end
end

function test_esm_soil_temperature_preserves_depth_time_orientation(testCase)
   % Each canonical soil-temperature column must preserve one native tsl row.
   source_dir = fullfile(testCase.TestData.tmp, 'esm-tsl-orientation');
   icemodel.helpers.ensureDirExists(source_dir);
   obsfile = fullfile(source_dir, 'obs_insitu_cdp_orientation.nc');
   [source_time, source_tsl] = writeEsmSoilTemperatureFixture(obsfile);

   % Exercise an interior time window so both depth and time indexing are
   % proven independently of the full-record dimensions.
   window = 2:4;
   observations = ...
      icemodel.verification.setup.buildEsmSnowmipObservations( ...
      "cdp", source_dir=source_dir, ...
      startdate=source_time(window(1)), enddate=source_time(window(end)));

   verifyEqual(testCase, observations.Time, source_time(window));
   for depth_index = 1:size(source_tsl, 1)
      variable_name = "soil_temp_" + string(depth_index) + "_C";
      verifyEqual(testCase, observations.(variable_name), ...
         source_tsl(depth_index, window)');
   end
end

function test_esm_default_import_uses_full_source_bounds(testCase)
   % Real imports without explicit bounds must stage the complete source record.
   src = icemodel.verification.setup.fetchEsmSnowmip( ...
      strict=false, silent=true);
   has_source = esmSnowmipCacheComplete(src);
   testCase.assumeTrue(has_source, ...
      sprintf('ESM-SnowMIP source cache not present at %s', src));

   % Derive the expected CDP bounds directly from the unwindowed builders so a
   % regression back to the one-year smoke window cannot satisfy this check.
   [forcing_tt, ~] = ...
      icemodel.verification.setup.buildEsmSnowmipForcing( ...
      "cdp", source_dir=src);
   [obs_tt, ~] = ...
      icemodel.verification.setup.buildEsmSnowmipObservations( ...
      "cdp", source_dir=src);
   source_times = [forcing_tt.Time(:); obs_tt.Time(:)];
   source_times = source_times(~isnat(source_times));
   source_start = min(source_times);
   source_end = max(source_times);

   % Import the same site without startdate/enddate and compare the manifest
   % period with the independently derived source bounds.
   eval_root = fullfile(testCase.TestData.tmp, 'source-bounds', 'eval');
   input_root = fullfile(testCase.TestData.tmp, 'source-bounds', 'input');
   manifest = icemodel.verification.setup.importEsmSnowmip(src, ...
      case_ids="cdp", evaluation_data_root=string(eval_root), ...
      input_data_root=string(input_root), overwrite=true, ...
      overwrite_family=true, dt_out="");
   verifyEqual(testCase, string(manifest.cases.native_timestep), "hourly");
   verifyEqual(testCase, string(manifest.cases.period.start), ...
      string(icemodel.verification.setup.formatManifestTime(source_start)));
   verifyEqual(testCase, string(manifest.cases.period.end), ...
      string(icemodel.verification.setup.formatManifestTime(source_end)));

   % Keep this test meaningful by proving the full source bounds extend beyond
   % the short metadata-only preview used by dry_run.
   [smoke_start, smoke_end] = ...
      icemodel.verification.helpers.default_smoke_window("cdp");
   verifyLessThan(testCase, source_start, smoke_start);
   verifyGreaterThan(testCase, source_end, smoke_end);

   met_files = dir(fullfile(input_root, 'met', 'esm_snowmip', ...
      'met_cdp_esm_snowmip_*_1hr.mat'));
   verifyNotEmpty(testCase, met_files);
   saved = load(fullfile(met_files(1).folder, met_files(1).name), 'met');
   verifyEqual(testCase, seconds(median(diff(saved.met.Time))), 3600);
end

function test_esm_dry_run_does_not_write_staging_tree(testCase)
   % ESM-SnowMIP dry runs should validate selected cases without source reads
   % or staged-file writes.
   src = fullfile(testCase.TestData.tmp, 'missing-esm-cache');

   output_root = fullfile(testCase.TestData.tmp, 'esm-dry-run');
   manifest = icemodel.verification.setup.importEsmSnowmip(src, ...
      case_ids="cdp", output_root=string(output_root), dry_run=true);

   verifyEqual(testCase, string(manifest.dataset_family), "esm_snowmip");
   verifyEqual(testCase, string(manifest.cases.case_id), "cdp");
   verifyEqual(testCase, string(manifest.cases.native_timestep), "15m");
   [smoke_start, smoke_end] = ...
      icemodel.verification.helpers.default_smoke_window("cdp");
   verifyEqual(testCase, string(manifest.cases.period.start), ...
      string(icemodel.verification.setup.formatManifestTime(smoke_start)));
   verifyEqual(testCase, string(manifest.cases.period.end), ...
      string(icemodel.verification.setup.formatManifestTime(smoke_end)));
   verifyFalse(testCase, isfolder(src));
   verifyFalse(testCase, isfolder(fullfile(output_root, 'eval')));
   verifyFalse(testCase, isfolder(fullfile(output_root, 'input')));
end

function test_esm_subset_refresh_preserves_other_cases(testCase)
   % Incremental ESM-SnowMIP refreshes should preserve unrelated manifest cases.
   src = icemodel.verification.setup.fetchEsmSnowmip( ...
      strict=false, silent=true);
   has_source = esmSnowmipCacheComplete(src);
   testCase.assumeTrue(has_source, ...
      sprintf('ESM-SnowMIP source cache not present at %s', src));

   eval_root = fullfile(testCase.TestData.tmp, 'subset-eval');
   input_root = fullfile(testCase.TestData.tmp, 'subset-input');
   icemodel.verification.setup.importEsmSnowmip(src, ...
      case_ids=["cdp", "wfj"], evaluation_data_root=string(eval_root), ...
      input_data_root=string(input_root), overwrite=true, overwrite_family=true);
   manifest = icemodel.verification.setup.importEsmSnowmip(src, ...
      case_ids="cdp", evaluation_data_root=string(eval_root), ...
      input_data_root=string(input_root), overwrite=true);

   ids = sort(string({manifest.cases.case_id}));
   verifyEqual(testCase, ids(:), ["cdp"; "wfj"]);
end

function test_esm_skip_missing_reaches_case_staging(testCase)
   % skip_missing=true should skip missing site files after non-strict cache probe.
   src = fullfile(testCase.TestData.tmp, 'empty-esm-cache');
   eval_root = fullfile(testCase.TestData.tmp, 'empty-esm-eval');
   input_root = fullfile(testCase.TestData.tmp, 'empty-esm-input');
   mkdir(src);

   manifest = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.importEsmSnowmip(src, ...
      case_ids="cdp", evaluation_data_root=string(eval_root), ...
      input_data_root=string(input_root), skip_missing=true), ...
      'icemodel:verification:importEsmSnowmip:caseSkipped');

   verifyEmpty(testCase, manifest.cases);
   verifyEqual(testCase, string(manifest.skipped.site), "cdp");
   verifyFalse(testCase, isfolder(fullfile(eval_root, ...
      'esm_snowmip', 'cdp')));
end

function test_esm_explicit_eval_root_pairs_input_tree(testCase)
   % ESM met staging should use the input tree beside an explicit eval root.
   src = icemodel.verification.setup.fetchEsmSnowmip( ...
      strict=false, silent=true);
   has_source = esmSnowmipCacheComplete(src);
   testCase.assumeTrue(has_source, ...
      sprintf('ESM-SnowMIP source cache not present at %s', src));

   eval_root = fullfile(testCase.TestData.tmp, 'paired-esm', 'eval');
   icemodel.verification.setup.importEsmSnowmip(src, ...
      case_ids="cdp", evaluation_data_root=string(eval_root), ...
      overwrite=true, overwrite_family=true);

   met_files = dir(fullfile(testCase.TestData.tmp, 'paired-esm', ...
      'input', 'met', 'esm_snowmip', 'met_cdp_esm_snowmip_*.mat'));
   verifyNotEmpty(testCase, met_files);
end

function test_esm_output_root_stages_eval_and_input(testCase)
   % output_root should isolate both ESM eval and input artifacts.
   src = icemodel.verification.setup.fetchEsmSnowmip( ...
      strict=false, silent=true);
   has_source = esmSnowmipCacheComplete(src);
   testCase.assumeTrue(has_source, ...
      sprintf('ESM-SnowMIP source cache not present at %s', src));

   output_root = fullfile(testCase.TestData.tmp, 'esm-output-root');
   icemodel.verification.setup.importEsmSnowmip(src, ...
      case_ids="cdp", output_root=string(output_root), ...
      overwrite=true, overwrite_family=true);

   verifyTrue(testCase, isfile(fullfile(output_root, 'eval', ...
      'esm_snowmip', 'manifest.json')));
   met_files = dir(fullfile(output_root, 'input', 'met', 'esm_snowmip', ...
      'met_cdp_esm_snowmip_*.mat'));
   verifyNotEmpty(testCase, met_files);

   c = icemodel.verification.loadmanifest("cdp", ...
      evaluation_data_root=fullfile(output_root, 'eval'), ...
      dataset_family="esm_snowmip");
   opts = icemodel.test.helpers.setModelOptsForCase(c);
   verifyEqual(testCase, string(opts.pathinput), ...
      string(fullfile(output_root, 'input')));
   verifyEqual(testCase, string(opts.forcings), "esm_snowmip");
   verifyTrue(testCase, all(contains(string(opts.metfname), ...
      string(fullfile(output_root, 'input')))));
end

function test_esm_default_missing_cache_errors(testCase)
   % Default imports should fail at cache validation instead of silently rewriting
   % requested cases into skipped manifest entries.
   src = fullfile(testCase.TestData.tmp, 'strict-empty-esm-cache');
   mkdir(src);

   verifyError(testCase, ...
      @() runSilently(@() icemodel.verification.setup.importEsmSnowmip( ...
      src, case_ids="cdp")), ...
      'icemodel:verification:fetchEsmSnowmip:missingSources');
end

function test_esm_blank_source_uses_default_cache(testCase)
   % Blank source_dir should use the same default cache validated by fetch.
   src = icemodel.verification.setup.fetchEsmSnowmip( ...
      strict=false, silent=true);
   has_source = esmSnowmipCacheComplete(src);
   testCase.assumeTrue(has_source, ...
      sprintf('ESM-SnowMIP source cache not present at %s', src));

   eval_root = fullfile(testCase.TestData.tmp, 'blank-esm-eval');
   input_root = fullfile(testCase.TestData.tmp, 'blank-esm-input');
   manifest = icemodel.verification.setup.importEsmSnowmip("", ...
      case_ids="cdp", evaluation_data_root=string(eval_root), ...
      input_data_root=string(input_root), overwrite=true, ...
      overwrite_family=true);

   verifyEqual(testCase, string(manifest.cases.case_id), "cdp");

   % Direct builders apply the same blank/default rule before file discovery.
   [forcing_default, forcing_default_meta] = ...
      icemodel.verification.setup.buildEsmSnowmipForcing("cdp");
   [forcing_blank, forcing_blank_meta] = ...
      icemodel.verification.setup.buildEsmSnowmipForcing( ...
      "cdp", source_dir="");
   verifyEqual(testCase, forcing_blank.Time, forcing_default.Time);
   verifyEqual(testCase, string(forcing_blank_meta.met_file), ...
      string(forcing_default_meta.met_file));
   verifyEqual(testCase, string(forcing_blank_meta.obs_file), ...
      string(forcing_default_meta.obs_file));

   [obs_default, obs_default_meta] = ...
      icemodel.verification.setup.buildEsmSnowmipObservations("cdp");
   [obs_blank, obs_blank_meta] = ...
      icemodel.verification.setup.buildEsmSnowmipObservations( ...
      "cdp", source_dir="");
   verifyEqual(testCase, obs_blank.Time, obs_default.Time);
   verifyEqual(testCase, string(obs_blank_meta.obs_file), ...
      string(obs_default_meta.obs_file));
end

% ----- helpers ---------------------------------------------------------

function [Time, tsl] = writeEsmSoilTemperatureFixture(pathname)
   %WRITEESMSOILTEMPERATUREFIXTURE Write a depth-by-time ESM observation file.
   Time = datetime(2001, 1, 1, 0:4, 0, 0, 'TimeZone', 'UTC')';
   sdepth = single([0.1; 0.5; 1.0]);
   tsl = [101 102 103 104 105; ...
          201 202 203 204 205; ...
          301 302 303 304 305];

   % Encode dimensions in the same order as the bundled ESM-SnowMIP files.
   nccreate(pathname, 'time', 'Dimensions', {'time', numel(Time)});
   ncwrite(pathname, 'time', hours(Time - Time(1)));
   ncwriteatt(pathname, 'time', 'units', 'hours since 2001-01-01 00:00:00.0');
   nccreate(pathname, 'sdepth', 'Dimensions', {'sdepth', numel(sdepth)}, ...
      'Datatype', 'single');
   ncwrite(pathname, 'sdepth', sdepth);
   nccreate(pathname, 'tsl', ...
      'Dimensions', {'sdepth', numel(sdepth), 'time', numel(Time)});
   ncwrite(pathname, 'tsl', tsl);
end

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

function runSilently(fcn)
   %RUNSILENTLY Capture expected command-window text while preserving outputs.
   captured_fcn = fcn; %#ok<NASGU>
   evalc('captured_fcn();');
end

function bytes = fileBytes(filename)
   %FILEBYTES Read a staged binary artifact for no-churn assertions.
   fid = fopen(filename, 'r');
   cleanup = onCleanup(@() fclose(fid));
   bytes = fread(fid, Inf, '*uint8');
end
