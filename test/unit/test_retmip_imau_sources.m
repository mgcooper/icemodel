function tests = test_retmip_imau_sources
   %TEST_RETMIP_IMAU_SOURCES Verify RetMIP/IMAU metadata and fetch contracts.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install the canonical verification path and allocate an isolated cache.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.cache = tempname;
end

function teardown(testCase)
   % Remove temporary source-cache artifacts.
   if isfolder(testCase.TestData.cache)
      rmdir(testCase.TestData.cache, 's')
   end
   clear testCase.TestData.cleanup
end

function test_fetch_retmip_reports_doi_urls(testCase)
   % fetchRetmip must expose the GEUS Dataverse DOI URLs without requiring the
   % large local products to exist.
   [returned, status] = icemodel.verification.setup.fetchRetmip( ...
      cache_dir=testCase.TestData.cache, strict=false, silent=true);

   expected = string(testCase.TestData.cache);
   testCase.verifyEqual(string(returned), expected);
   returned = string({status.product});
   expected = ["forcing", "outputs", "scripts"];
   testCase.verifyTrue(all(ismember(expected, returned)));
   testCase.verifyTrue(any(contains([status.landing_url], ...
      "10.22008/FK2/GZ3CSN")));
   testCase.verifyTrue(any(contains([status.landing_url], ...
      "10.22008/FK2/CVPUJL")));
   testCase.verifyFalse(all([status.present]));
end

function test_fetch_imau_reports_pangaea_urls(testCase)
   % fetchImau must keep the hourly case inventory and daily QA product distinct.
   [~, returned] = icemodel.verification.setup.fetchImau( ...
      cache_dir=testCase.TestData.cache, strict=false, silent=true);

   expected = ["hourly", "daily"];
   testCase.verifyEqual(sort(string({returned.product})), sort(expected));
   testCase.verifyTrue(any(contains([returned.landing_url], ...
      "10.1594/PANGAEA.971647")));
   testCase.verifyTrue(any(contains([returned.landing_url], ...
      "10.1594/PANGAEA.970127")));
end

function test_fetch_promice_reports_required_local_files(testCase)
   % fetchPromice must report the local L3 station product and AWS metadata
   % requirements without attempting a download.
   [returned, status] = icemodel.verification.setup.fetchPromice( ...
      cache_dir=testCase.TestData.cache, strict=false, silent=true);

   expected = string(testCase.TestData.cache);
   testCase.verifyEqual(string(returned), expected);
   testCase.verifyEqual(string({status.product}), ["metadata", "hour"]);
   testCase.verifyFalse(all([status.present]));
   metadata = status(string({status.product}) == "metadata");
   testCase.verifyTrue(ismember("AWS_sites_metadata.csv", ...
      metadata.missing_files));
end

function test_fetch_promice_station_filter_passes_full_cache(testCase)
   % A cache with the requested station NetCDF and required metadata satisfies
   % strict mode; station matching is case/underscore tolerant like the reader.
   cache = makePromiceCache(testCase.TestData.cache, "KAN_M");

   [~, status] = icemodel.verification.setup.fetchPromice( ...
      cache_dir=cache, stations="kanm", strict=true, silent=true);

   testCase.verifyTrue(all([status.present]));
   hour = status(string({status.product}) == "hour");
   testCase.verifyEmpty(hour.missing_stations);
end

function test_import_promice_accepts_positional_source_dir(testCase)
   % importPromiceSites(source_dir, sites=...) must parse the source positionally
   % while retaining skip-missing behavior for invalid tiny test NetCDFs.
   cache = makePromiceCache(testCase.TestData.cache, "KAN_M");
   eval_root = fullfile(testCase.TestData.cache, 'eval');
   input_root = fullfile(testCase.TestData.cache, 'input');
   mkdir(eval_root);
   mkdir(input_root);

   returned = icemodel.verification.setup.importPromiceSites(cache, ...
      sites="KAN_M", models="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root, icemodel_config_casename="", ...
      overwrite_family=true);

   testCase.verifyEmpty(returned.cases);
   testCase.verifyEqual(string(returned.skipped.site), "KAN_M");
end

function test_import_promice_keeps_promice_dir_alias(testCase)
   % The existing promice_dir name-value path remains an alias for source_dir.
   cache = makePromiceCache(testCase.TestData.cache, "KAN_M");
   eval_root = fullfile(testCase.TestData.cache, 'eval');
   input_root = fullfile(testCase.TestData.cache, 'input');
   mkdir(eval_root);
   mkdir(input_root);

   returned = icemodel.verification.setup.importPromiceSites( ...
      sites="KAN_M", promice_dir=cache, models="promice", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      icemodel_config_casename="", overwrite_family=true);

   testCase.verifyEmpty(returned.cases);
   testCase.verifyEqual(string(returned.skipped.site), "KAN_M");
end

function test_fetch_retmip_full_cache_passes(testCase)
   % A cache with one file per requested RetMIP product satisfies strict mode.
   for product = ["forcing", "outputs", "scripts"]
      folder = fullfile(testCase.TestData.cache, product);
      mkdir(folder);
      touch(fullfile(folder, product + ".txt"));
   end

   [~, returned] = icemodel.verification.setup.fetchRetmip( ...
      cache_dir=testCase.TestData.cache, strict=true, silent=true);

   testCase.verifyTrue(all([returned.present]));
end

function test_retmip_case_metadata_uses_protocol_table_and_aliases(testCase)
   % RetMIP metadata should use protocol Table 1 coordinates, windows, and ids.
   returned = icemodel.verification.setup.retmipCaseMetadata( ...
      ["Dye-2_16", "sum"]);

   ids = string({returned.case_id});
   testCase.verifyEqual(ids, ["dye2_2016", "summit"]);
   dye = returned(ids == "dye2_2016");
   testCase.verifyEqual(string(dye.retmip_station_id), "dye2_16");
   testCase.verifyEqual(string(dye.protocol_id), "Dye-2_16");
   testCase.verifyEqual(dye.site_location.lat_wgs84, 66.48001);
   testCase.verifyEqual(dye.site_location.lon_wgs84, -46.27889);
   testCase.verifyEqual(dye.site_location.elev_m, 2165);
   testCase.verifyEqual(string(dye.period.start), "2016-05-02 00:00:00");
   testCase.verifyEqual(string(dye.period.end), "2016-10-28 09:00:00");
   summit = returned(ids == "summit");
   testCase.verifyEqual(string(summit.protocol_id), "Summit");
   testCase.verifyEqual(summit.site_location.elev_m, 3254);
end

function test_retmip_import_dry_run_manifest(testCase)
   % importRetmip dry-run should build the five protocol cases and source links
   % without writing empty observations.
   returned = icemodel.verification.setup.importRetmip( ...
      testCase.TestData.cache, dry_run=true);

   ids = string({returned.cases.case_id});
   expected = ["kanu", "dye2_long", "dye2_2016", "summit", "fa"];
   testCase.verifyEqual(ids, expected);
   fa = returned.cases(ids == "fa");
   expected = "imau";
   testCase.verifyEqual(string(fa.colocation.source_association.family), ...
      expected);
   expected = "S21";
   testCase.verifyEqual(string(fa.colocation.source_association.source_id), ...
      expected);
   expected = "3hr";
   testCase.verifyEqual(string(fa.native_timestep), expected);
   kanu = returned.cases(ids == "kanu");
   testCase.verifyEqual(kanu.site_location.lat_wgs84, 67.0003);
   testCase.verifyEqual(string(kanu.period.end), "2016-12-31 06:00:00");
   testCase.verifyTrue(all(ismember(["tsfc", "melt", "snowf_subl", ...
      "density", "subsurface_temperature", "lwc"], ...
      string(kanu.comparison_variables))));
   summit = returned.cases(ids == "summit");
   expected = "pending_gcnet_import";
   testCase.verifyEqual(string(summit.colocation.retmip.native_met_status), ...
      expected);
   testCase.verifyEqual(string(summit.colocation.native_met.status), ...
      expected);
end

function test_retmip_import_stages_protocol_bundle(testCase)
   % A non-dry RetMIP import stages a real protocol observations.mat bundle
   % from local surface/profile files without creating a normal met timetable.
   cache = fullfile(testCase.TestData.cache, 'retmip-cache');
   forcing = fullfile(cache, 'forcing');
   outputs = fullfile(cache, 'outputs');
   scripts = fullfile(cache, 'scripts');
   mkdir(forcing);
   mkdir(outputs);
   mkdir(scripts);
   writeText(fullfile(forcing, 'RetMIP_surface_forcing_KANU.tab'), ...
      sprintf(['time;melt_mmweq;acc_subl_mmweq;Tsurf_K\n' ...
      '"01-May-2012 00:00:00;0;0;260"\n' ...
      '"01-May-2012 03:00:00;0;0;261"\n']));
   writeText(fullfile(forcing, 'RetMIP_initial_firn_density_KAN-U.tab'), ...
      sprintf(['depth_m;density_kgm3\n' ...
      '"0.1;350"\n']));
   nccreate(fullfile(outputs, 'RetMIP_GEUS_KANU_3hourly_values.nc'), ...
      'temp', 'Dimensions', {'time', 2});
   nccreate(fullfile(outputs, 'RetMIP_GEUS_KANU_3hourly_values.nc'), ...
      'rho', 'Dimensions', {'time', 2});
   touch(fullfile(scripts, 'README.md'));

   eval_root = fullfile(testCase.TestData.cache, 'eval');
   mkdir(eval_root);
   returned = icemodel.verification.setup.importRetmip(cache, ...
      case_ids="kanu", dry_run=false, evaluation_data_root=eval_root, ...
      overwrite_family=true);

   expected = "kanu/observations.mat";
   testCase.verifyEqual(string(returned.cases.evaluation_file), expected);
   testCase.verifyTrue(isfile(fullfile(eval_root, 'retmip', expected)));
   testCase.verifyTrue(logical(returned.cases.colocation.retmip.staged));
   testCase.verifyTrue(ismember("temp", ...
      string(returned.cases.comparison_variables)));
   testCase.verifyTrue(ismember("rho", ...
      string(returned.cases.comparison_variables)));
end

function test_imau_import_dry_run_manifest(testCase)
   % importImau dry-run should expose S21/S22/S23 and the S21->RetMIP FA link.
   returned = icemodel.verification.setup.importImau( ...
      testCase.TestData.cache, dry_run=true);

   ids = string({returned.cases.site_id});
   expected = ["S21", "S22", "S23"];
   testCase.verifyEqual(ids, expected);
   s21 = returned.cases(ids == "S21");
   expected = "retmip";
   testCase.verifyEqual(string(s21.colocation.source_association.family), ...
      expected);
   expected = "fa";
   testCase.verifyEqual(string(s21.colocation.source_association.source_id), ...
      expected);
   expected = "1hr";
   testCase.verifyEqual(string(s21.native_timestep), expected);
end

function test_sumup_colocation_records_mixed_anchor_metadata(testCase)
   % sumupColocation should return nearest family/source metadata from explicit
   % mixed anchors instead of assuming every anchor is PROMICE.
   anchors = [ ...
      struct('site', "FA", 'family', "retmip", 'source_id', "fa", ...
         'x_epsg3413', 0, 'y_epsg3413', 0)
      struct('site', "S21", 'family', "imau", 'source_id', "S21", ...
         'x_epsg3413', 50000, 'y_epsg3413', 0)];

   [tf, returned, distance_km] = ...
      icemodel.verification.helpers.sumupColocation(1000, 0, ...
      anchors=anchors, threshold_km=7.5);

   testCase.verifyTrue(tf);
   expected = 1;
   testCase.verifyEqual(distance_km, expected);
   expected = "retmip";
   testCase.verifyEqual(string(returned.family), expected);
   expected = "fa";
   testCase.verifyEqual(string(returned.source_id), expected);

   [tf2, returned2, distance2_km] = ...
      icemodel.verification.setup.anchorColocation(1000, 0, ...
      anchors=anchors, threshold_km=7.5);

   % The SUMup compatibility helper and the setup-level generic helper should
   % report the same nearest mixed anchor.
   testCase.verifyEqual(tf2, tf);
   testCase.verifyEqual(distance2_km, distance_km);
   testCase.verifyEqual(string(returned2.family), "retmip");
end

function test_sumup_colocation_keeps_nearest_when_outside_threshold(testCase)
   % The nearest anchor is still useful metadata even when the point is outside
   % the colocation threshold.
   anchors = struct('site', "HUMPHREY", 'family', "research_site", ...
      'source_id', "humphrey", 'x_epsg3413', 0, 'y_epsg3413', 0);

   [tf, returned, distance_km] = ...
      icemodel.verification.helpers.sumupColocation(20000, 0, ...
      anchors=anchors, threshold_km=7.5);

   testCase.verifyFalse(tf);
   expected = 20;
   testCase.verifyEqual(distance_km, expected);
   expected = "HUMPHREY";
   testCase.verifyEqual(string(returned.site), expected);
   expected = "research_site";
   testCase.verifyEqual(string(returned.family), expected);
end

function test_research_site_metadata_defines_humphrey_catchall(testCase)
   % Humphrey belongs to the generic research_site family, not a bespoke anchor
   % category or a SUMup subfamily.
   returned = icemodel.verification.setup.researchSiteMetadata("humphrey");

   expected = "research_site";
   testCase.verifyEqual(string(returned.family), expected);
   expected = "humphrey";
   testCase.verifyEqual(string(returned.source_id), expected);
   testCase.verifyGreaterThan(returned.lat_wgs84, 69);
   testCase.verifyLessThan(returned.lon_wgs84, -48);
end

function test_mixed_anchor_catalog_reads_staged_manifests(testCase)
   % mixedAnchorCatalog should tolerate partial eval trees and expose source
   % availability from the staged manifest rows it does find.
   eval_root = fullfile(testCase.TestData.cache, 'catalog-eval');
   writeTinyManifest(eval_root, "promice", "cp1", "CP1", 69.0, -49.0, ...
      1, 2, {'promice'}, {'promice'}, {'promice_obs'});
   writeTinyManifest(eval_root, "research_site", "humphrey", "humphrey", ...
      69.725714, -48.190512, 3, 4, {}, {}, {'sumup_obs'});

   anchors = icemodel.verification.setup.mixedAnchorCatalog( ...
      evaluation_data_root=eval_root);

   testCase.verifyEqual(numel(anchors), 2);
   families = sort(string({anchors.family}));
   testCase.verifyEqual(families, ["promice", "research_site"]);
   cp1 = anchors(string({anchors.case_id}) == "cp1");
   testCase.verifyEqual(string(cp1.site), "CP1");
   testCase.verifyEqual(string(cp1.met_sources), "promice");
   testCase.verifyEqual(string(cp1.userdata_sources), "promice");
end

function test_import_research_sites_dry_run_humphrey_manifest(testCase)
   % Dry-run research-site import should build the Humphrey manifest contract
   % without requiring local SUMup NetCDF products.
   eval_root = fullfile(testCase.TestData.cache, 'research-eval');
   writeTinyManifest(eval_root, "promice", "jar", "JAR", 69.7, -48.2, ...
      2000, 2000, {'promice'}, {'promice'}, {'promice_obs'});

   returned = icemodel.verification.setup.importResearchSites("", ...
      site_ids="humphrey", dry_run=true, evaluation_data_root=eval_root, ...
      input_data_root=fullfile(testCase.TestData.cache, 'input'));

   testCase.verifyEqual(string(returned.dataset_family), "research_site");
   testCase.verifyEqual(string(returned.cases.case_id), "humphrey");
   testCase.verifyEqual(string(returned.cases.forcing_sources), strings(0, 0));
   testCase.verifyEqual(string(returned.cases.eval_sources), "sumup_obs");
   testCase.verifyTrue(isfield(returned.cases.colocation, ...
      'nearest_noncolocated_promice'));
   nearest = returned.cases.colocation.nearest_noncolocated_promice;
   testCase.verifyTrue(any(string(nearest.nearest_anchor) == ["CP1", "JAR"]));
   testCase.verifyFalse(logical(returned.cases.colocation.research_site.staged));
end

function touch(filename)
   %TOUCH Create a tiny file for source-cache presence tests.
   fid = fopen(filename, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, 'placeholder\n');
   clear cleaner
end

function writeText(filename, text)
   %WRITETEXT Write one tiny source-cache fixture.
   fid = fopen(filename, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', text);
   clear cleaner
end

function writeTinyManifest(eval_root, family, case_id, site_id, lat, lon, x, y, ...
      forcing_sources, data_sources, eval_sources)
   %WRITETINYMANIFEST Write a minimal firn family manifest for catalog tests.
   family_root = fullfile(eval_root, family);
   mkdir(family_root);
   colocation = struct();
   for k = 1:numel(data_sources)
      src = char(data_sources{k});
      colocation.(src) = struct('kind', src, 'staged', true, ...
         'met_files', {{sprintf('%s-met.mat', src)}}, ...
         'data_files', {{sprintf('%s-data.mat', src)}});
   end
   values = { ...
      char(case_id)
      'firn_observational'
      char(site_id)
      char(site_id)
      'unknown'
      {'firn'}
      'unknown'
      struct('lat_wgs84', lat, 'lon_wgs84', lon, ...
         'x_epsg3413', x, 'y_epsg3413', y, 'elev_m', 1000)
      struct('start', '', 'end', '')
      ''
      forcing_sources
      eval_sources
      {}
      struct()
      colocation
      'irregular'
      'test manifest'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      family, "", "", "test", string(datetime('today')), entry);
   icemodel.verification.setup.writeManifest( ...
      fullfile(family_root, "manifest.json"), manifest);
end

function cache = makePromiceCache(root, station)
   %MAKEPROMICECACHE Create the smallest PROMICE-shaped cache for unit tests.
   cache = fullfile(root, 'promice-cache');
   hour = fullfile(cache, 'hour');
   mkdir(hour);
   writeText(fullfile(cache, "AWS_sites_metadata.csv"), ...
      sprintf('site_id,location_type\n%s,ice sheet\n', station));
   touch(fullfile(cache, "AWS_stations_metadata.csv"));
   touch(fullfile(cache, "AWS_variables.csv"));
   touch(fullfile(hour, station + "_hour.nc"));
end
