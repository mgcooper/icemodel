function tests = test_ktransect_sources
   %TEST_KTRANSECT_SOURCES Verify K-transect parser, builder, and staging contracts.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install the verification path and allocate a parser fixture folder.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.tmp = tempname;
   mkdir(fullfile(testCase.TestData.tmp, 'datasets'));
end

function teardown(testCase)
   % Remove temporary parser fixtures.
   if isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
   clear testCase.TestData.cleanup
end

%% Parser contracts

function test_parser_maps_native_channels(testCase)
   % Source channels map to icemodel-native names and units (K, Pa) with the
   % child DOI, series DOI, and event coordinates preserved.
   filename = fixtureFile(testCase, "AWS9", 2010);
   writeKtransectFixture(filename, "AWS9", 2010, "type0");

   [returned, meta] = icemodel.forcing.helpers.readKtransectTable(filename);

   expected = ["tair", "rh", "wspd", "wdir", "psfc", "swd", "swu", ...
      "lwd", "lwu", "height_rel", "aws_type"];
   testCase.verifyEqual(string(returned.Properties.VariableNames), expected);
   testCase.verifyEqual(string(returned.Properties.DimensionNames(1)), "time");
   testCase.verifyEqual(returned.tair(1), 268.15);
   testCase.verifyEqual(returned.psfc(1), 90000);
   testCase.verifyEqual(returned.rh(1), 80);
   testCase.verifyEqual(returned.wspd(1), 5);
   testCase.verifyEqual(returned.wdir(1), 180);
   testCase.verifyEqual(returned.aws_type(1), 0);
   testCase.verifyEqual(returned.time(1), ...
      datetime(2010, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'));
   testCase.verifyEqual(string(meta.site_id), "AWS9");
   testCase.verifyEqual(string(meta.doi), "10.1594/PANGAEA.950093");
   testCase.verifyEqual(string(meta.bundle_doi), "10.1594/PANGAEA.947483");
   testCase.verifyEqual(meta.event.lat_wgs84, 67.052460);
   testCase.verifyEqual(meta.event.elev_m, 1500.0);
   % The K-transect event line posts one bare DATE/TIME stamp; it maps to the
   % event start through the shared PANGAEA parser.
   testCase.verifyEqual(meta.event.start, ...
      datetime(2003, 8, 27, 20, 30, 0, 'TimeZone', 'UTC'));
end

function test_parser_handles_all_column_variants(testCase)
   % Type-0 AWS10 files omit "T tech", type-0 AWS6 files include it, and
   % type-1 files add "Ice melt"; the name-based mapping accepts all three.
   short_file = fixtureFile(testCase, "AWS10", 2012);
   writeKtransectFixture(short_file, "AWS10", 2012, "type0_short");
   full_file = fixtureFile(testCase, "AWS6", 2013);
   writeKtransectFixture(full_file, "AWS6", 2013, "type0");
   melt_file = fixtureFile(testCase, "AWS5", 2018);
   writeKtransectFixture(melt_file, "AWS5", 2018, "type1");

   short_data = icemodel.forcing.helpers.readKtransectTable(short_file);
   full_data = icemodel.forcing.helpers.readKtransectTable(full_file);
   melt_data = icemodel.forcing.helpers.readKtransectTable(melt_file);

   testCase.verifyFalse(ismember("ice_melt", ...
      string(short_data.Properties.VariableNames)));
   testCase.verifyFalse(ismember("ice_melt", ...
      string(full_data.Properties.VariableNames)));
   testCase.verifyTrue(ismember("ice_melt", ...
      string(melt_data.Properties.VariableNames)));
   testCase.verifyEqual(melt_data.ice_melt(1), 2.5);
   testCase.verifyEqual(melt_data.aws_type(1), 1);
end

function test_parser_extracts_maintenance_visits(testCase)
   % The once-yearly battery-voltage value of 100 marks the maintenance visit
   % and must surface as metadata rather than ride along as a channel.
   filename = fixtureFile(testCase, "AWS9", 2010);
   writeKtransectFixture(filename, "AWS9", 2010, "type0", visit_row=2);

   [~, meta] = icemodel.forcing.helpers.readKtransectTable(filename);

   expected = datetime(2010, 1, 1, 0, 30, 0, 'TimeZone', 'UTC');
   testCase.verifyEqual(meta.visits, expected);
end

function test_parser_accepts_seconds_bearing_timestamps(testCase)
   % Valid annual postings with explicit seconds follow the documented variant.
   filename = fixtureFile(testCase, "AWS9", 2010);
   writeKtransectFixture(filename, "AWS9", 2010, "type0");
   text = string(fileread(filename));
   text = regexprep(text, ...
      '(\d{4}-\d{2}-\d{2}T\d{2}:\d{2})(?=\t)', '$1:00');
   writeText(filename, text);

   returned = icemodel.forcing.helpers.readKtransectTable(filename);

   expected = datetime(2010, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + minutes(30) * (0:3)';
   testCase.verifyEqual(returned.Properties.RowTimes, expected);
end

function test_parser_preserves_blank_cells_as_nan(testCase)
   % Blank numeric PANGAEA cells are missing values, not zeros.
   filename = fixtureFile(testCase, "AWS9", 2010);
   writeKtransectFixture(filename, "AWS9", 2010, "type0", blank_met=true);

   returned = icemodel.forcing.helpers.readKtransectTable(filename);

   testCase.verifyTrue(isnan(returned.tair(1)));
   testCase.verifyTrue(isnan(returned.rh(1)));
   testCase.verifyTrue(isnan(returned.swd(1)));
end

function test_parser_masks_direction_below_speed_threshold(testCase)
   % Wind direction is undefined below 0.1 m/s even when the source cell
   % carries a finite compass value.
   filename = fixtureFile(testCase, "AWS9", 2010);
   writeKtransectFixture(filename, "AWS9", 2010, "type0", calm_row=1);

   returned = icemodel.forcing.helpers.readKtransectTable(filename);

   testCase.verifyEqual(returned.wspd(1), 0.05);
   testCase.verifyTrue(isnan(returned.wdir(1)));
   testCase.verifyEqual(returned.wdir(2), 180);
end

function test_parser_rejects_missing_child_doi(testCase)
   % Every annual child must carry its own pinned PANGAEA DOI.
   filename = fixtureFile(testCase, "AWS9", 2010);
   writeKtransectFixture(filename, "AWS9", 2010, "type0", ...
      child_doi="not-a-doi");

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.readKtransectTable(filename), ...
      'icemodel:forcing:readKtransectTable:missingChildDoi');
end

function test_parser_rejects_missing_tabular_section(testCase)
   % A cache artifact without the Date/Time table should fail clearly.
   filename = fullfile(testCase.TestData.tmp, 'datasets', 'bad.tab');
   writeText(filename, sprintf(['/* DATA DESCRIPTION:\n' ...
      'Citation:\tmissing table\n*/\nnot a data table\n']));

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.readKtransectTable(filename), ...
      'icemodel:forcing:readKtransectTable:missingHeader');
end

function test_parser_rejects_missing_required_columns(testCase)
   % Files without a required channel column must fail by name, not position.
   filename = fixtureFile(testCase, "AWS9", 2010);
   writeKtransectFixture(filename, "AWS9", 2010, "type0", drop_rh=true);

   testCase.verifyError(@() ...
      icemodel.forcing.helpers.readKtransectTable(filename), ...
      'icemodel:forcing:readKtransectTable:missingColumns');
end

%% Builder contracts

function test_builder_merges_annual_children(testCase)
   % Annual children merge onto one regular 30-minute axis with per-child DOI
   % provenance and source-faithful placeholder precipitation.
   writeKtransectFixture(fixtureFile(testCase, "AWS9", 2010), ...
      "AWS9", 2010, "type0");
   writeKtransectFixture(fixtureFile(testCase, "AWS9", 2011), ...
      "AWS9", 2011, "type0", start_time=datetime(2011, 1, 1, 0, 0, 0, ...
      'TimeZone', 'UTC'));

   [Data, metadata] = icemodel.forcing.buildKtransectData("AWS9", ...
      source_dir=testCase.TestData.tmp);

   % One regular half-hourly axis spanning both children.
   testCase.verifyEqual(unique(diff(Data.Time)), minutes(30));
   testCase.verifyEqual(Data.Time(1), ...
      datetime(2010, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'));
   testCase.verifyEqual(Data.Time(end), ...
      datetime(2011, 1, 1, 1, 30, 0, 'TimeZone', 'UTC'));
   % Placeholder precipitation channels stay all-missing by policy.
   testCase.verifyTrue(all(isnan(Data.rainf)))
   testCase.verifyTrue(all(isnan(Data.snowf)))
   % Child DOI pinning and family provenance.
   testCase.verifyEqual([metadata.children.year], [2010 2011]);
   testCase.verifyEqual(string(metadata.bundle_doi), ...
      "10.1594/PANGAEA.947483");
   testCase.verifyEqual(string(metadata.license), "CC-BY-4.0");
   testCase.verifyEqual(string(metadata.source_family), "ktransect");
   testCase.verifyEqual(string(metadata.station), "AWS9");
   testCase.verifyTrue(contains(metadata.doi, "10.1594/PANGAEA.950093"));
   testCase.verifyFalse(contains(string(jsonencode(metadata)), ...
      string(testCase.TestData.tmp)));
end

function test_builder_regularizes_interior_gaps_as_missing(testCase)
   % A missing half-hour posting becomes an explicit missing row, never an
   % interpolated value.
   writeKtransectFixture(fixtureFile(testCase, "AWS9", 2010), ...
      "AWS9", 2010, "type0", skip_third_row=true);

   Data = icemodel.forcing.buildKtransectData("AWS9", ...
      source_dir=testCase.TestData.tmp);

   gap_row = Data.Time == datetime(2010, 1, 1, 1, 0, 0, 'TimeZone', 'UTC');
   testCase.verifyTrue(any(gap_row));
   testCase.verifyTrue(isnan(Data.tair(gap_row)));
end

function test_builder_rejects_overlapping_annual_children(testCase)
   % Overlapping annual children would double-count samples; refuse loudly.
   % The second file deliberately repeats the first file's timestamps.
   writeKtransectFixture(fixtureFile(testCase, "AWS9", 2010), ...
      "AWS9", 2010, "type0");
   writeKtransectFixture(fixtureFile(testCase, "AWS9", 2011), ...
      "AWS9", 2011, "type0", start_time=datetime(2010, 1, 1, 0, 0, 0, ...
      'TimeZone', 'UTC'));

   testCase.verifyError(@() icemodel.forcing.buildKtransectData("AWS9", ...
      source_dir=testCase.TestData.tmp), ...
      'icemodel:forcing:buildKtransectData:overlappingAnnualFiles');
end

function test_builder_rejects_child_station_mismatch(testCase)
   % A filename match cannot override the station identity in the child
   % event metadata.
   filename = fixtureFile(testCase, "AWS9", 2010);
   writeKtransectFixture(filename, "AWS6", 2010, "type0");

   testCase.verifyError(@() ...
      icemodel.forcing.buildKtransectData("AWS9", ...
      source_dir=testCase.TestData.tmp), ...
      'icemodel:forcing:buildKtransectData:stationMismatch');
end

function test_builder_rejects_unknown_station(testCase)
   % Only the four published stations are valid.
   testCase.verifyError(@() icemodel.forcing.buildKtransectData("AWS7", ...
      source_dir=testCase.TestData.tmp), ...
      'icemodel:forcing:buildKtransectData:unknownStation');
end

function test_heights_reader_parses_station_blocks(testCase)
   % The sensor-height workbook parses into per-visit records keyed by the
   % station block, the year row, and the arrival/depart label, with the
   % per-type explanation strings preserved as notes.
   filename = writeHeightsFixture(testCase);

   [returned, notes] = icemodel.forcing.helpers.readKtransectHeights( ...
      filename);

   testCase.verifyEqual(string({returned.station}), ...
      ["AWS9", "AWS9", "AWS9", "AWS10"]);
   testCase.verifyEqual([returned.year], [2003 2004 2004 2010]);
   testCase.verifyEqual(string({returned.event}), ...
      ["arrival/leave", "arrival", "depart", "arrival"]);
   testCase.verifyEqual([returned.height_u], [602 602 610 480]);
   testCase.verifyEqual([returned.height_T], [550 550 560 450]);
   testCase.verifyEqual([returned.aws_type], [0 0 0 0]);
   testCase.verifyEqual(numel(notes), 2);
   testCase.verifyTrue(all(startsWith(notes, "For AWS type")));

   % The station filter returns only that block's rows.
   filtered = icemodel.forcing.helpers.readKtransectHeights( ...
      filename, station="AWS10");
   testCase.verifyEqual(string({filtered.station}), "AWS10");
end

function test_builder_attaches_sensor_heights(testCase)
   % With the workbook cached beside the annual files, the built Data carries
   % the station's sensor-height records; without it, the recorded reason
   % says so instead of erroring.
   writeKtransectFixture(fixtureFile(testCase, "AWS9", 2010), ...
      "AWS9", 2010, "type0");

   [~, bare] = icemodel.forcing.buildKtransectData("AWS9", ...
      source_dir=testCase.TestData.tmp);
   testCase.verifyFalse(bare.sensor_heights.present);
   testCase.verifyTrue(contains(bare.sensor_heights.reason, ...
      "not present"));

   writeHeightsFixture(testCase);
   [~, metadata] = icemodel.forcing.buildKtransectData("AWS9", ...
      source_dir=testCase.TestData.tmp);
   testCase.verifyTrue(metadata.sensor_heights.present);
   testCase.verifyEqual( ...
      string({metadata.sensor_heights.records.station}), ...
      ["AWS9", "AWS9", "AWS9"]);
end

%% Catalog, crosswalk, fetch, and namelist contracts

function test_site_catalog_contract(testCase)
   % The catalog publishes the four stations in canonical order and rejects
   % unknown ids before staging.
   returned = icemodel.verification.setup.ktransectSiteCatalog();

   testCase.verifyEqual(string({returned.site_id}), ...
      ["AWS5", "AWS6", "AWS9", "AWS10"]);
   testCase.verifyEqual(string({returned.case_id}), ...
      ["aws5", "aws6", "aws9", "aws10"]);
   selected = icemodel.verification.setup.ktransectSiteCatalog("AWS10");
   testCase.verifyEqual(string(selected.surface_zone), "accumulation");
   testCase.verifyError(@() ...
      icemodel.verification.setup.ktransectSiteCatalog("AWS1"), ...
      'icemodel:verification:ktransectSiteCatalog:unknownSite');
end

function test_alias_crosswalk_is_verified_with_evidence(testCase)
   % Crosswalk rows are human-approved for donor evaluation; each carries
   % measured evidence, and any row reverting to "hypothesis" re-blocks
   % donor use of that pairing.
   returned = icemodel.verification.setup.ktransectAliasCrosswalk();

   testCase.verifyEqual(string({returned.pangaea_id}), ...
      ["AWS5", "AWS6", "AWS9", "AWS10", ...
      "AWS5", "AWS6", "AWS6", "AWS9", "AWS10"]);
   testCase.verifyEqual(string({returned.alias_id}), ...
      ["S5", "S6", "S9", "S10", ...
      "KAN_L", "KAN_L", "KAN_M", "KAN_M", "KAN_U"]);
   testCase.verifyTrue(all(string({returned.status}) == "verified"));
   testCase.verifyTrue(all(strlength(string({returned.source})) > 0));
   testCase.verifyTrue(all(strlength(string({returned.evidence})) > 0));
end

function test_fetch_status_reports_cached_products(testCase)
   % The cache validator reports per-product presence without downloading.
   writeKtransectFixture(fixtureFile(testCase, "AWS9", 2010), ...
      "AWS9", 2010, "type0");

   [source_dir, status] = icemodel.verification.setup.fetchKtransect( ...
      cache_dir=testCase.TestData.tmp, strict=false, silent=true);

   testCase.verifyEqual(string(source_dir), string(testCase.TestData.tmp));
   names = string({status.product});
   testCase.verifyEqual(sort(names), sort(["datasets", "heights"]));
   datasets_row = status(names == "datasets");
   testCase.verifyTrue(datasets_row.present);
   testCase.verifyEqual(string(datasets_row.doi), "10.1594/PANGAEA.947483");
   heights_row = status(names == "heights");
   testCase.verifyFalse(heights_row.present);
end

function test_fetch_rejects_unknown_product(testCase)
   % Unknown product selectors fail before touching the cache.
   testCase.verifyError(@() icemodel.verification.setup.fetchKtransect( ...
      cache_dir=testCase.TestData.tmp, products="hourly"), ...
      'icemodel:verification:fetchKtransect:unknownProduct');
end

function test_family_registered_in_namelist(testCase)
   % ktransect is a manifest-discovered family: present in datasetfamily, and
   % deliberately absent from the static caseid dispatcher and the accepted
   % firn-preview list.
   testCase.verifyTrue(ismember("ktransect", ...
      icemodel.verification.namelists.datasetfamily()));
   testCase.verifyFalse(ismember("ktransect", ...
      icemodel.verification.namelists.firndatasetfamily()));
   % caseid raises an identifier-less error for non-static families, so pin
   % the rejection through its message rather than an error id.
   testCase.verifyThat(@() ...
      icemodel.verification.namelists.caseid("ktransect"), ...
      matlab.unittest.constraints.Throws(?MException));
end

function test_import_dry_run_manifest_shape(testCase)
   % A dry run returns the family manifest shape without reading caches or
   % writing artifacts, including the required PANGAEA attribution fields.
   manifest = icemodel.verification.setup.importKtransect("", dry_run=true);

   testCase.verifyEqual(string(manifest.dataset_family), "ktransect");
   testCase.verifyEqual(string(manifest.source_doi), ...
      "10.1594/PANGAEA.947483");
   testCase.verifyEqual(numel(manifest.cases), 4);
   testCase.verifyTrue(contains(string(manifest.citation), ...
      "10.1594/PANGAEA.947483"));
   testCase.verifyEqual(string(manifest.license), "CC-BY-4.0");
end

function test_import_stages_case_with_provenance(testCase)
   % A real (non-dry-run) import writes the observation artifact and a
   % manifest case whose native leg pins the child DOIs, license, and the
   % explicit no-forcing readiness reason.
   writeKtransectFixture(fixtureFile(testCase, "AWS9", 2010), ...
      "AWS9", 2010, "type0");
   writeKtransectFixture(fixtureFile(testCase, "AWS9", 2011), ...
      "AWS9", 2011, "type0", start_time=datetime(2011, 1, 1, 0, 0, 0, ...
      'TimeZone', 'UTC'));
   writeHeightsFixture(testCase);
   output_root = fullfile(testCase.TestData.tmp, 'staged');

   manifest = icemodel.verification.setup.importKtransect( ...
      testCase.TestData.tmp, case_ids="AWS9", output_root=output_root);

   testCase.verifyTrue(isfile(fullfile(output_root, 'eval', 'ktransect', ...
      'aws9', 'observations.mat')));
   testCase.verifyTrue(contains(string(manifest.citation), ...
      "10.1594/PANGAEA.947483"));
   testCase.verifyEqual(string(manifest.license), "CC-BY-4.0");
   testCase.verifyEqual(numel(manifest.cases), 1);
   entry = manifest.cases(1);
   testCase.verifyEqual(string(entry.case_id), "aws9");
   leg = entry.colocation.ktransect;
   testCase.verifyEqual(string(leg.kind), "annual_aws_met_and_eval");
    testCase.verifyEqual(string(leg.license), "CC-BY-4.0");
    testCase.verifyEqual([leg.children.year].', [2010; 2011]);
    testCase.verifyFalse(isfield(leg, 'source_dir'));
    testCase.verifyFalse(isfield(leg, 'source_file'));
    evaluation_file = fullfile(output_root, 'eval', 'ktransect', ...
       leg.evaluation_file);
    evaluation_info = dir(evaluation_file);
    testCase.verifyEqual(double(leg.evaluation_size_bytes), ...
       evaluation_info.bytes);
    testCase.verifyEqual(string(leg.evaluation_sha256), ...
       icemodel.verification.setup.fileSha256(evaluation_file));
    testCase.verifyFalse(contains(string(jsonencode(manifest)), ...
       string(testCase.TestData.tmp)));
    staged = load(fullfile(output_root, 'eval', 'ktransect', ...
       'aws9', 'observations.mat'), 'targets');
    testCase.verifyFalse(contains(string(jsonencode( ...
       staged.targets.metadata)), string(testCase.TestData.tmp)));
    testCase.verifyFalse(logical(leg.forcing_ready));
   testCase.verifyTrue(contains(string(leg.forcing_ready_reason), ...
      "build_forcing=false"));
   testCase.verifyEqual(string(entry.evaluation_file), ...
      string(fullfile('aws9', 'observations.mat')));
end

function test_import_stages_missing_height_attachment_with_warning(testCase)
   % POLICY A3: height provenance never blocks staging. A tab-only cache
   % stages with the absence recorded and a warning raised, so donor
   % assembly can annotate and the admission gates judge quality.
   writeKtransectFixture(fixtureFile(testCase, "AWS9", 2010), ...
      "AWS9", 2010, "type0");

   output_root = fullfile(testCase.TestData.tmp, 'staged');
   testCase.verifyWarning(@() ...
      icemodel.verification.setup.importKtransect( ...
      testCase.TestData.tmp, case_ids="AWS9", ...
      output_root=output_root), ...
      'icemodel:verification:importKtransect:missingSensorHeights');
   testCase.verifyTrue(isfile(fullfile(output_root, 'eval', ...
      'ktransect', 'manifest.json')));
end

function test_import_flags_native_identity_conflict(testCase)
   % When a prior case's native product identity no longer matches the fresh
   % source, the refreshed leg must demand an explicit native rebuild instead
   % of silently reusing stale forcing references.
   writeKtransectFixture(fixtureFile(testCase, "AWS9", 2010), ...
      "AWS9", 2010, "type0");
   writeHeightsFixture(testCase);
   output_root = fullfile(testCase.TestData.tmp, 'staged');
   icemodel.verification.setup.importKtransect( ...
      testCase.TestData.tmp, case_ids="AWS9", output_root=output_root);

   % Corrupt the staged product identity the way a changed upstream would.
   manifest_file = fullfile(output_root, 'eval', 'ktransect', ...
      'manifest.json');
   staged = jsondecode(fileread(manifest_file));
   staged.cases(1).colocation.ktransect.doi = '10.1594/PANGAEA.000000';
   fid = fopen(manifest_file, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', jsonencode(staged));
   clear cleaner

   manifest = icemodel.verification.setup.importKtransect( ...
      testCase.TestData.tmp, case_ids="AWS9", output_root=output_root, ...
      overwrite=true);

   leg = manifest.cases(1).colocation.ktransect;
   testCase.verifyEqual(string(leg.kind), "annual_aws_eval");
   testCase.verifyFalse(logical(leg.forcing_ready));
   testCase.verifyTrue(contains(string(leg.forcing_ready_reason), ...
      "identity changed"));
end

%% Fixture writers

function filename = fixtureFile(testCase, site_id, year)
   %FIXTUREFILE Canonical fixture path for one annual child.
   filename = fullfile(testCase.TestData.tmp, 'datasets', ...
      sprintf('K-transect_%s_%d.tab', site_id, year));
end

function writeKtransectFixture(filename, site_id, year, variant, kwargs)
   %WRITEKTRANSECTFIXTURE Write a tiny PANGAEA-shaped K-transect annual table.
   arguments
      filename (1, 1) string
      site_id (1, 1) string
      year (1, 1) double
      variant (1, 1) string ...
         {mustBeMember(variant, ["type0_short", "type0", "type1"])}
      kwargs.start_time (1, 1) datetime = NaT
      kwargs.visit_row (1, 1) double = 0
      kwargs.blank_met (1, 1) logical = false
      kwargs.skip_third_row (1, 1) logical = false
      kwargs.drop_rh (1, 1) logical = false
      kwargs.calm_row (1, 1) double = 0
      kwargs.child_doi (1, 1) string = "auto"
   end

   % Default four half-hourly postings starting at midnight on Jan 1.
   start_time = kwargs.start_time;
   if ismissing(start_time)
      start_time = datetime(year, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
   end
   times = start_time + minutes(30) * (0:3);
   if kwargs.skip_third_row
      times(3) = [];
   end

   child_doi = kwargs.child_doi;
   if child_doi == "auto"
      child_doi = childDoi(site_id, year);
   end
   metadata = sprintf(['/* DATA DESCRIPTION:\n' ...
      'Citation:\tSmeets, Paul C J P; et al. (2022): Automatic weather ' ...
      'station data from %s collected during %d at the Greenland ice ' ...
      'sheet along the K-transect, West-Greenland [dataset]. PANGAEA, ' ...
      'https://doi.org/%s,\n' ...
      '\tIn: Smeets, PCJP et al. (2022): Automatic weather station data ' ...
      'collected from 2003 to 2021 at the Greenland ice sheet along the ' ...
      'K-transect, West-Greenland [dataset publication series]. PANGAEA, ' ...
      'https://doi.org/10.1594/PANGAEA.947483\n' ...
      'Event(s):\tK-transect_%s * LATITUDE: 67.052460 * ' ...
      'LONGITUDE: -48.251800 * DATE/TIME: 2003-08-27T20:30:00 * ' ...
      'ELEVATION: 1500.0 m * LOCATION: Greenland icesheet\n' ...
      '*/\n'], site_id, year, child_doi, site_id);

   % Header variants mirror the published column sets: type0_short (AWS10,
   % no T tech), type0 (T tech present), type1 (Ice melt draw-wire added).
   % Assemble by concatenating fixed per-variant segments so no array grows.
   [melt_header, melt_value] = variantSegment(variant == "type1", ...
      "Ice melt [m]", "2.5");
   [tech_header, tech_value] = variantSegment(variant ~= "type0_short", ...
      "T tech [degC]", "-4.5");
   headers = [ ...
      "Date/Time", "dd [deg]", "ff [m/s]", "SWD [W/m**2]", ...
      "SWU [W/m**2]", "LWD [W/m**2]", "LWU [W/m**2]", "T body [degC]", ...
      "TTT [degC]", "RH [%]", "PPPP [hPa]", "Height rel [m]", ...
      melt_header, tech_header, "Vlog [V]", "ID"];
   if kwargs.drop_rh
      headers(headers == "RH [%]") = [];
   end

   % One value row per posting; the visit row carries the yearly Vlog==100
   % station-visit marker.
   lines = strings(1, numel(times));
   for k = 1:numel(times)
      vlog = "13.9";
      if kwargs.visit_row == k
         vlog = "100";
      end
      wspd = "5";
      if kwargs.calm_row == k
         wspd = "0.05";
      end
      row = [ ...
         string(datetime(times(k), 'Format', "yyyy-MM-dd'T'HH:mm")), ...
         "180", wspd, "100", "50", "250", "280", "-5.5", "-5", "80", ...
         "900", "1.0", melt_value, tech_value, vlog, awsType(variant)];
      if kwargs.blank_met
         row([4 9 10]) = "";
      end
      if kwargs.drop_rh
         row(10) = [];
      end
      lines(k) = strjoin(row, sprintf('\t'));
   end

   text = metadata + strjoin(headers, sprintf('\t')) + newline ...
      + strjoin(lines, newline) + newline;
   writeText(filename, text);
end

function [header, value] = variantSegment(include, header, value)
   %VARIANTSEGMENT Return one optional header/value pair or empty segments.
   if ~include
      header = strings(1, 0);
      value = strings(1, 0);
   end
end

function filename = writeHeightsFixture(testCase)
   %WRITEHEIGHTSFIXTURE Write a tiny sensor-height workbook fixture.
   % Mirrors the published sheet's block structure: header row, a numeric
   % station id opening each block, interleaved per-type explanation rows,
   % and per-year arrival/depart height rows.
   filename = fullfile(testCase.TestData.tmp, ...
      'metadata_GRL_AWS_56910_2003_2021.xlsx');
   % writecell rejects `missing` cells; NaN round-trips through the
   % spreadsheet as a blank for both the text and numeric predicates.
   m = NaN;
   rows = { ...
      'station', m, 'year', m, 'height u', 'height T', 'station type', 'comment'
      9, m, m, m, m, m, m, m
      m, m, m, 'For AWS type 0 the acoustic height ranger record is surface melt', m, m, m, m
      m, m, m, 'For AWS type 1 the acoustic height is the sensor height', m, m, m, m
      m, m, 2003, 'arrival/leave', 602, 550, 0, m
      m, m, 2004, 'arrival', 602, 550, 0, m
      m, m, m, 'depart', 610, 560, 0, m
      10, m, m, m, m, m, m, m
      m, m, 2010, 'arrival', 480, 450, 0, 'new station'};
   writecell(rows, filename, 'Sheet', 'sensor heights');
end

function value = awsType(variant)
   %AWSTYPE Return the generation flag for one fixture variant.
   value = "0";
   if variant == "type1"
      value = "1";
   end
end

function doi = childDoi(site_id, year)
   %CHILDDOI Return a stable per-child PANGAEA-shaped DOI for fixtures.
   % AWS9 2010 uses the real published child DOI so parser tests pin the
   % header-parsing contract against the genuine citation shape.
   if site_id == "AWS9" && year == 2010
      doi = "10.1594/PANGAEA.950093";
      return
   end
   doi = sprintf("10.1594/PANGAEA.9%04d%01d", mod(year, 10000), ...
      mod(sum(double(char(site_id))), 10));
end

function writeText(filename, text)
   %WRITETEXT Write one small parser fixture.
   fid = fopen(filename, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', text);
   clear cleaner
end
