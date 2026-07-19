function tests = test_mar_density_profiles
   %TEST_MAR_DENSITY_PROFILES Verify exact-date MAR profile source contracts.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Bootstrap the package path and build compact yearly source fixtures once.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.bootstrap_cleanup = cleanup;
   testCase.TestData.root = string(tempname);
   mkdir(testCase.TestData.root)
   writeMarFixture(testCase.TestData.root, 2019, true, true, ...
      true, true, false)
   writeMarFixture(testCase.TestData.root, 2020, false, false, ...
      true, true, false)
   writeMarFixture(testCase.TestData.root, 2022, true, false, ...
      true, true, false)
   writeMarFixture(testCase.TestData.root, 2023, true, false, ...
      false, true, false)
   writeMarFixture(testCase.TestData.root, 2024, true, false, ...
      true, false, false)
   writeMarFixture(testCase.TestData.root, 2025, true, false, ...
      false, false, false)
   writeMarFixture(testCase.TestData.root, 2026, true, true, ...
      true, true, true)
end

function teardownOnce(testCase)
   % Remove generated NetCDF sources after every test in this file completes.
   if isfolder(testCase.TestData.root)
      rmdir(testCase.TestData.root, 's')
   end
   clear testCase.TestData.bootstrap_cleanup
end

function test_reader_selects_exact_calendar_dates_and_preserves_schema(testCase)
   % Non-midnight requests select the same UTC date and remain separate groups.
   requests = [ ...
      datetime(2019, 1, 2, 18, 0, 0, TimeZone="UTC"); ...
      datetime(2019, 1, 1, 6, 0, 0, TimeZone="UTC"); ...
      datetime(2019, 1, 2, 2, 0, 0, TimeZone="UTC")];
   [profiles, status, dynamic_qa] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      requests, [1 1], source_dir=testCase.TestData.root, ...
      site_id="CEN", sample_method="nearest", ...
      requested_location=[77.0 -33.0]);

   expected_depth = marOutlay();
   expected_bounds = marOutlayBounds();
   testCase.verifyEqual(height(status), 2);
   testCase.verifyEqual(status.status, ["selected"; "selected"]);
   testCase.verifyEqual(status.selected_rows, [18; 18]);
   testCase.verifyEqual(height(profiles), 36);
   testCase.verifyEqual(numel(unique(profiles.profile_id)), 2);
   testCase.verifyEqual(unique(profiles.datetime), ...
      datetime(2019, 1, 1:2, 12, 0, 0, TimeZone="UTC")');
   for date = datetime(2019, 1, 1:2, TimeZone="UTC")
      % Every profile retains the complete source coordinate independently.
      rows = dateshift(profiles.datetime, 'start', 'day') == date;
      testCase.verifyEqual(profiles.depth(rows), expected_depth);
      testCase.verifyEqual(profiles.depth_lower_bound(rows), ...
         expected_bounds(:, 1), AbsTol=1e-12);
      testCase.verifyEqual(profiles.depth_upper_bound(rows), ...
         expected_bounds(:, 2), AbsTol=1e-12);
   end
   testCase.verifyEqual(profiles.source_variable, repmat("RO1", 36, 1));
   testCase.verifyEqual(profiles.product_role, ...
      repmat("model_output", 36, 1));
   testCase.verifyEqual(profiles.mar_version, repmat("MARv3.11", 36, 1));
   testCase.verifyTrue(all(contains( ...
      profiles.source_time_variables, "checked against YYYY/MM/DD/HH/MIN")));
   testCase.verifyEqual(profiles.depth_positive, repmat("down", 36, 1));
   testCase.verifyEqual(profiles.sample_method, repmat("nearest", 36, 1));
   testCase.verifyEqual(unique(profiles.grid_i), 1);
   testCase.verifyEqual(unique(profiles.grid_j), 1);
   testCase.verifyEqual(unique(profiles.requested_lat_wgs84), 77);
   testCase.verifyEqual(unique(profiles.requested_lon_wgs84), -33);
   testCase.verifyEqual(unique(profiles.sampled_lat_wgs84), 70);
   testCase.verifyEqual(unique(profiles.sampled_lon_wgs84), -40);
   density_index = strcmp(profiles.Properties.VariableNames, 'density');
   testCase.verifyEqual(profiles.Properties.VariableUnits{density_index}, ...
      'kg m-3');
   testCase.verifyEqual(profiles.Properties.UserData.public_variable, 'RO1');
   testCase.verifyTrue(contains( ...
      profiles.Properties.UserData.time_selection, 'no extrapolation'));
   testCase.verifyEqual(numel(dynamic_qa), 2);
   testCase.verifyEqual(dynamic_qa(1).diagnostic.status, "ok");
   testCase.verifyFalse(dynamic_qa(1).diagnostic.publishable);
end

function test_reader_never_extrapolates_and_accepts_reduced_sources(testCase)
   % Absent dates, reduced profile variables, and absent years remain distinct.
   requests = [datetime(2019, 1, 3, TimeZone="UTC"); ...
      datetime(2020, 1, 1, TimeZone="UTC"); ...
      datetime(2021, 1, 1, TimeZone="UTC")];
   [profiles, status, dynamic_qa] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      requests, [1 1], source_dir=testCase.TestData.root);

   testCase.verifyEmpty(profiles);
   testCase.verifyEmpty(dynamic_qa);
   testCase.verifyEqual(status.status, ...
      ["out_of_window"; "profile_not_available"; "source_missing"]);
   testCase.verifyTrue(contains(status.detail(1), "absent"));
   testCase.verifyTrue(contains(status.detail(2), "RO1"));
   testCase.verifyEqual(status.selected_rows, zeros(3, 1));
end

function test_reader_omits_non_ice_and_invalid_ice_profiles(testCase)
   % SRF identifies tundra zeros, while all-zero permanent ice is invalid data.
   day = datetime(2019, 1, 1, TimeZone="UTC");
   [nonice, nonice_status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      day, [2 2], source_dir=testCase.TestData.root);
   [zero_ice, zero_status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      day, [1 2], source_dir=testCase.TestData.root);

   testCase.verifyEmpty(nonice);
   testCase.verifyEmpty(zero_ice);
   testCase.verifyEqual(nonice_status.status, "non_ice");
   testCase.verifyEqual(zero_status.status, "invalid_profile");
   testCase.verifyTrue(contains(zero_status.detail, "all zero"));
end

function test_reader_rejects_density_range_and_grid_index_failures(testCase)
   % Physical range and native grid extent are both fail-closed profile checks.
   day = datetime(2019, 1, 1, TimeZone="UTC");
   [bad_density, density_status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      day, [2 1], source_dir=testCase.TestData.root);
   [bad_grid, grid_status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      day, [3 1], source_dir=testCase.TestData.root);

   testCase.verifyEmpty(bad_density);
   testCase.verifyEmpty(bad_grid);
   testCase.verifyEqual(density_status.status, "invalid_profile");
   testCase.verifyTrue(contains(density_status.detail, "250-1000"));
   testCase.verifyEqual(grid_status.status, "grid_out_of_range");
end

function test_reader_reports_missing_root_ambiguity_and_empty_requests(testCase)
   % Optional-source failures are statuses; an empty request is a typed no-op.
   day = datetime(2019, 1, 1, TimeZone="UTC");
   [missing, missing_status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      day, [1 1], source_dir=fullfile(testCase.TestData.root, "missing"));
   [empty, empty_status, empty_qa] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      NaT(0, 1), [1 1], source_dir=testCase.TestData.root);

   ambiguous_root = fullfile(testCase.TestData.root, "ambiguous");
   mkdir(ambiguous_root)
   source = fullfile(testCase.TestData.root, ...
      "MARv3.11-test-2019.nc");
   copyfile(source, fullfile(ambiguous_root, "MARv3.11-a-2019.nc"))
   copyfile(source, fullfile(ambiguous_root, "MARv3.11-b-2019.nc"))
   [ambiguous, ambiguous_status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      day, [1 1], source_dir=ambiguous_root);

   testCase.verifyEmpty(missing);
   testCase.verifyEqual(missing_status.status, "source_missing");
   testCase.verifyEmpty(empty);
   testCase.verifyEmpty(empty_status);
   testCase.verifyEmpty(empty_qa);
   testCase.verifyEqual(empty.Properties.UserData.public_variable, 'RO1');
   testCase.verifyEmpty(ambiguous);
   testCase.verifyEqual(ambiguous_status.status, "ambiguous_source");
end

function test_reader_requires_nonblank_sampling_provenance(testCase)
   % A public sampled profile cannot carry an empty sampling identity.
   day = datetime(2019, 1, 1, TimeZone="UTC");
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      day, [1 1], source_dir=testCase.TestData.root, sample_method=" "), ...
      'icemodel:forcing:readMarDensitySnapshots:blankSampleMethod');
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      day, [1 1], source_dir=testCase.TestData.root, ...
      sample_method="natural"), ...
      'icemodel:forcing:readMarDensitySnapshots:unsupportedSampleMethod');
end

function test_reader_rejects_missing_request_dates(testCase)
   % NaT cannot silently become an absent-year or nearest-date selection.
   testCase.verifyError(@() ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      NaT(1, 1), [1 1], source_dir=testCase.TestData.root), ...
      'icemodel:forcing:readMarDensitySnapshots:missingRequestedDate');
end

function test_reader_keeps_ro1_when_dynamic_fields_are_absent(testCase)
   % A full RO1/OUTLAY source remains usable when QA-only dynamic fields omit.
   [profiles, status, dynamic_qa] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      datetime(2022, 1, 1, TimeZone="UTC"), [1 1], ...
      source_dir=testCase.TestData.root);

   testCase.verifyEqual(height(profiles), 18);
   testCase.verifyEqual(status.status, "selected");
   testCase.verifyEqual(numel(dynamic_qa), 1);
   testCase.verifyEqual(dynamic_qa.diagnostic.status, "not_available");
   testCase.verifyFalse(dynamic_qa.diagnostic.available);
end

function test_reader_isolates_dynamic_source_shape_errors(testCase)
   % A malformed QA-only layer axis cannot invalidate authoritative RO1 rows.
   [profiles, status, dynamic_qa] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      datetime(2026, 1, 1, TimeZone="UTC"), [1 1], ...
      source_dir=testCase.TestData.root);

   testCase.verifyEqual(height(profiles), 18);
   testCase.verifyEqual(status.status, "selected");
   testCase.verifyEqual(dynamic_qa.diagnostic.status, "source_read_error");
   testCase.verifyTrue(contains(dynamic_qa.diagnostic.detail, ...
      "invalidDynamicShape"));
end

function test_reader_rejects_exact_public_axis_and_unit_faults(testCase)
   % Every public coordinate/RO1 axis and native unit fails before data read.
   faults = ["partial_ro1_missing", "partial_outlay_missing", ...
      "lon_missing", "lat_missing", "srf_missing", ...
      "grid_bad_x", "grid_bad_y", "time_bad_axis", ...
      "time_units_missing", "time_units_wrong", "outlay_bad_axis", ...
      "bounds_missing", "bounds_bad_axis", "bounds_bad_length", ...
      "ro1_missing_x", "ro1_missing_y", ...
      "ro1_missing_outlay", "ro1_missing_time", "ro1_extra_axis", ...
      "outlay_units_missing", "outlay_units_wrong", ...
      "ro1_units_missing", "ro1_units_wrong"];
   identifiers = [repmat("invalidPublicSchema", 1, 5), ...
      "invalidGridAxes", "invalidGridAxes", ...
      "invalidTimeAxis", "missingTimeUnits", "invalidTimeUnits", ...
      "invalidOutlayAxis", repmat("invalidOutlayBounds", 1, 3), ...
      repmat("invalidRo1Shape", 1, 5), ...
      "invalidOutlayUnits", "invalidOutlayUnits", ...
      "invalidRo1Units", "invalidRo1Units"];
   for k = 1:numel(faults)
      % Isolate one malformed yearly source so location remains unambiguous.
      root = fullfile(testCase.TestData.root, "public_" + faults(k));
      mkdir(root)
      writeMarFixture(root, 2019, true, true, true, true, false, faults(k))
      [profiles, status, dynamic_qa] = ...
         icemodel.forcing.helpers.readMarDensitySnapshots( ...
         datetime(2019, 1, 1, TimeZone="UTC"), [1 1], source_dir=root);

      testCase.verifyEmpty(profiles, faults(k));
      testCase.verifyEmpty(dynamic_qa, faults(k));
      testCase.verifyEqual(status.status, "invalid_source", faults(k));
      testCase.verifyTrue(contains(status.detail, identifiers(k)), faults(k));
   end
end

function test_reader_nests_exact_dynamic_axis_sector_and_unit_faults(testCase)
   % QA-only schema failures stay nested while authoritative RO1 remains valid.
   faults = ["dz_missing_x", "dz_missing_y", "dz_missing_time", ...
      "dz_missing_snolay", "dz_extra_axis", ...
      "dz_duplicate_axis", "rho_wrong_snolay", ...
      "shsn_missing_sector", "shsn_extra_axis", "sector_empty", ...
      "dz_units_missing", "dz_units_wrong", ...
      "rho_units_missing", "rho_units_wrong", ...
      "shsn_units_missing", "shsn_units_wrong"];
   for fault = faults
      % A fresh source exercises one exact validation branch without masking it.
      root = fullfile(testCase.TestData.root, "dynamic_" + fault);
      mkdir(root)
      writeMarFixture(root, 2019, true, true, true, true, false, fault)
      [profiles, status, dynamic_qa] = ...
         icemodel.forcing.helpers.readMarDensitySnapshots( ...
         datetime(2019, 1, 1, TimeZone="UTC"), [1 1], source_dir=root);

      testCase.verifyEqual(height(profiles), 18, fault);
      testCase.verifyEqual(status.status, "selected", fault);
      testCase.verifyEqual(dynamic_qa.diagnostic.status, ...
         "source_read_error", fault);
      expected = "invalidDynamicShape";
      if contains(fault, "units")
         % Missing and wrong attributes share the exact native-unit guard.
         expected = "invalidDynamicUnits";
      end
      diagnostic = fault + ": " + string(dynamic_qa.diagnostic.detail);
      testCase.verifyTrue(contains(dynamic_qa.diagnostic.detail, expected), ...
         diagnostic);
   end
end

function test_reader_rejects_duplicate_native_dates(testCase)
   % Duplicate native DATE values are ambiguous even when a request matches.
   duplicate_root = fullfile(testCase.TestData.root, "duplicate_dates");
   mkdir(duplicate_root)
   source = fullfile(testCase.TestData.root, ...
      "MARv3.11-test-2019.nc");
   duplicate = fullfile(duplicate_root, "MARv3.11-test-2019.nc");
   copyfile(source, duplicate)
   ncwrite(duplicate, 'TIME', [0.5; 0.5]);
   ncwrite(duplicate, 'DATE', single([2019010112; 2019010112]));
   ncwrite(duplicate, 'DD', single([1; 1]));

   [profiles, status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      datetime(2019, 1, 1, TimeZone="UTC"), [1 1], ...
      source_dir=duplicate_root);

   testCase.verifyEmpty(profiles);
   testCase.verifyEqual(status.status, "invalid_source_date");
end

function test_reader_ignores_quantized_float_date(testCase)
   % Float32 packed DATE duplicates distinct days and must never label profiles.
   packed = single([2019010512; 2019010612]);
   testCase.verifyEqual(packed(1), packed(2));
   quantized_root = fullfile(testCase.TestData.root, "quantized_date");
   mkdir(quantized_root)
   source = fullfile(testCase.TestData.root, ...
      "MARv3.11-test-2019.nc");
   quantized = fullfile(quantized_root, "MARv3.11-test-2019.nc");
   copyfile(source, quantized)
   ncwrite(quantized, 'DATE', packed);

   [profiles, status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      datetime(2019, 1, 1:2, TimeZone="UTC"), [1 1], ...
      source_dir=quantized_root);

   testCase.verifyEqual(status.status, ["selected"; "selected"]);
   testCase.verifyEqual(unique(profiles.datetime), ...
      datetime(2019, 1, 1:2, 12, 0, 0, TimeZone="UTC")');
   testCase.verifyEqual( ...
      profiles.Properties.UserData.ignored_time_variable, ...
      'DATE float32 YYYYMMDDHH is quantized and is never decoded');
end

function test_reader_requires_native_time_and_accepts_time_only(testCase)
   % TIME is authoritative; components validate it but cannot replace it.
   [component_only, component_status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      datetime(2023, 1, 1, TimeZone="UTC"), [1 1], ...
      source_dir=testCase.TestData.root);
   [time_profile, time_status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      datetime(2024, 1, 1, TimeZone="UTC"), [1 1], ...
      source_dir=testCase.TestData.root);
   [date_only, date_only_status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      datetime(2025, 1, 1, TimeZone="UTC"), [1 1], ...
      source_dir=testCase.TestData.root);

   testCase.verifyEmpty(component_only);
   testCase.verifyEqual(component_status.status, "invalid_source");
   testCase.verifyTrue(contains(component_status.detail, ...
      "invalidTimeAxis"));
   testCase.verifyEqual(time_status.status, "selected");
   testCase.verifyEqual(time_profile.datetime(1), ...
      datetime(2024, 1, 1, 12, 0, 0, TimeZone="UTC"));
   testCase.verifyTrue(startsWith( ...
      unique(time_profile.source_time_variables), "TIME (days since"));
   testCase.verifyEmpty(date_only);
   testCase.verifyEqual(date_only_status.status, "invalid_source");
   testCase.verifyTrue(contains(date_only_status.detail, ...
      "invalidTimeAxis"));
end

function test_reader_rejects_time_component_mismatch_and_partial_set(testCase)
   % Mixed or stale components cannot silently relabel the native TIME axis.
   mismatch_root = fullfile(testCase.TestData.root, "time_mismatch");
   mkdir(mismatch_root)
   full_source = fullfile(testCase.TestData.root, ...
      "MARv3.11-test-2019.nc");
   mismatch = fullfile(mismatch_root, "MARv3.11-test-2019.nc");
   copyfile(full_source, mismatch)
   ncwrite(mismatch, 'DD', single([1; 3]));

   partial_root = fullfile(testCase.TestData.root, "partial_time");
   mkdir(partial_root)
   time_source = fullfile(testCase.TestData.root, ...
      "MARv3.11-test-2024.nc");
   partial = fullfile(partial_root, "MARv3.11-test-2024.nc");
   copyfile(time_source, partial)
   nccreate(partial, 'YYYY', 'Dimensions', {'TIME', 2}, ...
      'Datatype', 'single');
   ncwrite(partial, 'YYYY', single([2024; 2024]));

   [mismatch_profile, mismatch_status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      datetime(2019, 1, 1, TimeZone="UTC"), [1 1], ...
      source_dir=mismatch_root);
   [partial_profile, partial_status] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      datetime(2024, 1, 1, TimeZone="UTC"), [1 1], ...
      source_dir=partial_root);

   testCase.verifyEmpty(mismatch_profile);
   testCase.verifyEqual(mismatch_status.status, "invalid_source");
   testCase.verifyTrue(contains(mismatch_status.detail, "timeMismatch"));
   testCase.verifyEmpty(partial_profile);
   testCase.verifyEqual(partial_status.status, "invalid_source");
   testCase.verifyTrue(contains(partial_status.detail, ...
      "partialTimeComponents"));
end

function test_dynamic_qa_reverses_layers_and_reconstructs(testCase)
   % Synthetic source order is bottom-to-surface and reverses to positive down.
   dz = 0.1 * ones(30, 1);
   rho = (301:330)';
   midpoint = (0.05:0.1:2.95)';
   ro1 = flipud(rho);
   qa = icemodel.forcing.helpers.marDynamicProfileQa( ...
      dz, rho, 3, ro1, midpoint);

   testCase.verifyEqual(qa.status, "ok");
   testCase.verifyTrue(qa.available);
   testCase.verifyTrue(qa.activity_mask_match);
   testCase.verifyEqual(qa.active_layer_count, 30);
   testCase.verifyEqual(qa.active_thickness_sum_m, 3, AbsTol=1e-12);
   testCase.verifyEqual(qa.thickness_residual_m, 0, AbsTol=1e-12);
   testCase.verifyEqual(qa.surface_down_source_indices, 30:-1:1);
   testCase.verifyEqual(qa.dynamic_depth_midpoint_m, midpoint, AbsTol=1e-12);
   testCase.verifyEqual(qa.dynamic_density_kg_m3, ro1);
   testCase.verifyEqual(qa.reconstruction_overlap_count, 30);
   testCase.verifyEqual(qa.reconstruction_rmse_kg_m3, 0, AbsTol=1e-12);
   testCase.verifyFalse(qa.publishable);
   testCase.verifyFalse(qa.correction_applied);
end

function test_dynamic_qa_classifies_unavailable_shape_mask_and_thickness(testCase)
   % Status precedence distinguishes absence, shape, source, mask, and thickness.
   unavailable = icemodel.forcing.helpers.marDynamicProfileQa( ...
      [], [], [], [], []);
   shape = icemodel.forcing.helpers.marDynamicProfileQa( ...
      ones(2, 1), ones(3, 1), 2, [], []);
   mismatch_density = [0; 300];
   mask = icemodel.forcing.helpers.marDynamicProfileQa( ...
      [1; 1], mismatch_density, 1, [], []);
   thickness = icemodel.forcing.helpers.marDynamicProfileQa( ...
      [1; 1], [300; 400], 3, [], []);
   invalid = icemodel.forcing.helpers.marDynamicProfileQa( ...
      [-1; 1], [300; 400], 1, [], []);

   testCase.verifyEqual(unavailable.status, "not_available");
   testCase.verifyEqual(shape.status, "shape_mismatch");
   testCase.verifyEqual(mask.status, "activity_mask_mismatch");
   testCase.verifyEqual(mask.activity_pair_mismatch_count, 1);
   testCase.verifyEqual(thickness.status, "thickness_mismatch");
   testCase.verifyEqual(thickness.thickness_residual_m, -1);
   testCase.verifyEqual(invalid.status, "invalid_values");
end

function test_full_mar_archive_exposes_profile_and_dynamic_headers(testCase)
   % Guarded source-backed proof covers every full RUH2 yearly file, not fixture.
   source_dir = string(getenv("ICEMODEL_MAR_DIR"));
   if source_dir == ""
      source_dir = "/Volumes/S03/DATA/greenland/mar3p11/RUH2";
   end
   testCase.assumeTrue(isfolder(source_dir), ...
      "full MAR RUH2 source archive is unavailable");
   files = dir(fullfile(source_dir, 'MARv3.11*.nc'));
   testCase.assumeEqual(numel(files), 40, ...
      "source-backed contract requires all 40 RUH2 yearly files");

   expected_outlay = marOutlay();
   for k = 1:numel(files)
      % Every header must retain the same exact public and QA-only contracts.
      filename = fullfile(files(k).folder, files(k).name);
      info = ncinfo(filename);
      names = string({info.Variables.Name});
      testCase.verifyTrue(all(ismember( ...
         ["TIME", "YYYY", "MM", "DD", "HH", "MIN", "RO1", "OUTLAY", ...
          "OUTLAY_bnds", "DZSN1", "ROSN1", "SHSN3"], ...
          names)), files(k).name);
      lon_info = ncinfo(filename, 'LON');
      grid_dims = string({lon_info.Dimensions.Name});
      testCase.verifyEqual(numel(grid_dims), 2, files(k).name);
      testCase.verifyEqual(sum(startsWith(upper(grid_dims), "X")), 1, ...
         files(k).name);
      testCase.verifyEqual(sum(startsWith(upper(grid_dims), "Y")), 1, ...
         files(k).name);
      verifyExactVariable(testCase, filename, 'LON', grid_dims, [], ...
         files(k).name);
      verifyExactVariable(testCase, filename, 'LAT', grid_dims, [], ...
         files(k).name);
      verifyExactVariable(testCase, filename, 'SRF', grid_dims, [], ...
         files(k).name);
      verifyExactVariable(testCase, filename, 'TIME', "TIME", [], ...
         files(k).name);
      verifyExactVariable(testCase, filename, 'OUTLAY', "OUTLAY", "m", ...
         files(k).name);
      verifyExactVariable(testCase, filename, 'RO1', ...
         [grid_dims, "OUTLAY", "TIME"], ["kg/m3", "kg m-3"], ...
         files(k).name);
      verifyExactVariable(testCase, filename, 'DZSN1', ...
         [grid_dims, "SNOLAY", "TIME"], "m", files(k).name);
      verifyExactVariable(testCase, filename, 'ROSN1', ...
         [grid_dims, "SNOLAY", "TIME"], ["kg/m3", "kg m-3"], ...
         files(k).name);
      shsn_info = verifyExactVariable(testCase, filename, 'SHSN3', ...
         [grid_dims, "SECTOR", "TIME"], "m", files(k).name);
      shsn_names = string({shsn_info.Dimensions.Name});
      sector_index = find(strcmpi(shsn_names, "SECTOR"));
      testCase.verifyEqual(numel(sector_index), 1, files(k).name);
      if isscalar(sector_index)
         % Permanent-ice diagnostics require a readable sector-one element.
         testCase.verifyGreaterThanOrEqual( ...
            shsn_info.Dimensions(sector_index).Length, 1, files(k).name);
      end
      dz_info = ncinfo(filename, 'DZSN1');
      rho_info = ncinfo(filename, 'ROSN1');
      dz_names = string({dz_info.Dimensions.Name});
      rho_names = string({rho_info.Dimensions.Name});
      testCase.verifyEqual(dz_info.Dimensions( ...
         strcmpi(dz_names, "SNOLAY")).Length, 30, files(k).name);
      testCase.verifyEqual(rho_info.Dimensions( ...
         strcmpi(rho_names, "SNOLAY")).Length, 30, files(k).name);
      outlay = double(ncread(filename, 'OUTLAY'));
      bounds = double(ncread(filename, 'OUTLAY_bnds'));
      testCase.verifyEqual(outlay(:), expected_outlay, AbsTol=1e-6);
      testCase.verifyEqual(sort(size(squeeze(bounds))), [2 18]);
      verifyNativeTime(testCase, filename, files(k).name)
   end
end

function test_reader_proves_asymmetric_real_2019_orientation_and_sector(testCase)
   % Guarded public-reader call proves native orientation, noon, and sector 1.
   source_dir = string(getenv("ICEMODEL_MAR_DIR"));
   if source_dir == ""
      source_dir = "/Volumes/S03/DATA/greenland/mar3p11/RUH2";
   end
   testCase.assumeTrue(isfolder(source_dir), ...
      "full MAR RUH2 source archive is unavailable");
   files = dir(fullfile(source_dir, '*-2019.nc'));
   testCase.assumeEqual(numel(files), 1, ...
      "source-backed reader proof requires one 2019 yearly file");

   [profiles, status, dynamic_qa] = ...
      icemodel.forcing.helpers.readMarDensitySnapshots( ...
      datetime(2019, 1, 1, TimeZone="UTC"), [49 91], ...
      source_dir=source_dir, site_id="real_2019");

   testCase.verifyEqual(status.status, "selected");
   testCase.verifyEqual(height(profiles), 18);
   testCase.verifyEqual(unique(profiles.datetime), ...
      datetime(2019, 1, 1, 12, 0, 0, TimeZone="UTC"));
   testCase.verifyEqual(profiles.depth, marOutlay(), AbsTol=1e-6);
   bounds = marOutlayBounds();
   testCase.verifyEqual(profiles.depth_lower_bound, bounds(:, 1), ...
      AbsTol=1e-6);
   testCase.verifyEqual(profiles.depth_upper_bound, bounds(:, 2), ...
      AbsTol=1e-6);
   testCase.verifyEqual(unique(profiles.grid_i), 49);
   testCase.verifyEqual(unique(profiles.grid_j), 91);
   testCase.verifyEqual(unique(profiles.sampled_lat_wgs84), 71.98178, ...
      AbsTol=1e-4);
   testCase.verifyEqual(unique(profiles.sampled_lon_wgs84), -40.87207, ...
      AbsTol=1e-4);
   diagnostic = dynamic_qa.diagnostic;
   testCase.verifyEqual(diagnostic.status, "ok");
   testCase.verifyEqual(diagnostic.active_layer_count, 29);
   testCase.verifyEqual(diagnostic.active_thickness_sum_m, 25.39769788, ...
      AbsTol=5e-6);
   testCase.verifyEqual(diagnostic.shsn3_permanent_ice_m, 25.3977, ...
      AbsTol=5e-6);
   testCase.verifyLessThanOrEqual(abs(diagnostic.thickness_residual_m), 2e-5);
end

function info = verifyExactVariable(testCase, filename, variable, ...
      expected_dimensions, accepted_units, diagnostic)
   %VERIFYEXACTVARIABLE Check one complete axis set and optional native units.
   info = ncinfo(filename, variable);
   actual = string({info.Dimensions.Name});
   expected_dimensions = string(expected_dimensions);
   testCase.verifyEqual(numel(actual), numel(expected_dimensions), diagnostic);
   testCase.verifyEqual(numel(unique(lower(actual))), numel(actual), diagnostic);
   testCase.verifyTrue(all(ismember( ...
      lower(expected_dimensions), lower(actual))), diagnostic);
   if ~isempty(accepted_units)
      % Unit spellings are source-faithful, not inferred from stamped output.
      units = strtrim(string(ncreadatt(filename, variable, 'units')));
      testCase.verifyTrue(any(units == string(accepted_units)), diagnostic);
   end
end

function verifyNativeTime(testCase, filename, diagnostic)
   %VERIFYNATIVETIME Prove TIME and exact components agree for one source year.
   units = string(ncreadatt(filename, 'TIME', 'units'));
   match = regexp(units, '^days since (.+)$', 'tokens', 'once');
   testCase.verifyNotEmpty(match, diagnostic);
   origin = datetime(string(match{1}), ...
      'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
   native = origin + days(double(ncread(filename, 'TIME')));
   components = datetime( ...
      double(ncread(filename, 'YYYY')), ...
      double(ncread(filename, 'MM')), ...
      double(ncread(filename, 'DD')), ...
      double(ncread(filename, 'HH')), ...
      double(ncread(filename, 'MIN')), ...
      zeros(size(native)), 'TimeZone', 'UTC');

   testCase.verifyEqual(native(:), components(:), diagnostic);
   testCase.verifyEqual(unique(hour(native)), 12, diagnostic);
   testCase.verifyEqual(unique(minute(native)), 0, diagnostic);
   testCase.verifyTrue(all(diff(native) == days(1)), diagnostic);
end

function test_stage_helper_is_additive_and_preserves_failed_refresh(testCase)
   % Sequential exact-date requests add profile groups without changing the
   % MAR forcing-ready contract; a later unavailable date preserves both.
   stage_root = fullfile(testCase.TestData.root, 'stage-additive');
   family_root = fullfile(stage_root, 'eval', 'sumup');
   userdata_root = fullfile(stage_root, 'input', 'userdata');
   mkdir(fullfile(family_root, 'case1'))
   mkdir(fullfile(userdata_root, 'mar3.11'))
   state = profileStageState(family_root, userdata_root);

   writeSumupTarget(family_root, datetime(2019, 1, 1, TimeZone="UTC"))
   state = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   writeSumupTarget(family_root, datetime(2019, 1, 2, TimeZone="UTC"))
   state = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);

   output_file = fullfile(userdata_root, ...
      state.colocation.mar.model_output_files);
   loaded = load(output_file, 'reference');
   groups = icemodel.verification.helpers.profileGroups( ...
      loaded.reference.data.density);
   testCase.verifyEqual(numel(groups), 2)
   testCase.verifyEqual(string(state.colocation.mar.model_output_format), ...
      "subsurface_profile_bundle")
   testCase.verifyEqual(string(state.colocation.mar.model_output_variables), ...
      "density")
   testCase.verifyTrue(state.colocation.mar.forcing_ready)

   % A source-missing requested year leaves the inherited sidecar and its
   % metadata atomic rather than relabelling only the manifest status.
   prior_leg = state.colocation.mar;
   writeSumupTarget(family_root, datetime(2021, 1, 1, TimeZone="UTC"))
   preserved = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   reloaded = load(output_file, 'reference');
   testCase.verifyEqual(numel( ...
      icemodel.verification.helpers.profileGroups( ...
      reloaded.reference.data.density)), 2)
   testCase.verifyEqual(preserved.colocation.mar, prior_leg)
end

function test_stage_helper_replay_preserves_inherited_sidecar(testCase)
   % A forcing-only replay without staged observations must not rewrite the
   % inherited sidecar group on either the first or a repeated callback.
   stage_root = fullfile(testCase.TestData.root, 'stage-replay');
   family_root = fullfile(stage_root, 'eval', 'sumup');
   userdata_root = fullfile(stage_root, 'input', 'userdata');
   mkdir(fullfile(family_root, 'case1'))
   mkdir(fullfile(userdata_root, 'mar3.11'))
   state = profileStageState(family_root, userdata_root);
   state.colocation.mar.model_output_files = ...
      'mar3.11/case1_mar3.11_density_profiles.mat';
   state.colocation.mar.model_output_format = 'subsurface_profile_bundle';
   state.colocation.mar.model_output_variables = {'density'};
   state.colocation.mar.model_output_status = 'staged';
   state.colocation.mar.model_output_note = ...
      '4 MAR RO1 profile group(s), exact SUMup UTC dates.';

   first = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   second = icemodel.verification.setup.stageMarDensityProfiles( ...
      first, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);

   testCase.verifyEqual(first.colocation.mar, state.colocation.mar)
   testCase.verifyEqual(second.colocation.mar, state.colocation.mar)
end

function test_stage_helper_keeps_optional_guard_failures_non_destructive(testCase)
   % Non-MAR callbacks, absent legs, unstaged legs, non-nearest sampling,
   % missing observations, and missing cell provenance all leave primary data.
   stage_root = fullfile(testCase.TestData.root, 'stage-guards');
   family_root = fullfile(stage_root, 'eval', 'sumup');
   userdata_root = fullfile(stage_root, 'input', 'userdata');
   mkdir(fullfile(family_root, 'case1'))
   mkdir(fullfile(userdata_root, 'mar3.11'))
   state = profileStageState(family_root, userdata_root);
   original = state;

   other_source = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "racmo", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   no_live = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, zeros(1, 0), "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   testCase.verifyEqual(other_source, original)
   testCase.verifyEqual(no_live, original)

   no_leg = rmfield(state, 'colocation');
   no_leg.colocation = struct();
   returned = icemodel.verification.setup.stageMarDensityProfiles( ...
      no_leg, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   testCase.verifyEqual(returned, no_leg)

   state.colocation.mar.staged = false;
   unstaged = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   testCase.verifyFalse(unstaged.colocation.mar.staged)
   state = original;

   state.colocation.mar.sample_method = 'natural';
   natural = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   testCase.verifyEqual(string(natural.colocation.mar.model_output_status), ...
      "profile_not_available")
   state = original;

   % Once an observation exists, missing data-artifact provenance remains an
   % explicit optional status rather than a forcing failure.
   writeSumupTarget(family_root, datetime(2019, 1, 1, TimeZone="UTC"))
   delete(fullfile(userdata_root, state.colocation.mar.data_files))
   missing_provenance = ...
      icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   testCase.verifyEqual(string( ...
      missing_provenance.colocation.mar.model_output_status), ...
      "missing_grid_provenance")
end

function test_stage_helper_preserves_or_replaces_incompatible_sidecar(testCase)
   % An incompatible prior table is never partially merged: additive mode keeps
   % it byte-for-byte, while explicit overwrite replaces it with public RO1.
   stage_root = fullfile(testCase.TestData.root, 'stage-incompatible');
   family_root = fullfile(stage_root, 'eval', 'sumup');
   userdata_root = fullfile(stage_root, 'input', 'userdata');
   mkdir(fullfile(family_root, 'case1'))
   mkdir(fullfile(userdata_root, 'mar3.11'))
   state = profileStageState(family_root, userdata_root);
   writeSumupTarget(family_root, datetime(2019, 1, 1, TimeZone="UTC"))

   output_file = fullfile(userdata_root, 'mar3.11', ...
      'case1_mar3.11_density_profiles.mat');
   reference = struct('format', 'subsurface_profile_bundle', ...
      'data', struct('density', table([1; 2], ...
      'VariableNames', {'wrong'})));
   save(output_file, 'reference')
   original_hash = icemodel.verification.setup.fileSha256(output_file);

   preserved = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   testCase.verifyEqual(string( ...
      preserved.colocation.mar.model_output_status), ...
      "incompatible_existing_sidecar")
   testCase.verifyEqual( ...
      icemodel.verification.setup.fileSha256(output_file), original_hash)

   replaced = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root, ...
      overwrite=true);
   loaded = load(output_file, 'reference');
   testCase.verifyEqual(string(replaced.colocation.mar.model_output_status), ...
      "staged")
   testCase.verifyTrue(ismember('profile_id', ...
      loaded.reference.data.density.Properties.VariableNames))
end

function test_stage_helper_classifies_optional_input_guards(testCase)
   % Missing source products and unusable observation dates remain optional;
   % exact saved grid_start provenance must still support a valid stage.
   stage_root = fullfile(testCase.TestData.root, 'stage-input-guards');
   family_root = fullfile(stage_root, 'eval', 'sumup');
   userdata_root = fullfile(stage_root, 'input', 'userdata');
   empty_mar_root = fullfile(stage_root, 'empty-mar');
   mkdir(fullfile(family_root, 'case1'))
   mkdir(fullfile(userdata_root, 'mar3.11'))
   mkdir(empty_mar_root)
   state = profileStageState(family_root, userdata_root);

   source_missing = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=empty_mar_root);
   testCase.verifyEqual(string( ...
      source_missing.colocation.mar.model_output_status), ...
      "profile_not_available")
   testCase.verifyTrue(source_missing.colocation.mar.forcing_ready)
   original_state = state;
   inherited_state = state;
   inherited_state.colocation.mar.model_output_files = ...
      'mar3.11/case1_mar3.11_density_profiles.mat';
   inherited_state.colocation.mar.model_output_format = ...
      'subsurface_profile_bundle';
   inherited_state.colocation.mar.model_output_variables = {'density'};
   inherited_state.colocation.mar.model_output_status = 'staged';
   inherited_state.colocation.mar.model_output_note = ...
      'Existing exact-date MAR density profiles.';
   preserved_source_missing = ...
      icemodel.verification.setup.stageMarDensityProfiles( ...
      inherited_state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=empty_mar_root);
   testCase.verifyEqual(preserved_source_missing.colocation.mar, ...
      inherited_state.colocation.mar)
   state = original_state;

   observation_missing = ...
      icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   testCase.verifyEqual(string( ...
      observation_missing.colocation.mar.model_output_status), ...
      "profile_not_available")

   density = table([0.5; 1.5], [350; 500], ...
      'VariableNames', {'midpoint', 'density'});
   targets = struct('format', 'subsurface_profile_bundle', ...
      'data', struct('density', density));
   observation_file = fullfile(family_root, 'case1', 'observations.mat');
   save(observation_file, 'targets')
   no_date = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   testCase.verifyEqual(string(no_date.colocation.mar.model_output_status), ...
      "profile_not_available")

   writeSumupTarget(family_root, NaT(1, 1, TimeZone="UTC"))
   no_finite_date = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   testCase.verifyEqual(string( ...
      no_finite_date.colocation.mar.model_output_status), ...
      "profile_not_available")

   % Replace coordinate fallback provenance with the exact native cell index.
   artifact_file = fullfile(userdata_root, state.colocation.mar.data_files);
   artifact_metadata = struct('grid_start', [1 1], ...
      'grid_count', [1 1]);
   save(artifact_file, 'artifact_metadata')
   writeSumupTarget(family_root, datetime(2019, 1, 1, TimeZone="UTC"))
   staged = icemodel.verification.setup.stageMarDensityProfiles( ...
      state, 1, "mar", family_root=family_root, ...
      userdata_outdir=userdata_root, mar_dir=testCase.TestData.root);
   testCase.verifyEqual(string(staged.colocation.mar.model_output_status), ...
      "staged")
end

function state = profileStageState(~, userdata_root)
   %PROFILESTAGESTATE Create one nearest-cell staged MAR test record.
   data_file = fullfile('mar3.11', 'case1_mar3.11_data.mat');
   artifact_metadata = struct('sample_method', 'nearest', ...
      'source_lat_wgs84', 70, 'source_lon_wgs84', -40);
   save(fullfile(userdata_root, data_file), 'artifact_metadata')
   state = struct('case_id', "case1", 'site_id', "CEN", ...
      'point', [77 -33], ...
      'evaluation_file_rel', fullfile('case1', 'observations.mat'), ...
      'colocation', struct('mar', struct('staged', true, ...
      'sample_method', 'nearest', 'data_files', data_file, ...
      'forcing_ready', true)));
end

function writeSumupTarget(family_root, profile_date)
   %WRITESUMUPTARGET Write one compact dated density observation profile.
   density = table(repmat("core", 2, 1), repmat(profile_date, 2, 1), ...
      [0.5; 1.5], [350; 500], ...
      'VariableNames', {'name', 'datetime', 'midpoint', 'density'});
   targets = struct('format', 'subsurface_profile_bundle', ...
      'data', struct('format', 'subsurface_profile_bundle', ...
      'density', density));
   save(fullfile(family_root, 'case1', 'observations.mat'), 'targets')
end

function writeTimeComponent(filename, name, values, time_dimension)
   %WRITETIMECOMPONENT Write one exact integer-valued source time channel.
   nccreate(filename, name, 'Dimensions', {time_dimension, 2}, ...
      'Datatype', 'single');
   ncwrite(filename, name, values);
end

function writeMarFixture(root, yyyy, include_profile, include_dynamic, ...
      include_time, include_components, malformed_dynamic, fault)
   %WRITEMARFIXTURE Create a two-day, two-by-two dimension-rich MAR source.
   if nargin < 8
      % Existing behavioral fixtures use the complete valid native schema.
      fault = "";
   end
   fault = string(fault);
   filename = fullfile(root, compose("MARv3.11-test-%d.nc", yyyy));
   time_dimension = 'TIME';
   lon_dimensions = {'X', 2, 'Y', 2};
   if fault == "grid_bad_x"
      % A non-X first axis proves grid axes are validated by name, not count.
      lon_dimensions = {'BADX', 2, 'Y', 2};
   elseif fault == "grid_bad_y"
      % A non-Y second axis covers the other required native grid coordinate.
      lon_dimensions = {'X', 2, 'BADY', 2};
   end
   if fault ~= "lon_missing"
      nccreate(filename, 'LON', 'Dimensions', lon_dimensions);
      ncwrite(filename, 'LON', [-40 -39; -38 -37]);
   end
   if fault ~= "lat_missing"
      nccreate(filename, 'LAT', 'Dimensions', {'X', 2, 'Y', 2});
      ncwrite(filename, 'LAT', [70 71; 72 73]);
   end
   if fault ~= "srf_missing"
      nccreate(filename, 'SRF', 'Dimensions', {'X', 2, 'Y', 2});
      ncwrite(filename, 'SRF', [4 4; 4 1]);
   end
   if include_time
      % Native elapsed TIME preserves noon support exactly.
      if fault == "time_bad_axis"
         % The coordinate variable name cannot legitimize an unknown axis.
         time_dimension = 'BADTIME';
      end
      nccreate(filename, 'TIME', 'Dimensions', {time_dimension, 2});
      ncwrite(filename, 'TIME', [0.5; 1.5]);
      if fault ~= "time_units_missing"
         time_units = char(compose( ...
            'days since %d-01-01 00:00:00', yyyy));
         if fault == "time_units_wrong"
            % A present but unsupported elapsed unit must fail closed too.
            time_units = char(compose( ...
               'fortnights since %d-01-01 00:00:00', yyyy));
         end
         ncwriteatt(filename, 'TIME', 'units', time_units);
      end
   end
   nccreate(filename, 'DATE', 'Dimensions', {time_dimension, 2}, ...
      'Datatype', 'single');
   packed_date = single([yyyy * 1e6 + 10112; yyyy * 1e6 + 10212]);
   ncwrite(filename, 'DATE', packed_date);
   if include_components
      % Small integer components remain exact even though source storage is float.
      writeTimeComponent(filename, 'YYYY', single([yyyy; yyyy]), time_dimension)
      writeTimeComponent(filename, 'MM', single([1; 1]), time_dimension)
      writeTimeComponent(filename, 'DD', single([1; 2]), time_dimension)
      writeTimeComponent(filename, 'HH', single([12; 12]), time_dimension)
      writeTimeComponent(filename, 'MIN', single([0; 0]), time_dimension)
   end
   if ~include_profile
      return
   end

   % Fixed source depths and bounds intentionally preserve the negative surface
   % half-cell edge while the coordinate itself remains positive down.
   outlay = marOutlay();
   bounds = marOutlayBounds();
   outlay_dimension = 'OUTLAY';
   if fault == "outlay_bad_axis"
      % OUTLAY values on a differently named axis are not the public product.
      outlay_dimension = 'BADOUTLAY';
   end
   if fault ~= "partial_outlay_missing"
      nccreate(filename, 'OUTLAY', 'Dimensions', {outlay_dimension, 18});
      ncwrite(filename, 'OUTLAY', outlay);
      if fault ~= "outlay_units_missing"
         % Wrong-but-present units exercise the same required source attribute.
         outlay_units = 'm';
         if fault == "outlay_units_wrong"
            outlay_units = 'cm';
         end
         ncwriteatt(filename, 'OUTLAY', 'units', outlay_units);
      end
   end
   if fault ~= "bounds_missing" && fault ~= "partial_outlay_missing"
      bounds_dimensions = {'BND', 2, outlay_dimension, 18};
      bounds_data = bounds';
      if fault == "bounds_bad_axis"
         % A bounds array must use the authoritative OUTLAY dimension itself.
         bounds_dimensions = {'BND', 2, 'BADOUTLAY', 18};
      elseif fault == "bounds_bad_length"
         % Exactly two endpoints are required for every retained depth level.
         bounds_dimensions = {'BND', 3, 'OUTLAY', 18};
         bounds_data = zeros(3, 18);
      end
      nccreate(filename, 'OUTLAY_bnds', ...
         'Dimensions', bounds_dimensions);
      ncwrite(filename, 'OUTLAY_bnds', bounds_data);
      ncwriteatt(filename, 'OUTLAY', 'bounds', 'OUTLAY_bnds');
   end

   % Four cells exercise valid, out-of-range, all-zero ice, and non-ice status.
   ro1 = zeros(2, 2, 18, 2);
   for t = 1:2
      valid = (300:10:470)' + 5 * (t - 1);
      invalid = valid;
      invalid(10) = 1001;
      ro1(1, 1, :, t) = reshape(valid, 1, 1, 18, 1);
      ro1(2, 1, :, t) = reshape(invalid, 1, 1, 18, 1);
   end
   ro1_dimensions = { ...
      'X', 2, 'Y', 2, outlay_dimension, 18, time_dimension, 2};
   ro1_data = ro1;
   switch fault
      case "ro1_missing_x"
         % Each missing required public axis receives an isolated fixture.
         ro1_dimensions = {'Y', 2, 'OUTLAY', 18, 'TIME', 2};
         ro1_data = squeeze(ro1(1, :, :, :));
      case "ro1_missing_y"
         ro1_dimensions = {'X', 2, 'OUTLAY', 18, 'TIME', 2};
         ro1_data = squeeze(ro1(:, 1, :, :));
      case "ro1_missing_outlay"
         ro1_dimensions = {'X', 2, 'Y', 2, 'TIME', 2};
         ro1_data = squeeze(ro1(:, :, 1, :));
      case "ro1_missing_time"
         ro1_dimensions = {'X', 2, 'Y', 2, 'OUTLAY', 18};
         ro1_data = ro1(:, :, :, 1);
      case "ro1_extra_axis"
         % A populated extra axis proves unknown dimensions never default.
         ro1_dimensions = [ro1_dimensions, {'EXTRA', 2}];
         ro1_data = repmat(ro1, [1 1 1 1 2]);
   end
   if fault ~= "partial_ro1_missing"
      nccreate(filename, 'RO1', 'Dimensions', ro1_dimensions);
      ncwrite(filename, 'RO1', ro1_data);
      if fault ~= "ro1_units_missing"
         ro1_units = 'kg m-3';
         if fault == "ro1_units_wrong"
            ro1_units = 'kg m-2';
         end
         ncwriteatt(filename, 'RO1', 'units', ro1_units);
      end
   end
   if ~include_dynamic
      return
   end

   % Native dynamic layers are stored bottom-to-surface; sector 2 is a sentinel
   % proving the reader selects permanent-ice SHSN3 sector 1 only.
   dz = zeros(2, 2, 30, 2);
   rho = zeros(2, 2, 30, 2);
   shsn3 = zeros(2, 2, 2, 2);
   for t = 1:2
      dz(1, 1, :, t) = 0.1;
      rho(1, 1, :, t) = reshape((301:330)' + t - 1, 1, 1, 30, 1);
      shsn3(1, 1, 1, t) = 3;
      shsn3(1, 1, 2, t) = 99;
   end
   dz_dimensions = {'X', 2, 'Y', 2, 'SNOLAY', 30, time_dimension, 2};
   dz_data = dz;
   switch fault
      case "dz_missing_x"
         % Each required dynamic axis receives an isolated malformed source.
         dz_dimensions = {'Y', 2, 'SNOLAY', 30, 'TIME', 2};
         dz_data = squeeze(dz(1, :, :, :));
      case "dz_missing_y"
         dz_dimensions = {'X', 2, 'SNOLAY', 30, 'TIME', 2};
         dz_data = squeeze(dz(:, 1, :, :));
      case "dz_missing_time"
         dz_dimensions = {'X', 2, 'Y', 2, 'SNOLAY', 30};
         dz_data = dz(:, :, :, 1);
      case "dz_missing_snolay"
         % Missing, extra, and duplicated axes exercise distinct exact-set gaps.
         dz_dimensions = {'X', 2, 'Y', 2, 'TIME', 2};
         dz_data = squeeze(dz(:, :, 1, :));
      case "dz_extra_axis"
         dz_dimensions = [dz_dimensions, {'EXTRA', 2}];
         dz_data = repmat(dz, [1 1 1 1 2]);
      case "dz_duplicate_axis"
         dz_dimensions = {'X', 2, 'Y', 2, 'SNOLAY', 30, ...
            'SNOLAY', 30, 'TIME', 2};
         dz_data = repmat(reshape(dz, [2 2 30 1 2]), [1 1 1 30 1]);
   end
   nccreate(filename, 'DZSN1', 'Dimensions', dz_dimensions);
   ncwrite(filename, 'DZSN1', dz_data);
   if fault ~= "dz_units_missing"
      dz_units = 'm';
      if fault == "dz_units_wrong"
         dz_units = 'cm';
      end
      ncwriteatt(filename, 'DZSN1', 'units', dz_units);
   end
   density_layer_name = 'SNOLAY';
   if malformed_dynamic || fault == "rho_wrong_snolay"
      density_layer_name = 'OTHERLAY';
   end
   nccreate(filename, 'ROSN1', 'Dimensions', ...
      {'X', 2, 'Y', 2, density_layer_name, 30, time_dimension, 2});
   ncwrite(filename, 'ROSN1', rho);
   if fault ~= "rho_units_missing"
      rho_units = 'kg m-3';
      if fault == "rho_units_wrong"
         rho_units = 'kg m-2';
      end
      ncwriteatt(filename, 'ROSN1', 'units', rho_units);
   end
   shsn_dimensions = {'X', 2, 'Y', 2, 'SECTOR', 2, time_dimension, 2};
   shsn_data = shsn3;
   write_shsn = true;
   switch fault
      case "shsn_missing_sector"
         % Missing/extra sector shapes must not return the first scalar silently.
         shsn_dimensions = {'X', 2, 'Y', 2, time_dimension, 2};
         shsn_data = squeeze(shsn3(:, :, 1, :));
      case "shsn_extra_axis"
         shsn_dimensions = [shsn_dimensions, {'EXTRA', 2}];
         shsn_data = repmat(shsn3, [1 1 1 1 2]);
      case "sector_empty"
         % An empty unlimited dimension has no permanent-ice sector-one entry.
         % Put the sole unlimited axis last for the default netcdf4-classic file.
         shsn_dimensions = {'X', 2, 'Y', 2, time_dimension, 2, 'SECTOR', Inf};
         write_shsn = false;
   end
   nccreate(filename, 'SHSN3', 'Dimensions', shsn_dimensions);
   if write_shsn
      ncwrite(filename, 'SHSN3', shsn_data);
   end
   if fault ~= "shsn_units_missing"
      shsn_units = 'm';
      if fault == "shsn_units_wrong"
         shsn_units = 'cm';
      end
      ncwriteatt(filename, 'SHSN3', 'units', shsn_units);
   end
end

function outlay = marOutlay()
   %MAROUTLAY Return the exact fixed-depth coordinate in full RUH2 sources.
   outlay = [0; 0.05; 0.1; 0.2; 0.3; 0.4; 0.5; 0.65; 0.8; 1; ...
      1.5; 2; 3; 5; 7.5; 10; 15; 20];
end

function bounds = marOutlayBounds()
   %MAROUTLAYBOUNDS Return the exact full-RUH2 fixed-depth interval edges.
   bounds = [-0.025 0.025; 0.025 0.075; 0.075 0.15; 0.15 0.25; ...
      0.25 0.35; 0.35 0.45; 0.45 0.575; 0.575 0.725; 0.725 0.9; ...
      0.9 1.25; 1.25 1.75; 1.75 2.5; 2.5 4; 4 6.25; 6.25 8.75; ...
      8.75 12.5; 12.5 17.5; 17.5 22.5];
end
