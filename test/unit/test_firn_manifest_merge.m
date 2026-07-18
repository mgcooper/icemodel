function tests = test_firn_manifest_merge
   %TEST_FIRN_MANIFEST_MERGE Verify incremental family-manifest staging.
   %
   % Exercises icemodel.verification.setup.writeFamilyManifestMerge directly
   % (no raw forcing sources needed, so it never data-gates): staging a NEW
   % site into a family manifest that already holds a committed site must ADD
   % the new case and PRESERVE the existing case + its sidecar files byte for
   % byte. This is the KAN no-churn guarantee the firn importers rely on.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   tmp = tempname;
   mkdir(tmp)
   testCase.TestData.tmp = tmp;
   testCase.TestData.manifest_file = fullfile(tmp, 'manifest.json');
end

function teardown(testCase)
   if isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
end

function test_new_site_preserves_existing_case_bytes(testCase)
   % Stage a first site, snapshot its JSON, then merge a second site and assert
   % the first case re-encodes IDENTICALLY (no churn) and the new case is added.

   mf = testCase.TestData.manifest_file;

   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   before = jsondecode(fileread(mf));
   kanl_before = encodeCase(before.cases(1));

   second = familyWith({caseEntry("kanm", "KAN_M", "ablation")});
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, second, ...
      requested_ids="kanm");

   after = jsondecode(fileread(mf));
   ids = arrayfun(@(c) string(c.case_id), after.cases);

   % Both cases present; first stays first (existing order preserved).
   testCase.verifyEqual(numel(after.cases), 2);
   testCase.verifyEqual(ids(:), ["kanl"; "kanm"]);

   % The untouched first case re-encodes byte for byte.
   kanl_after = encodeCase(after.cases(string({after.cases.case_id}) == "kanl"));
   testCase.verifyEqual(kanl_after, kanl_before, ...
      'existing case churned when a new site was staged');

   % Returned manifest matches the file.
   testCase.verifyEqual(numel(merged.cases), 2);
end

function test_new_site_preserves_existing_case_null_coordinates(testCase)
   % MATLAB decodes JSON null as numeric []; rewriting another case must keep
   % nullable site-location scalars as null rather than changing them to [].
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   fields = ["lat_wgs84", "lon_wgs84", "x_epsg3413", ...
      "y_epsg3413", "elev_m"];
   for field = fields
      first.cases.site_location.(char(field)) = NaN;
   end
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("kanm", "KAN_M", "ablation")}), ...
      requested_ids="kanm");

   raw = fileread(mf);
   for field = fields
      name = regexptranslate('escape', char(field));
      testCase.verifyEqual(numel(regexp(raw, ...
         ['"' name '"\s*:\s*null'], 'match')), 1);
      testCase.verifyEmpty(regexp(raw, ...
         ['"' name '"\s*:\s*\[\]'], 'once'));
   end
end

function test_write_repairs_existing_site_location_arrays_to_null(testCase)
   % A prior writer may have changed a nullable scalar from null to []. The
   % next canonical write must repair that raw JSON type despite equal decoding.
   mf = testCase.TestData.manifest_file;
   manifest = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   manifest.cases.site_location.x_epsg3413 = NaN;
   icemodel.verification.setup.writeManifest(mf, manifest);

   canonical = fileread(mf);
   drifted = strrep(canonical, ...
      '"x_epsg3413": null', '"x_epsg3413": []');
   testCase.verifyNotEqual(drifted, canonical);
   fid = fopen(mf, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fwrite(fid, drifted, 'char');
   clear cleaner

   icemodel.verification.setup.writeManifest(mf, manifest);

   repaired = fileread(mf);
   testCase.verifyTrue(contains(repaired, '"x_epsg3413": null'));
   testCase.verifyFalse(contains(repaired, '"x_epsg3413": []'));
end

function test_unknown_family_field_is_preserved(testCase)
   % A hand-added family field (e.g. "schema") on the committed manifest must
   % survive an incremental stage (the importers never emit it).

   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.schema = "metadata_only";
   icemodel.verification.setup.writeManifest(mf, first);

   second = familyWith({caseEntry("kanm", "KAN_M", "ablation")});
   icemodel.verification.setup.writeFamilyManifestMerge(mf, second, ...
      requested_ids="kanm");

   decoded = jsondecode(fileread(mf));
   testCase.verifyTrue(isfield(decoded, 'schema'));
   testCase.verifyEqual(string(decoded.schema), "metadata_only");
end

function test_single_case_manifest_writes_json_arrays(testCase)
   % JSON files must keep array shape even when a family has one case or one
   % skipped record, because external reviewers and schema checks read the file
   % before MATLAB's jsondecode scalarizes one-element arrays.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanm", "KAN_M", "ablation")});
   first.skipped = struct('site', "missing", 'reason', "source absent");

   icemodel.verification.setup.writeManifest(mf, first);

   raw = fileread(mf);
   testCase.verifyTrue(contains(raw, '"cases": ['));
   testCase.verifyTrue(contains(raw, '"skipped": ['));
end

function test_formatManifestTime_includes_midnight_clock(testCase)
   % Manifest period strings should not depend on MATLAB's default datetime
   % display, which omits the clock for midnight values.
   returned = icemodel.verification.setup.formatManifestTime( ...
      datetime(1998, 10, 1, 0, 0, 0, 'TimeZone', 'UTC'));

   testCase.verifyEqual(returned, '1998-10-01 00:00:00');
   testCase.verifyEqual( ...
      icemodel.verification.setup.formatManifestTime(""), '');
end

function test_manifestWindow_uses_canonical_time_serialization(testCase)
   % Every importer and RCM stage must serialize finite and open bounds alike.
   window_start = datetime(2015, 1, 2, 0, 0, 0, TimeZone="UTC");
   open_end = NaT(TimeZone="UTC");

   window = icemodel.verification.setup.manifestWindow( ...
      window_start, open_end);

   testCase.verifyEqual(string(window.start), "2015-01-02 00:00:00");
   testCase.verifyEqual(string(window.end), "");
end

function test_restage_same_site_is_idempotent(testCase)
   % Re-staging the SAME site updates exactly its entry; the OTHER site is
   % untouched and no duplicate case is created.

   mf = testCase.TestData.manifest_file;
   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("kanl", "KAN_L", "ablation")}), ...
      requested_ids="kanl", overwrite_family=true);
   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("kanm", "KAN_M", "ablation")}), ...
      requested_ids="kanm");

   before = jsondecode(fileread(mf));
   kanl_before = encodeCase(before.cases(string({before.cases.case_id}) == "kanl"));

   % Re-stage kanm with a changed note; only kanm updates.
   updated = caseEntry("kanm", "KAN_M", "ablation");
   updated.notes = "restaged";
   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({updated}), requested_ids="kanm");

   after = jsondecode(fileread(mf));
   testCase.verifyEqual(numel(after.cases), 2);

   kanl_after = encodeCase(after.cases(string({after.cases.case_id}) == "kanl"));
   testCase.verifyEqual(kanl_after, kanl_before);

   kanm = after.cases(string({after.cases.case_id}) == "kanm");
   testCase.verifyEqual(string(kanm.notes), "restaged");
end

function test_stale_duplicate_case_ids_are_collapsed(testCase)
   % A stale manifest can have duplicate case_id rows from an older merge bug.
   % The next write must restore one stable entry per staged case folder.

   mf = testCase.TestData.manifest_file;
   old = familyWith({caseEntry("kanu", "KAN_U", "percolation")
      caseEntry("kanu", "KAN_U", "percolation")
      caseEntry("dy2", "DYE-2", "percolation")});
   icemodel.verification.setup.writeManifest(mf, old);

   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("kanm", "KAN_M", "ablation")}), ...
      requested_ids="kanm");

   decoded = jsondecode(fileread(mf));
   ids = string({decoded.cases.case_id});
   testCase.verifyEqual(ids(:), ["kanu"; "dy2"; "kanm"]);
end

function test_new_duplicate_case_ids_are_rejected(testCase)
   % Importers must not introduce duplicate case identities in a staged family.

   mf = testCase.TestData.manifest_file;
   duplicate = familyWith({caseEntry("kanu", "KAN_U", "percolation")
      caseEntry("kanu", "KAN_U", "percolation")});

   testCase.verifyError(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, duplicate, ...
      requested_ids="kanu", overwrite_family=true), ...
      'icemodel:verification:writeFamilyManifestMerge:duplicateCaseId');
end

function test_skipped_for_other_site_is_preserved(testCase)
   % A skip recorded for site A survives staging site B; staging A again clears
   % A's stale skip.

   mf = testCase.TestData.manifest_file;
   m = familyWith({});
   m.skipped = struct('site', "zzz", 'reason', "missing");
   icemodel.verification.setup.writeFamilyManifestMerge(mf, m, ...
      requested_ids="zzz", overwrite_family=true);

   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("kanm", "KAN_M", "ablation")}), ...
      requested_ids="kanm");

   decoded = jsondecode(fileread(mf));
   testCase.verifyEqual(string(decoded.skipped.site), "zzz");

   % Re-stage zzz successfully; its skip clears.
   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("zzz", "ZZZ", "ablation")}), ...
      requested_ids="zzz");
   decoded = jsondecode(fileread(mf));
   testCase.verifyEmpty(decoded.skipped);
end

function test_promice_skip_id_clears_on_normalized_request(testCase)
   % PROMICE skipped.site can store raw station ids. A successful
   % restage keyed by the normalized manifest case id must clear that stale skip.
   mf = testCase.TestData.manifest_file;
   m = familyWith({});
   m.skipped = struct('site', "KAN_M", 'reason', "missing");
   icemodel.verification.setup.writeFamilyManifestMerge(mf, m, ...
      requested_ids="kanm", overwrite_family=true);

   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("kanm", "KAN_M", "ablation")}), ...
      requested_ids="kanm");

   decoded = jsondecode(fileread(mf));
   testCase.verifyEmpty(decoded.skipped);
end

function test_distinct_skip_id_is_preserved(testCase)
   % Compact matching should only clear records equivalent to the requested id.
   mf = testCase.TestData.manifest_file;
   m = familyWith({});
   m.skipped = struct('site', "KAN_U", 'reason', "missing");
   icemodel.verification.setup.writeFamilyManifestMerge(mf, m, ...
      requested_ids="kanu", overwrite_family=true);

   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("kanm", "KAN_M", "ablation")}), ...
      requested_ids="kanm");

   decoded = jsondecode(fileread(mf));
   testCase.verifyEqual(string(decoded.skipped.site), "KAN_U");
end

function test_fresh_skip_preserves_existing_case_entry(testCase)
   % A transient requested refresh failure must not delete committed case state.
   mf = testCase.TestData.manifest_file;
   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("kanm", "KAN_M", "ablation")}), ...
      requested_ids="kanm", overwrite_family=true);

   missing = familyWith({});
   missing.skipped = struct('site', "kanm", 'reason', "missing source");
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      missing, requested_ids="kanm");

   decoded = jsondecode(fileread(mf));
   testCase.verifyEqual(string(decoded.cases.case_id), "kanm");
   testCase.verifyEqual(string(decoded.skipped.site), "kanm");
   testCase.verifyEqual(string(merged.cases.case_id), "kanm");
end

function test_same_case_patch_preserves_omitted_colocation_leg(testCase)
   % Source selectors patch one case without declaring its full colocation graph.
   mf = testCase.TestData.manifest_file;
   existing = caseEntry("kanm", "KAN_M", "ablation");
   existing.colocation.mar = struct('kind', 'point_met_and_data', ...
      'staged', true, 'source_id', 'mar3.11', ...
      'met_files', "mar3.11/met_kanm_mar3.11_2012_15m.mat", ...
      'data_files', "mar3.11/kanm_mar3.11_2012.mat");
   existing.forcing_sources = {'promice'; 'mar3.11'};
   existing.eval_sources = {'promice_obs'; 'mar3.11'};
   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({existing}), requested_ids="kanm", overwrite_family=true);

   incoming = caseEntry("kanm", "KAN_M", "ablation");
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({incoming}), requested_ids="kanm");

   testCase.verifyTrue(isfield(merged.cases.colocation, 'mar'));
   testCase.verifyEqual(jsonencode(merged.cases.colocation.mar), ...
      jsonencode(existing.colocation.mar));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(merged.cases.forcing_sources)));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(merged.cases.eval_sources)));
end

function test_additive_shorter_refresh_preserves_enclosing_windows(testCase)
   % A surgical shorter call patches metadata without contracting durable scope.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2010 00:00:00', '31-Dec-2021 00:00:00');
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   shorter = first;
   shorter.cases.period = struct('start', '01-Jan-2015 00:00:00', ...
      'end', '31-Dec-2018 00:00:00');
   shorter.cases.colocation.mar.window = struct( ...
      'start', '01-Jan-2016 00:00:00', 'end', '31-Dec-2017 00:00:00');
   shorter.cases.colocation.mar.met_files = {'met-short.mat'};
   shorter.cases.colocation.mar.data_files = {'data-short.mat'};
   shorter.cases.colocation.mar.staged = false;
   shorter.cases.forcing_sources = {'promice'};
   shorter.cases.eval_sources = {'promice_obs'};
   shorter.cases.notes = 'surgical metadata refresh';
   testCase.verifyWarningFree(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, shorter, ...
      requested_ids="kanl"));
   merged = jsondecode(fileread(mf));

   testCase.verifyEqual(jsonencode(merged.cases.period), ...
      jsonencode(first.cases.period));
   testCase.verifyEqual(jsonencode(merged.cases.colocation.mar.window), ...
      jsonencode(first.cases.colocation.mar.window));
   testCase.verifyEqual(string(merged.cases.colocation.mar.met_files), ...
      "met-a.mat");
   testCase.verifyEqual(string(merged.cases.colocation.mar.data_files), ...
      "data-a.mat");
   testCase.verifyTrue(merged.cases.colocation.mar.staged);
   testCase.verifyTrue(ismember("mar3.11", ...
      string(merged.cases.forcing_sources)));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(merged.cases.eval_sources)));
   testCase.verifyEqual(string(merged.cases.notes), ...
      "surgical metadata refresh");
end

function test_additive_wider_refresh_expands_windows(testCase)
   % A successful wider call expands both case and staged-leg coverage.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.period = struct('start', '01-Jan-2012 00:00:00', ...
      'end', '31-Dec-2018 00:00:00');
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2013 00:00:00', '31-Dec-2017 00:00:00');
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   wider = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   wider.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2010 00:00:00', '31-Dec-2021 00:00:00');
   wider.cases.colocation.mar.met_files = {'met-wide.mat'};
   wider.cases.colocation.mar.data_files = {'data-wide.mat'};
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, wider, ...
      requested_ids="kanl");

   testCase.verifyEqual(jsonencode(merged.cases.period), ...
      jsonencode(wider.cases.period));
   testCase.verifyEqual(jsonencode(merged.cases.colocation.mar.window), ...
      jsonencode(wider.cases.colocation.mar.window));
   testCase.verifyEqual(string(merged.cases.colocation.mar.met_files), ...
      "met-wide.mat");
   testCase.verifyEqual(string(merged.cases.colocation.mar.data_files), ...
      "data-wide.mat");
end

function test_wider_forcing_refresh_preserves_existing_model_output(testCase)
   % Widening met/Data without rebuilding a profile keeps its artifact group.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("fa", "firn aquifer (fa)", "accumulation")});
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2012 00:00:00', '31-Dec-2015 00:00:00', "met-old.mat");
   first.cases.colocation.mar.model_output_files = ...
      {'input/userdata/mar3.11/firn aquifer (fa)_mar3.11_density_profiles.mat'};
   first.cases.colocation.mar.model_output_format = 'density_profile_timeseries';
   first.cases.colocation.mar.model_output_variables = {'density'};
   first.cases.colocation.mar.model_output_status = 'staged';
   first.cases.colocation.mar.model_output_note = 'exact-date profiles';
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="fa", overwrite_family=true);

   wider = first;
   wider.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2012 00:00:00', '31-Dec-2018 00:00:00', "met-wide.mat");
   wider.cases.colocation.mar.model_output_status = 'profile_not_available';
   wider.cases.colocation.mar.model_output_note = ...
      'Observation source was unavailable during forcing-only refresh.';
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, wider, ...
      requested_ids="fa");
   leg = merged.cases.colocation.mar;

   testCase.verifyEqual(string(leg.met_files), "met-wide.mat");
   testCase.verifyEqual(string(leg.data_files), "data-wide.mat");
   testCase.verifyEqual(string(leg.model_output_files), ...
      "input/userdata/mar3.11/firn aquifer (fa)_mar3.11_density_profiles.mat");
   testCase.verifyEqual(string(leg.model_output_format), ...
      "density_profile_timeseries");
   testCase.verifyEqual(string(leg.model_output_variables), "density");
   testCase.verifyEqual(string(leg.model_output_status), "staged");
   testCase.verifyEqual(string(leg.model_output_note), "exact-date profiles");
end

function test_new_model_output_replaces_existing_artifact_group(testCase)
   % A rebuilt sidecar replaces, rather than unions with, stale profile state.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("fa", "firn aquifer (fa)", "accumulation")});
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2012 00:00:00', '31-Dec-2015 00:00:00', "met-old.mat");
   first.cases.colocation.mar.model_output_files = {'profile-old.mat'};
   first.cases.colocation.mar.model_output_format = 'legacy_profile';
   first.cases.colocation.mar.model_output_variables = {'old_density'};
   first.cases.colocation.mar.model_output_status = 'staged';
   first.cases.colocation.mar.model_output_note = 'old profile';
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="fa", overwrite_family=true);

   rebuilt = first;
   rebuilt.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2012 00:00:00', '31-Dec-2018 00:00:00', "met-wide.mat");
   rebuilt.cases.colocation.mar.model_output_files = {'profile-new.mat'};
   rebuilt.cases.colocation.mar.model_output_format = ...
      'density_profile_timeseries';
   rebuilt.cases.colocation.mar.model_output_variables = {'density'};
   rebuilt.cases.colocation.mar.model_output_status = 'staged';
   rebuilt.cases.colocation.mar.model_output_note = 'new exact-date profile';
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, rebuilt, ...
      requested_ids="fa");
   leg = merged.cases.colocation.mar;

   testCase.verifyEqual(string(leg.model_output_files), "profile-new.mat");
   testCase.verifyEqual(string(leg.model_output_format), ...
      "density_profile_timeseries");
   testCase.verifyEqual(string(leg.model_output_variables), "density");
   testCase.verifyEqual(string(leg.model_output_status), "staged");
   testCase.verifyEqual(string(leg.model_output_note), ...
      "new exact-date profile");
end

function test_additive_partial_extension_unions_windows_and_references(testCase)
   % Partial coverage extensions retain both file sets supporting the union.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2010 00:00:00', '31-Dec-2015 00:00:00', "met-old.mat");
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   extension = first;
   extension.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2014 00:00:00', '31-Dec-2021 00:00:00', "met-new.mat");
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      extension, requested_ids="kanl");

   testCase.verifyEqual(string(merged.cases.colocation.mar.window.start), ...
      "01-Jan-2010 00:00:00");
   testCase.verifyEqual(string(merged.cases.colocation.mar.window.end), ...
      "31-Dec-2021 00:00:00");
   testCase.verifyEqual(string(merged.cases.colocation.mar.met_files), ...
      ["met-old.mat"; "met-new.mat"]);
   testCase.verifyEqual(string(merged.cases.colocation.mar.data_files), ...
      ["data-old.mat"; "data-new.mat"]);
   testCase.verifyTrue(merged.cases.colocation.mar.staged);
end

function test_additive_disjoint_refresh_does_not_invent_gap_coverage(testCase)
   % Separate artifacts cannot justify one continuous scalar window across a gap.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.period = struct('start', '01-Jan-2010 00:00:00', ...
      'end', '31-Dec-2012 00:00:00');
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2010 00:00:00', '31-Dec-2012 00:00:00', "met-old.mat");
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   disjoint = first;
   disjoint.cases.period = struct('start', '01-Jan-2014 00:00:00', ...
      'end', '31-Dec-2016 00:00:00');
   disjoint.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2014 00:00:00', '31-Dec-2016 00:00:00', "met-new.mat");
   testCase.verifyWarning(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, disjoint, ...
      requested_ids="kanl"), ...
      'icemodel:verification:writeFamilyManifestMerge:disjointCoverage');
   merged = jsondecode(fileread(mf));

   testCase.verifyEqual(jsonencode(merged.cases.period), ...
      jsonencode(first.cases.period));
   testCase.verifyEqual(jsonencode(merged.cases.colocation.mar.window), ...
      jsonencode(first.cases.colocation.mar.window));
   testCase.verifyEqual(string(merged.cases.colocation.mar.met_files), ...
      "met-old.mat");
end

function test_additive_bounded_refresh_preserves_unbounded_prior_state(testCase)
   % Blank all-available coverage is durable and cannot become a bounded subset.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.period = struct('start', '', 'end', '');
   first.cases.colocation.mar = windowedArtifactLeg('', '', "met-all.mat");
   first.cases.colocation.mar = rmfield( ...
      first.cases.colocation.mar, 'staged');
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   bounded = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   bounded.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2016 00:00:00', '31-Dec-2017 00:00:00', "met-bounded.mat");
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, bounded, ...
      requested_ids="kanl");

   testCase.verifyEqual(string(merged.cases.period.start), "");
   testCase.verifyEqual(string(merged.cases.period.end), "");
   testCase.verifyEqual(string(merged.cases.colocation.mar.window.start), "");
   testCase.verifyEqual(string(merged.cases.colocation.mar.window.end), "");
   testCase.verifyEqual(string(merged.cases.colocation.mar.met_files), ...
      "met-all.mat");
   testCase.verifyFalse(isfield(merged.cases.colocation.mar, 'staged'));
end

function test_additive_identity_conflict_replaces_instead_of_union(testCase)
   % A different concrete sampling method denotes a replacement artifact leg.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2010 00:00:00', '31-Dec-2015 00:00:00', "met-nearest.mat");
   first.cases.colocation.mar.method = 'nearest';
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   replacement = first;
   replacement.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2014 00:00:00', '31-Dec-2021 00:00:00', "met-natural.mat");
   replacement.cases.colocation.mar.sample_method = 'natural';
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      replacement, requested_ids="kanl");

   testCase.verifyEqual(jsonencode(merged.cases.colocation.mar.window), ...
      jsonencode(replacement.cases.colocation.mar.window));
   testCase.verifyEqual(string(merged.cases.colocation.mar.met_files), ...
      "met-natural.mat");
end

function test_requested_invalidation_consumes_replacement_signal(testCase)
   % A validated requested-leg invalidation replaces prior artifact state, and
   % its orchestration-only signal must not enter the durable manifest.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2010 00:00:00', '31-Dec-2015 00:00:00', "met-old.mat");
   first.cases.colocation.mar.model_output_files = {'profile-old.mat'};
   first.cases.colocation.mar.model_output_format = 'legacy_profile';
   first.cases.colocation.mar.model_output_variables = {'old_density'};
   first.cases.colocation.mar.model_output_status = 'staged';
   first.cases.colocation.mar.model_output_note = 'old profile';
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   invalidated = first;
   invalidated.cases.forcing_sources = {'promice'};
   invalidated.cases.colocation.mar = struct( ...
      'kind', first.cases.colocation.mar.kind, ...
      'staged', false, 'source_id', 'mar3.11', ...
      'window', first.cases.colocation.mar.window, ...
      'reason', 'prior staged artifact file is missing', ...
      'replace_prior_artifacts', true);
   merged = icemodel.verification.setup.writeFamilyManifestMerge( ...
      mf, invalidated, requested_ids="kanl");

   testCase.verifyFalse(logical(merged.cases.colocation.mar.staged));
   testCase.verifyFalse(isfield(merged.cases.colocation.mar, 'met_files'));
   testCase.verifyFalse(isfield(merged.cases.colocation.mar, 'data_files'));
   testCase.verifyFalse(isfield( ...
      merged.cases.colocation.mar, 'replace_prior_artifacts'));
   testCase.verifyFalse(any(startsWith(string(fieldnames( ...
      merged.cases.colocation.mar)), "model_output_")));
end

function test_additive_method_alias_equivalence_unions_matching_cadence(testCase)
   % Production method and manifest sample_method are equivalent aliases. A
   % valid mixed 15-minute-met/hourly-Data leg can extend and union references.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2010 00:00:00', '31-Dec-2015 00:00:00', ...
      "met_kanl_mar3.11_20100101_20151231_15m.mat");
   first.cases.colocation.mar.data_files = ...
      {'kanl_mar3.11_20100101_20151231.mat'};
   first.cases.colocation.mar.method = 'nearest';
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   extension = first;
   extension.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2014 00:00:00', '31-Dec-2021 00:00:00', ...
      "met_kanl_mar3.11_20140101_20211231_15m.mat");
   extension.cases.colocation.mar.data_files = ...
      {'kanl_mar3.11_20140101_20211231.mat'};
   extension.cases.colocation.mar.sample_method = 'nearest';
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      extension, requested_ids="kanl");

   testCase.verifyEqual(string(merged.cases.colocation.mar.window.start), ...
      "01-Jan-2010 00:00:00");
   testCase.verifyEqual(string(merged.cases.colocation.mar.window.end), ...
      "31-Dec-2021 00:00:00");
   testCase.verifyEqual(string(merged.cases.colocation.mar.met_files), [ ...
      "met_kanl_mar3.11_20100101_20151231_15m.mat"; ...
      "met_kanl_mar3.11_20140101_20211231_15m.mat"]);
   testCase.verifyEqual(string(merged.cases.colocation.mar.data_files), [ ...
      "kanl_mar3.11_20100101_20151231.mat"; ...
      "kanl_mar3.11_20140101_20211231.mat"]);
end

function test_additive_met_cadence_conflict_replaces_instead_of_union(testCase)
   % Known 15-minute versus 30-minute met cadence is a replacement boundary,
   % even when the paired hourly Data cadence and scalar identity agree.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2010 00:00:00', '31-Dec-2015 00:00:00', ...
      "met_kanl_mar3.11_20100101_20151231_15m.mat");
   first.cases.colocation.mar.data_files = ...
      {'kanl_mar3.11_20100101_20151231.mat'};
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   replacement = first;
   replacement.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2014 00:00:00', '31-Dec-2021 00:00:00', ...
      "met_kanl_mar3.11_20140101_20211231_30m.mat");
   replacement.cases.colocation.mar.data_files = ...
      {'kanl_mar3.11_20140101_20211231.mat'};
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      replacement, requested_ids="kanl");

   testCase.verifyEqual(jsonencode(merged.cases.colocation.mar.window), ...
      jsonencode(replacement.cases.colocation.mar.window));
   testCase.verifyEqual(string(merged.cases.colocation.mar.met_files), ...
      "met_kanl_mar3.11_20140101_20211231_30m.mat");
   testCase.verifyEqual(string(merged.cases.colocation.mar.data_files), ...
      "kanl_mar3.11_20140101_20211231.mat");
end

function test_additive_userdata_cadence_conflict_replaces_instead_of_union(testCase)
   % Suffix-free hourly userdata cannot union with a 30-minute Data variant,
   % while matching 15-minute model-met references do not hide that conflict.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2010 00:00:00', '31-Dec-2015 00:00:00', ...
      "met_kanl_mar3.11_20100101_20151231_15m.mat");
   first.cases.colocation.mar.data_files = ...
      {'kanl_mar3.11_20100101_20151231.mat'};
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   replacement = first;
   replacement.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2014 00:00:00', '31-Dec-2021 00:00:00', ...
      "met_kanl_mar3.11_20140101_20211231_15m.mat");
   replacement.cases.colocation.mar.data_files = ...
      {'kanl_mar3.11_20140101_20211231_30m.mat'};
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      replacement, requested_ids="kanl");

   testCase.verifyEqual(jsonencode(merged.cases.colocation.mar.window), ...
      jsonencode(replacement.cases.colocation.mar.window));
   testCase.verifyEqual(string(merged.cases.colocation.mar.met_files), ...
      "met_kanl_mar3.11_20140101_20211231_15m.mat");
   testCase.verifyEqual(string(merged.cases.colocation.mar.data_files), ...
      "kanl_mar3.11_20140101_20211231_30m.mat");
end

function test_additive_explicit_single_class_cadence_conflict_replaces(testCase)
   % A generic explicit cadence is authoritative for a met-only leg even when
   % legacy reference names cannot encode cadence themselves.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2010 00:00:00', '31-Dec-2015 00:00:00', "met-old.mat");
   first.cases.colocation.mar = rmfield( ...
      first.cases.colocation.mar, 'data_files');
   first.cases.colocation.mar.artifact_cadence_seconds = 900;
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   replacement = first;
   replacement.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2014 00:00:00', '31-Dec-2021 00:00:00', "met-new.mat");
   replacement.cases.colocation.mar = rmfield( ...
      replacement.cases.colocation.mar, 'data_files');
   replacement.cases.colocation.mar.artifact_metadata = ...
      struct('artifact_cadence_seconds', 1800);
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      replacement, requested_ids="kanl");

   testCase.verifyEqual(jsonencode(merged.cases.colocation.mar.window), ...
      jsonencode(replacement.cases.colocation.mar.window));
   testCase.verifyEqual(string(merged.cases.colocation.mar.met_files), ...
      "met-new.mat");
end

function test_additive_point_conflict_replaces_instead_of_union(testCase)
   % A relocated case must not retain artifacts sampled at its prior point.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2010 00:00:00', '31-Dec-2015 00:00:00', "met-old-point.mat");
   first.cases.colocation.merra = windowedArtifactLeg( ...
      '01-Jan-2011 00:00:00', '31-Dec-2014 00:00:00', "met-merra.mat");
   first.cases.colocation.merra.source_id = 'merra2';
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   replacement = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   replacement.cases.site_location.lat_wgs84 = 68;
   replacement.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2014 00:00:00', '31-Dec-2021 00:00:00', "met-new-point.mat");
   merged = icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      replacement, requested_ids="kanl");

   testCase.verifyEqual(jsonencode(merged.cases.colocation.mar.window), ...
      jsonencode(replacement.cases.colocation.mar.window));
   testCase.verifyEqual(string(merged.cases.colocation.mar.met_files), ...
      "met-new-point.mat");
   testCase.verifyEqual(jsonencode(merged.cases.colocation.merra), ...
      jsonencode(first.cases.colocation.merra));
   testCase.verifyTrue(ismember("merra2", ...
      string(merged.cases.forcing_sources)));
end

function test_overwrite_family_explicitly_narrows_windows_with_warning(testCase)
   % The destructive family boundary may contract coverage, but must warn.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = windowedArtifactLeg( ...
      '01-Jan-2010 00:00:00', '31-Dec-2021 00:00:00');
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   narrower = first;
   narrower.cases.period = struct('start', '01-Jan-2015 00:00:00', ...
      'end', '31-Dec-2018 00:00:00');
   narrower.cases.colocation.mar.window = struct( ...
      'start', '01-Jan-2016 00:00:00', 'end', '31-Dec-2017 00:00:00');
   testCase.verifyWarning(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, narrower, ...
      requested_ids="kanl", overwrite_family=true), ...
      'icemodel:verification:writeFamilyManifestMerge:overwriteFamily');
   replaced = jsondecode(fileread(mf));

   testCase.verifyEqual(jsonencode(replaced.cases.period), ...
      jsonencode(narrower.cases.period));
   testCase.verifyEqual(jsonencode(replaced.cases.colocation.mar.window), ...
      jsonencode(narrower.cases.colocation.mar.window));
end

function test_overwrite_family_warns_when_unbounded_period_becomes_bounded(testCase)
   % All-available coverage is wider than a bounded replacement and must warn.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.period = struct('start', '', 'end', '');
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   bounded = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   testCase.verifyWarning(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, bounded, ...
      requested_ids="kanl", overwrite_family=true), ...
      'icemodel:verification:writeFamilyManifestMerge:overwriteFamily');
end

function test_overwrite_family_drops_prior_cases(testCase)
   % overwrite_family=true forces a full rewrite from the requested cases.

   mf = testCase.TestData.manifest_file;
   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("kanl", "KAN_L", "ablation")}), ...
      requested_ids="kanl", overwrite_family=true);
   testCase.verifyWarning(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
         familyWith({caseEntry("kanm", "KAN_M", "ablation")}), ...
         requested_ids="kanm", overwrite_family=true), ...
      'icemodel:verification:writeFamilyManifestMerge:overwriteFamily');

   decoded = jsondecode(fileread(mf));
   testCase.verifyEqual(numel(decoded.cases), 1);
   testCase.verifyEqual(string(decoded.cases.case_id), "kanm");
end

function test_overwrite_family_warns_when_source_removed(testCase)
   % Replacing a retained case still warns when a prior source disappears.
   mf = testCase.TestData.manifest_file;
   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("kanl", "KAN_L", "ablation")}), ...
      requested_ids="kanl", overwrite_family=true);
   replacement = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   replacement.cases.forcing_sources = {'promice'};

   testCase.verifyWarning(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, replacement, ...
         requested_ids="kanl", overwrite_family=true), ...
      'icemodel:verification:writeFamilyManifestMerge:overwriteFamily');
end

function test_overwrite_family_source_addition_is_warning_free(testCase)
   % Kill-safe replacement may add source state without reporting a removal.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.forcing_sources = {'promice'};
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   testCase.verifyWarningFree(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
         familyWith({caseEntry("kanl", "KAN_L", "ablation")}), ...
         requested_ids="kanl", overwrite_family=true));
end

function test_overwrite_family_warns_when_skipped_record_removed(testCase)
   % A skipped record is durable family state, so dropping it during a full
   % replacement must surface the same destructive-overwrite warning.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.skipped = struct('site', "zzz", 'reason', "missing source");
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids=["kanl", "zzz"], overwrite_family=true);

   testCase.verifyWarning(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
         familyWith({caseEntry("kanl", "KAN_L", "ablation")}), ...
         requested_ids="kanl", overwrite_family=true), ...
      'icemodel:verification:writeFamilyManifestMerge:overwriteFamily');

   decoded = jsondecode(fileread(mf));
   testCase.verifyEmpty(decoded.skipped);
end

function test_overwrite_family_warns_when_retained_case_period_narrows(testCase)
   % A same-id replacement that clips either prior endpoint removes coverage.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   replacement = first;
   replacement.cases.period.start = '01-Jan-2010 00:00:00';
   testCase.verifyWarning(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, replacement, ...
      requested_ids="kanl", overwrite_family=true), ...
      'icemodel:verification:writeFamilyManifestMerge:overwriteFamily');
end

function test_overwrite_family_warns_when_case_reference_file_is_lost(testCase)
   % Retaining a case id does not make a dropped reference artifact safe.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.reference_file = 'kanl/reference.mat';
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   replacement = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   testCase.verifyWarning(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, replacement, ...
      requested_ids="kanl", overwrite_family=true), ...
      'icemodel:verification:writeFamilyManifestMerge:overwriteFamily');
end

function test_overwrite_family_warns_when_evaluation_file_is_cleared(testCase)
   % Clearing the retained case's primary evaluation ref is destructive state.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   replacement = first;
   replacement.cases.evaluation_file = '';
   testCase.verifyWarning(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, replacement, ...
      requested_ids="kanl", overwrite_family=true), ...
      'icemodel:verification:writeFamilyManifestMerge:overwriteFamily');
end

function test_overwrite_family_warns_when_retained_leg_loses_staged_state(testCase)
   % A retained leg name still loses durable state when staged flips true to false.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = artifactLeg(true, ["met-a.mat", "met-b.mat"]);
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   replacement = first;
   replacement.cases.colocation.mar.staged = false;
   testCase.verifyWarning(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, replacement, ...
      requested_ids="kanl", overwrite_family=true), ...
      'icemodel:verification:writeFamilyManifestMerge:overwriteFamily');
end

function test_overwrite_family_warns_when_retained_leg_loses_files(testCase)
   % Fewer met or Data refs are destructive even when the leg remains staged.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = artifactLeg(true, ["met-a.mat", "met-b.mat"]);
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   replacement = first;
   replacement.cases.colocation.mar.met_files = {'met-a.mat'};
   replacement.cases.colocation.mar.data_files = {'data-a.mat'};
   testCase.verifyWarning(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, replacement, ...
      requested_ids="kanl", overwrite_family=true), ...
      'icemodel:verification:writeFamilyManifestMerge:overwriteFamily');
end

function test_overwrite_family_coverage_and_artifact_addition_is_warning_free(testCase)
   % A broader period with only added sources/files is an additive replacement.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   first.cases.colocation.mar = artifactLeg(true, "met-a.mat");
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   replacement = first;
   replacement.cases.period.start = '01-Jan-2008 00:00:00';
   replacement.cases.period.end = '31-Dec-2023 00:00:00';
   replacement.cases.forcing_sources{end + 1} = 'merra2';
   replacement.cases.colocation.mar.met_files = {'met-a.mat'; 'met-b.mat'};
   replacement.cases.colocation.mar.data_files = {'data-a.mat'; 'data-b.mat'};
   testCase.verifyWarningFree(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, replacement, ...
      requested_ids="kanl", overwrite_family=true));
end

function test_identical_additive_merge_does_not_rewrite_manifest(testCase)
   % Repeating an ordinary additive merge with an identical payload must leave
   % both manifest bytes and filesystem modification time unchanged.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);
   before_text = fileread(mf);
   before_info = dir(mf);

   % Advance beyond coarse one-second filesystem timestamp resolution so an
   % accidental rewrite cannot alias the original modification time.
   pause(1.1)
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl");

   after_info = dir(mf);
   testCase.verifyEqual(fileread(mf), before_text);
   testCase.verifyEqual(after_info.datenum, before_info.datenum);
end

function test_field_order_only_merge_does_not_rewrite_manifest(testCase)
   % Recursive JSON-object field order is not a semantic manifest change.
   mf = testCase.TestData.manifest_file;
   canonical = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   reordered = orderfields(canonical, flipud(fieldnames(canonical)));
   reordered.cases = orderfields( ...
      reordered.cases, flipud(fieldnames(reordered.cases)));
   reordered.cases.period = orderfields( ...
      reordered.cases.period, flipud(fieldnames(reordered.cases.period)));
   icemodel.verification.setup.writeManifest(mf, reordered);
   before_text = fileread(mf);
   before_info = dir(mf);

   % Advance past coarse timestamp resolution so a field-order rewrite cannot
   % alias the original mtime.
   pause(1.1)
   icemodel.verification.setup.writeFamilyManifestMerge(mf, canonical, ...
      requested_ids="kanl");

   after_info = dir(mf);
   testCase.verifyEqual(fileread(mf), before_text);
   testCase.verifyEqual(after_info.datenum, before_info.datenum);
end

function test_transient_reuse_note_does_not_rewrite_manifest(testCase)
   % Exact full-coverage replay notes are runtime diagnostics, not provenance.
   % A real source note remains semantic and must still rewrite the manifest.
   mf = testCase.TestData.manifest_file;
   first_case = caseEntry("kanl", "KAN_L", "ablation");
   second_case = caseEntry("kanm", "KAN_M", "ablation");
   first_leg = windowedArtifactLeg( ...
      '01-Jan-2012 00:00:00', '31-Dec-2012 23:00:00');
   first_leg.source = 'mar';
   first_leg.sample_method = 'nearest';
   first_leg.model_output_files = {'mar3.11/kanl/profile.mat'};
   second_leg = windowedArtifactLeg( ...
      '01-Jan-2012 00:00:00', '31-Dec-2012 23:00:00');
   second_leg.source = 'mar';
   second_leg.sample_method = 'nearest';
   first_case.colocation.mar = first_leg;
   second_case.colocation.mar = second_leg;
   initial = familyWith({first_case, second_case});
   icemodel.verification.setup.writeManifest(mf, initial);
   before_text = fileread(mf);
   before_info = dir(mf);

   % Both complete-reuse messages are nonpersistent even after filesystem
   % timestamp resolution has advanced.
   pause(1.1)
   transient = [ ...
      "Existing staged met/Data files fully cover requested window."
      "Existing staged Data file fully covers requested window."];
   for note = reshape(transient, 1, [])
      replay = initial;
      replay.cases(1).colocation.mar.note = char(note);
      icemodel.verification.setup.writeManifest(mf, replay);
      testCase.verifyEqual(fileread(mf), before_text);
      testCase.verifyEqual(dir(mf).datenum, before_info.datenum);
   end

   % A scientific or manually curated note is durable manifest content.
   pause(1.1)
   changed = initial;
   changed.cases(1).colocation.mar.note = ...
      'Source rows were manually reviewed.';
   icemodel.verification.setup.writeManifest(mf, changed);
   decoded = jsondecode(fileread(mf));
   testCase.verifyNotEqual(fileread(mf), before_text);
   testCase.verifyGreaterThan(dir(mf).datenum, before_info.datenum);
   testCase.verifyEqual(string(decoded.cases(1).colocation.mar.note), ...
      "Source rows were manually reviewed.");
end

function test_semantic_compare_writes_real_case_order_change(testCase)
   % JSON array order remains semantic even when every case value is unchanged.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation"), ...
      caseEntry("kanm", "KAN_M", "ablation")});
   icemodel.verification.setup.writeManifest(mf, first);
   before_text = fileread(mf);
   before_info = dir(mf);

   % Reverse the case array through the low-level writer because additive merge
   % intentionally preserves established case order.
   reordered = first;
   reordered.cases = reordered.cases([2, 1]);
   pause(1.1)
   icemodel.verification.setup.writeManifest(mf, reordered);

   after_info = dir(mf);
   decoded = jsondecode(fileread(mf));
   testCase.verifyNotEqual(fileread(mf), before_text);
   testCase.verifyGreaterThan(after_info.datenum, before_info.datenum);
   testCase.verifyEqual(string({decoded.cases.case_id}), ["kanm", "kanl"]);
end

function test_semantic_compare_replaces_malformed_existing_json(testCase)
   % Invalid JSON cannot be treated as a semantic no-op by the shared writer.
   mf = testCase.TestData.manifest_file;
   fid = fopen(mf, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fwrite(fid, '{malformed', 'char');
   clear cleaner

   manifest = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   icemodel.verification.setup.writeManifest(mf, manifest);

   decoded = jsondecode(fileread(mf));
   testCase.verifyEqual(string(decoded.cases.case_id), "kanl");
end

function test_case_field_mismatch_preserves_error_identifier(testCase)
   % Schema drift between preserved and touched case entries must stay catchable.
   mf = testCase.TestData.manifest_file;
   first = familyWith({caseEntry("kanl", "KAN_L", "ablation")});
   icemodel.verification.setup.writeFamilyManifestMerge(mf, first, ...
      requested_ids="kanl", overwrite_family=true);

   second = familyWith({caseEntry("kanm", "KAN_M", "ablation")});
   second.cases.extra_field = "unexpected";

   testCase.verifyError(@() ...
      icemodel.verification.setup.writeFamilyManifestMerge(mf, second, ...
      requested_ids="kanm"), ...
      'icemodel:verification:writeFamilyManifestMerge:fieldMismatch');
end

%% Local helpers
function m = familyWith(case_cell)
   %FAMILYWITH Build a family manifest holding the given case entries.
   if isempty(case_cell)
      cases = struct([]);
   else
      cases = vertcat(case_cell{:});
   end
   m = icemodel.verification.setup.makeFamilyManifest("promice", "", ...
      "https://promice.org", "colocated[test]", "01-Jan-2026", cases, ...
      struct('site', {}, 'reason', {}));
end

function entry = caseEntry(case_id, site_id, zone)
   %CASEENTRY Minimal schema-valid firn case entry for the merge tests.
   values = { ...
      char(case_id)
      'firn_observational'
      char(site_id)
      char(replace(site_id, "_", "-"))
      char(zone)
      cellstr(["seasonal_snow"; "bare_ice"])
      'none'
      struct('lat_wgs84', 67, 'lon_wgs84', -50, ...
      'x_epsg3413', -2e5, 'y_epsg3413', -2.5e6, 'elev_m', 600)
      struct('start', '01-Jan-2009 00:00:00', ...
      'end', '31-Dec-2022 00:00:00')
      char(fullfile(char(case_id), 'observations.mat'))
      cellstr(["promice"; "mar3.11"])
      cellstr(["promice_obs"; "racmo2.3p3"])
      cellstr(["tice1"; "tice2"])
      struct('subsurface_temperature', 'tice1..tice8')
      struct('promice', struct('kind', 'station_met_and_eval'))
      'daily'
      'merge test fixture'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
end

function leg = artifactLeg(staged, met_files)
   %ARTIFACTLEG Build one retained source leg with paired met/Data references.
   met_files = cellstr(reshape(string(met_files), [], 1));
   data_files = replace(string(met_files), "met-", "data-");
   leg = struct('kind', 'point_met_and_data', 'staged', staged, ...
      'source_id', 'mar3.11', 'met_files', {met_files}, ...
      'data_files', {cellstr(data_files)});
end

function leg = windowedArtifactLeg(startdate, enddate, met_file)
   %WINDOWEDARTIFACTLEG Build a retained source leg with bounded coverage.
   if nargin < 3
      met_file = "met-a.mat";
   end
   leg = artifactLeg(true, met_file);
   leg.window = struct('start', startdate, 'end', enddate);
end

function s = encodeCase(case_struct)
   %ENCODECASE Stable JSON encoding of one case for byte comparison.
   s = jsonencode(case_struct, PrettyPrint=true);
end
