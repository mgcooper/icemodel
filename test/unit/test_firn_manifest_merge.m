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

function test_overwrite_family_drops_prior_cases(testCase)
   % overwrite_family=true forces a full rewrite from the requested cases.

   mf = testCase.TestData.manifest_file;
   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("kanl", "KAN_L", "ablation")}), ...
      requested_ids="kanl", overwrite_family=true);
   icemodel.verification.setup.writeFamilyManifestMerge(mf, ...
      familyWith({caseEntry("kanm", "KAN_M", "ablation")}), ...
      requested_ids="kanm", overwrite_family=true);

   decoded = jsondecode(fileread(mf));
   testCase.verifyEqual(numel(decoded.cases), 1);
   testCase.verifyEqual(string(decoded.cases.case_id), "kanm");
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
      "https://promice.org", "colocated[test]", "01-Jan-2026", cases);
   m.skipped = struct('site', {}, 'reason', {});
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
      struct('start', '01-Jan-2009', 'end', '31-Dec-2022')
      char(fullfile(char(case_id), 'observations.mat'))
      cellstr(["promice"; "mar"])
      cellstr(["promice_obs"; "racmo"])
      cellstr(["tice1"; "tice2"])
      struct('subsurface_temperature', 'tice1..tice8')
      struct('promice', struct('kind', 'station_met_and_eval'))
      'daily'
      'merge test fixture'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
end

function s = encodeCase(case_struct)
   %ENCODECASE Stable JSON encoding of one case for byte comparison.
   s = jsonencode(case_struct, PrettyPrint=true);
end
