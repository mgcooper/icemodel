function tests = test_firn_verification_contracts
   %TEST_FIRN_VERIFICATION_CONTRACTS Verify the firn-evaluation contracts.
   %
   % Exercises the firn_observational lane end to end against the committed
   % co-located PROMICE firn cases under demo/data/eval/promice/:
   %   - listcases enumerates the firn family alongside the snow families;
   %   - each firn manifest's case_type validates against the namelist;
   %   - the candidate adapter resolves the declared comparison variables;
   %   - comparecase reports SOFT (diagnostic) metrics with no hard PASS/FAIL;
   %   - the co-located multi-model bundle (promice/mar/merra met + racmo Data)
   %     the manifest declares resolves on disk.
   %
   % This suite must NOT vacuously pass: setupOnce asserts the staged firn
   % cases are actually present so a missing firn tree fails the suite rather
   % than silently skipping it.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Install the canonical test/demo config (resolves the demo/data eval root)
   % and require the committed firn cases to be present.

   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;

   testCase.TestData.tmpdir = tempname(fullfile( ...
      icemodel.getpath('test'), 'artifacts', 'tmp'));
   icemodel.helpers.ensureDirExists(testCase.TestData.tmpdir);

   % The staged PROMICE firn anchors. These are committed (M2a); the suite is
   % only meaningful when they enumerate, so require them explicitly. The KAN
   % transect (kanl/kanm/kanu) is joined by the firn-accumulation demo fixture
   % egp (accumulation/dry-firn, 2015-2016) staged for CI firn-zone coverage
   % beyond the KAN ablation/percolation anchors.
   testCase.TestData.expected_firn_ids = ["kanl"; "kanm"; "kanu"];
   testCase.TestData.expected_accum_firn_ids = "egp";

   firn_cases = icemodel.verification.listcases(dataset_family="promice");
   testCase.assertNotEmpty(firn_cases, ...
      'no staged promice firn cases enumerated; M2a firn tree missing');
   testCase.TestData.firn_cases = firn_cases;

   % The committed SUMup minimal fixture: 3 KAN-co-located firn_observational
   % cases. The suite is only meaningful when they enumerate, so require them.
   % SUMup cases share the compact KAN ids with PROMICE (the family folder
   % is the namespace; no redundant sumup prefix). Resolve them with the
   % dataset_family filter to disambiguate from the PROMICE kanl/kanm/kanu.
   testCase.TestData.expected_sumup_ids = ["kanl"; "kanm"; "kanu"];
   testCase.TestData.expected_accum_sumup_ids = "egp";
   sumup_cases = icemodel.verification.listcases(dataset_family="sumup");
   testCase.assertNotEmpty(sumup_cases, ...
      'no staged sumup firn cases enumerated; committed fixture missing');
   testCase.TestData.sumup_cases = sumup_cases;
end

function teardownOnce(testCase)
   if isfield(testCase.TestData, 'tmpdir') && ...
         exist(testCase.TestData.tmpdir, 'dir') == 7
      rmdir(testCase.TestData.tmpdir, 's')
   end
   clear testCase.TestData.cleanup
end

function test_namelist_includes_firn_observational(testCase)
   % The case-type namelist must publish firn_observational (M2c). It must NOT
   % publish firn_analytical (deferred with the Meyer-Hewitt namespace).

   case_types = icemodel.verification.namelists.casetype();
   testCase.verifyTrue(ismember("firn_observational", case_types));
   testCase.verifyFalse(ismember("firn_analytical", case_types));

   % The firn families are valid dataset-family filter selectors.
   families = icemodel.verification.namelists.datasetfamily();
   testCase.verifyTrue(ismember("promice", families));
   testCase.verifyTrue(ismember("sumup", families));
end

function test_listcases_enumerates_firn_family_alongside_snow(testCase)
   % listcases (unfiltered) must surface the firn promice cases together with
   % the snow cases. This is the M2a gap M2c closes: familyManifestFiles now
   % globs the firn root as well as the snow root.

   all_cases = icemodel.verification.listcases();
   ids = [all_cases.case_id];

   % Firn cases enumerate alongside the snow lane.
   testCase.verifyTrue(all(ismember(testCase.TestData.expected_firn_ids, ids)), ...
      'firn promice cases not enumerated by unfiltered listcases');

   % The SUMup firn family enumerates alongside promice and the snow lane.
   testCase.verifyTrue(all(ismember(testCase.TestData.expected_sumup_ids, ids)), ...
      'sumup firn cases not enumerated by unfiltered listcases');

   % The snow lane is still present (regression guard for the harmonized
   % heterogeneous-schema vertcat).
   testCase.verifyTrue(all(ismember(["cdp"; "colbeck1976"], ids)), ...
      'snow cases dropped when firn families were added');
end

function test_each_firn_case_type_validates(testCase)
   % Every enumerated firn case must carry case_type=firn_observational and
   % that value must be in the canonical namelist.

   case_types = icemodel.verification.namelists.casetype();
   for id = testCase.TestData.expected_firn_ids'
      manifest = icemodel.verification.loadmanifest(id);
      testCase.verifyEqual(string(manifest.case_type), "firn_observational", ...
         sprintf('%s is not firn_observational', id));
      testCase.verifyTrue(ismember(string(manifest.case_type), case_types), ...
         sprintf('%s case_type not in namelist', id));
      testCase.verifyEqual(string(manifest.dataset_family), "promice");
   end
end

function test_each_firn_case_carries_valid_surface_zone(testCase)
   % Every promice firn case must carry a surface_zone single-sourced from
   % promicesiteinfo and validating against the canonical (glaciological-zone-
   % only) namelist. The KAN transect zones are pinned: kanl=ablation,
   % kanm=ablation, kanu=percolation. bare_ice / seasonal_snow are NOT zones.

   allowed = icemodel.verification.namelists.surfacezone();
   expected = struct('kanl', "ablation", ...
      'kanm', "ablation", 'kanu', "percolation");

   % The zone namelist must NOT carry capability descriptors.
   testCase.verifyFalse(ismember("bare_ice", allowed));
   testCase.verifyFalse(ismember("seasonal_snow", allowed));

   for id = testCase.TestData.expected_firn_ids'
      manifest = icemodel.verification.loadmanifest(id);
      zone = string(manifest.surface_zone);
      testCase.verifyTrue(ismember(zone, allowed), ...
         sprintf('%s surface_zone "%s" not in namelist', id, zone));
      testCase.verifyEqual(zone, expected.(char(id)), ...
         sprintf('%s surface_zone mismatch', id));
      % The manifest must agree with the single source of truth.
      testCase.verifyEqual(zone, ...
         string(icemodel.verification.helpers.promicesiteinfo( ...
         manifest.site_id).surface_zone));
   end
end

function test_each_firn_case_carries_valid_permafrost_zone(testCase)
   % Every promice firn case carries a permafrost_zone (ORTHOGONAL to
   % surface_zone) single-sourced from promicesiteinfo and validating against the
   % canonical namelist. The KAN anchors sit on the ice sheet (not permafrost
   % ground) so all three are "none".

   allowed = icemodel.verification.namelists.permafrostzone();
   testCase.verifyTrue(ismember("none", allowed));
   testCase.verifyTrue(ismember("continuous", allowed));
   % permafrost_zone is NOT a surface_zone value (orthogonal vocabularies).
   testCase.verifyFalse(ismember("continuous", ...
      icemodel.verification.namelists.surfacezone()));

   for id = testCase.TestData.expected_firn_ids'
      manifest = icemodel.verification.loadmanifest(id);
      testCase.verifyTrue(isfield(manifest, 'permafrost_zone'), ...
         sprintf('%s missing permafrost_zone field', id));
      pfz = string(manifest.permafrost_zone);
      testCase.verifyTrue(ismember(pfz, allowed), ...
         sprintf('%s permafrost_zone "%s" not in namelist', id, pfz));
      testCase.verifyEqual(pfz, "none", ...
         sprintf('%s (KAN ice-sheet anchor) permafrost_zone should be none', id));
      % The manifest must agree with the single source of truth.
      testCase.verifyEqual(pfz, string( ...
         icemodel.verification.helpers.promicesiteinfo( ...
         manifest.site_id).permafrost_zone));
   end
end

function test_each_firn_case_carries_valid_eval_target(testCase)
   % Every promice firn case must carry an eval_target string array single-
   % sourced from promicesiteinfo and validating against the eval_target
   % namelist. KAN pins: kanl/kanm=["seasonal_snow","bare_ice"],
   % kanu=["seasonal_snow","firn"].

   allowed = icemodel.verification.namelists.evaltarget();
   expected = struct( ...
      'kanl', ["seasonal_snow"; "bare_ice"], ...
      'kanm', ["seasonal_snow"; "bare_ice"], ...
      'kanu', ["seasonal_snow"; "firn"]);

   for id = testCase.TestData.expected_firn_ids'
      manifest = icemodel.verification.loadmanifest(id);
      target = string(manifest.eval_target);
      testCase.verifyTrue(all(ismember(target, allowed)), ...
         sprintf('%s eval_target not in namelist', id));
      testCase.verifyEqual(sort(target(:)), sort(expected.(char(id))(:)), ...
         sprintf('%s eval_target mismatch', id));
      % The manifest must agree with the single source of truth.
      testCase.verifyEqual(sort(target(:)), sort(string( ...
         icemodel.verification.helpers.promicesiteinfo( ...
         manifest.site_id).eval_target(:))));
   end
end

function test_manifest_is_metadata_only(testCase)
   % The firn manifest is METADATA-ONLY: it records WHICH forcing/eval sources
   % are available (by id) and the colocation regime, NOT a bundled
   % evaluation.mat/reference.mat data copy. No per-case bundle file is staged.

   needed = icemodel.verification.setup.firnCaseManifestFieldNames();
   family_root = fileparts(icemodel.verification.loadmanifest("kanl").manifest_path);

   for id = testCase.TestData.expected_firn_ids'
      manifest = icemodel.verification.loadmanifest(id);

      % Metadata-only schema fields present.
      for f = ["period", "forcing_sources", "eval_sources", "colocation"]
         testCase.verifyTrue(isfield(manifest, f), ...
            sprintf('%s missing metadata-only field %s', id, f));
      end
      testCase.verifyTrue(all(ismember(needed, string(fieldnames(manifest)))), ...
         sprintf('%s manifest missing a canonical firn field', id));

      % Source ids recorded by name (not bundled data).
      testCase.verifyTrue(all(ismember(["promice", "mar", "merra"], ...
         string(manifest.forcing_sources))), ...
         sprintf('%s forcing_sources incomplete', id));
      testCase.verifyTrue(all(ismember(["promice_obs", "racmo"], ...
         string(manifest.eval_sources))), ...
         sprintf('%s eval_sources incomplete', id));

      % NO bundled evaluation.mat / reference.mat colocation copy on disk.
      testCase.verifyFalse(isfile(fullfile(family_root, id, "evaluation.mat")), ...
         sprintf('%s still carries a bundled evaluation.mat', id));
      testCase.verifyFalse(isfile(fullfile(family_root, id, "reference.mat")), ...
         sprintf('%s still carries a bundled reference.mat', id));
   end
end

function test_candidate_adapter_resolves_declared_firn_variables(testCase)
   % The firn candidate adapter must map the firn comparison variables the
   % staged manifests declare (ablation, snow_depth, tsfc, tice1..tice8) from
   % synthetic icemodel-shaped output, including the thermistor profile
   % T(z,t) sampled from ice2.T at staged thermistor depths.

   manifest = icemodel.verification.loadmanifest("kanm");
   vars = string(manifest.comparison_variables);

   % Build a small synthetic icemodel output covering the declared axes.
   time = (datetime(2013, 6, 1, 0, 0, 0) + hours(0:5))';
   n = numel(time);
   ntice = nnz(startsWith(vars, "tice"));

   ice1 = struct("Time", time, ...
      "ablation", linspace(0, 0.1, n)', ...
      "snow_depth_m", linspace(0.5, 0.2, n)', ...
      "Tsfc", 263.15 + (1:n)' * 0.1);

   % A column T(z,t): rows are depth nodes (dz_thermal spacing), cols are time.
   nz = 40;
   T_column = 263.15 + repmat((0:nz - 1)' * 0.05, 1, n);
   ice2 = struct("T", T_column);
   dz = 0.04;

   % Provide the thermistor depths the firn adapter samples T(z,t) at.
   obsvars = manifest.observation_variables;
   if ~isstruct(obsvars)
      obsvars = struct();
   end
   obsvars.thermistor_depths_m = (1:max(ntice, 8))' * 0.5;

   opts = struct("smbmodel", "icemodel", "sitename", "kanm", ...
      "simyears", 2013, "dz_thermal", dz);
   adapter_manifest = struct( ...
      "case_type", "firn_observational", ...
      "comparison_variables", vars, ...
      "observation_variables", obsvars);

   candidate = icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, ice2, opts, adapter_manifest);

   testCase.verifyEqual(string(candidate.format), "timeseries");
   resolved = string(candidate.data.Properties.VariableNames);

   Tf = icemodel.physicalConstant('Tf');

   % Surface / ablation / snow-depth axes map and convert correctly.
   testCase.verifyTrue(ismember("ablation", resolved));
   testCase.verifyEqual(candidate.data.ablation, ice1.ablation);
   testCase.verifyTrue(ismember("snow_depth", resolved));
   testCase.verifyEqual(candidate.data.snow_depth, ice1.snow_depth_m);
   testCase.verifyTrue(ismember("tsfc", resolved));
   testCase.verifyEqual(candidate.data.tsfc, ice1.Tsfc - Tf);

   % The thermistor profile T(z,t) resolves for every staged tice axis.
   tice_vars = vars(startsWith(vars, "tice"));
   testCase.verifyNotEmpty(tice_vars);
   for tv = tice_vars'
      testCase.verifyTrue(ismember(tv, resolved), ...
         sprintf('thermistor axis %s not resolved by firn adapter', tv));
   end
end

function test_comparecase_soft_gate_no_hard_fail(testCase)
   % comparecase on a firn case must report SOFT diagnostic metrics with no
   % hard PASS/FAIL gate. The staged smoke reference (RACMO Data) does not
   % carry the firn observation axes, so the firn obs variables are reported
   % as missing_candidate_variable - a per-variable diagnostic, never a fatal
   % case error.

   result = icemodel.verification.comparecase("kanl", ...
      "artifact_dir", testCase.TestData.tmpdir, ...
      "make_plot", false);

   % Soft gate classification.
   testCase.verifyEqual(string(result.gate_mode), "soft");

   % Metrics are produced (one row per comparison variable) and the case does
   % not throw. status is a per-variable diagnostic, not a hard gate.
   testCase.verifyGreaterThan(height(result.metrics), 0);
   testCase.verifyTrue(exist(result.metrics_path, 'file') == 2);

   % No row carries a hard PASS/FAIL marker; the soft lane only emits the
   % diagnostic status vocabulary.
   allowed = ["ok"; "missing_target_variable"; ...
      "missing_candidate_variable"; "no_overlap"];
   testCase.verifyTrue(all(ismember(string(result.metrics.status), allowed)), ...
      'firn soft gate produced an unexpected (hard) status');
end

function test_colocated_files_resolve_on_disk(testCase)
   % The individual forcing/Data files each metadata-only firn manifest declares
   % (promice/mar/merra met + promice/racmo Data) must resolve on disk under the
   % standard icemodel input layout. Colocation is recorded as metadata
   % (manifest.colocation) pointing at these individual files, NOT a bundle.

   input_root = icemodel.verification.helpers.inputDataRoot();
   met_dir = fullfile(input_root, 'met');
   ud_dir = fullfile(input_root, 'userdata');

   for id = testCase.TestData.expected_firn_ids'
      manifest = icemodel.verification.loadmanifest(id);
      cf = manifest.colocation;

      % All four model legs are recorded as colocation metadata.
      testCase.verifyTrue(all(isfield(cf, ...
         {'promice', 'mar', 'merra', 'racmo'})), ...
         sprintf('%s colocation missing a model leg', id));

      % Met files (promice/mar/merra) resolve under met/.
      verifyFilesResolve(testCase, met_dir, cf.promice.met_files, id, "promice met");
      verifyFilesResolve(testCase, met_dir, cf.mar.met_files, id, "mar met");
      verifyFilesResolve(testCase, met_dir, cf.merra.met_files, id, "merra met");

      % Data files (promice/racmo) resolve under userdata/.
      verifyFilesResolve(testCase, ud_dir, cf.promice.data_files, id, "promice data");
      verifyFilesResolve(testCase, ud_dir, cf.racmo.data_files, id, "racmo data");
   end
end

function test_sumup_cases_are_firn_observational(testCase)
   % Every committed SUMup case is a metadata-only firn_observational case
   % (soft gate), in the sumup family, with a valid case_type.

   case_types = icemodel.verification.namelists.casetype();
   for id = testCase.TestData.expected_sumup_ids'
      manifest = icemodel.verification.loadmanifest(id, dataset_family="sumup");
      testCase.verifyEqual(string(manifest.case_type), "firn_observational", ...
         sprintf('%s is not firn_observational', id));
      testCase.verifyTrue(ismember(string(manifest.case_type), case_types), ...
         sprintf('%s case_type not in namelist', id));
      testCase.verifyEqual(string(manifest.dataset_family), "sumup");
   end
end

function test_sumup_cases_inherit_kan_zone_and_target(testCase)
   % The KAN-co-located SUMup cases inherit the anchor classification from
   % promicesiteinfo (the single source of truth): kanl/kanm=ablation +
   % [seasonal_snow,bare_ice]; kanu=percolation + [seasonal_snow,firn]. All
   % three carry permafrost_zone=none (KAN ice-sheet anchors).

   zone_ok = icemodel.verification.namelists.surfacezone();
   target_ok = icemodel.verification.namelists.evaltarget();
   pfz_ok = icemodel.verification.namelists.permafrostzone();

   expected = struct( ...
      'kanl', struct('site', "KAN_L", 'zone', "ablation"), ...
      'kanm', struct('site', "KAN_M", 'zone', "ablation"), ...
      'kanu', struct('site', "KAN_U", 'zone', "percolation"));

   for id = testCase.TestData.expected_sumup_ids'
      manifest = icemodel.verification.loadmanifest(id, dataset_family="sumup");
      exp = expected.(char(id));
      anchor = icemodel.verification.helpers.promicesiteinfo(exp.site);

      zone = string(manifest.surface_zone);
      testCase.verifyTrue(ismember(zone, zone_ok), ...
         sprintf('%s surface_zone "%s" not in namelist', id, zone));
      testCase.verifyEqual(zone, exp.zone, ...
         sprintf('%s surface_zone mismatch', id));
      testCase.verifyEqual(zone, string(anchor.surface_zone), ...
         sprintf('%s surface_zone does not match anchor %s', id, exp.site));

      target = string(manifest.eval_target);
      testCase.verifyTrue(all(ismember(target, target_ok)), ...
         sprintf('%s eval_target not in namelist', id));
      testCase.verifyEqual(sort(target(:)), sort(string(anchor.eval_target(:))), ...
         sprintf('%s eval_target does not match anchor %s', id, exp.site));

      pfz = string(manifest.permafrost_zone);
      testCase.verifyTrue(ismember(pfz, pfz_ok), ...
         sprintf('%s permafrost_zone "%s" not in namelist', id, pfz));
      testCase.verifyEqual(pfz, "none", ...
         sprintf('%s (KAN ice-sheet anchor) permafrost_zone should be none', id));
   end
end

function test_bundled_eval_target_is_forcing_agnostic_data_only(testCase)
   % The re-architecture central contract: the eval target is a data-only
   % observations.mat bundle (a `targets` struct with format + data/metadata,
   % NO forcing), referenced via evaluation_path. Verify it across the two
   % observational families that ship the bundle on disk (esm_snowmip cdp +
   % the SUMup KAN cases): the bundle loads through the standard loadArtifact
   % path, carries a recognized format and a data payload, and bundles no
   % forcing field. This is the forcing-agnostic eval contract; the forcing is
   % discovered separately at runtime, never inside observations.mat.

   forcing_fields = ["forcing", "met", "forcings", "tair", "swd", "lwd"];

   % esm_snowmip cdp: timeseries obs bundle.
   esm = icemodel.verification.loadmanifest("cdp");
   testCase.assertTrue(isfile(esm.evaluation_path), ...
      'cdp observations.mat missing on disk');
   esm_targets = icemodel.verification.helpers.loadArtifact( ...
      esm.evaluation_path, "targets");
   testCase.verifyEqual(string(esm_targets.format), "timeseries");
   testCase.verifyTrue(isfield(esm_targets, 'data'), ...
      'cdp eval target carries no data payload');
   testCase.verifyTrue(istimetable(esm_targets.data), ...
      'cdp eval target data is not a timetable');
   testCase.verifyFalse(any(isfield(esm_targets, cellstr(forcing_fields))), ...
      'cdp observations.mat bundles forcing (must be forcing-agnostic)');

   % SUMup KAN cases: subsurface_profile_bundle obs bundles.
   for id = testCase.TestData.expected_sumup_ids'
      su = icemodel.verification.loadmanifest(id, dataset_family="sumup");
      testCase.assertTrue(isfile(su.evaluation_path), ...
         sprintf('%s sumup observations.mat missing on disk', id));
      su_targets = icemodel.verification.helpers.loadArtifact( ...
         su.evaluation_path, "targets");
      testCase.verifyTrue(strlength(string(su_targets.format)) > 0, ...
         sprintf('%s sumup eval target has no format', id));
      testCase.verifyTrue(isfield(su_targets, 'data'), ...
         sprintf('%s sumup eval target carries no data payload', id));
      testCase.verifyFalse(any(isfield(su_targets, cellstr(forcing_fields))), ...
         sprintf('%s sumup observations.mat bundles forcing', id));
   end
end

function test_forcing_sources_is_informational_not_load_bearing(testCase)
   % forcing_sources is INFORMATIONAL only: the eval (observations.mat /
   % evaluation_path) and the comparison contract (comparison_variables) drive
   % comparecase, never forcing_sources. Verify each firn case carries
   % forcing_sources as a metadata string list AND that the comparison driver
   % resolves the case without consulting it (comparecase reports soft metrics
   % regardless of what forcing_sources lists).

   for case_entry = reshape(testCase.TestData.firn_cases, 1, [])
      id = string(case_entry.case_id);
      manifest = icemodel.verification.loadmanifest(id, dataset_family="promice");
      testCase.verifyTrue(isfield(manifest, 'forcing_sources'), ...
         sprintf('%s missing forcing_sources', id));
      % It is a recorded id list (strings), not a data handle.
      testCase.verifyTrue(isstring(string(manifest.forcing_sources)), ...
         sprintf('%s forcing_sources is not a string id list', id));
   end
end

function test_sumup_obs_files_resolve_on_disk(testCase)
   % Each committed SUMup case must reference an observation profile bundle
   % (colocation.sumup.obs_file) that resolves on disk, and listcases must
   % surface it as the case evaluation_path.

   for case_entry = reshape(testCase.TestData.sumup_cases, 1, [])
      id = case_entry.case_id;
      testCase.verifyTrue(isfield(case_entry.colocation, 'sumup'), ...
         sprintf('%s missing sumup colocation leg', id));
      testCase.verifyNotEmpty(case_entry.colocation.sumup.obs_file, ...
         sprintf('%s declares no obs_file', id));
      testCase.verifyTrue(strlength(case_entry.evaluation_path) > 0, ...
         sprintf('%s evaluation_path not resolved from obs_file', id));
      testCase.verifyTrue(isfile(case_entry.evaluation_path), ...
         sprintf('%s obs bundle missing on disk: %s', id, ...
         case_entry.evaluation_path));
   end
end

function test_sumup_candidate_adapter_maps_profile_variables(testCase)
   % The firn candidate adapter must map the SUMup profile-bundle comparison
   % variables (density rho(z), subsurface_temperature T(z,t)) from synthetic
   % icemodel column state, only as far as the staged cases declare. SUMup
   % cases declare profile-bundle axes, NOT the PROMICE thermistor series, so
   % the candidate is a subsurface_profile_bundle.

   manifest = icemodel.verification.loadmanifest("kanu", dataset_family="sumup");
   vars = string(manifest.comparison_variables);
   testCase.assertTrue(ismember("density", vars));
   testCase.assertTrue(ismember("subsurface_temperature", vars));

   nz = 40;
   n = 6;
   time = (datetime(2013, 6, 1) + hours(0:n - 1))';
   ice1 = struct("Time", time, "accumulation", linspace(0, 0.2, n)');
   ice2 = struct( ...
      "T", 263.15 + repmat((0:nz - 1)' * 0.05, 1, n), ...
      "ro_sno", 350 + (0:nz - 1)' * 5);
   opts = struct("smbmodel", "icemodel", "sitename", "kanu", ...
      "simyears", 2013, "dz_thermal", 0.04);

   adapter_manifest = struct( ...
      "case_type", "firn_observational", ...
      "comparison_variables", vars, ...
      "observation_variables", manifest.observation_variables);

   candidate = icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, ice2, opts, adapter_manifest);

   % The SUMup case yields a profile bundle, not a timeseries.
   testCase.verifyEqual(string(candidate.format), "subsurface_profile_bundle");
   Tf = icemodel.physicalConstant('Tf');

   % Density profile rho(z) resolves with a depth axis.
   testCase.verifyTrue(isfield(candidate.data, 'density'));
   testCase.verifyEqual(candidate.data.density.value, ice2.ro_sno);

   % Subsurface temperature profile T(z,t) resolves in degrees C.
   testCase.verifyTrue(isfield(candidate.data, 'subsurface_temperature'));
   testCase.verifyEqual(candidate.data.subsurface_temperature, ice2.T - Tf);
end

function test_committed_accum_firn_fixtures_enumerate(testCase)
   % The firn-accumulation demo fixture (egp, dry-firn) must enumerate alongside
   % the KAN transect in BOTH the promice and sumup families, so CI carries
   % firn-zone dry-firn coverage beyond the KAN ablation/percolation anchors.

   promice_ids = [icemodel.verification.listcases(dataset_family="promice").case_id];
   testCase.verifyTrue(all(ismember( ...
      testCase.TestData.expected_accum_firn_ids, promice_ids)), ...
      'firn-accumulation promice fixture (egp) not enumerated');

   sumup_ids = [icemodel.verification.listcases(dataset_family="sumup").case_id];
   testCase.verifyTrue(all(ismember( ...
      testCase.TestData.expected_accum_sumup_ids, sumup_ids)), ...
      'firn-accumulation sumup fixture (egp) not enumerated');
end

function test_all_committed_promice_fixtures_validate(testCase)
   % DATA-DRIVEN contract: EVERY committed promice firn fixture (KAN transect +
   % the egp accumulation fixture) must validate against the canonical
   % namelists AND agree with promicesiteinfo (the single source of truth) for
   % case_type / surface_zone / eval_target / permafrost_zone. This enumerates
   % the fixtures rather than hard-coding KAN-specific values, so new committed
   % firn-accumulation cases are validated automatically.

   case_types = icemodel.verification.namelists.casetype();
   zone_ok = icemodel.verification.namelists.surfacezone();
   target_ok = icemodel.verification.namelists.evaltarget();
   pfz_ok = icemodel.verification.namelists.permafrostzone();

   for case_entry = reshape(testCase.TestData.firn_cases, 1, [])
      id = string(case_entry.case_id);
      manifest = icemodel.verification.loadmanifest(id, dataset_family="promice");
      anchor = icemodel.verification.helpers.promicesiteinfo(manifest.site_id);

      testCase.verifyEqual(string(manifest.case_type), "firn_observational", ...
         sprintf('%s case_type', id));
      testCase.verifyTrue(ismember(string(manifest.case_type), case_types), ...
         sprintf('%s case_type not in namelist', id));

      zone = string(manifest.surface_zone);
      testCase.verifyTrue(ismember(zone, zone_ok), ...
         sprintf('%s surface_zone "%s" not in namelist', id, zone));
      testCase.verifyEqual(zone, string(anchor.surface_zone), ...
         sprintf('%s surface_zone disagrees with promicesiteinfo', id));

      target = string(manifest.eval_target);
      testCase.verifyTrue(all(ismember(target, target_ok)), ...
         sprintf('%s eval_target not in namelist', id));
      testCase.verifyEqual(sort(target(:)), sort(string(anchor.eval_target(:))), ...
         sprintf('%s eval_target disagrees with promicesiteinfo', id));

      pfz = string(manifest.permafrost_zone);
      testCase.verifyTrue(ismember(pfz, pfz_ok), ...
         sprintf('%s permafrost_zone "%s" not in namelist', id, pfz));
      testCase.verifyEqual(pfz, string(anchor.permafrost_zone), ...
         sprintf('%s permafrost_zone disagrees with promicesiteinfo', id));
   end
end

function test_accum_firn_fixtures_soft_gate_no_hard_fail(testCase)
   % comparecase on the firn-accumulation fixture (egp) must report a SOFT
   % diagnostic gate with no hard PASS/FAIL, mirroring the KAN soft-gate
   % contract. This guards that the new committed cases run end to end.

   allowed = ["ok"; "missing_target_variable"; ...
      "missing_candidate_variable"; "no_overlap"];
   for id = testCase.TestData.expected_accum_firn_ids'
      result = icemodel.verification.comparecase(id, ...
         "artifact_dir", testCase.TestData.tmpdir, "make_plot", false);
      testCase.verifyEqual(string(result.gate_mode), "soft", ...
         sprintf('%s is not soft-gated', id));
      testCase.verifyGreaterThan(height(result.metrics), 0, ...
         sprintf('%s produced no metrics', id));
      testCase.verifyTrue(all(ismember(string(result.metrics.status), allowed)), ...
         sprintf('%s soft gate produced an unexpected (hard) status', id));
   end
end

function test_accum_firn_colocation_files_resolve_on_disk(testCase)
   % The individual forcing/Data files each firn-accumulation fixture declares
   % (promice/mar/merra met + promice/racmo Data) must resolve on disk under the
   % standard icemodel input layout, exactly as the KAN fixtures do.

   input_root = icemodel.verification.helpers.inputDataRoot();
   met_dir = fullfile(input_root, 'met');
   ud_dir = fullfile(input_root, 'userdata');

   for id = testCase.TestData.expected_accum_firn_ids'
      manifest = icemodel.verification.loadmanifest(id, dataset_family="promice");
      cf = manifest.colocation;
      verifyFilesResolve(testCase, met_dir, cf.promice.met_files, id, "promice met");
      verifyFilesResolve(testCase, met_dir, cf.mar.met_files, id, "mar met");
      verifyFilesResolve(testCase, met_dir, cf.merra.met_files, id, "merra met");
      verifyFilesResolve(testCase, ud_dir, cf.promice.data_files, id, "promice data");
      verifyFilesResolve(testCase, ud_dir, cf.racmo.data_files, id, "racmo data");
   end
end

function test_all_committed_sumup_fixtures_resolve(testCase)
   % DATA-DRIVEN contract: every committed sumup fixture (KAN + egp) must be
   % a firn_observational case whose obs bundle resolves on disk and whose
   % anchor classification agrees with promicesiteinfo.

   case_types = icemodel.verification.namelists.casetype();
   for case_entry = reshape(testCase.TestData.sumup_cases, 1, [])
      id = string(case_entry.case_id);
      manifest = icemodel.verification.loadmanifest(id, dataset_family="sumup");
      testCase.verifyEqual(string(manifest.case_type), "firn_observational", ...
         sprintf('%s sumup case_type', id));
      testCase.verifyTrue(ismember(string(manifest.case_type), case_types), ...
         sprintf('%s sumup case_type not in namelist', id));

      % obs bundle resolves (evaluation_path comes from colocation.sumup.obs_file).
      testCase.verifyTrue(strlength(case_entry.evaluation_path) > 0, ...
         sprintf('%s sumup evaluation_path empty', id));
      testCase.verifyTrue(isfile(case_entry.evaluation_path), ...
         sprintf('%s sumup obs bundle missing: %s', id, ...
         case_entry.evaluation_path));

      % Anchor classification (when co-located) matches promicesiteinfo.
      anchor = icemodel.verification.helpers.promicesiteinfo(manifest.site_id);
      if strlength(string(manifest.surface_zone)) > 0
         testCase.verifyEqual(string(manifest.surface_zone), ...
            string(anchor.surface_zone), ...
            sprintf('%s sumup surface_zone disagrees with anchor', id));
      end
   end
end

%% Local helpers
function verifyFilesResolve(testCase, base, names, id, label)
   %VERIFYFILESRESOLVE Assert each relative filename resolves under base.
   names = string(names);
   names = names(strlength(names) > 0);
   testCase.verifyNotEmpty(names, ...
      sprintf('%s declared no %s files', id, label));
   for nm = reshape(names, 1, [])
      testCase.verifyTrue(exist(fullfile(base, nm), 'file') == 2, ...
         sprintf('%s %s file missing on disk: %s', id, label, nm));
   end
end
