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
   % only meaningful when they enumerate, so require them explicitly.
   testCase.TestData.expected_firn_ids = ["kanl"; "kanm"; "kanu"];

   firn_cases = icemodel.verification.listcases(dataset_family="promice");
   testCase.assertNotEmpty(firn_cases, ...
      'no staged promice firn cases enumerated; M2a firn tree missing');
   testCase.TestData.firn_cases = firn_cases;
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
