function tests = test_firn_verification_contracts
   %TEST_FIRN_VERIFICATION_CONTRACTS Verify the firn-evaluation contracts.
   %
   % Exercises the firn_observational lane end to end against the committed
   % co-located PROMICE firn cases under demo/data/eval/firn/promice/:
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

function test_colocated_bundle_resolves_on_disk(testCase)
   % The co-located multi-model bundle each firn manifest declares
   % (promice/mar/merra met + promice/racmo Data) must resolve on disk under
   % the standard icemodel input layout.

   input_root = icemodel.verification.helpers.inputDataRoot();
   met_dir = fullfile(input_root, 'met');
   ud_dir = fullfile(input_root, 'userdata');

   for id = testCase.TestData.expected_firn_ids'
      manifest = icemodel.verification.loadmanifest(id);
      cf = manifest.colocated_forcing;

      % Evaluation + reference artifacts resolve.
      testCase.verifyTrue(exist(manifest.evaluation_path, 'file') == 2, ...
         sprintf('%s evaluation.mat missing', id));
      testCase.verifyTrue(exist(manifest.reference_path, 'file') == 2, ...
         sprintf('%s reference.mat missing', id));

      % All four model legs are recorded.
      testCase.verifyTrue(all(isfield(cf, ...
         {'promice', 'mar', 'merra', 'racmo'})), ...
         sprintf('%s colocated_forcing missing a model leg', id));

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
