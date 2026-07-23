function tests = test_snow_verification_contracts
   %TEST_SNOW_VERIFICATION_CONTRACTS Verify the snow-validation framework.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install the full verification config used by this top-level data suite.

   [~, ~, ~, evaluation_root, cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment( ...
      icemodel_config_casename="verification");
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.tmpdir = tempname(fullfile( ...
      icemodel.getpath('test'), 'artifacts', 'tmp'));
   icemodel.helpers.ensureDirExists(testCase.TestData.tmpdir);

   % The full ESM-SnowMIP archive is local scientific data, not a release
   % fixture. A clean checkout must skip this archive contract explicitly.
   testCase.assumeTrue(isfile(fullfile( ...
      evaluation_root, 'esm_snowmip', 'manifest.json')), ...
      ['Full snow-verification archive is not installed under the ' ...
      'verification data root.'])
end

function teardown(testCase)
   % Remove temporary artifacts and restore the caller's config.

   if exist(testCase.TestData.tmpdir, 'dir') == 7
      rmdir(testCase.TestData.tmpdir, 's')
   end
   clear testCase.TestData.cleanup
end

function test_listcases_returns_expected_ids(testCase)
   % LISTCASES should expose the curated smoke cases from top-level data/eval.

   cases = icemodel.verification.listcases();
   ids = [cases.case_id];

   testCase.verifyTrue(all(ismember(expectedCaseIds(), ids)));
end

function test_runner_discovers_only_the_selected_root_inventory(testCase)
   % Initial discovery and later comparison must use one resolved root pair.

   selected_root = fullfile(testCase.TestData.tmpdir, 'selected-inventory');
   stageSelectedSnowInventory(selected_root, "cdp");

   % The repository verification tree also contains WFJ. Requesting both ids
   % proves a default-root discovery would return a second, unintended case.
   results = run_snow_verification_suite( ...
      cases=["cdp", "wfj"], data_root=selected_root, ...
      candidate_provider=@selectedTargetCandidate, make_plots=false);

   testCase.verifyEqual(numel(results.cases), 1)
   testCase.verifyEqual(string(results.cases.case_id), "cdp")
   selected_eval = string(fullfile(selected_root, 'eval'));
   testCase.verifyTrue(startsWith( ...
      string(results.cases.manifest_path), selected_eval + filesep))
   testCase.verifyFalse(any(string(results.summary.case_id) == "wfj"))
end

function test_verification_namelists_expose_curated_selectors(testCase)
   % The verification namespace should publish one canonical selector catalog.

   testCase.verifyEqual(icemodel.verification.namelists.datasetfamily(), ...
      expectedDatasetFamilies());
   testCase.verifyEqual(icemodel.verification.namelists.casetype(), ...
      expectedCaseTypes());
   testCase.verifyEqual(icemodel.verification.namelists.caseid(), ...
      expectedCaseIds());
   testCase.verifyEqual( ...
      icemodel.verification.namelists.caseid("esm_snowmip"), ...
      expectedEsmCaseIds());

   % Family-specific namelists used by case-id resolution and arg validation.
   testCase.verifyEqual(icemodel.verification.namelists.snowmipsite(), ...
      expectedEsmCaseIds());
   testCase.verifyEqual(icemodel.verification.namelists.laughtests(), ...
      "colbeck1976");

   % The richer site catalog is keyed by the same site-name namelist.
   catalog = icemodel.verification.setup.esmSnowmipSiteCatalog();
   testCase.verifyEqual(string({catalog.sitename})', expectedEsmCaseIds());
   info_cdp = icemodel.verification.setup.esmSnowmipSiteCatalog("cdp");
   testCase.verifyEqual(info_cdp.long_name, "Col de Porte");
end

function test_each_site_stages_expected_comparison_variables(testCase)
   % Per-site staged comparison_variables must include the canonical snow
   % column for every ESM-SnowMIP site, and surface_temp_C for sites whose
   % upstream obs files contain a usable surface-temperature channel.
   % This is the self-verification that catches regressions in the obs
   % builder / manifest schema when one site silently loses a variable.

   for sitename = expectedEsmCaseIds()'
      manifest = icemodel.verification.loadmanifest(sitename);
      vars = string(manifest.comparison_variables);

      % Every site stages snow depth; this is the most basic sanity check.
      testCase.verifyTrue(ismember("snow_depth_m", vars), ...
         sprintf('snow_depth_m missing for %s', sitename));
   end

   % Surface temperature is observed at WFJ (Weissfluhjoch). Pin it explicitly
   % here; sites without a usable 'ts' channel
   % (e.g. sod) are expected to drop the variable and are excluded.
   wfj = icemodel.verification.loadmanifest("wfj");
   testCase.verifyTrue(ismember("surface_temp_C", ...
      string(wfj.comparison_variables)), ...
      'WFJ must stage surface_temp_C as a comparison variable');
end

function test_each_esm_case_carries_land_zone_seasonal_snow_target(testCase)
   % Every ESM-SnowMIP case is an off-ice seasonal snowpack at a land site:
   % surface_zone="land" (glaciological zone) and eval_target=["seasonal_snow"]
   % (the capability the case exercises). seasonal_snow is NOT a zone - it moved
   % to the eval_target descriptor.

   zones = icemodel.verification.namelists.surfacezone();
   targets = icemodel.verification.namelists.evaltarget();
   testCase.verifyTrue(ismember("land", zones));
   testCase.verifyFalse(ismember("seasonal_snow", zones));
   testCase.verifyTrue(ismember("seasonal_snow", targets));

   for sitename = expectedEsmCaseIds()'
      manifest = icemodel.verification.loadmanifest(sitename);
      testCase.verifyEqual(string(manifest.surface_zone), "land", ...
         sprintf('%s surface_zone is not land', sitename));
      testCase.verifyEqual(string(manifest.eval_target), "seasonal_snow", ...
         sprintf('%s eval_target is not seasonal_snow', sitename));
   end
end

function test_each_esm_case_carries_valid_permafrost_zone(testCase)
   % Every ESM-SnowMIP case is an off-ice land site and carries a permafrost_zone
   % (ORTHOGONAL to surface_zone) sampled by point-in-polygon from the Obu et al.
   % (2019) permafrost-zone map: Sodankyla (boreal Lapland) is continuous, Senator
   % Beck (Colorado alpine) discontinuous, Swamp Angel / Weissfluhjoch sporadic,
   % Old Jack Pine isolated, and the remaining lower/mid-latitude sites fall
   % outside any permafrost polygon -> none. The value must validate against the
   % canonical namelist.

   allowed = icemodel.verification.namelists.permafrostzone();
   testCase.verifyTrue(ismember("none", allowed));
   testCase.verifyTrue(ismember("isolated", allowed));

   expectedByName = struct( ...
      'sod', "continuous", 'snb', "discontinuous", ...
      'swa', "sporadic",   'wfj', "sporadic", ...
      'ojp', "isolated");

   for sitename = expectedEsmCaseIds()'
      manifest = icemodel.verification.loadmanifest(sitename);
      testCase.verifyTrue(isfield(manifest, 'permafrost_zone'), ...
         sprintf('%s missing permafrost_zone field', sitename));
      pfz = string(manifest.permafrost_zone);
      testCase.verifyTrue(ismember(pfz, allowed), ...
         sprintf('%s permafrost_zone "%s" not in namelist', sitename, pfz));
      key = char(lower(string(sitename)));
      if isfield(expectedByName, key)
         expected = expectedByName.(key);
      else
         expected = "none";
      end
      testCase.verifyEqual(pfz, expected, ...
         sprintf('%s permafrost_zone mismatch', sitename));
   end
end

function test_forcing_includes_rainf_snowf_passthrough(testCase)
   % Site forcing files must include rainf/snowf channels so future
   % rain/snow-aware downstream consumers can use them directly.

   % Resolve the staged met file via the standard chain (createMetFileNames)
   % using the same opts a runIcemodelSnowCandidate run would build.
   manifest = icemodel.verification.loadmanifest("cdp");
   opts = icemodel.test.helpers.setModelOptsForCase(manifest);
   met_files = opts.metfname;
   testCase.assertNotEmpty(met_files);

   loaded = load(met_files{1}, 'met');
   names = string(loaded.met.Properties.VariableNames);
   testCase.verifyTrue(ismember("rainf", names));
   testCase.verifyTrue(ismember("snowf", names));
   testCase.verifyTrue(ismember("ppt", names));
end

function test_manifest_schema_helpers_are_setup_owned(testCase)
   % Manifest construction schemas should live with setup/update tooling.

   family_names = icemodel.verification.setup.familyManifestFieldNames();
   case_names = icemodel.verification.setup.caseManifestFieldNames();

   testCase.verifyEqual(family_names(:), expectedFamilyManifestFields());
   testCase.verifyTrue(ismember("case_id", case_names));
   testCase.verifyTrue(ismember("case_type", case_names));
end

function test_loadmanifest_resolves_repo_data_paths(testCase)
   % LOADMANIFEST should resolve full verification data from the top-level root.
   % Forcing is no longer in the manifest - it is staged under data/input/met/ via
   % the standard icemodel naming convention, so checks below cover evaluation
   % paths only.

   manifest = icemodel.verification.loadmanifest("cdp");
   verification_root = string(icemodel.internal.fullpath("data"));

   testCase.verifyEqual(manifest.dataset_family, "esm_snowmip");
   testCase.verifyEqual(manifest.case_type, "esm_site");
   testCase.verifyTrue(contains(manifest.evaluation_path, ...
      fullfile("eval", "esm_snowmip", "cdp")));
   testCase.verifyTrue(startsWith(string(manifest.evaluation_path), ...
      verification_root + filesep));
   testCase.verifyTrue(exist(manifest.evaluation_path, 'file') == 2);
end

function test_comparecase_metadata_only_writes_soft_metrics(testCase)
   % COMPARECASE on a metadata-only esm_site case must run from the committed
   % observations.mat obs bundle and report SOFT diagnostic metrics. With no
   % bundled reference.mat smoke copy and no model candidate supplied, the
   % default candidate is empty, so every comparison variable is reported as
   % missing_candidate_variable - a per-variable diagnostic, never a hard fail.

   result = icemodel.verification.comparecase("cdp", ...
      "artifact_dir", testCase.TestData.tmpdir, ...
      "make_plot", false);

   testCase.verifyEqual(string(result.gate_mode), "soft");
   testCase.verifyTrue(exist(result.metrics_path, 'file') == 2);
   testCase.verifyGreaterThan(height(result.metrics), 0);

   % The soft lane only emits the diagnostic status vocabulary; no hard marker.
   allowed = ["ok"; "missing_target_variable"; ...
      "missing_candidate_variable"; "no_overlap"];
   testCase.verifyTrue(all(ismember(string(result.metrics.status), allowed)), ...
      'esm soft gate produced an unexpected (hard) status');
end

function test_comparecase_handles_colbeck_bundle(testCase)
   % The synthetic Colbeck bundle should compare all three staged experiments.

   result = icemodel.verification.comparecase("colbeck1976", ...
      "artifact_dir", testCase.TestData.tmpdir, ...
      "make_plot", false);

   testCase.verifyEqual(sort(unique(result.metrics.experiment))', ...
      expectedColbeckExperimentIds());
   testCase.verifyTrue(all(result.metrics.n > 0));
end

function test_icemodel_candidate_provider_runs_model_entry_point(testCase)
   % The verification lane must accept candidates produced by icemodel(opts).

   manifest = icemodel.verification.loadmanifest("cdp");
   candidate = icemodel.verification.runIcemodelSnowCandidate(manifest);
   result = icemodel.verification.comparecase("cdp", ...
      "artifact_dir", testCase.TestData.tmpdir, ...
      "candidate", candidate, ...
      "make_plot", false);

   testCase.verifyEqual(candidate.format, "timeseries");
   testCase.verifyTrue(ismember("snow_depth_m", ...
      string(candidate.data.Properties.VariableNames)));
   testCase.verifyTrue(all(result.metrics.status == "ok"));
   testCase.verifyTrue(any(result.metrics.rmse > 0));
end

function test_icemodel_candidate_preserves_missing_targets(testCase)
   % Synthetic candidates should not turn missing observations into zeros.

   manifest = icemodel.verification.loadmanifest("wfj");
   targets = icemodel.verification.helpers.loadArtifact( ...
      manifest.evaluation_path, "targets");
   candidate = icemodel.verification.runIcemodelSnowCandidate(manifest);
   result = icemodel.verification.comparecase("wfj", ...
      "artifact_dir", testCase.TestData.tmpdir, ...
      "candidate", candidate, ...
      "make_plot", false);

   missing_target = ~isfinite(targets.data.snow_depth_m);
   paired = isfinite(targets.data.snow_depth_m) ...
      & isfinite(candidate.data.snow_depth_m);
   row = result.metrics(result.metrics.variable == "snow_depth_m", :);

   testCase.verifyTrue(any(missing_target));
   testCase.verifyTrue(all(~isfinite(candidate.data.snow_depth_m(missing_target))));
   testCase.verifyEqual(double(row.n), nnz(paired));
   testCase.verifyLessThan(double(row.n), height(targets.data));
end

function test_candidate_adapter_derives_swe_from_depth_and_density(testCase)
   % Snow-model outputs can use core-adjacent names and derive SWE in adapter.

   % These named fixture values make the adapter contract explicit: depth and
   % density are native model-like outputs, while SWE and Celsius surface
   % temperature are derived verification variables.
   fixture_start_time = datetime(2000, 1, 1, 0, 0, 0);
   fixture_sample_hours = 0:2;
   snow_depth_m = [0.1; 0.2; 0.3];
   snow_density_kg_m3 = [250; 300; 350];
   surface_temp_K = [263.15; 264.15; 265.15];

   time = fixture_start_time + hours(fixture_sample_hours);
   ice1 = struct( ...
      "Time", time(:), ...
      "snow_depth", snow_depth_m, ...
      "snow_density_kg_m3", snow_density_kg_m3, ...
      "Tsfc", surface_temp_K);
   ice2 = struct();
   opts = struct("smbmodel", "icemodel", "sitename", "verification", ...
      "simyears", year(fixture_start_time));
   manifest = struct( ...
      "case_type", "esm_site", ...
      "comparison_variables", ...
      expectedCoreSiteVariables());

   candidate = icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, ice2, opts, manifest);

   Tf = icemodel.physicalConstant('Tf');
   testCase.verifyEqual(candidate.data.snow_depth_m, snow_depth_m);
   testCase.verifyEqual(candidate.data.swe_kg_m2, ...
      snow_depth_m .* snow_density_kg_m3);
   testCase.verifyEqual(candidate.data.surface_temp_C, surface_temp_K - Tf);
end

function test_candidate_adapter_samples_soil_temp_from_ice2(testCase)
   % soil_temp_<k>_C must be sampled from ice2.T at the manifested depth.

   time = datetime(2000, 1, 1, 0, 0, 0) + hours(0:2);
   T_column = [273.15 273.05 273.25; ...
               272.15 272.10 272.20; ...
               271.15 271.20 271.30];
   ice1 = struct("Time", time(:), "Tsfc", T_column(1, :)');
   ice2 = struct("T", T_column);
   opts = struct("smbmodel", "icemodel", "sitename", "wfj", ...
      "simyears", 2000, "dz_thermal", 0.04);
   manifest = struct( ...
      "case_type", "esm_site", ...
      "comparison_variables", "soil_temp_1_C", ...
      "observation_variables", struct("soil_depths_m", 0.01));

   candidate = icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, ice2, opts, manifest);

   Tf = icemodel.physicalConstant('Tf');
   testCase.verifyTrue(ismember("soil_temp_1_C", ...
      candidate.data.Properties.VariableNames));
   testCase.verifyEqual(candidate.data.soil_temp_1_C, T_column(1, :)' - Tf);
end

function test_plotcase_writes_figure_without_candidate(testCase)
   % Users should be able to visualize staged targets directly.

   outfile = fullfile(testCase.TestData.tmpdir, 'cdp_plot.png');
   f = icemodel.verification.plotcase("cdp", ...
      "visible", "off", ...
      "output_file", outfile);

   testCase.verifyTrue(exist(outfile, 'file') == 2);
   lines = findall(f, 'Type', 'Line');
   testCase.verifyTrue(all(string(get(lines, 'Marker')) == "none"));
   testCase.verifyEmpty(findall(f, 'Type', 'Line', 'LineStyle', 'none'));
   close(f)
end

function test_colbeck_compare_solutions_is_line_only(testCase)
   % The dedicated four-way Colbeck driver shares the no-marker contract.
   old_close = get(groot, 'DefaultFigureCloseRequestFcn');
   restore_close = onCleanup(@() set(groot, ...
      'DefaultFigureCloseRequestFcn', old_close));
   set(groot, 'DefaultFigureCloseRequestFcn', @(~, ~) []);
   before = findall(groot, 'Type', 'Figure');

   result = icemodel.verification.colbeck.compareSolutions( ...
      experiment_names="exp1", make_plot=true, save_plot=false, ...
      plot_visible="off");
   clear restore_close

   after = findall(groot, 'Type', 'Figure');
   is_new = arrayfun(@(candidate) ~any(candidate == before), after);
   fig = after(is_new);
   testCase.addTeardown(@() delete(fig(isgraphics(fig))));
   testCase.verifyNumElements(fig, 1);
   lines = findall(fig, 'Type', 'Line');
   testCase.verifyNotEmpty(lines);
   testCase.verifyTrue(all(string(get(lines, 'Marker')) == "none"));
   testCase.verifyEmpty(findall(fig, 'Type', 'Line', 'LineStyle', 'none'));
   expected_root = string(fullfile(icemodel.getpath('test'), 'data'));
   testCase.verifyEqual(string(fileparts(result.evaluation_data_root)), ...
      expected_root)
end

function test_comparecase_writes_separate_scatter_for_site_cases(testCase)
   % ESM site comparisons should keep time-series and scatter figures separate.

   manifest = icemodel.verification.loadmanifest("wfj");
   candidate = icemodel.verification.runIcemodelSnowCandidate(manifest);
   result = icemodel.verification.comparecase("wfj", ...
      "artifact_dir", testCase.TestData.tmpdir, ...
      "candidate", candidate, ...
      "make_plot", true, ...
      "save_plot", true, ...
      "plot_visible", "off");

   testCase.verifyTrue(exist(result.figure_path, 'file') == 2);
   testCase.verifyTrue(exist(result.scatter_figure_path, 'file') == 2);
end

function test_comparecase_skips_scatter_for_colbeck_bundle(testCase)
   % Current Colbeck bundle comparisons should not emit scatter figures.

   result = icemodel.verification.comparecase("colbeck1976", ...
      "artifact_dir", testCase.TestData.tmpdir, ...
      "make_plot", true, ...
      "save_plot", true, ...
      "plot_visible", "off");

   testCase.verifyTrue(exist(result.figure_path, 'file') == 2);
   testCase.verifyEqual(result.scatter_figure_path, "");
end

function test_colbeck_evaluation_carries_two_target_sources(testCase)
   % The unified Colbeck evaluation.mat carries both numerical_summa and
   % analytical_clark2017 target bundles.

   manifest = icemodel.verification.loadmanifest("colbeck1976");
   targets = icemodel.verification.helpers.loadArtifact( ...
      manifest.evaluation_path, "targets");

   testCase.verifyTrue(isfield(targets, 'numerical_summa'));
   testCase.verifyTrue(isfield(targets, 'analytical_clark2017'));
   testCase.verifyTrue(isfield(targets.numerical_summa, 'experiments'));
   testCase.verifyTrue(isfield(targets.analytical_clark2017, 'experiments'));
end

function test_plot_timeseries_shows_sparse_points_with_markers(testCase)
   % Generic sparse observation series should remain visible on a dense axis.

   % A short dense hourly axis with only three finite observations reproduces
   % the sparse SWE target pattern that originally rendered as a blank line.
   fixture_start_time = datetime(2000, 1, 1, 0, 0, 0);
   sample_hours = 0:9;
   observed_indices = [2 6 10];
   observed_values = [1 2 3];

   time = fixture_start_time + hours(sample_hours);
   values = nan(numel(sample_hours), 1);
   values(observed_indices) = observed_values;

   f = figure('Visible', 'off');
   cleaner = onCleanup(@() close(f));
   ax = axes(f);
   icemodel.plot.timeseries(time(:), values, axes=ax);

   lines = findall(ax, 'Type', 'line');

   testCase.verifyTrue(any(string(get(lines, 'Marker')) ~= "none"));
end

function ids = expectedCaseIds()
   %EXPECTEDCASEIDS Canonical case catalog: 10 ESM-SnowMIP sites + Colbeck.

   ids = ["cdp"; "oas"; "obs"; "ojp"; "rme"; "sap"; "snb"; "sod"; ...
      "swa"; "wfj"; "colbeck1976"];
end

function ids = expectedEsmCaseIds()
   %EXPECTEDESMCASEIDS All 10 ESM-SnowMIP reference sites.

   ids = ["cdp"; "oas"; "obs"; "ojp"; "rme"; "sap"; "snb"; "sod"; ...
      "swa"; "wfj"];
end

function names = expectedDatasetFamilies()
   %EXPECTEDDATASETFAMILIES Dataset families currently staged for verification.
   % snow/: esm_snowmip, laugh_tests. firn/source families include promice,
   % sumup, retmip, imau, and research_site.

   names = ["esm_snowmip"; "laugh_tests"; "promice"; "sumup"; ...
      "retmip"; "imau"; "research_site"];
end

function names = expectedCaseTypes()
   %EXPECTEDCASETYPES Verification case types used by manifests and runners.
   % firn_observational is the soft-gated firn evaluation case type; the
   % firn_analytical type is deferred with the Meyer-Hewitt namespace.

   names = ["esm_site"; "synthetic_process"; "firn_observational"];
end

function fields = expectedFamilyManifestFields()
   %EXPECTEDFAMILYMANIFESTFIELDS Required top-level family-manifest fields.

   fields = {'dataset_family'; 'source_doi'; 'source_url'; ...
      'source_version'; 'retrieval_date'; 'cases'; 'skipped'};
end

function ids = expectedColbeckExperimentIds()
   %EXPECTEDCOLBECKEXPERIMENTIDS Colbeck experiments staged from Laugh-Tests.

   ids = ["exp1", "exp2", "exp3"];
end

function names = expectedCoreSiteVariables()
   %EXPECTEDCORESITEVARIABLES Minimal site variables derived from model output.

   names = ["snow_depth_m"; "swe_kg_m2"; "surface_temp_C"];
end

function stageSelectedSnowInventory(data_root, case_id)
   %STAGESELECTEDSNOWINVENTORY Copy one case into an isolated manifest tree.

   source_family = fullfile(icemodel.internal.fullpath('data'), ...
      'eval', 'esm_snowmip');
   target_family = fullfile(data_root, 'eval', 'esm_snowmip');
   icemodel.helpers.ensureDirExists(target_family)
   icemodel.helpers.ensureDirExists(fullfile(data_root, 'input'))

   % Filter the family manifest before copying its one selected case folder.
   manifest = jsondecode(fileread(fullfile(source_family, 'manifest.json')));
   entries = reshape(manifest.cases, [], 1);
   keep = string({entries.case_id}) == case_id;
   manifest.cases = entries(keep);
   assert(isscalar(manifest.cases), 'selected snow fixture case is missing')
   copyfile(fullfile(source_family, case_id), ...
      fullfile(target_family, case_id))
   writeJson(fullfile(target_family, 'manifest.json'), manifest)
end

function candidate = selectedTargetCandidate(case_manifest)
   %SELECTEDTARGETCANDIDATE Reuse the isolated observation bundle as candidate.

   candidate = icemodel.verification.helpers.loadArtifact( ...
      case_manifest.evaluation_path, "targets");
end

function writeJson(filename, value)
   %WRITEJSON Write one temporary manifest with deterministic indentation.

   text = jsonencode(value, PrettyPrint=true);
   fid = fopen(filename, 'w');
   cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '%s\n', text);
   clear cleanup
end
