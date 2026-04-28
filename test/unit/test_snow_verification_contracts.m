function tests = test_snow_verification_contracts
   %TEST_SNOW_VERIFICATION_CONTRACTS Verify the snow-validation framework.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install the canonical test/demo config for the verification tests.

   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.tmpdir = tempname(fullfile( ...
      icemodel.getpath('test'), 'artifacts', 'tmp'));
   icemodel.helpers.ensureDirExists(testCase.TestData.tmpdir);
end

function teardown(testCase)
   % Remove temporary artifacts and restore the caller's config.

   if exist(testCase.TestData.tmpdir, 'dir') == 7
      rmdir(testCase.TestData.tmpdir, 's')
   end
   clear testCase.TestData.cleanup
end

function test_listcases_returns_expected_ids(testCase)
   % LISTCASES should expose the curated smoke cases from committed demo/data.

   cases = icemodel.verification.listcases();
   ids = [cases.case_id];

   testCase.verifyTrue(all(ismember(expectedCaseIds(), ids)));
end

function test_verification_namelists_expose_curated_selectors(testCase)
   % The verification namespace should publish one canonical selector catalog.

   testCase.verifyEqual(icemodel.verification.namelists.datasetfamily(), ...
      expectedDatasetFamilies());
   testCase.verifyEqual(icemodel.verification.namelists.casetype(), ...
      expectedCaseTypes());
   testCase.verifyEqual(icemodel.verification.namelists.tier(), ...
      expectedTiers());
   testCase.verifyEqual(icemodel.verification.namelists.caseid(), ...
      expectedCaseIds());
   testCase.verifyEqual( ...
      icemodel.verification.namelists.caseid("esm_snowmip"), ...
      expectedEsmCaseIds());
end

function test_manifest_schema_helpers_are_setup_owned(testCase)
   % Manifest construction schemas should live with setup/update tooling.

   family_names = icemodel.verification.setup.familyManifestFieldNames();
   case_names = icemodel.verification.setup.caseManifestFieldNames();

   testCase.verifyEqual(family_names(:), expectedFamilyManifestFields());
   testCase.verifyTrue(ismember("case_id", case_names));
   testCase.verifyTrue(ismember("case_type", case_names));
end

function test_loadmanifest_resolves_demo_data_paths(testCase)
   % LOADMANIFEST should resolve case files under the demo/data tree.

   manifest = icemodel.verification.loadmanifest("cdp");

   testCase.verifyEqual(manifest.dataset_family, "esm_snowmip");
   testCase.verifyEqual(manifest.case_type, "esm_site");
   testCase.verifyTrue(contains(manifest.forcing_path, ...
      fullfile("demo", "data", "eval", "snow", "esm_snowmip", "cdp")));
   testCase.verifyTrue(exist(manifest.evaluation_path, 'file') == 2);
end

function test_comparecase_smoke_reference_writes_metrics(testCase)
   % COMPARECASE should run from the committed smoke reference artifacts.

   result = icemodel.verification.comparecase("cdp", ...
      "artifact_dir", testCase.TestData.tmpdir, ...
      "make_plot", false);

   testCase.verifyEqual(unique(result.metrics.status), "ok");
   testCase.verifyTrue(exist(result.metrics_path, 'file') == 2);
   testCase.verifyGreaterThan(height(result.metrics), 0);
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

function test_plotcase_writes_figure_without_candidate(testCase)
   % Users should be able to visualize staged targets directly.

   outfile = fullfile(testCase.TestData.tmpdir, 'cdp_plot.png');
   f = icemodel.verification.plotcase("cdp", ...
      "visible", "off", ...
      "output_file", outfile);

   testCase.verifyTrue(exist(outfile, 'file') == 2);
   close(f)
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

function test_plot_timeseries_shows_sparse_points_with_markers(testCase)
   % Sparse observation series should remain visible on a dense time axis.

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
   markers = string(get(lines, 'Marker'));

   testCase.verifyTrue(any(markers ~= "none"));
end

function ids = expectedCaseIds()
   %EXPECTEDCASEIDS Canonical curated smoke-case catalog.

   ids = ["cdp"; "wfj"; "colbeck1976"];
end

function ids = expectedEsmCaseIds()
   %EXPECTEDESMCASEIDS ESM-SnowMIP cases in the committed smoke subset.

   ids = ["cdp"; "wfj"];
end

function names = expectedDatasetFamilies()
   %EXPECTEDDATASETFAMILIES Dataset families currently staged for verification.

   names = ["esm_snowmip"; "laugh_tests"];
end

function names = expectedCaseTypes()
   %EXPECTEDCASETYPES Verification case types used by manifests and runners.

   names = ["esm_site"; "synthetic_process"];
end

function tiers = expectedTiers()
   %EXPECTEDTIERS Public runner tiers supported by the smoke harness.

   tiers = ["smoke"; "full"];
end

function fields = expectedFamilyManifestFields()
   %EXPECTEDFAMILYMANIFESTFIELDS Required top-level family-manifest fields.

   fields = {'dataset_family'; 'source_doi'; 'source_url'; ...
      'source_version'; 'retrieval_date'; 'cases'};
end

function ids = expectedColbeckExperimentIds()
   %EXPECTEDCOLBECKEXPERIMENTIDS Colbeck experiments staged from Laugh-Tests.

   ids = ["exp1", "exp2", "exp3"];
end

function names = expectedCoreSiteVariables()
   %EXPECTEDCORESITEVARIABLES Minimal site variables derived from model output.

   names = ["snow_depth_m"; "swe_kg_m2"; "surface_temp_C"];
end
