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

   testCase.verifyTrue(all(ismember(["cdp", "wfj", "colbeck1976"], ids)));
end

function test_verification_namelists_expose_curated_selectors(testCase)
   % The verification namespace should publish one canonical selector catalog.

   testCase.verifyEqual(icemodel.verification.namelists.datasetfamily(), ...
      ["esm_snowmip"; "laugh_tests"]);
   testCase.verifyEqual(icemodel.verification.namelists.casetype(), ...
      ["esm_site"; "synthetic_process"]);
   testCase.verifyEqual(icemodel.verification.namelists.tier(), ...
      ["smoke"; "full"]);
   testCase.verifyEqual(icemodel.verification.namelists.caseid(), ...
      ["cdp"; "wfj"; "colbeck1976"]);
   testCase.verifyEqual( ...
      icemodel.verification.namelists.caseid("esm_snowmip"), ...
      ["cdp"; "wfj"]);
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
      ["exp1", "exp2", "exp3"]);
   testCase.verifyTrue(all(result.metrics.n > 0));
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

function test_plot_timeseries_shows_sparse_points_with_markers(testCase)
   % Sparse observation series should remain visible on a dense time axis.

   time = datetime(2000, 1, 1, 0, 0, 0) + hours(0:9);
   values = nan(10, 1);
   values([2 6 10]) = [1 2 3];

   f = figure('Visible', 'off');
   cleaner = onCleanup(@() close(f));
   ax = axes(f);
   icemodel.plot.timeseries(time(:), values, axes=ax);

   lines = findall(ax, 'Type', 'line');
   markers = string(get(lines, 'Marker'));

   testCase.verifyTrue(any(markers ~= "none"));
end
