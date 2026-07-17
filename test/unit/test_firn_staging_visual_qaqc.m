function tests = test_firn_staging_visual_qaqc
   %TEST_FIRN_STAGING_VISUAL_QAQC Verify staged firn artifact plotting helpers.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Build a tiny staged tree that matches the real eval/input contract without
   % touching large source fixtures or external RCM products.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   root = tempname;
   mkdir(root);
   testCase.addTeardown(@() rmdir(root, 's'));
   testCase.TestData.root = root;
   writeTinyFirnTree(root);
   writeTinyEsmTree(root);
   writeTinyLaughTree(root);
   writeEmptyResearchTree(root);
end

function teardown(testCase)
   % Release the bootstrap cleanup handle.
   testCase.TestData.cleanup = [];
end

function test_plot_verification_artifacts_summarizes_staged_tree(testCase)
   % plotVerificationArtifacts should discover observations, met files, and
   % userdata from a manifest rooted at output_root without committed fixtures.
   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);

   testCase.verifyGreaterThanOrEqual(height(returned), 4);
   testCase.verifyTrue(all(returned.dataset_family == "promice"));
   testCase.verifyTrue(all(returned.case_id == "kanm"));
   testCase.verifyTrue(all(returned.target_format == "timeseries"));
   testCase.verifyTrue(any(contains(returned.target_variables, "tice10m")));
   testCase.verifyTrue(any(contains(returned.met_sources, "promice")));
   testCase.verifyTrue(any(contains(returned.userdata_sources, "promice")));
   testCase.verifyTrue(any(contains(returned.met_sources, "mar3.11")));
   testCase.verifyTrue(any(contains(returned.userdata_sources, "mar3.11")));

   plotted = strjoin(returned.plotted_variables, ",");
   expected = ["ablation", "tice10m", "tair", "swd", "lwd", "ppt"];
   for v = expected
      testCase.verifyTrue(contains(plotted, v), ...
         sprintf('visual QA did not plot fixture variable %s', v));
   end

   accounted = strjoin(returned.plotted_variables + "," ...
      + returned.unplotted_variables, ",");
   for v = expected
      testCase.verifyTrue(contains(accounted, v), ...
         sprintf('visual QA did not account for fixture variable %s', v));
   end
end

function test_plot_verification_artifacts_pairs_explicit_eval_root(testCase)
   % An explicit eval root should imply the sibling input root unless overridden.
   returned = icemodel.verification.plotVerificationArtifacts( ...
      evaluation_data_root=fullfile(testCase.TestData.root, 'eval'), ...
      dataset_family="promice", save_figs=false);

   testCase.verifyTrue(any(contains(returned.met_sources, "promice")));
   testCase.verifyTrue(any(contains(returned.userdata_sources, "promice")));
   testCase.verifyTrue(any(contains(returned.plotted_variables, "tair")));
end

function test_saved_observations_are_stamped(testCase)
   % Observation bundles should be self-describing on disk, not only after the
   % plotting loader stamps them in memory.

   loaded = load(fullfile(testCase.TestData.root, 'eval', 'promice', ...
      'kanm', 'observations.mat'), 'targets');

   testCase.verifyEqual( ...
      loaded.targets.data.Properties.DimensionNames{1}, 'Time');
   testCase.verifyEqual(string(loaded.targets.data.Properties.VariableUnits), ...
      ["m", "K", "K"]);
end

function test_plot_verification_artifacts_discovers_legacy_esm_met(testCase)
   % ESM-SnowMIP manifests predate colocation records, so the visual QA plotter
   % discovers their staged met files by the nested filename convention.

   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="esm_snowmip", save_figs=false);

   testCase.verifyTrue(any(returned.met_sources == "esm_snowmip"));
   testCase.verifyTrue(any(contains(returned.plotted_variables, "tair")));
end

function test_plot_verification_artifacts_discovers_flat_legacy_esm_met(testCase)
   % Older ESM-SnowMIP staged trees used input/met directly, without a source
   % subfolder, and should remain discoverable for visual QA.
   nested = fullfile(testCase.TestData.root, 'input', 'met', 'esm_snowmip', ...
      'met_cdp_esm_snowmip_20120101_20120102_1hr.mat');
   delete(nested)

   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="esm_snowmip", save_figs=false);

   testCase.verifyTrue(any(returned.met_sources == "esm_snowmip"));
   testCase.verifyTrue(any(contains(returned.plotted_variables, "tair")));
end

function test_plot_verification_artifacts_groups_esm_native_aliases(testCase)
   % Native ESM snow-depth/store and temperature aliases belong in the existing
   % domain groups instead of creating an oversized other-variables figure.
   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="esm_snowmip", save_figs=false);

   surface = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "surface_height_depth"));
   temperature = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "temperature_state"));
   expected_surface = ["snow_depth_m", "swe_kg_m2", ...
      "snd_auto_m", "snd_man_m"];
   expected_temperature = ["surface_temp_C", "soil_temp_1_C", ...
      "soil_temp_2_C", "soil_temp_3_C", "soil_temp_4_C"];
   testCase.verifyTrue(all(ismember(expected_surface, surface)));
   testCase.verifyTrue(all(ismember(expected_temperature, temperature)));
   testCase.verifyFalse(any(returned.figure_group == "other_variables"));

   % The data-driven group registry must uphold the seven-panel ceiling even
   % when a site exposes every representative native alias.
   for k = 1:height(returned)
      variables = listedVariables(returned.plotted_variables(k));
      testCase.verifyLessThanOrEqual(numel(variables), 7, ...
         sprintf('%s preview has too many stacked panels', ...
         returned.figure_group(k)));
   end
end

function test_sparse_evaluation_state_uses_available_daily_samples(testCase)
   % Sparse evaluation values on a complete source grid remain visible while
   % dense forcing keeps its stricter complete-day contract.
   observation_file = fullfile(testCase.TestData.root, 'eval', ...
      'esm_snowmip', 'cdp', 'observations.mat');
   loaded = load(observation_file, 'targets');
   targets = loaded.targets;
   sparse_variables = ["snow_depth_m", "swe_kg_m2", ...
      "snd_auto_m", "snd_man_m"];
   keep = [6, 30];
   for varname = sparse_variables
      values = nan(height(targets.data), 1);
      values(keep) = targets.data.(varname)(keep);
      targets.data.(varname) = values;
   end
   save(observation_file, 'targets')

   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="esm_snowmip", save_figs=false);

   surface = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "surface_height_depth"));
   testCase.verifyTrue(all(ismember(sparse_variables, surface)));
   testCase.verifyFalse(any(contains(returned.plotted_variables( ...
      returned.figure_group == "other_variables"), sparse_variables)));
end

function test_plot_firn_artifacts_writes_figure(testCase)
   % The compatibility wrapper must write figures and forward case cleanup.
   figure_root = fullfile(testCase.TestData.root, 'figures');
   obsolete = fullfile(figure_root, 'promice', 'kanm', 'energy_fluxes.png');
   icemodel.helpers.ensureDirExists(fileparts(obsolete));
   fclose(fopen(obsolete, 'w'));

   returned = icemodel.verification.plotFirnArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", figure_root=figure_root, save_figs=true, ...
      overwrite=true);

   testCase.verifyTrue(all(isfile(returned.figure_file)));
   testCase.verifyFalse(isfile(obsolete));
end

function test_plot_firn_artifacts_all_excludes_snow_families(testCase)
   % The firn-specific wrapper's default "all" selector should not forward the
   % snow-only family set used by the neutral verification plotter.
   returned = icemodel.verification.plotFirnArtifacts( ...
      output_root=string(testCase.TestData.root), save_figs=false);

   testCase.verifyTrue(all(ismember(returned.dataset_family, ...
      icemodel.verification.namelists.firndatasetfamily())));
   testCase.verifyFalse(any(returned.dataset_family == "esm_snowmip"));
end

function test_plot_firn_artifacts_rejects_snow_family(testCase)
   % The firn wrapper should not accept explicit snow-only family selectors.
   testCase.verifyError(@() icemodel.verification.plotFirnArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="esm_snowmip", save_figs=false), ...
      'icemodel:validators:mustBeFirnDatasetFamilySelection');
end

function test_plot_verification_artifacts_accepts_date_window(testCase)
   % The visual QA entry point must support one-window previews of a longer
   % staged artifact tree without requiring the artifacts to be restaged.

   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false, ...
      startdate=datetime(2012, 1, 2, 0, 0, 0, 'TimeZone', 'UTC'), ...
      enddate=datetime(2012, 1, 2, 23, 0, 0, 'TimeZone', 'UTC'));

   testCase.verifyGreaterThan(height(returned), 0);
end

function test_plot_verification_artifacts_treats_nat_as_unset(testCase)
   % NaT bounds come from runtime option defaults and should mean unbounded plots.
   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false, startdate=NaT, enddate=NaT);

   testCase.verifyGreaterThan(height(returned), 0);
   testCase.verifyTrue(any(contains(returned.plotted_variables, "tair")));
end

function test_plot_verification_artifacts_handles_empty_case(testCase)
   % Missing artifacts should produce a summary row instead of failing during
   % variable grouping.
   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="research_site", save_figs=false);

   testCase.verifyEqual(returned.figure_group, "no_plottable_data");
   testCase.verifyEqual(returned.case_id, "empty");
end

function test_thermistor_string_plotted_once(testCase)
   % The thermistor-string summary is one group and includes its 10 m primary.
   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);
   has_tice1 = arrayfun(@hasListedTice1, returned.plotted_variables);

   testCase.verifyEqual(sum(has_tice1), 1);
   testCase.verifyEqual(returned.figure_group(has_tice1), ...
      "subsurface_temperature_string");
   thermistor = returned(returned.figure_group == ...
      "subsurface_temperature_string", :);
   testCase.verifyNumElements(thermistor.figure_group, 1);
   testCase.verifyTrue(ismember("tice10m", ...
      listedVariables(thermistor.plotted_variables)));
end

function test_thermistor_string_uses_one_axes_and_emphasizes_primary(testCase)
   % Reuse the accepted PROMICE QA design: one axes, thin ticeN, thick tice10m.
   observation_file = fullfile(testCase.TestData.root, 'eval', 'promice', ...
      'kanm', 'observations.mat');
   loaded = load(observation_file, 'targets');
   targets = loaded.targets;
   for sensor = 2:8
      targets.data.("tice" + sensor) = targets.data.tice1 + sensor / 100;
   end
   save(observation_file, 'targets')

   % Retain the normally closed figure only long enough to inspect its objects.
   old_close = get(groot, 'DefaultFigureCloseRequestFcn');
   restore_close = onCleanup(@() set(groot, ...
      'DefaultFigureCloseRequestFcn', old_close));
   set(groot, 'DefaultFigureCloseRequestFcn', @(~, ~) []);
   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);
   clear restore_close

   figures = findall(groot, 'Type', 'Figure');
   names = string(get(figures, 'Name'));
   figures = figures(startsWith(names, "promice/kanm"));
   names = string(get(figures, 'Name'));
   testCase.addTeardown(@() delete(figures(isgraphics(figures))));
   thermistor_fig = figures(contains(names, "thermistor string"));
   testCase.verifyNumElements(thermistor_fig, 1);
   ax = findall(thermistor_fig, 'Type', 'Axes');
   testCase.verifyNumElements(ax, 1);

   % Every source thermistor shares that axes and the primary is drawn last,
   % thicker, and without introducing any marker-only graphics objects.
   lines = findall(ax, 'Type', 'Line');
   labels = string(get(lines, 'DisplayName'));
   named = strlength(labels) > 0;
   lines = lines(named);
   labels = labels(named);
   expected = ["tice" + string((1:8)'); "tice10m (PRIMARY)"];
   testCase.verifyEqual(sort(labels), sort(expected));
   primary = lines(labels == "tice10m (PRIMARY)");
   diagnostic = lines(startsWith(labels, "tice") ...
      & labels ~= "tice10m (PRIMARY)");
   testCase.verifyNumElements(primary, 1);
   testCase.verifyEqual(primary.LineWidth, 2.0);
   testCase.verifyTrue(all(abs([diagnostic.LineWidth] - 0.6) < 1e-12));
   testCase.verifyTrue(all(string(get(lines, 'Marker')) == "none"));
   testCase.verifyEmpty(findall(ax, 'Type', 'Line', 'LineStyle', 'none'));

   row = returned(returned.figure_group == ...
      "subsurface_temperature_string", :);
   testCase.verifyNumElements(row.figure_group, 1);
   testCase.verifyEqual(sort(listedVariables(row.plotted_variables)), ...
      sort(["tice" + string((1:8)'); "tice10m"]));
end

function test_surface_albedo_uses_shared_shortwave_weighted_daily_mean(testCase)
   % Dedicated and met panels share source-aware radiometer/model reductions.
   observation_file = fullfile(testCase.TestData.root, 'eval', 'promice', ...
      'kanm', 'observations.mat');
   loaded = load(observation_file, 'targets');
   targets = loaded.targets;
   targets.data.albedo = nan(height(targets.data), 1);
   targets.data.swd = zeros(height(targets.data), 1);
   first = 5:10;
   second = 29:34;
   targets.data.albedo(first) = [0.4; 0.4; 0.4; 0.9; 0.9; 0.9];
   targets.data.swd(first) = [10; 10; 10; 500; 500; 500];
   targets.data.albedo(second) = [0.5; 0.5; 0.5; 0.8; 0.8; 0.8];
   targets.data.swd(second) = [20; 20; 20; 400; 400; 400];
   save(observation_file, 'targets')

   met_file = fullfile(testCase.TestData.root, 'input', 'met', 'promice', ...
      'met_kanm_promice_20120101_20120102_1hr.mat');
   loaded = load(met_file, 'met');
   met = loaded.met;
   met.albedo = targets.data.albedo;
   met.swd = targets.data.swd;
   save(met_file, 'met')

   % A native MAR state remains defined through zero-shortwave periods. This
   % distinguishes its arithmetic daily state from radiometer energy weighting.
   mar_file = fullfile(testCase.TestData.root, 'input', 'userdata', ...
      'mar3.11', 'kanm_mar3.11_20120101_20120102.mat');
   loaded = load(mar_file, 'Data');
   Data = loaded.Data;
   model_day = days(dateshift(Data.Time, 'start', 'day') ...
      - dateshift(Data.Time(1), 'start', 'day')) + 1;
   Data.albedo = 0.7 + 0.1 .* model_day;
   Data.swd(:) = 0;
   save(mar_file, 'Data')

   mar_met_file = fullfile(testCase.TestData.root, 'input', 'met', ...
      'mar3.11', 'met_kanm_mar3.11_20120101_20120102_1hr.mat');
   loaded = load(mar_met_file, 'met');
   met = loaded.met;
   model_day = days(dateshift(met.Time, 'start', 'day') ...
      - dateshift(met.Time(1), 'start', 'day')) + 1;
   met.albedo = 0.7 + 0.1 .* model_day;
   met.swd(:) = 0;
   save(mar_met_file, 'met')

   % Keep the rendered albedo figure so the legend-to-line contract is testable.
   old_close = get(groot, 'DefaultFigureCloseRequestFcn');
   restore_close = onCleanup(@() set(groot, ...
      'DefaultFigureCloseRequestFcn', old_close));
   set(groot, 'DefaultFigureCloseRequestFcn', @(~, ~) []);
   icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);
   clear restore_close

   figures = findall(groot, 'Type', 'Figure');
   names = string(get(figures, 'Name'));
   figures = figures(startsWith(names, "promice/kanm"));
   names = string(get(figures, 'Name'));
   testCase.addTeardown(@() delete(figures(isgraphics(figures))));
   albedo_fig = figures(contains(names, "surface albedo"));
   testCase.verifyNumElements(albedo_fig, 1);
   lines = findall(albedo_fig, 'Type', 'Line');
   labels = string(get(lines, 'DisplayName'));
   named = strlength(labels) > 0;
   lines = lines(named);
   labels = labels(named);
   promice = lines(labels == "PROMICE observations");

   expected = [(0.4 * 10 + 0.9 * 500) / 510; ...
      (0.5 * 20 + 0.8 * 400) / 420];
   model_expected = [0.8; 0.9];
   testCase.verifyNumElements(promice, 1);
   testCase.verifyEqual(promice.YData(:), expected, 'AbsTol', 1e-12);
   testCase.verifyEqual(numel(promice.XData), 2);
   testCase.verifyTrue(all(arrayfun(@(line) any(isfinite(line.YData)), lines)));
   testCase.verifyTrue(all(string(get(lines, 'Marker')) == "none"));
   mar = lines(labels == "MAR 3.11");
   testCase.verifyNumElements(mar, 1);
   testCase.verifyEqual(mar.YData(:), model_expected, 'AbsTol', 1e-12);

   met_fig = figures(contains(names, "met forcing"));
   testCase.verifyNumElements(met_fig, 1);
   met_axes = findall(met_fig, 'Type', 'Axes');
   albedo_ax = met_axes(string(get(met_axes, 'Tag')) == "surface albedo");
   met_lines = findall(albedo_ax, 'Type', 'Line');
   met_labels = string(get(met_lines, 'DisplayName'));
   met_promice = met_lines(met_labels == "PROMICE forcing");
   met_mar = met_lines(met_labels == "MAR 3.11 forcing");
   testCase.verifyNumElements(met_promice, 1);
   testCase.verifyNumElements(met_mar, 1);
   testCase.verifyEqual(met_promice.YData(:), expected, 'AbsTol', 1e-12);
   testCase.verifyEqual(met_mar.YData(:), model_expected, 'AbsTol', 1e-12);
end

function test_radiation_native_partial_hours_break_gaps(testCase)
   % Partial PROMICE radiation stays hourly and an omitted interval is not joined.
   observation_file = fullfile(testCase.TestData.root, 'eval', 'promice', ...
      'kanm', 'observations.mat');
   loaded = load(observation_file, 'targets');
   targets = loaded.targets;
   n_rows = height(targets.data);
   swd = 200 + (1:n_rows)';
   swd([1:6 25:30]) = NaN;
   targets.data.swd = swd;
   targets.data.swu = swd .* 0.6;
   targets.data.lwd = 250 + zeros(n_rows, 1);
   targets.data.lwu = 300 + zeros(n_rows, 1);
   targets.data(10:15, :) = [];
   native_time = targets.data.Time;
   save(observation_file, 'targets')

   % Retain the radiation figure so its native source line can be inspected.
   old_close = get(groot, 'DefaultFigureCloseRequestFcn');
   restore_close = onCleanup(@() set(groot, ...
      'DefaultFigureCloseRequestFcn', old_close));
   set(groot, 'DefaultFigureCloseRequestFcn', @(~, ~) []);
   icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);
   clear restore_close

   figures = findall(groot, 'Type', 'Figure');
   names = string(get(figures, 'Name'));
   figures = figures(startsWith(names, "promice/kanm"));
   names = string(get(figures, 'Name'));
   testCase.addTeardown(@() delete(figures(isgraphics(figures))));
   radiation_fig = figures(contains(names, "radiation fluxes"));
   testCase.verifyNumElements(radiation_fig, 1);
   axes_handles = findall(radiation_fig, 'Type', 'Axes');
   tags = string(get(axes_handles, 'Tag'));
   swd_ax = axes_handles(tags == "swd");
   testCase.verifyNumElements(swd_ax, 1);
   lines = findall(swd_ax, 'Type', 'Line');
   labels = string(get(lines, 'DisplayName'));
   promice = lines(labels == "PROMICE observations");

   % One midpoint NaN is inserted for the omitted six-hour interval; all actual
   % native timestamps and finite samples remain, with no marker-only overlay.
   testCase.verifyNumElements(promice, 1);
   testCase.verifyEqual(numel(promice.XData), numel(native_time) + 1);
   testCase.verifyTrue(all(ismember(native_time, promice.XData)));
   inserted = ~ismember(promice.XData, native_time);
   testCase.verifyEqual(nnz(inserted), 1);
   testCase.verifyTrue(isnan(promice.YData(inserted)));
   testCase.verifyTrue(any(isfinite(promice.YData)));
   testCase.verifyEqual(string(promice.Marker), "none");
   testCase.verifyEmpty(findall(swd_ax, 'Type', 'Line', 'LineStyle', 'none'));
end

function test_snow_depth_aliases_share_one_canonical_panel(testCase)
   % MAR snowd and PROMICE snow_depth overlay without renaming staged artifacts.
   observation_file = fullfile(testCase.TestData.root, 'eval', 'promice', ...
      'kanm', 'observations.mat');
   loaded = load(observation_file, 'targets');
   targets = loaded.targets;
   targets.data.snow_depth = linspace(0.1, 0.3, height(targets.data))';
   targets.data.Properties.VariableUnits{end} = 'm';
   targets.data.surface_height = linspace(-0.2, 0.1, height(targets.data))';
   targets.data.Properties.VariableUnits{end} = 'm';
   save(observation_file, 'targets')

   mar_file = fullfile(testCase.TestData.root, 'input', 'userdata', ...
      'mar3.11', 'kanm_mar3.11_20120101_20120102.mat');
   loaded = load(mar_file, 'Data');
   Data = loaded.Data;
   Data.snowd = linspace(0.2, 0.4, height(Data))';
   Data.Properties.VariableUnits{end} = 'm';
   save(mar_file, 'Data')

   % Retain the surface figure long enough to prove one canonical plot panel.
   old_close = get(groot, 'DefaultFigureCloseRequestFcn');
   restore_close = onCleanup(@() set(groot, ...
      'DefaultFigureCloseRequestFcn', old_close));
   set(groot, 'DefaultFigureCloseRequestFcn', @(~, ~) []);
   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);
   clear restore_close

   figures = findall(groot, 'Type', 'Figure');
   names = string(get(figures, 'Name'));
   figures = figures(startsWith(names, "promice/kanm"));
   names = string(get(figures, 'Name'));
   testCase.addTeardown(@() delete(figures(isgraphics(figures))));
   surface_fig = figures(contains(names, ...
      "surface height change, snow depth, and stores"));
   testCase.verifyNumElements(surface_fig, 1);
   axes_handles = findall(surface_fig, 'Type', 'Axes');
   tags = string(get(axes_handles, 'Tag'));
   snow_ax = axes_handles(tags == "snow_depth");
   surface_ax = axes_handles(tags == "surface_height");
   testCase.verifyNumElements(snow_ax, 1);
   testCase.verifyNumElements(surface_ax, 1);
   testCase.verifyEqual(string(snow_ax.Title.String), "");
   testCase.verifyEqual(string(surface_ax.Title.String), ...
      "surface height change (+up; not snow depth)");
   testCase.verifyFalse(any(tags == "snowd"));
   labels = string(get(findall(snow_ax, 'Type', 'Line'), 'DisplayName'));
   testCase.verifyTrue(all(ismember(["PROMICE observations", "MAR 3.11"], ...
      labels)));

   surface = returned(returned.figure_group == "surface_height_depth", :);
   variables = listedVariables(surface.plotted_variables);
   testCase.verifyEqual(nnz(variables == "snow_depth"), 1);
   testCase.verifyEqual(nnz(variables == "surface_height"), 1);
   testCase.verifyFalse(any(variables == "snowd"));

   % The plotting alias is in-memory only: staged MAR retains its native name.
   staged = load(mar_file, 'Data');
   staged_names = string(staged.Data.Properties.VariableNames);
   testCase.verifyTrue(ismember("snowd", staged_names));
   testCase.verifyFalse(ismember("snow_depth", staged_names));
end

function test_plot_verification_artifacts_uses_canonical_shared_groups(testCase)
   % The shared renderer should keep model forcing, balance terms, and related
   % components in stable domain-oriented groups without diagnostic spillover.
   observation_file = fullfile(testCase.TestData.root, 'eval', 'promice', ...
      'kanm', 'observations.mat');
   loaded = load(observation_file, 'targets');
   targets = loaded.targets;
   for sensor = 2:11
      targets.data.("tice" + sensor) = targets.data.tice1 + sensor / 100;
   end
   save(observation_file, 'targets')

   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);

   % Every shared group remains independently selectable by report tooling.
   expected = ["met_forcing", "radiation_fluxes", "energy_balance", ...
      "turbulent_fluxes", "mass_balance_terms", ...
      "mass_balance_components", "mar_combined_mass_diagnostics", ...
      "surface_height_depth", "surface_albedo", ...
      "meteorological_diagnostics", ...
      "quality_flags"];
   testCase.verifyTrue(all(ismember(expected, returned.figure_group)));
   for group = expected(expected ~= "met_forcing")
      row = returned(returned.figure_group == group, :);
      variables = listedVariables(row.plotted_variables);
      testCase.verifyLessThanOrEqual(numel(variables), 7, ...
         sprintf('%s preview has too many stacked panels', group));
   end

   % Group membership encodes the plotting contract and protects against the
   % former net-radiation and diagnostic-divergence catch-all panels.
   radiation = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "radiation_fluxes"));
   energy = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "energy_balance"));
   turbulent = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "turbulent_fluxes"));
   mass = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "mass_balance_terms"));
   components = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "mass_balance_components"));
   mar_combined = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "mar_combined_mass_diagnostics"));
   temperature = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "temperature_state"));
   surface = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "surface_height_depth"));
   diagnostics = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "meteorological_diagnostics"));
   flags = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "quality_flags"));
   thermistor_rows = startsWith(returned.figure_group, ...
      "subsurface_temperature_string");
   thermistors = listedVariables(strjoin( ...
      returned.plotted_variables(thermistor_rows), ','));
   testCase.verifyEqual(radiation, ["swd"; "swu"; "lwd"; "lwu"]);
   testCase.verifyFalse(ismember("netr", radiation));
   testCase.verifyEqual(energy, ["swn"; "lwn"; "thf"; "netr"]);
   testCase.verifyEqual(turbulent, ["shf"; "lhf"]);
   testCase.verifyFalse(ismember("sndiv", mass));
   testCase.verifyTrue(ismember("subl", mass));
   testCase.verifyTrue(all(ismember(["melt", "meltin"], components)));
   testCase.verifyEqual(mar_combined, ...
      ["subl_evap"; "refreeze_deposition"]);
   testCase.verifyTrue(ismember("tice10m", temperature));
   testCase.verifyEqual(surface, "ablation");
   testCase.verifyEqual(diagnostics, "wdir");
   testCase.verifyEqual(flags, "tice10m_qc_flag");
   testCase.verifyEqual(sort(thermistors), ...
      sort(["tice" + string((1:11)'); "tice10m"]));
   testCase.verifyEqual(nnz(thermistor_rows), 1);

   % All eight canonical forcing channels are represented in at most five
   % combined panels; the summary records channels rather than panel count.
   met_variables = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "met_forcing"));
   testCase.verifyEqual(sort(met_variables), sort(["swd"; "lwd"; "tair"; ...
      "rh"; "wspd"; "albedo"; "ppt"; "psfc"]));
   all_plotted = listedVariables(strjoin(returned.plotted_variables, ','));
   all_unplotted = listedVariables(strjoin(returned.unplotted_variables, ','));
   testCase.verifyFalse(ismember("sndiv", all_plotted));
   testCase.verifyTrue(ismember("sndiv", all_unplotted));
   testCase.verifyFalse(ismember("step_magnitude", all_plotted));
   testCase.verifyTrue(ismember("step_magnitude", all_unplotted));
   testCase.verifyFalse(ismember("tice10m_source", all_plotted));
   testCase.verifyTrue(ismember("tice10m_source", all_unplotted));
   testCase.verifyTrue(ismember("dtice1", all_unplotted));
   testCase.verifyTrue(ismember("modis", all_plotted));
   testCase.verifyFalse(ismember("modis", all_unplotted));
   modis_rows = contains(returned.plotted_variables, "modis");
   testCase.verifyEqual(nnz(modis_rows), 1);
   testCase.verifyEqual(returned.figure_group(modis_rows), "surface_albedo");
   testCase.verifyFalse(any(returned.figure_group == "other_variables"));
   testCase.verifyFalse(any(startsWith( ...
      returned.figure_group, "measurement_geometry")));
end

function test_plot_derives_mar_energy_balance_from_native_components(testCase)
   % MAR userdata stages native forcing components; the renderer should expose
   % processmet-consistent balances without mutating the staged artifact.
   mar_file = fullfile(testCase.TestData.root, 'input', 'userdata', ...
      'mar3.11', 'kanm_mar3.11_20120101_20120102.mat');
   loaded = load(mar_file, 'Data');
   Data = removevars(loaded.Data, {'swu', 'lwu', 'swn', 'lwn'});
   Data.albedo = repmat(0.8, height(Data), 1);
   Data.tsfc = repmat(263.15, height(Data), 1);
   Data = icemodel.forcing.helpers.stampMetadata(Data, strict=false);
   save(mar_file, 'Data')

   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);

   energy = listedVariables(returned.plotted_variables( ...
      returned.figure_group == "energy_balance"));
   testCase.verifyEqual(energy, ["swn"; "lwn"; "thf"; "netr"]);
   staged = load(mar_file, 'Data');
   staged_names = string(staged.Data.Properties.VariableNames);
   testCase.verifyFalse(any(ismember(["swn", "lwn", "netr"], staged_names)));
end

function test_profiles_exclude_support_metadata_and_blank_rows(testCase)
   % Profile support columns never become science panels, and a profile with no
   % finite science pair does not create a blank summary row.
   observation_file = fullfile(testCase.TestData.root, 'eval', 'promice', ...
      'kanm', 'observations.mat');
   loaded = load(observation_file, 'targets');
   targets = loaded.targets;
   profile = table([0.1; 0.2], [NaN; NaN], [-8; -9], [11; 12], ...
      'VariableNames', {'depth', 'error', 'subsurface_temperature', ...
      'measurement_id'});
   profile.Properties.VariableUnits = {'m', 'degC', 'degC', '1'};
   targets.subsurface_temperature = profile;
   save(observation_file, 'targets')

   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);
   profile_row = returned(returned.figure_group == "profiles", :);
   testCase.verifyEqual(height(profile_row), 1);
   testCase.verifyEqual(listedVariables(profile_row.plotted_variables), ...
      "subsurface_temperature");
   testCase.verifyFalse(any(contains(returned.plotted_variables, ...
      ["error", "measurement_id"])));

   % Removing the only finite profile values also removes the profile row rather
   % than promising a plotted variable backed by a white export.
   targets.subsurface_temperature.subsurface_temperature(:) = NaN;
   save(observation_file, 'targets')
   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);
   testCase.verifyFalse(any(returned.figure_group == "profiles"));
   testCase.verifyFalse(any(contains(returned.plotted_variables, ...
      ["error", "measurement_id"])));
end

function test_named_profiles_and_interval_totals_render_readably(testCase)
   % Named profiles get distinct colors, and interval-total observations remain
   % disjoint from both unobserved years and incremental daily model SMB.
   observation_file = fullfile(testCase.TestData.root, 'eval', 'promice', ...
      'kanm', 'observations.mat');
   loaded = load(observation_file, 'targets');
   targets = loaded.targets;
   profile = table([1; 2; 1; 2], [-10; -11; -12; -13], ...
      ["core A"; "core A"; "core B"; "core B"], ...
      'VariableNames', {'depth', 'subsurface_temperature', 'name'});
   profile.Properties.VariableUnits = {'m', 'degC', ''};
   targets.subsurface_temperature = profile;
   targets.data.snowf_subl = repmat(0.001, height(targets.data), 1);
   targets.data.Properties.VariableUnits{end} = 'mWE/h';
   targets.smb = table( ...
      [datetime(2012, 1, 1); datetime(2013, 1, 1)], ...
      [datetime(2012, 2, 1); datetime(2013, 2, 1)], ...
      [0.42; 0.31], 'VariableNames', {'start_date', 'end_date', 'smb'});
   targets.smb.Properties.VariableUnits = {'', '', 'm w.e.'};
   save(observation_file, 'targets')

   old_close = get(groot, 'DefaultFigureCloseRequestFcn');
   restore_close = onCleanup(@() set(groot, ...
      'DefaultFigureCloseRequestFcn', old_close));
   set(groot, 'DefaultFigureCloseRequestFcn', @(~, ~) []);
   icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);
   clear restore_close

   figures = findall(groot, 'Type', 'Figure');
   names = string(get(figures, 'Name'));
   figures = figures(startsWith(names, "promice/kanm"));
   names = string(get(figures, 'Name'));
   testCase.addTeardown(@() delete(figures(isgraphics(figures))));

   profile_fig = figures(contains(names, "profiles"));
   testCase.verifyNumElements(profile_fig, 1);
   profile_lines = findall(profile_fig, 'Type', 'Line');
   profile_labels = string(get(profile_lines, 'DisplayName'));
   testCase.verifyEqual(sort(profile_labels), sort( ...
      ["PROMICE observations - core A"; ...
       "PROMICE observations - core B"]));
   testCase.verifyEqual(size(unique(vertcat(profile_lines.Color), 'rows'), 1), 2);
   profile_legend = findall(profile_fig, 'Type', 'Legend');
   testCase.verifyNumElements(profile_legend, 1);
   testCase.verifyEqual(profile_legend.Color, [1 1 1]);
   testCase.verifyEqual(profile_legend.TextColor, [0 0 0]);

   mass_fig = figures(contains(names, "canonical surface-mass-balance terms"));
   testCase.verifyNumElements(mass_fig, 1);
   mass_axes = findall(mass_fig, 'Type', 'Axes');
   mass_tags = string(get(mass_axes, 'Tag'));
   smb_ax = mass_axes(mass_tags == "smb");
   smb_lines = findall(smb_ax, 'Type', 'Line');
   smb_labels = string(get(smb_lines, 'DisplayName'));
   testCase.verifyFalse(any(smb_labels == "PROMICE observations"));
   testCase.verifyTrue(any(smb_labels == "MAR 3.11"));

   interval_fig = figures(contains(names, ...
      "observed interval surface mass balance"));
   testCase.verifyNumElements(interval_fig, 1);
   interval_axes = findall(interval_fig, 'Type', 'Axes');
   interval_titles = arrayfun(@(ax) string(ax.Title.String), interval_axes);
   interval_ax = interval_axes(interval_titles == ...
      "surface mass balance accumulated over observation interval");
   interval_lines = findall(interval_ax, 'Type', 'Line');
   interval_labels = string(get(interval_lines, 'DisplayName'));
   testCase.verifyFalse(any(interval_labels == "MAR 3.11"));
   observation_line = interval_lines( ...
      interval_labels == "PROMICE observations");
   testCase.verifyNumElements(observation_line, 1);
   testCase.verifyEqual(observation_line.YData(:), ...
      [0.42; 0.42; NaN; 0.31; 0.31; NaN]);
   testCase.verifyEqual(string(interval_ax.YLabel.String), ...
      "interval SMB [mWE]");

   component_fig = figures(contains(names, ...
      "precipitation and meltwater components"));
   testCase.verifyNumElements(component_fig, 1);
   component_axes = findall(component_fig, 'Type', 'Axes');
   component_titles = arrayfun(@(ax) string(ax.Title.String), component_axes);
   accumulation_ax = component_axes(startsWith(component_titles, ...
      "net snow accumulation"));
   testCase.verifyNumElements(accumulation_ax, 1);
   testCase.verifyEqual(string(accumulation_ax.YLabel.String), ...
      "snowf_subl [mmWE]");
end

function test_group_names_preserve_reductions_while_tiles_omit_redundancy(testCase)
   % Standalone figure names preserve reduction semantics while stacked tiles
   % leave repeated variable/support prose to axes labels and report captions.
   old_close = get(groot, 'DefaultFigureCloseRequestFcn');
   restore_close = onCleanup(@() set(groot, ...
      'DefaultFigureCloseRequestFcn', old_close));
   set(groot, 'DefaultFigureCloseRequestFcn', @(~, ~) []);
   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);
   clear restore_close

   % The stable summary tokens remain machine-readable while exported figure
   % names carry the reduction/support meaning needed for visual review.
   expected_groups = ["energy_balance", "mass_balance_terms", ...
      "radiation_fluxes", "surface_albedo"];
   testCase.verifyTrue(all(ismember(expected_groups, ...
      returned.figure_group)));
   figures = findall(groot, 'Type', 'Figure');
   names = string(get(figures, 'Name'));
   figures = figures(startsWith(names, "promice/kanm"));
   names = string(get(figures, 'Name'));
   testCase.addTeardown(@() delete(figures(isgraphics(figures))));

   testCase.verifyEqual(nnz(contains(names, ...
      "energy-balance terms (daily means; complete days only)")), 1);
   testCase.verifyEqual(nnz(contains(names, ...
      "canonical surface-mass-balance terms " ...
      + "(daily totals; complete days only)")), 1);
   testCase.verifyEqual(nnz(contains(names, ...
      "radiation fluxes (native support)")), 1);
    testCase.verifyEqual(nnz(contains(names, ...
       "surface albedo (daily source-aware means; " ...
       + "source-specific support rules)")), 1);
     testCase.verifyEqual(nnz(contains(names, ...
        "surface height change, snow depth, and stores " ...
        + "(daily means; complete forcing " ...
        + "days, available observation samples)")), 1);

   % The exported figure title keeps case/group context but leaves reduction
   % details to the report caption instead of consuming tile space.
   energy_fig = figures(contains(names, "energy-balance terms"));
   layouts = energy_fig.Children(arrayfun(@(h) ...
      isa(h, 'matlab.graphics.layout.TiledChartLayout'), ...
      energy_fig.Children));
   testCase.verifyNumElements(layouts, 1);
   testCase.verifyEqual(string(layouts.Title.String), ...
      "promice / kanm - energy-balance terms");

   % Met tiles use stable nonvisual tags and no redundant visible titles.
   met_fig = figures(contains(names, "met forcing"));
   testCase.verifyNumElements(met_fig, 1);
   met_axes = findall(met_fig, 'Type', 'Axes');
   met_tags = string(get(met_axes, 'Tag'));
   met_titles = arrayfun(@(ax) string(ax.Title.String), met_axes);
   testCase.verifyTrue(all(met_titles == ""));
   testCase.verifyEqual(sort(met_tags), sort(["downwelling radiation"; ...
      "air temperature and humidity"; "wind speed and surface pressure"; ...
      "surface albedo"; "precipitation"]));

   % Each paired ruler contributes scalar calls to one legend. Both the source
   % and variable must therefore be explicit, with no duplicate display text.
   temperature_ax = met_axes(met_tags == "air temperature and humidity");
   wind_ax = met_axes(met_tags == "wind speed and surface pressure");
   temperature_labels = string(get(findall(temperature_ax, ...
      'Type', 'Line'), 'DisplayName'));
   wind_labels = string(get(findall(wind_ax, ...
      'Type', 'Line'), 'DisplayName'));
   temperature_labels(temperature_labels == "") = [];
   wind_labels(wind_labels == "") = [];
   expected_temperature = ["PROMICE forcing - tair"; ...
      "MAR 3.11 forcing - tair"; "PROMICE forcing - rh"; ...
      "MAR 3.11 forcing - rh"];
   expected_wind = ["PROMICE forcing - wspd"; ...
      "MAR 3.11 forcing - wspd"; "PROMICE forcing - psfc"; ...
      "MAR 3.11 forcing - psfc"];
   testCase.verifyEqual(sort(temperature_labels), ...
      sort(expected_temperature));
   testCase.verifyEqual(sort(wind_labels), sort(expected_wind));
   testCase.verifyEqual(numel(unique(temperature_labels)), ...
      numel(temperature_labels));
   testCase.verifyEqual(numel(unique(wind_labels)), numel(wind_labels));

   % Paired scalar variables use the same source color while line style carries
   % the variable identity: left-axis channels are solid and right-axis channels
   % are dotted. The radiation pair follows the same solid/dotted convention.
   temperature_lines = findall(temperature_ax, 'Type', 'Line');
   temperature_styles = string(get(temperature_lines, 'LineStyle'));
   testCase.verifyTrue(all(temperature_styles( ...
      endsWith(temperature_labels, " - tair")) == "-"));
   testCase.verifyTrue(all(temperature_styles( ...
      endsWith(temperature_labels, " - rh")) == ":"));
   wind_lines = findall(wind_ax, 'Type', 'Line');
   wind_styles = string(get(wind_lines, 'LineStyle'));
   testCase.verifyTrue(all(wind_styles( ...
      endsWith(wind_labels, " - wspd")) == "-"));
   testCase.verifyTrue(all(wind_styles( ...
      endsWith(wind_labels, " - psfc")) == ":"));
   radiation_ax = met_axes(met_tags == "downwelling radiation");
   radiation_lines = findall(radiation_ax, 'Type', 'Line');
   radiation_labels = string(get(radiation_lines, 'DisplayName'));
   radiation_styles = string(get(radiation_lines, 'LineStyle'));
   testCase.verifyTrue(all(radiation_styles( ...
      endsWith(radiation_labels, " - swd")) == "-"));
   testCase.verifyTrue(all(radiation_styles( ...
      endsWith(radiation_labels, " - lwd")) == ":"));
end

function test_shared_figures_use_one_column_and_deduplicated_sources(testCase)
   % Retain the normally closed figures long enough to inspect the rendered
   % layout and provenance contract without writing report artifacts.
   old_close = get(groot, 'DefaultFigureCloseRequestFcn');
   restore_close = onCleanup(@() set(groot, ...
      'DefaultFigureCloseRequestFcn', old_close));
   set(groot, 'DefaultFigureCloseRequestFcn', @(~, ~) []);
   icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);
   clear restore_close

   % Select only this fixture's figures and always delete them after inspection.
   figures = findall(groot, 'Type', 'Figure');
   names = string(get(figures, 'Name'));
   figures = figures(startsWith(names, "promice/kanm"));
   names = string(get(figures, 'Name'));
   testCase.addTeardown(@() delete(figures(isgraphics(figures))));

   % The canonical met view uses no more than five full-width stacked axes.
   met_fig = figures(contains(names, "met forcing"));
   testCase.verifyNumElements(met_fig, 1);
   met_axes = findall(met_fig, 'Type', 'Axes');
   testCase.verifyGreaterThanOrEqual(numel(met_axes), 1);
   testCase.verifyLessThanOrEqual(numel(met_axes), 5);
   left = arrayfun(@(ax) ax.Position(1), met_axes);
   width = arrayfun(@(ax) ax.Position(3), met_axes);
   testCase.verifyEqual(left, repmat(left(1), size(left)), 'AbsTol', 1e-12);
   testCase.verifyEqual(width, repmat(width(1), size(width)), 'AbsTol', 1e-12);
   for k = 2:numel(met_axes)
      testCase.verifyEqual(met_axes(k).XLim, met_axes(1).XLim);
   end
   met_legends = findall(met_fig, 'Type', 'Legend');
   testCase.verifyGreaterThanOrEqual(numel(met_legends), 1);
   % Final layout freezes the position first selected by BEST so later y-limit
   % clearance cannot make MATLAB relocate the legend across the data crown.
   testCase.verifyTrue(all(string(get(met_legends, 'Location')) == "none"));
   testCase.verifyTrue(all(arrayfun(@(lgd) lgd.NumColumns == ...
      min(4, numel(lgd.String)), met_legends)));

   % The resolved legend geometry reserves only the occupied visual side of
   % the downwelling-radiation tile. Its peak crown must remain below the
   % legend border instead of being hidden by the compact in-axes legend.
   downwelling_ax = met_axes(string(get(met_axes, 'Tag')) == ...
      "downwelling radiation");
   testCase.verifyNumElements(downwelling_ax, 1);
   drawnow limitrate nocallbacks
   axes_pixels = getpixelposition(downwelling_ax, true);
   downwelling_legend = gobjects(0, 1);
   for lgd = reshape(met_legends, 1, [])
      legend_pixels = getpixelposition(lgd, true);
      legend_center = legend_pixels(1:2) + legend_pixels(3:4) / 2;
      if legend_center(1) >= axes_pixels(1) ...
            && legend_center(1) <= axes_pixels(1) + axes_pixels(3) ...
            && legend_center(2) >= axes_pixels(2) ...
            && legend_center(2) <= axes_pixels(2) + axes_pixels(4)
         downwelling_legend(end + 1, 1) = lgd; %#ok<AGROW>
      end
   end
   testCase.verifyNumElements(downwelling_legend, 1);
   legend_pixels = getpixelposition(downwelling_legend, true);
   legend_bottom = (legend_pixels(2) - axes_pixels(2)) / axes_pixels(4);
   testCase.verifyGreaterThan(legend_bottom, 2 / 3);
   lines = findall(downwelling_ax, 'Type', 'Line');
   ydata = get(lines, 'YData');
   if ~iscell(ydata)
      ydata = {ydata};
   end
   finite_max = max(cellfun(@(v) max(v, [], 'omitnan'), ydata));
   data_top = (finite_max - downwelling_ax.YLim(1)) ...
      / diff(downwelling_ax.YLim);
   testCase.verifyLessThan(data_top, legend_bottom);

   % The generic stacked renderer follows the same compact in-axes policy and
   % leaves a small manual y margin for best-position legend placement.
   radiation_fig = figures(contains(names, "radiation fluxes"));
   radiation_legends = findall(radiation_fig, 'Type', 'Legend');
   testCase.verifyGreaterThanOrEqual(numel(radiation_legends), 1);
   testCase.verifyTrue(all(string(get(radiation_legends, 'Location')) == "none"));
   testCase.verifyTrue(all(arrayfun(@(lgd) lgd.NumColumns == ...
      min(4, numel(lgd.String)), radiation_legends)));
   swd_ax = findall(radiation_fig, 'Type', 'Axes', 'Tag', 'swd');
   testCase.verifyNumElements(swd_ax, 1);
   testCase.verifyEqual(string(swd_ax.YLimMode), "manual");

   % Only the two physically paired panels create yyaxis rulers. Radiation,
   % albedo, and precipitation retain one ruler without an empty themed peer.
   met_tags = string(get(met_axes, 'Tag'));
   met_titles = arrayfun(@(ax) string(ax.Title.String), met_axes);
   single = ismember(met_tags, ...
      ["downwelling radiation", "surface albedo", "precipitation"]);
   paired = ismember(met_tags, ...
      ["air temperature and humidity", "wind speed and surface pressure"]);
     testCase.verifyEqual(nnz(single), 3);
     testCase.verifyEqual(nnz(paired), 2);
     testCase.verifyTrue(all(met_titles == ""));
     testCase.verifyTrue(all(arrayfun(@(ax) isscalar(ax.YAxis), ...
        met_axes(single))));
   testCase.verifyTrue(all(arrayfun(@(ax) numel(ax.YAxis) == 2, ...
      met_axes(paired))));

   % Prefer native PROMICE observations over their repeated userdata columns,
   % and collapse MODIS columns repeated by MAR and RACMO into one overlay.
   height_fig = figures(contains(names, ...
      "surface height change, snow depth, and stores"));
   albedo_fig = figures(contains(names, "surface albedo"));
   testCase.verifyNumElements(height_fig, 1);
   testCase.verifyNumElements(albedo_fig, 1);
   height_labels = figureLegendLabels(height_fig);
   albedo_labels = figureLegendLabels(albedo_fig);
    % Separate axes may each repeat the same source label, but native PROMICE
    % evaluation channels must never fall back to the ambiguous bare source.
    testCase.verifyGreaterThanOrEqual( ...
       nnz(height_labels == "PROMICE observations"), 1);
   testCase.verifyFalse(any(height_labels == "PROMICE"));
   testCase.verifyEqual(nnz(albedo_labels == "MODIS (GEUS)"), 1);

   % Every renderer emits human-facing labels and explicit light legend colors.
   legends = findall(figures, 'Type', 'Legend');
   labels = figureLegendLabels(figures);
   testCase.verifyFalse(any(contains(lower(labels), "userdata")));
   testCase.verifyFalse(any(contains(labels, ":")));
   testCase.verifyTrue(all(arrayfun(@(lgd) isequal(lgd.Color, [1 1 1]), ...
      legends)));
   testCase.verifyTrue(all(arrayfun(@(lgd) isequal(lgd.TextColor, [0 0 0]), ...
      legends)));
   lines = findall(figures, 'Type', 'Line');
   testCase.verifyTrue(all(string(get(lines, 'Marker')) == "none"));
   testCase.verifyEmpty(findall(figures, 'Type', 'Line', 'LineStyle', 'none'));
end

function test_timeseries_tiles_share_union_of_finite_time_support(testCase)
   % Different variable edge support must not give adjacent diagnostic panels
   % different time ranges or hide which event aligns across energy channels.
   mar_file = fullfile(testCase.TestData.root, 'input', 'userdata', ...
      'mar3.11', 'kanm_mar3.11_20120101_20120102.mat');
   loaded = load(mar_file, 'Data');
   Data = loaded.Data;
   Data.swn(1:24) = NaN;
   Data.lwn(25:48) = NaN;
   save(mar_file, 'Data')

   old_close = get(groot, 'DefaultFigureCloseRequestFcn');
   restore_close = onCleanup(@() set(groot, ...
      'DefaultFigureCloseRequestFcn', old_close));
   set(groot, 'DefaultFigureCloseRequestFcn', @(~, ~) []);
   icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="promice", save_figs=false);
   clear restore_close

   % Inspect the normally closed energy figure and clean it up after the test.
   figures = findall(groot, 'Type', 'Figure');
   names = string(get(figures, 'Name'));
   figures = figures(startsWith(names, "promice/kanm"));
   names = string(get(figures, 'Name'));
   testCase.addTeardown(@() delete(figures(isgraphics(figures))));
   energy_fig = figures(contains(names, "energy-balance terms"));
   testCase.verifyNumElements(energy_fig, 1);
   energy_axes = findall(energy_fig, 'Type', 'Axes');
   tags = string(get(energy_axes, 'Tag'));
   time_axes = energy_axes(ismember(tags, ["swn", "lwn", "thf"]));
   testCase.verifyNumElements(time_axes, 3);

   expected = [datetime(2012, 1, 1, 'TimeZone', 'UTC'), ...
      datetime(2012, 1, 2, 'TimeZone', 'UTC')];
   for ax = reshape(time_axes, 1, [])
      returned = ax.XLim;
      testCase.verifyEqual(returned, expected);
   end
end

function test_plot_verification_artifacts_overwrite_clears_case_pngs(testCase)
   % Overwrite should remove obsolete group PNGs only from selected cases.
   figure_root = fullfile(testCase.TestData.root, 'overwrite_figures');
   case_dir = fullfile(figure_root, 'promice', 'kanm');
   unselected_dir = fullfile(figure_root, 'esm_snowmip', 'cdp');
   icemodel.helpers.ensureDirExists(case_dir);
   icemodel.helpers.ensureDirExists(unselected_dir);

   % Seed obsolete and interrupted-run names, an invalid current-name PNG, and
   % one unselected-family sentinel. Overwrite must replace only the chosen case.
   obsolete = fullfile(case_dir, ["energy_fluxes.png", "mass_fluxes.png", ...
      "partial_only.png", "measurement_geometry.png", ...
      "measurement_geometry_2.png", "step_magnitude.png", ...
      "tice10m_source.png"]);
   invalid_current = fullfile(case_dir, 'met_forcing.png');
   unselected = fullfile(unselected_dir, 'energy_fluxes.png');
   for file = [obsolete, string(invalid_current), string(unselected)]
      fclose(fopen(file, 'w'));
   end

   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), dataset_family="promice", ...
      figure_root=figure_root, save_figs=true, overwrite=true);
   current = fullfile(case_dir, 'radiation_fluxes.png');
   current_met = fullfile(case_dir, 'met_forcing.png');
   current_surface = fullfile(case_dir, 'surface_height_depth.png');

   testCase.verifyFalse(any(isfile(obsolete)));
   testCase.verifyTrue(isfile(current));
   testCase.verifyTrue(isfile(current_met));
   testCase.verifyTrue(isfile(current_surface));
   met_info = imfinfo(current_met);
   testCase.verifyEqual(string(met_info.Format), "png");
   testCase.verifyGreaterThan(met_info.Width, 0);
   testCase.verifyGreaterThan(met_info.Height, 0);
   testCase.verifyTrue(any(returned.figure_file == string(current)));
   surface = returned(returned.figure_group == "surface_height_depth", :);
   testCase.verifyEqual(surface.figure_file, string(current_surface));
   testCase.verifyFalse(contains(surface.plotted_variables, "step_magnitude"));
   diagnostics = returned( ...
      returned.figure_group == "meteorological_diagnostics", :);
   testCase.verifyFalse(contains( ...
      diagnostics.plotted_variables, "tice10m_source"));
   unplotted = strjoin(returned.unplotted_variables, ',');
   testCase.verifyTrue(contains(unplotted, "step_magnitude"));
   testCase.verifyTrue(contains(unplotted, "tice10m_source"));
   testCase.verifyFalse(any(startsWith( ...
      returned.figure_group, "measurement_geometry")));
   testCase.verifyTrue(isfile(unselected));
end

function test_plot_verification_artifacts_uses_colbeck_grid(testCase)
   % Sparse Laugh-Test experiments should bypass the generic legend-heavy view.
   figure_root = fullfile(testCase.TestData.root, 'laugh_figures');
   returned = icemodel.verification.plotVerificationArtifacts( ...
      output_root=string(testCase.TestData.root), ...
      dataset_family="laugh_tests", figure_root=figure_root, save_figs=true);

   testCase.verifyEqual(height(returned), 1);
   testCase.verifyEqual(returned.figure_group, "colbeck_comparison");
   testCase.verifyTrue(isfile(returned.figure_file));
   testCase.verifyTrue(contains(returned.plotted_variables, ...
      "snow_liquid_water_storage_m"));
   testCase.verifyTrue(contains(returned.plotted_variables, ...
      "bottom_outflow_mps"));
end

function test_colbeck_grid_keeps_legend_outside_data(testCase)
   % The reused experiment-grid view should reserve the plot area for data.
   eval_root = fullfile(testCase.TestData.root, 'eval');
   loaded = load(fullfile(eval_root, 'laugh_tests', 'colbeck1976', ...
      'evaluation.mat'), 'targets');
   fig = icemodel.verification.plotcase("colbeck1976", ...
      evaluation_data_root=string(eval_root), ...
      input_data_root=fullfile(testCase.TestData.root, 'input'), ...
      dataset_family="laugh_tests", source="compare", ...
      candidate=loaded.targets.analytical_clark2017, visible="off");
   testCase.addTeardown(@() close(fig));

   legends = findobj(fig, 'Type', 'Legend');
   testCase.verifyNumElements(legends, 1);
   testCase.verifyEqual(string(legends.Location), "eastoutside");
   lines = findall(fig, 'Type', 'Line');
   testCase.verifyTrue(all(string(get(lines, 'Marker')) == "none"));
   testCase.verifyEmpty(findall(fig, 'Type', 'Line', 'LineStyle', 'none'));
end

function tf = hasListedTice1(value)
   %HASLISTEDTICE1 True when tice1 appears as a complete plotted variable token.
   parts = strip(split(string(value), ","));
   tf = any(parts == "tice1");
end

function variables = listedVariables(value)
   %LISTEDVARIABLES Parse one or more comma-separated summary cells.
   variables = strip(split(string(value), ","));
   variables(variables == "") = [];
end

function labels = figureLegendLabels(figures)
   %FIGURELEGENDLABELS Collect every rendered legend label from figure handles.
   legends = findall(figures, 'Type', 'Legend');
   counts = arrayfun(@(lgd) numel(lgd.String), legends);
   labels = strings(sum(counts), 1);
   first = 1;
   for k = 1:numel(legends)
      rows = first:(first + counts(k) - 1);
      labels(rows) = string(legends(k).String(:));
      first = first + counts(k);
   end
end

function test_plot_forcing_filters_to_date_window(testCase)
   % The first-class forcing plot supports hourly inspection of a short window
   % inside a multi-day or multi-year met file.

   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   met = makePlotMet();

   returned = icemodel.plot.forcing(met, axes=ax, names="met", ...
      frequency="hourly", ...
      startdate=datetime(2012, 1, 2, 0, 0, 0, 'TimeZone', 'UTC'), ...
      enddate=datetime(2012, 1, 2, 23, 0, 0, 'TimeZone', 'UTC'));

   x = returned.energy(1).XData;
   testCase.verifyGreaterThanOrEqual(min(x), ...
      datetime(2012, 1, 2, 0, 0, 0, 'TimeZone', 'UTC'));
   testCase.verifyLessThanOrEqual(max(x), ...
      datetime(2012, 1, 2, 23, 0, 0, 'TimeZone', 'UTC'));
end

function test_plot_forcing_uses_compact_readable_legend(testCase)
   % Shared legends stay compact and explicitly readable under dark themes.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);

   icemodel.plot.forcing(makePlotMet(), axes=ax, names="met");

   legends = findobj(fig, 'Type', 'Legend');
   testCase.verifyNumElements(legends, 1);
   testCase.verifyEqual(string(legends.Location), "best");
   testCase.verifyEqual(legends.Color, [1 1 1]);
   testCase.verifyEqual(legends.TextColor, [0 0 0]);
   testCase.verifyLessThanOrEqual(legends.NumColumns, 3);
end

function test_compare_timeseries_filters_to_date_window(testCase)
   % Grouped time-series panels use the same inclusive window controls as the
   % full artifact plotter.

   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   met = makePlotMet();

   returned = icemodel.plot.compareTimeseries(met, "tair", axes=ax, ...
      frequency="hourly", ...
      startdate=datetime(2012, 1, 2, 0, 0, 0, 'TimeZone', 'UTC'), ...
      enddate=datetime(2012, 1, 2, 23, 0, 0, 'TimeZone', 'UTC'));

   x = returned.lines(1).XData;
   testCase.verifyGreaterThanOrEqual(min(x), ...
      datetime(2012, 1, 2, 0, 0, 0, 'TimeZone', 'UTC'));
   testCase.verifyLessThanOrEqual(max(x), ...
      datetime(2012, 1, 2, 23, 0, 0, 'TimeZone', 'UTC'));
end

function test_timeseries_breaks_safely_inferred_omitted_time_gap(testCase)
   % Repeated hourly support makes one later outage unambiguous to the renderer.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours([0; 1; 2; 10; 11]);
   values = (1:numel(time))';

   line = icemodel.plot.timeseries(time, values, axes=ax);

   testCase.verifyEqual(numel(line.XData), numel(time) + 1);
   testCase.verifyEqual(nnz(isnan(line.YData)), 1);
   testCase.verifyTrue(all(ismember(time, line.XData)));
   testCase.verifyEqual(line.YData(isfinite(line.YData))', values);
end

function test_timeseries_uses_small_sparse_markers(testCase)
   % Generic sparse series remain visible without an oversized marker carpet.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:4)');
   values = [1; NaN; NaN; NaN; 2];

   icemodel.plot.timeseries(time, values, axes=ax);

   markers = findall(ax, 'Type', 'Line', 'LineStyle', 'none');
   testCase.verifyNumElements(markers, 1);
   testCase.verifyEqual(markers.MarkerSize, 6);
end

function test_source_colors_match_runoff_palette_and_ignore_order(testCase)
   % Canonical and role colors are keyed by labels, never the active subset.
   names = ["MAR 3.11", "RACMO 2.3p3", "MERRA-2", ...
      "PROMICE forcing", "PROMICE observations", ...
      "observations:data observations", "ESM-SnowMIP forcing", ...
      "Laugh Tests observations", "MODIS (GEUS)", ...
      "SUMup observations", "future source"];
   expected_models = [0.866 0.329 0; 0.25 0.80 0.54; ...
      0 0 0.172413793103448];
   colors = icemodel.plot.sourceColor(names);
   order = [11 3 1 7 2 8 5 4 6 9 10];

   testCase.verifyEqual(colors(1:3, :), expected_models, 'AbsTol', 1e-15);
   testCase.verifyEqual(icemodel.plot.sourceColor(names(order)), ...
      colors(order, :), 'AbsTol', 1e-15);
   testCase.verifyEqual(colors(4, :), colors(5, :));
   testCase.verifyNotEqual(colors(5, :), colors(6, :));
   testCase.verifyEqual(colors(9, :), [0.55 0.55 0.55], ...
      'AbsTol', 1e-15);
   testCase.verifyNotEqual(colors(9, :), colors(3, :));
   testCase.verifyEqual(colors(10, :), [0.0000 0.4470 0.7410], ...
      'AbsTol', 1e-15);
   testCase.verifyNotEqual(colors(10, :), colors(1, :));
   testCase.verifyNotEqual(colors(10, :), colors(2, :));
   testCase.verifyNotEqual(colors(10, :), colors(3, :));
   testCase.verifyNotEqual(colors(10, :), colors(9, :));

   % SUMup keeps one explicit color in every plotting role and panel.
   sumup = icemodel.plot.sourceColor(["SUMup observations", ...
      "SUMup forcing - smb"]);
   testCase.verifyEqual(sumup, repmat(colors(10, :), 2, 1), ...
      'AbsTol', 1e-15);

   % Unknown native families retain one fallback color across plot roles and
   % variable-qualified legend labels.
   retmip = icemodel.plot.sourceColor(["RetMIP forcing - swd", ...
      "RetMIP forcing - lwd", "RetMIP observations"]);
   imau = icemodel.plot.sourceColor(["IMAU forcing - tair", ...
      "IMAU forcing - psfc", "IMAU observations"]);
   testCase.verifyEqual(retmip, repmat(retmip(1, :), 3, 1));
   testCase.verifyEqual(imau, repmat(imau(1, :), 3, 1));
end

function test_compare_timeseries_uses_source_colors(testCase)
   % The grouped overlay forwards canonical source labels to the shared palette.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   met = makePlotMet();
   names = ["RACMO 2.3p3", "MAR 3.11"];

   returned = icemodel.plot.compareTimeseries({met, met}, "tair", ...
      axes=ax, names=names, frequency="hourly");

   testCase.verifyEqual(vertcat(returned.lines.Color), ...
      icemodel.plot.sourceColor(names), 'AbsTol', 1e-15);
   legends = findobj(fig, 'Type', 'Legend');
   testCase.verifyNumElements(legends, 1);
   testCase.verifyEqual(legends.Color, [1 1 1]);
   testCase.verifyEqual(legends.TextColor, [0 0 0]);
end

function test_compare_timeseries_uses_per_source_daily_completeness(testCase)
   % Sparse observations use finite samples without weakening dense forcing.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:47))';
   observation_values = nan(size(Time));
   observation_values([6, 30]) = [0.2, 0.4];
   forcing_values = ones(size(Time));
   forcing_values(4) = NaN;
   observations = timetable(Time, observation_values, ...
      'VariableNames', {'snow_depth_m'});
   forcing = timetable(Time, forcing_values, ...
      'VariableNames', {'snow_depth_m'});

   returned = icemodel.plot.compareTimeseries( ...
      {observations, forcing}, "snow_depth_m", axes=ax, ...
      names=["ESM-SnowMIP observations", "MAR 3.11 forcing"], ...
      frequency="daily", aggregation=["mean_omitmissing", "mean"]);

   testCase.verifyEqual(returned.plotted, [true; true]);
   testCase.verifyEqual(returned.lines(1).YData(:), [0.2; 0.4], ...
      'AbsTol', 1e-15);
   forcing_daily = returned.lines(2).YData(:);
   testCase.verifyTrue(isnan(forcing_daily(1)));
   testCase.verifyEqual(forcing_daily(2), 1, 'AbsTol', 1e-15);
   testCase.verifyError(@() icemodel.plot.compareTimeseries( ...
      {observations, forcing}, "snow_depth_m", axes=ax, ...
      aggregation=["mean", "mean", "mean"]), ...
      'icemodel:plot:compareTimeseries:badAggregation');
end

function test_compare_timeseries_preserves_interval_observation_spans(testCase)
   % Accumulated interval observations should not be daily-integrated again.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   time = [datetime(2012, 1, 1, 6, 0, 0, 'TimeZone', 'UTC'); ...
      datetime(2012, 1, 31, 18, 0, 0, 'TimeZone', 'UTC')];
   obs = timetable([0.42; 0.42], 'RowTimes', time, ...
      'VariableNames', {'smb'});
   obs.Properties.VariableUnits = {'m w.e.'};
   obs.Properties.UserData = struct('interval_observation', true, ...
      'interval_variables', "smb");

   returned = icemodel.plot.compareTimeseries(obs, "smb", axes=ax, ...
      aggregation="daily_total");

   x = returned.lines(1).XData;
   y = returned.lines(1).YData;
   testCase.verifyEqual(x(1), time(1));
   testCase.verifyEqual(x(end), time(end));
   testCase.verifyEqual(y(:), [0.42; 0.42]);
   testCase.verifyEqual(returned.unit, "mWE");
   testCase.verifyEqual(string(ax.YLabel.String), "smb [mWE]");
end

function test_compare_timeseries_displays_daily_mass_totals_in_mmwe(testCase)
   % Rates from source-specific metre conventions should overlay as daily mmWE.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2012, 1, 2, 23, 0, 0, 'TimeZone', 'UTC'))';
   retmip = timetable(repmat(0.001, numel(Time), 1), 'RowTimes', Time, ...
      'VariableNames', {'melt'});
   retmip.Properties.VariableUnits = {'mWE/h'};
   model = timetable(repmat(0.001 / 3600, numel(Time), 1), ...
      'RowTimes', Time, 'VariableNames', {'melt'});
   model.Properties.VariableUnits = {'m s-1'};

   returned = icemodel.plot.compareTimeseries({retmip, model}, "melt", ...
      axes=ax, names=["RetMIP observations", "MAR 3.11"], ...
      aggregation="daily_total", frequency="daily", display_unit="mmWE");

   testCase.verifyEqual(returned.unit, "mmWE");
   testCase.verifyEqual(returned.lines(1).YData(:), [24; 24], ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(returned.lines(2).YData(:), [24; 24], ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(string(ax.YLabel.String), "melt [mmWE]");
end

function test_compare_timeseries_masks_promice_partial_daily_means(testCase)
   % Preserve EGP's 24/24 extreme while masking the observed KAN_M 18/24 and
   % KAN_U 16/24 coverage patterns as incomplete dense daily summaries.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:71)');
   values = ones(size(time));
   values(25:30) = NaN;
   values(49:56) = NaN;
   T = timetable(values, 'RowTimes', time, 'VariableNames', {'tair'});
   T.Properties.VariableUnits = {'degC'};

   returned = icemodel.plot.compareTimeseries(T, "tair", axes=ax, ...
      frequency="daily", aggregation="mean");

   testCase.verifyEqual(returned.lines(1).YData(1), 1);
   testCase.verifyTrue(all(isnan(returned.lines(1).YData(2:3))));
end

function test_compare_timeseries_integrates_native_without_gap_inflation(testCase)
   % Daily rate integration accepts a complete day and rejects even one omitted
   % posting instead of presenting a partial accumulation as a daily total.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   full_time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:47)');
   keep = true(size(full_time));
   keep(35) = false;
   time = full_time(keep);
   rate = 1e-5 * ones(size(time));
   T = timetable(rate, 'RowTimes', time, 'VariableNames', {'ppt'});
   T.Properties.VariableUnits = {'m s-1'};

   returned = icemodel.plot.compareTimeseries(T, "ppt", axes=ax, ...
      frequency="daily", aggregation="daily_total");

   testCase.verifyEqual(returned.lines(1).YData(1), 24 * 3600e-5, ...
      'AbsTol', 1e-12);
   testCase.verifyTrue(isnan(returned.lines(1).YData(2)));
end

function test_compare_timeseries_preserves_missing_daily_totals(testCase)
   % Sparse totals retain native timestamps rather than acquiring a daily grid.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   time = [datetime(2012, 1, 1, 6, 0, 0, TimeZone="UTC"); ...
      datetime(2012, 1, 3, 6, 0, 0, TimeZone="UTC")];
   T = timetable([0.1; 0.2], 'RowTimes', time, ...
      'VariableNames', {'smb'});
   T.Properties.VariableUnits = {'m w.e.'};

   returned = icemodel.plot.compareTimeseries(T, "smb", axes=ax, ...
      frequency="daily", aggregation="daily_total");

   testCase.verifyEqual(returned.lines(1).XData(:), time);
   testCase.verifyEqual(returned.lines(1).YData(:), [0.1; 0.2]);

   % A lone rate posting has no inferable represented duration. It must not
   % become either an invented total or an all-NaN ghost legend entry.
   single = timetable(1e-5, 'RowTimes', time(1), ...
      'VariableNames', {'ppt'});
   single.Properties.VariableUnits = {'m s-1'};
   fig_single = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_single));
   returned = icemodel.plot.compareTimeseries(single, "ppt", ...
      axes=axes(fig_single), frequency="daily", aggregation="daily_total");
   testCase.verifyFalse(any(returned.plotted));
   testCase.verifyEmpty(returned.lines);
end

function test_compare_timeseries_accepts_exact_dense_15m_days(testCase)
   % Exact 15-minute support proves both daily mean and rate-total behavior.
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + minutes((0:191)' .* 15);
   T = timetable(2 .* ones(size(time)), 1e-5 .* ones(size(time)), ...
      0.1 .* ones(size(time)), 0.01 .* ones(size(time)), 'RowTimes', time, ...
      'VariableNames', {'tair', 'ppt', 'runoff', 'smb'});
   T.Properties.VariableUnits = {'degC', 'm s-1', 'mWE/h', 'm w.e.'};

   fig_mean = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_mean));
   mean_out = icemodel.plot.compareTimeseries(T, "tair", ...
      axes=axes(fig_mean), frequency="daily", aggregation="mean");
   testCase.verifyEqual(mean_out.lines(1).YData(:), [2; 2]);

   % A rate is integrated over all 96 expected postings in each UTC day.
   fig_total = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_total));
   total_out = icemodel.plot.compareTimeseries(T, "ppt", ...
      axes=axes(fig_total), frequency="daily", aggregation="daily_total");
   testCase.verifyEqual(total_out.lines(1).YData(:), ...
      [86400e-5; 86400e-5], 'AbsTol', 1e-12);

   % The second supported rate unit and already-discrete amounts share the same
   % exact-grid gate but retain their distinct integration semantics.
   fig_mwe = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_mwe));
   mwe_out = icemodel.plot.compareTimeseries(T, "runoff", ...
      axes=axes(fig_mwe), frequency="daily", aggregation="daily_total");
   testCase.verifyEqual(mwe_out.lines(1).YData(:), [2.4; 2.4], ...
      'AbsTol', 1e-12);
   fig_amount = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_amount));
   amount_out = icemodel.plot.compareTimeseries(T, "smb", ...
      axes=axes(fig_amount), frequency="daily", aggregation="daily_total");
   testCase.verifyEqual(amount_out.lines(1).YData(:), [0.96; 0.96], ...
      'AbsTol', 1e-12);

   % A non-UTC display zone still aggregates by the same UTC-day instants.
   zoned = T;
   % Timetable row times require whole-variable assignment; nested assignment
   % dispatches through table expansion rather than datetime's TimeZone setter.
   time_dimension = zoned.Properties.DimensionNames{1};
   zoned_time = zoned.(time_dimension);
   zoned_time.TimeZone = 'America/Los_Angeles';
   zoned.(time_dimension) = zoned_time;
   fig_zoned = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_zoned));
   zoned_out = icemodel.plot.compareTimeseries(zoned, "tair", ...
      axes=axes(fig_zoned), frequency="daily", aggregation="mean");
   expected_day = datetime(2012, 1, (1:2)', TimeZone="UTC");
   testCase.verifyEqual(zoned_out.lines(1).XData(:), expected_day);
   testCase.verifyEqual(zoned_out.lines(1).YData(:), [2; 2]);
end

function test_compare_timeseries_vectorized_days_match_grid(testCase)
   % One vectorized pass must keep finite-value and timestamp defects local to
   % their day while producing the same complete-day means and rate totals.
   n_slots = 96;
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + minutes((0:(5 * n_slots - 1))' .* 15);
   tair = reshape(repmat(1:5, n_slots, 1), [], 1);
   ppt = 1e-5 .* ones(size(time));
   tair(n_slots + 8) = NaN;
   time(2 * n_slots + 8) = time(2 * n_slots + 7);
   time(3 * n_slots + 8) = time(3 * n_slots + 8) + minutes(5);
   T = timetable(tair, ppt, 'RowTimes', time, ...
      'VariableNames', {'tair', 'ppt'});
   T.Properties.VariableUnits = {'degC', 'm s-1'};

   fig_mean = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_mean));
   mean_out = icemodel.plot.compareTimeseries(T, "tair", ...
      axes=axes(fig_mean), frequency="daily", aggregation="mean");
   testCase.verifyEqual(mean_out.lines(1).YData(:), ...
      [1; NaN; NaN; NaN; 5]);

   fig_total = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_total));
   total_out = icemodel.plot.compareTimeseries(T, "ppt", ...
      axes=axes(fig_total), frequency="daily", aggregation="daily_total");
   expected_total = 86400e-5;
   testCase.verifyEqual(total_out.lines(1).YData(:), ...
      [expected_total; expected_total; NaN; NaN; expected_total], ...
      'AbsTol', 1e-12);
end

function test_compare_timeseries_dense_daily_runtime_is_linear(testCase)
   % Twenty-five years of 15-minute data must stay below a broad runtime guard
   % that the former O(n_days*n_rows) mask loop cannot satisfy.
   n_days = 25 * 365;
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + minutes((0:(n_days * 96 - 1))' .* 15);
   T = timetable(ones(size(time)), 'RowTimes', time, ...
      'VariableNames', {'tair'});
   T.Properties.VariableUnits = {'degC'};
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   timer = tic;
   returned = icemodel.plot.compareTimeseries(T, "tair", axes=axes(fig), ...
      frequency="daily", aggregation="mean");
   elapsed = toc(timer);

   testCase.verifyEqual(numel(returned.lines(1).YData), n_days);
   testCase.verifyTrue(all(returned.lines(1).YData == 1));
   testCase.verifyLessThan(elapsed, 20);
end

function test_compare_timeseries_accepts_centered_native_phase(testCase)
   % MERRA-style interval-center timestamps retain their repeated native phase.
   time = datetime(2012, 1, 1, 1, 30, 0, TimeZone="UTC") ...
      + hours((0:15)' .* 3);
   values = [ones(8, 1); 2 .* ones(8, 1)];
   T = timetable(values, 'RowTimes', time, 'VariableNames', {'tair'});
   T.Properties.VariableUnits = {'degC'};
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries(T, "tair", axes=axes(fig), ...
      frequency="daily", aggregation="mean");

   testCase.verifyEqual(returned.lines(1).YData(:), [1; 2]);
end

function test_compare_timeseries_omitmissing_requires_full_grid(testCase)
   % Structural value gaps are allowed, but an omitted native timestamp is not.
   full_time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:47)');
   keep = true(size(full_time));
   keep(30) = false;
   time = full_time(keep);
   albedo = nan(size(time));
   albedo(5:7) = [0.7; 0.71; 0.72];
   albedo(29:31) = [0.65; 0.66; 0.67];
   T = timetable(albedo, 'RowTimes', time, 'VariableNames', {'albedo'});
   T.Properties.VariableUnits = {'-'};
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries(T, "albedo", axes=axes(fig), ...
      frequency="daily", aggregation="mean_omitmissing");

   testCase.verifyEqual(returned.lines(1).YData(1), 0.71, 'AbsTol', 1e-12);
   testCase.verifyTrue(isnan(returned.lines(1).YData(2)));
end

function test_compare_timeseries_daily_albedo_mean_is_shortwave_weighted(testCase)
   % Every shared daily-albedo call inherits reflected/incoming energy semantics.
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:47)');
   albedo = nan(size(time));
   swd = zeros(size(time));
   first = 7:12;
   second = 31:36;
   albedo(first) = [0.4; 0.4; 0.4; 0.9; 0.9; 0.9];
   swd(first) = [10; 10; 10; 500; 500; 500];
   albedo(second) = [0.5; 0.5; 0.5; 0.8; 0.8; 0.8];
   swd(second) = [20; 20; 20; 400; 400; 400];
   T = timetable(albedo, swd, 'RowTimes', time);
   T.Properties.VariableUnits = {'-', 'W m-2'};
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries(T, "albedo", axes=axes(fig), ...
      frequency="daily", aggregation="mean");

   expected = [(0.4 * 10 + 0.9 * 500) / 510; ...
      (0.5 * 20 + 0.8 * 400) / 420];
   testCase.verifyEqual(returned.lines(1).YData(:), expected, ...
      'AbsTol', 1e-12);
end

function test_compare_timeseries_daily_albedo_requires_sunlit_support(testCase)
   % A complete hourly axis with only a few valid polar-shoulder observations
   % must not present those samples as a representative daily albedo.
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:47)');
   albedo = nan(size(time));
   swd = zeros(size(time));
   albedo(7:10) = 0.4;
   swd(7:10) = 100;
   albedo(31:36) = 0.8;
   swd(31:36) = 100;
   T = timetable(albedo, swd, 'RowTimes', time);
   T.Properties.VariableUnits = {'-', 'W m-2'};
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries(T, "albedo", axes=axes(fig), ...
      frequency="daily", aggregation="mean");

   testCase.verifyTrue(isnan(returned.lines(1).YData(1)));
   testCase.verifyEqual(returned.lines(1).YData(2), 0.8, ...
      'AbsTol', 1e-12);
end

function test_compare_timeseries_albedo_weighting_has_sparse_fallback(testCase)
   % Albedo products without SWD retain the finite-sample daily reduction.
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:23)');
   albedo = nan(size(time));
   albedo([7, 13]) = [0.4; 0.9];
   T = timetable(albedo, 'RowTimes', time);
   T.Properties.VariableUnits = {'-'};
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries(T, "albedo", axes=axes(fig), ...
      frequency="daily", aggregation="shortwave_weighted_mean");

   testCase.verifyEqual(returned.lines(1).YData, 0.65, ...
      'AbsTol', 1e-12);
end

function test_compare_timeseries_source_aware_model_albedo(testCase)
   % Native model states remain valid at night, model-derived ratios need only
   % positive modeled shortwave, and radiometers retain the six-hour gate.
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:47)');
   albedo = [0.82 .* ones(24, 1); 0.7 .* ones(12, 1); 0.9 .* ones(12, 1)];
   swd = [zeros(24, 1); 10 .* ones(12, 1); 500 .* ones(12, 1)];
   radiometer = timetable(albedo, swd, 'RowTimes', time);
   radiometer.Properties.VariableUnits = {'-', 'W m-2'};

   % Native model-state support is allowed to be structurally missing, but the
   % exact timestamp grid must still be complete before a finite daily value is
   % presented.
   state = radiometer;
   state_missing = radiometer;
   state_missing.albedo(1:24) = NaN;
   state_missing.albedo(25:44) = NaN;

   % Fewer than six positive-SWD model samples still define a deterministic
   % daily energy ratio; zero-SWD polar night remains honestly undefined.
   ratio = radiometer;
   ratio.swd(25:48) = 0;
   ratio.swd(31:34) = [10; 10; 500; 500];
   ratio.albedo(31:34) = [0.7; 0.7; 0.9; 0.9];
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries( ...
      {radiometer, state, state_missing, ratio}, "albedo", ...
      axes=axes(fig), names=["PROMICE", "MAR 3.11", "MERRA-2", ...
      "RACMO 2.3p3"], frequency="daily", ...
      aggregation="shortwave_weighted_mean", ...
      albedo_source=["radiometer", "model_state", "model_state", ...
      "model_ratio"]);

   radiometer_daily = returned.lines(1).YData(:);
   testCase.verifyTrue(isnan(radiometer_daily(1)));
   testCase.verifyEqual(radiometer_daily(2), ...
      (0.7 * 10 + 0.9 * 500) / 510, ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(returned.lines(2).YData(:), [0.82; 0.8], ...
      'AbsTol', 1e-12);
   testCase.verifyTrue(isnan(returned.lines(3).YData(1)));
   testCase.verifyEqual(returned.lines(3).YData(2), 0.9, ...
      'AbsTol', 1e-12);
   testCase.verifyTrue(isnan(returned.lines(4).YData(1)));
   testCase.verifyEqual(returned.lines(4).YData(2), ...
      (0.7 * 10 + 0.9 * 500) / 510, 'AbsTol', 1e-12);
end

function test_compare_timeseries_rejects_mismatched_albedo_sources(testCase)
   % Source-role vectors must align one-for-one with their input timetables.
   Time = datetime(2012, 1, 1, TimeZone="UTC");
   T = timetable(0.8, 'RowTimes', Time, 'VariableNames', {'albedo'});
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   testCase.verifyError(@() icemodel.plot.compareTimeseries({T, T}, ...
      "albedo", axes=axes(fig), ...
      albedo_source=["radiometer", "model_state", "model_ratio"]), ...
      'icemodel:plot:compareTimeseries:badAlbedoSource');
end

function test_compare_timeseries_uses_circular_daily_wind_direction(testCase)
   % Daily wind direction must average across the 360/0 wrap.
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:23)');
   wdir = repmat([350; 10], 12, 1);
   T = timetable(wdir, 'RowTimes', time);
   T.Properties.VariableUnits = {'degrees'};
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries(T, "wdir", axes=axes(fig), ...
      frequency="daily", aggregation="mean");

   testCase.verifyEqual(returned.lines(1).YData, 360, 'AbsTol', 1e-12);
end

function test_compare_timeseries_uses_circular_omitmissing_wind_direction(testCase)
   % Omitting missing directions must not revert the circular variable to a
   % scalar arithmetic mean.
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:23)');
   wdir = repmat([350; 10], 12, 1);
   wdir(3:4) = NaN;
   T = timetable(wdir, 'RowTimes', time);
   T.Properties.VariableUnits = {'degrees'};
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries(T, "wdir", axes=axes(fig), ...
      frequency="daily", aggregation="mean_omitmissing");

   testCase.verifyEqual(returned.lines(1).YData, 360, 'AbsTol', 1e-12);
end

function test_compare_timeseries_masks_undefined_circular_mean(testCase)
   % Opposing direction vectors have no physically defined daily mean.
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:23)');
   wdir = repmat([90; 270], 12, 1);
   T = timetable(wdir, 'RowTimes', time);
   T.Properties.VariableUnits = {'degrees'};
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries(T, "wdir", axes=axes(fig), ...
      frequency="daily", aggregation="mean");

   testCase.verifyFalse(returned.plotted);
end

function test_compare_timeseries_speed_weights_daily_wind_direction(testCase)
   % Daily mean direction follows the mean wind vector when speed is present.
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:23)');
   wdir = repmat([360; 180], 12, 1);
   wspd = repmat([10; 1], 12, 1);
   T = timetable(wdir, wspd, 'RowTimes', time);
   T.Properties.VariableUnits = {'degrees', 'm s-1'};
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries(T, "wdir", axes=axes(fig), ...
      frequency="daily", aggregation="mean");

   testCase.verifyEqual(returned.lines(1).YData, 360, 'AbsTol', 1e-12);
end

function test_compare_timeseries_masks_partial_boundary_days(testCase)
   % A complete interior day remains usable while both clipped boundaries mask.
   time = datetime(2012, 1, 1, 12, 0, 0, TimeZone="UTC") ...
      + hours((0:47)');
   T = timetable(ones(size(time)), 'RowTimes', time, ...
      'VariableNames', {'tair'});
   T.Properties.VariableUnits = {'degC'};
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries(T, "tair", axes=axes(fig), ...
      frequency="daily", aggregation="mean");

   testCase.verifyEqual(numel(returned.lines(1).YData), 3);
   testCase.verifyTrue(isnan(returned.lines(1).YData(1)));
   testCase.verifyEqual(returned.lines(1).YData(2), 1);
   testCase.verifyTrue(isnan(returned.lines(1).YData(3)));
end

function test_compare_timeseries_rejects_bad_dense_days(testCase)
   % Finite counts cannot hide a duplicate-plus-omitted or shifted posting.
   origin = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC");
   regular = origin + hours((0:47)');

   duplicate = regular;
   duplicate(6) = duplicate(5);
   fig_duplicate = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_duplicate));
   duplicate_table = timetable(ones(size(duplicate)), ...
      'RowTimes', duplicate, 'VariableNames', {'tair'});
   % Isolate time-axis validation from metadata-driven Kelvin conversion.
   duplicate_table.Properties.VariableUnits = {'degC'};
   duplicate_out = icemodel.plot.compareTimeseries( ...
      duplicate_table, "tair", axes=axes(fig_duplicate), ...
      frequency="daily", aggregation="mean");
   testCase.verifyTrue(isnan(duplicate_out.lines(1).YData(1)));
   testCase.verifyEqual(duplicate_out.lines(1).YData(2), 1);

   % A half-hour shift preserves the row count but violates the inferred grid.
   offgrid = regular;
   offgrid(6) = offgrid(6) + minutes(30);
   fig_offgrid = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_offgrid));
   offgrid_table = timetable(ones(size(offgrid)), ...
      'RowTimes', offgrid, 'VariableNames', {'tair'});
   % Use the same explicit temperature unit on both malformed-axis fixtures.
   offgrid_table.Properties.VariableUnits = {'degC'};
   offgrid_out = icemodel.plot.compareTimeseries( ...
      offgrid_table, "tair", axes=axes(fig_offgrid), ...
      frequency="daily", aggregation="mean");
   testCase.verifyTrue(isnan(offgrid_out.lines(1).YData(1)));
   testCase.verifyEqual(offgrid_out.lines(1).YData(2), 1);
end

function test_compare_timeseries_reversed_dense_axis_fails_closed(testCase)
   % Reversing an otherwise complete dense axis cannot yield plausible means.
   time = flip(datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:23)'));
   T = timetable(ones(size(time)), 'RowTimes', time, ...
      'VariableNames', {'tair'});
   T.Properties.VariableUnits = {'degC'};
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries(T, "tair", axes=axes(fig), ...
      frequency="daily", aggregation="mean");

   testCase.verifyFalse(any(returned.plotted));
   testCase.verifyEmpty(returned.lines);
end

function test_compare_timeseries_infers_dense_support_before_filter(testCase)
   % Filtering a known dense source to one posting does not make a complete day.
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:47)');
   T = timetable(ones(size(time)), 'RowTimes', time, ...
      'VariableNames', {'tair'});
   T.Properties.VariableUnits = {'degC'};
   selected = time(31);
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));

   returned = icemodel.plot.compareTimeseries(T, "tair", axes=axes(fig), ...
      frequency="daily", aggregation="mean", ...
      startdate=selected, enddate=selected);

   testCase.verifyFalse(any(returned.plotted));
   testCase.verifyEmpty(returned.lines);
end

function test_compare_timeseries_empty_accumulators_fail_closed(testCase)
   % Empty finite/slot index sets must yield no line rather than error.
   time = datetime(2012, 1, 1, 0, 0, 0, TimeZone="UTC") ...
      + hours((0:47)');
   values = [ones(24, 1); nan(24, 1)];
   T = timetable(values, 'RowTimes', time, 'VariableNames', {'tair'});
   testCase.verifyEqual(string(T.Properties.DimensionNames{1}), "Time");
   T.Properties.VariableUnits = {'degC'};
   day_two_start = datetime(2012, 1, 2, 0, 0, 0, TimeZone="UTC");
   day_two_end = datetime(2012, 1, 2, 23, 59, 59, TimeZone="UTC");
   fig_missing = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_missing));

   missing = icemodel.plot.compareTimeseries(T, "tair", ...
      axes=axes(fig_missing), frequency="daily", aggregation="mean", ...
      startdate=day_two_start, enddate=day_two_end);
   testCase.verifyFalse(any(missing.plotted));
   testCase.verifyEmpty(missing.lines);

   % A filtered day whose every posting misses the dominant native phase leaves
   % the slot accumulator empty and must fail closed by the same contract.
   shifted = T;
   shifted_time = shifted.Time;
   shifted_time(25:48) = shifted_time(25:48) + minutes(10);
   shifted.Time = shifted_time;
   shifted.tair(:) = 1;
   fig_shifted = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_shifted));
   offgrid = icemodel.plot.compareTimeseries(shifted, "tair", ...
      axes=axes(fig_shifted), frequency="daily", aggregation="mean", ...
      startdate=day_two_start, enddate=day_two_end);
   testCase.verifyFalse(any(offgrid.plotted));
   testCase.verifyEmpty(offgrid.lines);
end

function test_compare_timeseries_preserves_sparse_and_skips_empty(testCase)
   % Sparse singleton amounts/means stay native; empty/all-missing inputs skip.
   time = datetime(2012, 1, 1, 6, 0, 0, TimeZone="UTC");
   single = timetable(-5, 0.2, 'RowTimes', time, ...
      'VariableNames', {'tair', 'smb'});
   single.Properties.VariableUnits = {'degC', 'm w.e.'};

   fig_mean = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_mean));
   mean_out = icemodel.plot.compareTimeseries(single, "tair", ...
      axes=axes(fig_mean), frequency="daily", aggregation="mean");
   testCase.verifyEqual(mean_out.lines(1).XData, time);
   testCase.verifyEqual(mean_out.lines(1).YData, -5);

   % A sparse amount is already an interval value, even without span metadata.
   fig_total = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_total));
   total_out = icemodel.plot.compareTimeseries(single, "smb", ...
      axes=axes(fig_total), frequency="daily", aggregation="daily_total");
   testCase.verifyEqual(total_out.lines(1).XData, time);
   testCase.verifyEqual(total_out.lines(1).YData, 0.2);

   % Neither an empty table nor an all-missing channel creates a plot handle.
   fig_empty = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig_empty));
   empty_out = icemodel.plot.compareTimeseries(single([], :), "tair", ...
      axes=axes(fig_empty), frequency="daily");
   testCase.verifyFalse(any(empty_out.plotted));
   missing = single;
   missing.tair(:) = NaN;
   missing_out = icemodel.plot.compareTimeseries(missing, "tair", ...
      axes=axes(fig_empty), frequency="daily");
   testCase.verifyFalse(any(missing_out.plotted));
end

function test_profile_filters_to_date_window(testCase)
   % Profile plots filter dated profile rows so users can inspect one year or
   % one campaign inside a longer staged observation table.

   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   T = table([1; 2], [320; 330], ...
      [datetime(2012, 1, 1, 'TimeZone', 'UTC'); ...
      datetime(2012, 1, 2, 'TimeZone', 'UTC')], ...
      'VariableNames', {'depth', 'density', 'datetime'});

   returned = icemodel.plot.profile(T, "density", axes=ax, names="density", ...
      startdate=datetime(2012, 1, 2, 'TimeZone', 'UTC'), ...
      enddate=datetime(2012, 1, 2, 'TimeZone', 'UTC'));

   testCase.verifyEqual(numel(returned.lines(1).XData), 1);
   testCase.verifyEqual(returned.lines(1).XData, 330);
end

function test_profile_handles_empty_dated_table(testCase)
   % Empty profile tables should be skipped before date-window indexing.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   T = table(zeros(0, 1), zeros(0, 1), NaT(0, 1, 'TimeZone', 'UTC'), ...
      'VariableNames', {'depth', 'density', 'datetime'});

   returned = icemodel.plot.profile(T, "density", axes=ax, names="density", ...
      startdate=datetime(2012, 1, 1, 'TimeZone', 'UTC'), ...
      enddate=datetime(2012, 1, 2, 'TimeZone', 'UTC'));

   testCase.verifyFalse(any(returned.plotted));
end

function test_profile_preserves_hold_state(testCase)
   % Profile plotting should not change a caller-composed axes hold state.
   fig = figure('Visible', 'off');
   testCase.addTeardown(@() close(fig));
   ax = axes(fig);
   hold(ax, 'on')
   T = table(1, 320, datetime(2012, 1, 1, 'TimeZone', 'UTC'), ...
      'VariableNames', {'depth', 'density', 'datetime'});

   icemodel.plot.profile(T, "density", axes=ax, names="density");

   testCase.verifyTrue(ishold(ax));
end

function writeTinyFirnTree(root)
   %WRITETINYFIRNTREE Create a minimal manifest + MAT artifact tree.
   eval_case = fullfile(root, 'eval', 'promice', 'kanm');
   met_dir = fullfile(root, 'input', 'met', 'promice');
   ud_dir = fullfile(root, 'input', 'userdata', 'promice');
   mar_met_dir = fullfile(root, 'input', 'met', 'mar3.11');
   mar_ud_dir = fullfile(root, 'input', 'userdata', 'mar3.11');
   racmo_ud_dir = fullfile(root, 'input', 'userdata', 'racmo2.3p3');
   icemodel.helpers.ensureDirExists(eval_case);
   icemodel.helpers.ensureDirExists(met_dir);
   icemodel.helpers.ensureDirExists(ud_dir);
   icemodel.helpers.ensureDirExists(mar_met_dir);
   icemodel.helpers.ensureDirExists(mar_ud_dir);
   icemodel.helpers.ensureDirExists(racmo_ud_dir);

   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') + hours(0:47))';
   n = numel(Time);
   % One sparse step proves screening internals stay staged but unplotted.
   step_magnitude = zeros(n, 1);
   step_magnitude(12) = 0.05;
   obs = timetable(Time, linspace(0, 0.3, n)', ...
      263.15 + linspace(0, 2, n)', 260.15 + linspace(0, 2, n)', ...
      'VariableNames', {'ablation', 'tice10m', 'tice1'});
   targets = struct('format', 'timeseries', 'data', obs, ...
      'metadata', struct('source', 'tiny fixture'));
   targets = icemodel.verification.setup.stampArtifactMetadata(targets);
   save(fullfile(eval_case, 'observations.mat'), 'targets');

   met = timetable(Time, 263.15 + (0:n - 1)', 100 + (0:n - 1)', ...
      200 + (0:n - 1)', 80 + zeros(n, 1), 5 + zeros(n, 1), ...
      0.7 + zeros(n, 1), zeros(n, 1), 70000 + zeros(n, 1), ...
      'VariableNames', {'tair', 'swd', 'lwd', 'rh', 'wspd', 'albedo', ...
      'ppt', 'psfc'});
   met = icemodel.forcing.helpers.stampMetadata(met);
   save(fullfile(met_dir, 'met_kanm_promice_20120101_20120102_1hr.mat'), ...
      'met');

   Data = timetable(Time, linspace(0, 0.3, n)', step_magnitude, ...
      263.15 + linspace(0, 2, n)', 260.15 + linspace(0, 2, n)', ...
      263.15 + linspace(0, 2, n)', zeros(n, 1), ones(n, 1), ...
      linspace(0, 0.01, n)', linspace(0, 0.005, n)', ...
      linspace(0, 0.002, n)', 0.7 + zeros(n, 1), 180 + zeros(n, 1), ...
      'VariableNames', {'ablation', 'step_magnitude', 'tice10m', 'tice1', ...
      'tice10m_source', 'tice10m_qc_flag', 'dtice1', 'melt', ...
      'runoff', 'refreeze', 'albedo', 'wdir'});
   Data = icemodel.forcing.helpers.stampMetadata(Data, strict=false);
   save(fullfile(ud_dir, 'kanm_promice_20120101_20120102.mat'), 'Data');

   met = timetable(Time, 260.15 + (0:n - 1)', 90 + (0:n - 1)', ...
      190 + (0:n - 1)', 75 + zeros(n, 1), 4 + zeros(n, 1), ...
      0.72 + zeros(n, 1), zeros(n, 1), 70500 + zeros(n, 1), ...
      'VariableNames', {'tair', 'swd', 'lwd', 'rh', 'wspd', 'albedo', ...
      'ppt', 'psfc'});
   met = icemodel.forcing.helpers.stampMetadata(met);
   save(fullfile(mar_met_dir, ...
      'met_kanm_mar3.11_20120101_20120102_1hr.mat'), 'met');

   Data = timetable(Time, linspace(0, 0.1, n)', ...
      linspace(-2e-4, 2e-4, n)', ...
      [repmat(-1e-4, 24, 1); repmat(2e-4, 24, 1)], ...
      [repmat(1e-4, 24, 1); repmat(2e-4, 24, 1)], ...
      100 + zeros(n, 1), -20 + zeros(n, 1), ...
      200 + zeros(n, 1), -180 + zeros(n, 1), ...
      20 + zeros(n, 1), -40 + zeros(n, 1), 10 + zeros(n, 1), ...
      5 + zeros(n, 1), linspace(-0.01, 0.01, n)', ...
      0.71 + zeros(n, 1), ...
      'VariableNames', {'smb', 'subl', 'subl_evap', ...
      'refreeze_deposition', 'swd', 'swu', 'lwd', 'lwu', 'swn', ...
      'lwn', 'shf', 'lhf', 'sndiv', 'modis'});
   Data = icemodel.forcing.helpers.stampMetadata(Data, strict=false);
   save(fullfile(mar_ud_dir, 'kanm_mar3.11_20120101_20120102.mat'), ...
      'Data');

   % RACMO carries both distinct meltwater component channels and a second
   % replicated MODIS column so the visual-QA fixture exercises deduplication.
   Data = timetable(Time, linspace(0, 0.01, n)', ...
      linspace(0, 0.006, n)', linspace(0, 0.002, n)', ...
      linspace(-0.01, 0.01, n)', 0.71 + zeros(n, 1), ...
      'VariableNames', {'melt', 'meltin', 'subl', 'sndiv', 'modis'});
   Data = icemodel.forcing.helpers.stampMetadata(Data, strict=false);
   save(fullfile(racmo_ud_dir, ...
      'kanm_racmo2.3p3_20120101_20120102.mat'), 'Data');

   colocation = struct();
   colocation.promice = struct('staged', true, ...
      'met_files', "promice/met_kanm_promice_20120101_20120102_1hr.mat", ...
      'data_files', "promice/kanm_promice_20120101_20120102.mat", ...
      'forcing_ready', true);
   colocation.mar = struct('staged', true, 'source_id', "mar3.11", ...
      'met_files', "mar3.11/met_kanm_mar3.11_20120101_20120102_1hr.mat", ...
      'data_files', "mar3.11/kanm_mar3.11_20120101_20120102.mat");
   colocation.racmo = struct('staged', true, 'source_id', "racmo2.3p3", ...
      'data_files', ...
      "racmo2.3p3/kanm_racmo2.3p3_20120101_20120102.mat");
   values = { ...
      'kanm'
      'firn_observational'
      'KAN_M'
      'KAN_M'
      'ablation'
      {'firn'}
      'none'
      struct('lat_wgs84', 67, 'lon_wgs84', -50, ...
         'x_epsg3413', 0, 'y_epsg3413', 0, 'elev_m', 1000)
      struct('start', '2012-01-01 00:00:00', ...
         'end', '2012-01-02 23:00:00')
      'kanm/observations.mat'
      {'promice', 'mar3.11', 'racmo2.3p3'}
      {'promice_obs', 'mar3.11', 'racmo2.3p3'}
      {'ablation', 'tice10m'}
      struct('present_groups', "tiny")
      colocation
      '1hr'
      'tiny visual QA fixture'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "promice", "", "", "tiny", "today", entry);
   icemodel.verification.setup.writeManifest( ...
      fullfile(root, 'eval', 'promice', 'manifest.json'), manifest);
end

function writeTinyEsmTree(root)
   %WRITETINYESMTREE Create a minimal legacy-schema ESM-SnowMIP fixture.
   eval_case = fullfile(root, 'eval', 'esm_snowmip', 'cdp');
   met_dir = fullfile(root, 'input', 'met', 'esm_snowmip');
   icemodel.helpers.ensureDirExists(eval_case);
   icemodel.helpers.ensureDirExists(met_dir);

   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') + hours(0:47))';
   n = numel(Time);
   depth = linspace(0.1, 0.2, n)';
   obs = timetable(Time, depth, depth + 0.01, depth + 0.02, ...
      100 + (0:n - 1)', -5 + zeros(n, 1), -4 + zeros(n, 1), ...
      -3 + zeros(n, 1), -2 + zeros(n, 1), -1 + zeros(n, 1), ...
      'VariableNames', {'snow_depth_m', 'snd_auto_m', 'snd_man_m', ...
      'swe_kg_m2', 'surface_temp_C', 'soil_temp_1_C', 'soil_temp_2_C', ...
      'soil_temp_3_C', 'soil_temp_4_C'});
   targets = struct('format', 'timeseries', 'data', obs, ...
      'metadata', struct('source', 'tiny esm fixture'));
   targets = icemodel.verification.setup.stampArtifactMetadata(targets);
   save(fullfile(eval_case, 'observations.mat'), 'targets');

   met = timetable(Time, 263.15 + (0:n - 1)', 100 + (0:n - 1)', ...
      200 + (0:n - 1)', zeros(n, 1), ...
      'VariableNames', {'tair', 'swd', 'lwd', 'ppt'});
   met = icemodel.forcing.helpers.stampMetadata(met);
   save(fullfile(met_dir, 'met_cdp_esm_snowmip_20120101_20120102_1hr.mat'), ...
      'met');
   save(fullfile(root, 'input', 'met', ...
      'met_cdp_esm_snowmip_20120101_20120102_1hr.mat'), 'met');

   values = { ...
      'cdp'
      'esm_site'
      'cdp'
      'Col de Porte'
      'land'
      {'seasonal_snow'}
      'none'
      char(fullfile('cdp', 'observations.mat'))
      ''
      'hourly'
      struct('start', '2012-01-01 00:00:00', ...
      'end', '2012-01-02 23:00:00')
      {'snow_depth_m', 'snd_auto_m', 'snd_man_m', 'swe_kg_m2', ...
      'surface_temp_C', 'soil_temp_1_C', 'soil_temp_2_C', ...
      'soil_temp_3_C', 'soil_temp_4_C'}
      struct('obs_file', char(fullfile('cdp', 'observations.mat')))
      'tiny visual QA ESM fixture'};
   entry = icemodel.verification.setup.makeCaseManifestEntry(values);
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "esm_snowmip", "", "", "tiny", "today", entry);
   icemodel.verification.setup.writeManifest( ...
      fullfile(root, 'eval', 'esm_snowmip', 'manifest.json'), manifest);
end

function writeTinyLaughTree(root)
   %WRITETINYLAUGHTREE Create a sparse dual-source Colbeck fixture.
   eval_case = fullfile(root, 'eval', 'laugh_tests', 'colbeck1976');
   icemodel.helpers.ensureDirExists(eval_case);

   % Three experiments with only a few samples reproduce the layout problem
   % without loading or computing the production Colbeck benchmark.
   Time = datetime(1990, 1, 1, 0, [1; 30; 59], 0, 'TimeZone', 'UTC');
   experiment_names = ["exp1", "exp2", "exp3"];
   numerical_experiments = struct();
   analytical_experiments = struct();
   for k = 1:numel(experiment_names)
      name = char(experiment_names(k));
      numerical_experiments.(name) = timetable(Time, ...
         k * [0.001; 0.003; 0.005], k * [0; 1e-6; 2e-6], ...
         'VariableNames', {'snow_liquid_water_storage_m', ...
         'bottom_outflow_mps'});
      analytical_experiments.(name) = timetable(Time, ...
         k * [0.0011; 0.0031; 0.0051], k * [0; 0.9e-6; 1.9e-6], ...
         'VariableNames', {'snow_liquid_water_storage_m', ...
         'bottom_outflow_mps'});
   end

   % Keep the production multi-source target shape so plotcase exercises its
   % existing Colbeck target selection and experiment-grid comparison path.
   numerical_summa = struct('format', 'experiment_bundle', ...
      'experiments', numerical_experiments);
   analytical_clark2017 = struct('format', 'experiment_bundle', ...
      'experiments', analytical_experiments);
   targets = struct('numerical_summa', numerical_summa, ...
      'analytical_clark2017', analytical_clark2017);
   reference = analytical_clark2017;
   save(fullfile(eval_case, 'evaluation.mat'), 'targets');
   save(fullfile(eval_case, 'reference.mat'), 'reference');

   % The canonical Laugh manifest fields let the visual-QA entry point use the
   % same listcases/loadmanifest path as production staging.
   values = { ...
      'colbeck1976'
      'synthetic_process'
      'colbeck1976'
      'Colbeck 1976 tiny visual-QA fixture'
      ''
      strings(0, 1)
      ''
      'colbeck1976/evaluation.mat'
      'colbeck1976/reference.mat'
      'sparse test output'
      struct('start', '1990-01-01 00:01:00', ...
      'end', '1990-01-01 00:59:00')
      {'snow_liquid_water_storage_m', 'bottom_outflow_mps'}
      struct('snow_liquid_water_storage_m', 'storage', ...
      'bottom_outflow_mps', 'outflow')
      'tiny Laugh-Test visual QA fixture'};
   entry = icemodel.verification.setup.makeCaseManifestEntry(values);
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "laugh_tests", "", "", "tiny", "today", entry);
   icemodel.verification.setup.writeManifest( ...
      fullfile(root, 'eval', 'laugh_tests', 'manifest.json'), manifest);
end

function writeEmptyResearchTree(root)
   %WRITEEMPTYRESEARCHTREE Create a manifest row with no loadable artifacts.
   icemodel.helpers.ensureDirExists(fullfile(root, 'eval', 'research_site'));
   values = { ...
      'empty'
      'firn_observational'
      'EMPTY'
      'EMPTY'
      'unknown'
      {'firn'}
      'none'
      struct('lat_wgs84', NaN, 'lon_wgs84', NaN, ...
         'x_epsg3413', NaN, 'y_epsg3413', NaN, 'elev_m', NaN)
      struct('start', '', 'end', '')
      'empty/missing_observations.mat'
      {}
      {}
      {}
      struct()
      struct()
      'none'
      'empty visual QA fixture'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "research_site", "", "", "tiny", "today", entry);
   icemodel.verification.setup.writeManifest( ...
      fullfile(root, 'eval', 'research_site', 'manifest.json'), manifest);
end

function met = makePlotMet()
   %MAKEPLOTMET Build a two-day met timetable for plotting-helper tests.
   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') + hours(0:47))';
   n = numel(Time);
   met = timetable(Time, 263.15 + (0:n - 1)', 100 + (0:n - 1)', ...
      200 + (0:n - 1)', zeros(n, 1), ...
      'VariableNames', {'tair', 'swd', 'lwd', 'ppt'});
end
