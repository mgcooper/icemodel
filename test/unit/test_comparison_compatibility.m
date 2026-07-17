function tests = test_comparison_compatibility
   %TEST_COMPARISON_COMPATIBILITY Verify additive compatibility discovery.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install the canonical test path and allocate an isolated eval/input tree.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.tmp = tempname;
   mkdir(testCase.TestData.tmp)
end

function teardown(testCase)
   % Remove synthetic manifests and tiny artifacts.
   if isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
   clear testCase.TestData.cleanup
end

function test_reports_artifact_pair_classes(testCase)
   % A mixed-source manifest should report data-data, model-data, and
   % model-model pairs from staged artifact/header variables.
   eval_root = writeCompatibilityTree(testCase.TestData.tmp);

   returned = icemodel.verification.comparisonCompatibility("all", ...
      evaluation_data_root=eval_root, dataset_family="research_site");

   report = returned.cases;
   testCase.verifyEqual(string(report.case_id), "all");
   expected_sources = ["promice_obs", "reference", "promice", "retmip", ...
      "imau", "sumup", "research_site", "mar3.11", "merra2", ...
      "racmo2.3p3"];
   sources = unique(string({report.artifacts.source}), 'stable');
   for source = expected_sources
      testCase.verifyTrue(any(sources == source), source);
   end

   classes = unique(string({returned.pairs.class}), 'stable');
   testCase.verifyTrue(any(classes == "data-data"));
   testCase.verifyTrue(any(classes == "model-data"));
   testCase.verifyTrue(any(classes == "model-model"));
   testCase.verifyTrue(hasPair(returned.pairs, "promice", "imau", ...
      "data-data", "smb"));
   testCase.verifyTrue(hasPair(returned.pairs, "promice", "mar3.11", ...
      "model-data", "smb"));
   testCase.verifyTrue(hasPair(returned.pairs, "promice_obs", "reference", ...
      "model-data", "tsfc"));
   testCase.verifyTrue(hasPair(returned.pairs, "retmip", "mar3.11", ...
      "model-model", "subsurface_temperature"));
   testCase.verifyTrue(hasPair(returned.pairs, "promice", "retmip", ...
      "model-data", "subsurface_temperature"));
   testCase.verifyTrue(hasPair(returned.pairs, "mar3.11", "merra2", ...
      "model-model", "smb"));
end

function test_keeps_manifest_declarations_additive(testCase)
   % Manifest fields remain complete declarations even when artifact headers
   % provide a narrower set of discovered variables.
   eval_root = writeCompatibilityTree(testCase.TestData.tmp);

   returned = icemodel.verification.comparisonCompatibility("all", ...
      evaluation_data_root=eval_root, dataset_family="research_site");
   report = returned.cases;

   testCase.verifyEqual(report.declared.comparison_variables, ...
      ["smb"; "density"; "tice1"]);
   testCase.verifyEqual(report.declared.forcing_sources, ...
      ["mar3.11"; "merra2"]);
   testCase.verifyTrue(any(report.declared.eval_sources == "retmip_protocol"));

   promice = firstArtifact(report.artifacts, "promice");
   testCase.verifyEqual(string(promice.evidence), "mat_payload");
   testCase.verifyTrue(any(promice.artifact_variables == "smb"));
   testCase.verifyTrue(any(promice.artifact_variables ...
      == "subsurface_temperature"));
   testCase.verifyTrue(any(promice.declared_variables == "density"));
   testCase.verifyEqual(string(promice.variables(:)), ...
      ["smb"; "subsurface_temperature"; "density"]);

   retmip = firstArtifact(report.artifacts, "retmip");
   testCase.verifyEqual(string(retmip.evidence), "netcdf_header");
   testCase.verifyTrue(any(retmip.artifact_variables ...
      == "subsurface_temperature"));
end

function test_resolves_explicit_input_root(testCase)
   % Separate eval/input roots should still resolve relative met/userdata paths.
   input_root = fullfile(testCase.TestData.tmp, 'standalone-input');
   eval_root = writeCompatibilityTree( ...
      fullfile(testCase.TestData.tmp, 'standalone-eval'), input_root);

   returned = icemodel.verification.comparisonCompatibility("all", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      dataset_family="research_site");

   testCase.verifyTrue(hasPair(returned.pairs, "promice", "mar3.11", ...
      "model-data", "smb"));
end

function test_loadmanifest_carries_explicit_input_root_to_colocated_loader(testCase)
   % Runtime helpers should reuse the root resolved by loadmanifest.
   input_root = fullfile(testCase.TestData.tmp, 'loader-input');
   eval_root = writeCompatibilityTree( ...
      fullfile(testCase.TestData.tmp, 'loader-eval'), input_root);

   manifest = icemodel.verification.loadmanifest("all", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      dataset_family="research_site");
   bundle = icemodel.verification.helpers.loadColocatedData(manifest, ...
      "promice");

   testCase.verifyEqual(string(manifest.input_data_root), string(input_root));
   testCase.verifyEqual(string(bundle.format), "timeseries");
   testCase.verifyGreaterThan(height(bundle.data), 0);
end

function test_colocated_loader_filters_to_leg_window(testCase)
   % Broad cached RCM/Data artifacts should expose only the staged overlap.
   input_root = fullfile(testCase.TestData.tmp, 'window-input');
   eval_root = writeCompatibilityTree( ...
      fullfile(testCase.TestData.tmp, 'window-eval'), input_root);
   manifest_file = fullfile(eval_root, 'research_site', 'manifest.json');
   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.colocation.promice.window = ...
      struct('start', '2000-01-01 01:00:00', ...
      'end', '2000-01-01 01:00:00');
   writeJson(manifest_file, manifest);

   manifest = icemodel.verification.loadmanifest("all", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      dataset_family="research_site");
   bundle = icemodel.verification.helpers.loadColocatedData(manifest, ...
      "promice");

   testCase.verifyEqual(height(bundle.data), 1);
   testCase.verifyEqual(bundle.data.Time, datetime(2000, 1, 1, 1, 0, 0));
end

function test_loadmanifest_normalizes_manifest_lists_to_columns(testCase)
   % JSON-decoded arrays should be column-oriented for visual inspection.
   eval_root = writeCompatibilityTree(testCase.TestData.tmp);

   manifest = icemodel.verification.loadmanifest("all", ...
      evaluation_data_root=eval_root, dataset_family="research_site");

   testCase.verifySize(manifest.comparison_variables, [3 1]);
   testCase.verifySize(manifest.forcing_sources, [2 1]);
   testCase.verifySize(manifest.observation_variables.native_variables, [2 1]);
   testCase.verifyGreaterThan(size(manifest.eval_sources, 1), 1);
   testCase.verifyEqual(size(manifest.eval_sources, 2), 1);
end

function test_comparecase_uses_explicit_input_root_for_default_candidate(testCase)
   % The no-candidate smoke path should find RCM userdata under input_data_root.
   input_root = fullfile(testCase.TestData.tmp, 'compare-input');
   eval_root = writeCompatibilityTree( ...
      fullfile(testCase.TestData.tmp, 'compare-eval'), input_root);
   writeTargetsPayload(fullfile(eval_root, 'research_site', 'all', ...
      'observations.mat'), "smb");
   manifest_file = fullfile(eval_root, 'research_site', 'manifest.json');
   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.reference_file = '';
   writeJson(manifest_file, manifest);

   result = icemodel.verification.comparecase("all", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      dataset_family="research_site", make_plot=false);

   testCase.verifyTrue(any(string(result.metrics.status) == "ok"));
end

function test_comparecase_uses_available_rcm_candidate_when_racmo_absent(testCase)
   % A MAR-only case should not default to an empty RACMO candidate.
   input_root = fullfile(testCase.TestData.tmp, 'mar-only-input');
   eval_root = writeCompatibilityTree( ...
      fullfile(testCase.TestData.tmp, 'mar-only-eval'), input_root);
   writeTargetsPayload(fullfile(eval_root, 'research_site', 'all', ...
      'observations.mat'), "smb");
   manifest_file = fullfile(eval_root, 'research_site', 'manifest.json');
   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.reference_file = '';
   manifest.cases.colocation = rmfield(manifest.cases.colocation, 'racmo');
   manifest.cases.eval_sources = {'promice_obs', 'mar3.11'};
   manifest.cases.forcing_sources = {'mar3.11'};
   writeJson(manifest_file, manifest);

   result = icemodel.verification.comparecase("all", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      dataset_family="research_site", make_plot=false);

   ok = string(result.metrics.variable) == "smb" ...
      & string(result.metrics.status) == "ok";
   testCase.verifyTrue(any(ok));
end

function test_colocated_loader_combines_profile_sidecar_with_userdata(testCase)
   % A MAR profile sidecar becomes the primary candidate while exact scalar
   % userdata channels remain available for soft comparisons.
   input_root = fullfile(testCase.TestData.tmp, 'profile-candidate-input');
   userdata_root = fullfile(input_root, 'userdata', 'mar3.11');
   mkdir(userdata_root)
   Time = datetime(2020, 1, 1:2, TimeZone="UTC")';
   Data = timetable(Time, [0.1; 0.2], 'VariableNames', {'smb'});
   save(fullfile(userdata_root, 'case_data.mat'), 'Data')

   profile_time = datetime(2020, 1, 1, 12, 0, 0, TimeZone="UTC");
   density = table(repmat("mar_case_20200101", 2, 1), ...
      repmat(profile_time, 2, 1), [0; 1], [300; 500], ...
      'VariableNames', {'profile_id', 'datetime', 'depth', 'density'});
   reference = struct('format', 'subsurface_profile_bundle', ...
      'data', struct('density', density), 'metadata', struct());
   save(fullfile(userdata_root, 'case_profiles.mat'), 'reference')

   manifest = struct('input_data_root', input_root, ...
      'eval_sources', {{'mar3.11'}}, ...
      'colocation', struct('mar', struct('staged', true, ...
      'source_id', 'mar3.11', 'data_files', 'mar3.11/case_data.mat', ...
      'model_output_files', 'mar3.11/case_profiles.mat')));
   bundle = icemodel.verification.helpers.loadColocatedData(manifest, "mar3.11");
   resolved = icemodel.verification.helpers.resolveCandidateBundle(manifest);

   testCase.verifyEqual(string(bundle.format), "subsurface_profile_bundle")
   testCase.verifyEqual(bundle.data.density, density)
   testCase.verifyEqual(bundle.data.smb, Data(:, 'smb'))
   testCase.verifyEqual(string(resolved.format), "subsurface_profile_bundle")
   testCase.verifyEqual(resolved.data.density, density)
end

function test_met_only_rcm_leg_is_not_comparison_source(testCase)
   % RCM met files are forcing artifacts until met-channel comparison exists.
   eval_root = writeCompatibilityTree(testCase.TestData.tmp);
   manifest_file = fullfile(eval_root, 'research_site', 'manifest.json');
   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.colocation.mar = stagedLeg('met_files', {'mar_met.mat'});
   manifest.cases.forcing_sources = {'mar3.11'};
   manifest.cases.eval_sources = {'promice_obs'};
   writeJson(manifest_file, manifest);

   returned = icemodel.verification.comparisonCompatibility("all", ...
      evaluation_data_root=eval_root, dataset_family="research_site");
   sources = string({returned.cases.artifacts.source});

   testCase.verifyFalse(any(sources == "mar3.11"));
   testCase.verifyFalse(hasPair(returned.pairs, "promice_obs", ...
      "mar3.11", "model-data", "smb"));
end

function test_empty_match_returns_empty(testCase)
   % No-match filters should not return a scalar blank compatibility record.
   eval_root = writeCompatibilityTree(testCase.TestData.tmp);

   returned = icemodel.verification.comparisonCompatibility("missing", ...
      evaluation_data_root=eval_root, dataset_family="research_site");

   testCase.verifyEmpty(returned.cases);
   testCase.verifyEmpty(returned.pairs);
end

function test_missing_artifact_does_not_derive_pairs(testCase)
   % Missing declared artifacts should stay diagnostic-only and not create
   % runnable comparison pairs from manifest declarations.
   eval_root = writeCompatibilityTree(testCase.TestData.tmp);
   delete(fullfile(eval_root, 'research_site', 'all', 'reference.mat'));

   returned = icemodel.verification.comparisonCompatibility("all", ...
      evaluation_data_root=eval_root, dataset_family="research_site");
   report = returned.cases;
   reference = firstArtifact(report.artifacts, "reference");

   testCase.verifyFalse(reference.exists);
   testCase.verifyEmpty(reference.variables);
   testCase.verifyFalse(hasPair(returned.pairs, "promice_obs", ...
      "reference", "model-data", "tsfc"));
end

function test_duplicate_observation_artifact_is_collected_once(testCase)
   % A colocation obs_file that resolves to the top-level evaluation artifact
   % should not create a second source record or self-comparison pair.
   eval_root = writeCompatibilityTree(testCase.TestData.tmp);
   manifest_file = fullfile(eval_root, 'research_site', 'manifest.json');
   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.colocation.sumup.obs_file = {'all/observations.mat'};
   writeJson(manifest_file, manifest);

   returned = icemodel.verification.comparisonCompatibility("all", ...
      evaluation_data_root=eval_root, dataset_family="research_site");
   report = returned.cases;
   obs_path = string(fullfile(eval_root, 'research_site', 'all', ...
      'observations.mat'));
   artifact_paths = string({report.artifacts.path});
   artifact_kinds = string({report.artifacts.kind});

   testCase.verifyEqual(sum(artifact_paths == obs_path ...
      & ismember(artifact_kinds, ["evaluation", "observation"])), 1);
   testCase.verifyFalse(hasPair(returned.pairs, "promice_obs", ...
      "sumup", "data-data", "smb"));
end

function test_reference_payload_variables_are_discovered(testCase)
   % Reference bundles saved as a top-level `reference` struct should expose
   % their experiment timetable channels for pair derivation.
   eval_root = writeCompatibilityTree(testCase.TestData.tmp);
   case_root = fullfile(eval_root, 'research_site', 'all');
   writeTargetsPayload(fullfile(case_root, 'observations.mat'), ...
      ["snow_liquid_water_storage_m", "soil_temp_1_C"]);
   writeReferencePayload(fullfile(case_root, 'reference.mat'), ...
      ["snow_liquid_water_storage_m", "soil_temp_1_C"]);

   returned = icemodel.verification.comparisonCompatibility("all", ...
      evaluation_data_root=eval_root, dataset_family="research_site");
   report = returned.cases;
   reference = firstArtifact(report.artifacts, "reference");

   testCase.verifyEqual(string(reference.evidence), "mat_payload");
   testCase.verifyTrue(any(reference.artifact_variables ...
      == "snow_liquid_water_storage_m"));
   testCase.verifyTrue(any(reference.artifact_variables == "soil_temp_1_C"));
   testCase.verifyTrue(hasPair(returned.pairs, "promice_obs", ...
      "reference", "model-data", "snow_liquid_water_storage_m"));
   testCase.verifyTrue(hasPair(returned.pairs, "promice_obs", ...
      "reference", "model-data", "soil_temp_1_C"));
end

function test_file_bearing_colocation_leg_without_staged_is_ignored(testCase)
   % Current manifests must mark staged legs explicitly before paths are usable.
   eval_root = writeCompatibilityTree(testCase.TestData.tmp);
   manifest_file = fullfile(eval_root, 'research_site', 'manifest.json');
   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.colocation.promice = ...
      rmfield(manifest.cases.colocation.promice, 'staged');
   writeJson(manifest_file, manifest);

   returned = icemodel.verification.comparisonCompatibility("all", ...
      evaluation_data_root=eval_root, dataset_family="research_site");

   sources = string({returned.cases.artifacts.source});
   testCase.verifyFalse(any(sources == "promice"));
   testCase.verifyFalse(hasPair(returned.pairs, "promice", "imau", ...
      "data-data", "smb"));
end

%% Local fixture helpers
function eval_root = writeCompatibilityTree(root, input_root)
   %WRITECOMPATIBILITYTREE Build one mixed-source manifest and tiny artifacts.
   if nargin < 2
      input_root = fullfile(root, 'input');
   end
   eval_root = fullfile(root, 'eval');
   family_root = fullfile(eval_root, 'research_site');
   case_root = fullfile(family_root, 'all');
   mkdir(case_root);
   mkdir(fullfile(input_root, 'met'));
   mkdir(fullfile(input_root, 'userdata'));

   obs_file = fullfile(case_root, 'observations.mat');
   writeMatVars(obs_file, ["smb", "surface_temp_C"]);
   reference_file = fullfile(case_root, 'reference.mat');
   writeMatVars(reference_file, "tsfc");
   promice_file = fullfile(input_root, 'userdata', 'promice_data.mat');
   writeDataPayload(promice_file, ["smb", "tice1"]);
   imau_file = fullfile(input_root, 'userdata', 'imau_data.mat');
   writeDataPayload(imau_file, "smb");
   sumup_file = fullfile(input_root, 'userdata', 'sumup_obs.mat');
   writeDataPayload(sumup_file, "smb");
   research_file = fullfile(input_root, 'userdata', 'research_obs.mat');
   writeDataPayload(research_file, "smb");
   retmip_file = fullfile(input_root, 'userdata', 'retmip_output.nc');
   writeNcVars(retmip_file, ["T_ice", "Depth"]);
   mar_file = fullfile(input_root, 'userdata', 'mar_data.mat');
   writeDataPayload(mar_file, ["smb", "subsurface_temperature"]);
   merra_file = fullfile(input_root, 'userdata', 'merra_data.mat');
   writeDataPayload(merra_file, "smb");
   racmo_file = fullfile(input_root, 'userdata', 'racmo_data.mat');
   writeDataPayload(racmo_file, "smb");

   colocation = struct();
   colocation.promice = stagedLeg('data_files', {'promice_data.mat'});
   colocation.retmip = stagedLeg('model_output_files', {'retmip_output.nc'});
   colocation.imau = stagedLeg('data_files', {'imau_data.mat'});
   colocation.sumup = stagedLeg('obs_file', {'sumup_obs.mat'});
   colocation.research_site = stagedLeg('obs_file', {'research_obs.mat'});
   colocation.mar = stagedLeg('data_files', {'mar_data.mat'});
   colocation.merra = stagedLeg('data_files', {'merra_data.mat'});
   colocation.racmo = stagedLeg('data_files', {'racmo_data.mat'});
   colocation.anchor = struct('family', 'promice', 'source_id', 'KAN_M');

   c = struct();
   c.case_id = 'all';
   c.case_type = 'firn_observational';
   c.site_id = 'all';
   c.site_name = 'all';
   c.surface_zone = 'accumulation';
   c.eval_target = {'firn'};
   c.permafrost_zone = 'none';
   c.site_location = struct('lat_wgs84', 0, 'lon_wgs84', 0, ...
      'x_epsg3413', 0, 'y_epsg3413', 0, 'elev_m', 0);
   c.period = struct('start', '2000-01-01', 'end', '2000-01-02');
   c.evaluation_file = 'all/observations.mat';
   c.reference_file = 'all/reference.mat';
   c.forcing_sources = {'mar3.11', 'merra2'};
   c.eval_sources = {'promice_obs', 'sumup_obs', 'retmip_protocol', ...
      'imau_obs', 'research_site_obs', 'mar3.11', 'merra2', ...
      'racmo2.3p3'};
   c.comparison_variables = {'smb', 'density', 'tice1'};
   c.observation_variables = struct('present_groups', {{'smb'}}, ...
      'native_variables', {{'raw_smb', 'raw_density'}});
   c.colocation = colocation;
   c.native_timestep = '1hr';
   c.notes = 'synthetic compatibility fixture';

   manifest = struct( ...
      'dataset_family', 'research_site', ...
      'source_doi', '', ...
      'source_url', '', ...
      'source_version', 'synthetic', ...
      'retrieval_date', '2026-06-28', ...
      'cases', c, ...
      'skipped', []);
   writeJson(fullfile(family_root, 'manifest.json'), manifest);
end

function leg = stagedLeg(fieldname, values)
   %STAGEDLEG Build one staged colocation leg with a file field.
   leg = struct('staged', true);
   leg.(fieldname) = values;
end

function writeMatVars(filename, names)
   %WRITEMATVARS Save tiny scalar variables with canonical names.
   data = struct();
   for name = reshape(string(names), 1, [])
      data.(char(name)) = 1;
   end
   save(filename, '-struct', 'data');
end

function writeDataPayload(filename, names)
   %WRITEDATAPAYLOAD Save a tiny production-shaped userdata MAT payload.
   Time = (datetime(2000, 1, 1, 0, 0, 0) + hours(0:1))';
   values = repmat({ones(numel(Time), 1)}, 1, numel(names));
   Data = timetable(Time, values{:}, 'VariableNames', cellstr(names));
   save(filename, 'Data');
end

function writeTargetsPayload(filename, names)
   %WRITETARGETSPAYLOAD Save a tiny observation bundle with timetable channels.
   targets = struct('format', 'timeseries', ...
      'data', channelTimetable(names), 'metadata', struct());
   save(filename, 'targets');
end

function writeReferencePayload(filename, names)
   %WRITEREFERENCEPAYLOAD Save a tiny reference bundle with experiment channels.
   reference = struct('format', 'experiment_bundle', ...
      'experiments', struct('exp1', channelTimetable(names)), ...
      'metadata', struct());
   save(filename, 'reference');
end

function tt = channelTimetable(names)
   %CHANNELTIMETABLE Build a two-row timetable for arbitrary channel names.
   Time = (datetime(2000, 1, 1, 0, 0, 0) + hours(0:1))';
   values = repmat({ones(numel(Time), 1)}, 1, numel(names));
   tt = timetable(Time, values{:}, 'VariableNames', cellstr(names));
end

function writeNcVars(filename, names)
   %WRITENCVARS Write tiny one-dimensional NetCDF header/data variables.
   for name = reshape(string(names), 1, [])
      nccreate(filename, name, 'Dimensions', {'time', 2});
      ncwrite(filename, name, [1 2]);
   end
end

function writeJson(filename, value)
   %WRITEJSON Persist a small manifest fixture.
   fid = fopen(filename, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fwrite(fid, jsonencode(value));
   clear cleaner
end

function tf = hasPair(pairs, source_a, source_b, class, variable)
   %HASPAIR True when a pair matches sources, class, and common variable.
   tf = false;
   for k = 1:numel(pairs)
      sources = sort([string(pairs(k).source_a), string(pairs(k).source_b)]);
      expected_sources = sort([string(source_a), string(source_b)]);
      if isequal(sources, expected_sources) ...
            && string(pairs(k).class) == class ...
            && any(string(pairs(k).variables) == variable)
         tf = true;
         return
      end
   end
end

function artifact = firstArtifact(artifacts, source)
   %FIRSTARTIFACT Return the first artifact from a source.
   idx = find(string({artifacts.source}) == source, 1);
   artifact = artifacts(idx);
end
