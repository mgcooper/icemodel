function tests = test_firn_staging_promice
   %TEST_FIRN_STAGING_PROMICE Verify the co-located firn staging driver.
   %
   % Exercises icemodel.verification.setup.importPromiceSites end to end:
   % the co-located PROMICE/MAR/MERRA/RACMO bundle is staged for a PROMICE
   % anchor site and the per-site manifest entry resolves. Reads PROMICE from
   % the verification cache and RCMs from the small test/data/forcing fixtures.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   % Resolve every source the bundle needs; skip the whole file when any is
   % absent (the bundle is co-located, so a partial fixture tree cannot stage
   % a faithful bundle).
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   cfg = icemodel.config('getenv', true);
   forcing_root = string(fullfile(cfg.ICEMODEL_DATA_PATH, 'forcing'));

   promice = firstWithData( ...
      string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'promice')), ...
      @(p) ~isempty(dir(fullfile(p, "hour", "*_hour.nc"))));
   mar = firstWithData( ...
      fullfile(forcing_root, 'mar'), ...
      @(p) ~isempty(dir(fullfile(p, "MARv3.11*.nc"))));
   merra = firstWithData( ...
      fullfile(forcing_root, 'merra2'), ...
      @(p) ~isempty(dir(fullfile(p, "slv", "*_Nx.*.nc4*"))));
   racmo = firstWithData( ...
      fullfile(forcing_root, 'racmo'), ...
      @(p) ~isempty(dir(fullfile(p, "*.RACMO23p3_*.nc"))));

   testCase.assumeTrue(all([ ...
      strlength(promice) > 0, strlength(mar) > 0, ...
      strlength(merra) > 0, strlength(racmo) > 0]), ...
      'co-located firn fixture sources not all available');

   testCase.TestData.promice = promice;
   testCase.TestData.mar = mar;
   testCase.TestData.merra = merra;
   testCase.TestData.racmo = racmo;

   % Stage into a private eval / input root so the committed demo tree is
   % untouched by the test.
   tmp = tempname;
   mkdir(tmp)
   testCase.TestData.tmp = tmp;
   testCase.TestData.eval_root = fullfile(tmp, 'eval');
   testCase.TestData.input_root = fullfile(tmp, 'input');
   mkdir(testCase.TestData.eval_root)
   mkdir(fullfile(testCase.TestData.input_root, 'met'))
   mkdir(fullfile(testCase.TestData.input_root, 'userdata'))
end

function teardownOnce(testCase)
   if isfield(testCase.TestData, 'tmp') && isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
end

function test_dry_run_uncataloged_site_uses_unknown_anchor(testCase)
   % Uncataloged PROMICE ids should reach the site metadata fallback.
   manifest = icemodel.verification.setup.importPromiceSites( ...
      case_ids="UNCAT_999", dry_run=true, ...
      promice_dir=testCase.TestData.promice, ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      icemodel_config_casename="", skip_missing=true);

   testCase.verifyEqual(numel(manifest.cases), 1);
   c = manifest.cases(1);
   leg_fields = ["kind", "staged", "eval_staged", ...
      "met_files", "data_files", "window"];
   testCase.verifyTrue(all(isfield(c.colocation.promice, leg_fields)));
   testCase.verifyEqual(string(c.case_id), "uncat999");
   testCase.verifyEqual(string(c.surface_zone), "unknown");
   testCase.verifyEmpty(string(c.eval_target(:)));
   testCase.verifyEmpty(manifest.skipped);
   testCase.verifyFalse(c.colocation.promice.staged);
   testCase.verifyFalse(c.colocation.promice.eval_staged);
   testCase.verifyEmpty(c.colocation.promice.met_files);
   testCase.verifyEmpty(c.colocation.promice.data_files);
   testCase.verifyEqual(string(c.eval_sources), "promice_obs");
end

function test_native_promice_stage_uses_hourly_userdata_and_15m_met(testCase)
   % Native PROMICE userdata must stay hourly for met swapping; met can be 15m.
   manifest = icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", ...
      promice_dir=testCase.TestData.promice, ...
      build_forcing=true, ...
      startdate="2012-01-01", enddate="2012-01-02", ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      icemodel_config_casename="", overwrite=true);

   c = manifest.cases(string({manifest.cases.case_id}) == "kanm");
   leg_fields = ["kind", "staged", "eval_staged", ...
      "met_files", "data_files", "window"];
   testCase.verifyTrue(all(isfield(c.colocation.promice, leg_fields)));
   testCase.verifyTrue(c.colocation.promice.staged);
   testCase.verifyTrue(c.colocation.promice.eval_staged);
   testCase.verifyEqual(string(manifest.source_version), ...
      "pypromice-L3-hour");
   testCase.verifyEqual(string(c.native_timestep), "hourly");
   comparison_variables = string(c.comparison_variables);
   testCase.verifyTrue(ismember("tice10m", comparison_variables));
   testCase.verifyFalse(any(ismember( ...
      ["tice10m_source", "tice10m_qc_flag"], comparison_variables)));
   met_file = fullfile(testCase.TestData.input_root, 'met', ...
      c.colocation.promice.met_files{1});
   data_file = fullfile(testCase.TestData.input_root, 'userdata', ...
      c.colocation.promice.data_files{1});
   met_bundle = load(met_file, 'met');
   data_bundle = load(data_file, 'Data');

   testCase.verifyTrue(endsWith(string(met_file), "_15m.mat"));
   testCase.verifyEqual(seconds(median(diff(met_bundle.met.Time))), 900);
   testCase.verifyEqual(seconds(median(diff(data_bundle.Data.Time))), 3600);

   % A wider same-case request is not an identical no-op: the fixed observation
   % bundle and window-named native artifacts must widen together visibly.
   observation_file = fullfile(testCase.TestData.eval_root, 'promice', ...
      c.evaluation_file);
   initial = load(observation_file, 'targets');

   % Observations and userdata intentionally persist the same native PROMICE
   % table, so corrected radiation/cloud fraction and the canonical thermistor
   % source-mask contract cannot diverge by sink.
   native_channels = ["swd", "swu", "cfrac", "tice10m", ...
      "tice10m_source", "tice10m_qc_flag"];
   testCase.verifyTrue(all(ismember(native_channels, ...
      string(data_bundle.Data.Properties.VariableNames))));
   testCase.verifyEqual(initial.targets.data.Time, data_bundle.Data.Time);
   for name = native_channels
      testCase.verifyEqual(initial.targets.data.(name), ...
          data_bundle.Data.(name));
   end
   flagged = data_bundle.Data.tice10m_qc_flag > 0;
   testCase.verifyTrue(all(isnan(data_bundle.Data.tice10m(flagged))));
   testCase.verifyEqual(data_bundle.Data.tice10m(~flagged), ...
      data_bundle.Data.tice10m_source(~flagged));

   widened = icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", promice_dir=testCase.TestData.promice, ...
      startdate="2012-01-01", enddate="2012-01-03", ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      icemodel_config_casename="", build_forcing=true, overwrite=false);
   current = load(observation_file, 'targets');
   widened_case = widened.cases( ...
      string({widened.cases.case_id}) == "kanm");

   testCase.verifyGreaterThan(height(current.targets.data), ...
      height(initial.targets.data));
   testCase.verifyEqual(string(widened_case.period.end), ...
      "2012-01-03 00:00:00");
   testCase.verifyTrue(contains(string( ...
      widened_case.colocation.promice.met_files), "_20120103_15m.mat"));
end

function test_rcm_only_forcing_source_skips_promice_runtime_files(testCase)
   % forcing_sources controls runtime artifacts only: PROMICE observations still
   % define the case, but PROMICE met/userdata are skipped when not selected.
   [eval_root, input_root] = privateRoots(testCase);

   manifest = icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", forcing_sources="mar", build_forcing=true, ...
      promice_dir=testCase.TestData.promice, ...
      mar_dir=testCase.TestData.mar, ...
      startdate="2012-01-01", enddate="2012-01-31", ...
      evaluation_data_root=eval_root, ...
      input_data_root=input_root, ...
      icemodel_config_casename="", overwrite=true);

   c = manifest.cases(string({manifest.cases.case_id}) == "kanm");
   leg_fields = ["kind", "staged", "eval_staged", ...
      "met_files", "data_files", "window"];
   testCase.verifyTrue(all(isfield(c.colocation.promice, leg_fields)));
   testCase.verifyTrue(isfile(fullfile(eval_root, ...
      'promice', 'kanm', 'observations.mat')));
   testCase.verifyTrue(ismember("promice_obs", string(c.eval_sources)));
   testCase.verifyFalse(ismember("promice", string(c.forcing_sources)));
   testCase.verifyTrue(ismember("mar3.11", string(c.forcing_sources)));
   testCase.verifyFalse(c.colocation.promice.staged);
   testCase.verifyTrue(c.colocation.promice.eval_staged);
   testCase.verifyEmpty(c.colocation.promice.met_files);
   testCase.verifyEmpty(c.colocation.promice.data_files);
   testCase.verifyEmpty(dir(fullfile(input_root, ...
      'met', 'promice', 'met_kanm_promice_*.mat')));
   testCase.verifyEmpty(dir(fullfile(input_root, ...
      'userdata', 'promice', 'kanm_promice_*.mat')));

   % The artifact audit must not reinterpret the evaluation-only placeholder
   % as a missing native PROMICE runtime leg.
   audit = icemodel.verification.auditArtifacts( ...
      families="promice", evaluation_data_root=eval_root, ...
      input_data_root=input_root, report_dir="");
   audit_codes = string({audit.findings.code});
   testCase.verifyFalse(any(ismember(audit_codes, ...
      ["staged_leg_without_artifact", "missing_artifact"])));

end

function test_build_observations_false_reuses_existing_case(testCase)
   % Fast forcing refresh should reuse the existing observation contract and
   % keep already-staged sources that this call is not rebuilding.
   [eval_root, input_root] = privateRoots(testCase);

   prior = icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", ...
      promice_dir=testCase.TestData.promice, ...
      build_forcing=true, ...
      startdate="2012-01-01", enddate="2012-01-31", ...
      evaluation_data_root=eval_root, ...
      input_data_root=input_root, ...
      icemodel_config_casename="", overwrite=true);
   prior_case = prior.cases(string({prior.cases.case_id}) == "kanm");
   testCase.verifyNotEmpty(prior_case.colocation.promice.met_files);
   testCase.verifyNotEmpty(prior_case.colocation.promice.data_files);

   obs_file = fullfile(eval_root, 'promice', 'kanm', 'observations.mat');
   before = dir(obs_file);

   manifest = icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", build_observations=false, build_forcing=true, ...
      forcing_sources="mar", ...
      promice_dir=fullfile(testCase.TestData.tmp, 'absent-promice'), ...
      mar_dir=testCase.TestData.mar, ...
      evaluation_data_root=eval_root, ...
      input_data_root=input_root, ...
      icemodel_config_casename="", overwrite=true);

   after = dir(obs_file);
   c = manifest.cases(string({manifest.cases.case_id}) == "kanm");
   testCase.verifyEqual(after.bytes, before.bytes);
   testCase.verifyTrue(ismember("promice_obs", string(c.eval_sources)));
   testCase.verifyTrue(ismember("mar3.11", string(c.forcing_sources)));
   testCase.verifyTrue(ismember("promice", string(c.forcing_sources)));
   testCase.verifyNotEmpty(c.colocation.promice.met_files);
   testCase.verifyNotEmpty(c.colocation.promice.data_files);
   testCase.verifyFalse(isfolder(fullfile( ...
      testCase.TestData.tmp, 'absent-promice')));

end

function test_whole_family_forcing_replacement_strips_omitted_runtime(testCase)
   % Whole-family forcing replacement retains observations but removes every
   % omitted native/RCM runtime leg.
   [eval_root, input_root] = privateRoots(testCase);
   initial = icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", ...
      forcing_sources=["promice", "mar", "merra"], build_forcing=true, ...
      promice_dir=testCase.TestData.promice, ...
      mar_dir=testCase.TestData.mar, merra_dir=testCase.TestData.merra, ...
      startdate="2012-01-01", enddate="2012-01-31", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      icemodel_config_casename="", overwrite=true, overwrite_family=true);
   initial_case = initial.cases( ...
      string({initial.cases.case_id}) == "kanm");
   testCase.verifyTrue(logical(initial_case.colocation.promice.staged));
   testCase.verifyTrue(isfield(initial_case.colocation, 'merra'));
   absent_promice = fullfile(testCase.TestData.tmp, 'absent-promice');

   replaced = icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", build_observations=false, build_forcing=true, ...
      forcing_sources="mar", promice_dir=absent_promice, ...
      mar_dir=testCase.TestData.mar, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      icemodel_config_casename="", overwrite=true, overwrite_family=true);

   c = replaced.cases(string({replaced.cases.case_id}) == "kanm");
   testCase.verifyFalse(logical(c.colocation.promice.staged));
   testCase.verifyTrue(logical(c.colocation.promice.eval_staged));
   testCase.verifyEmpty(c.colocation.promice.met_files);
   testCase.verifyEmpty(c.colocation.promice.data_files);
   testCase.verifyFalse(isfield(c.colocation, 'merra'));
   testCase.verifyEqual(string(c.forcing_sources), "mar3.11");
   testCase.verifyTrue(ismember("promice_obs", string(c.eval_sources)));
   testCase.verifyTrue(ismember("mar3.11", string(c.eval_sources)));
   testCase.verifyFalse(ismember("merra2", string(c.eval_sources)));
   testCase.verifyFalse(isfolder(absent_promice));
end

function test_observation_rebuild_preserves_omitted_promice_runtime(testCase)
   % Rebuilding observations while adding MAR keeps the omitted PROMICE leg.
   [eval_root, input_root] = privateRoots(testCase);
   prior = icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", promice_dir=testCase.TestData.promice, ...
      build_forcing=true, ...
      startdate="2012-01-01", enddate="2012-01-31", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      icemodel_config_casename="", overwrite=true);
   prior_case = prior.cases(string({prior.cases.case_id}) == "kanm");
   prior_promice = jsonencode(prior_case.colocation.promice);

   manifest = icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", forcing_sources="mar", build_forcing=true, ...
      promice_dir=testCase.TestData.promice, mar_dir=testCase.TestData.mar, ...
      startdate="2012-01-01", enddate="2012-01-31", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      icemodel_config_casename="", overwrite=true);

   c = manifest.cases(string({manifest.cases.case_id}) == "kanm");
   testCase.verifyEqual(jsonencode(c.colocation.promice), prior_promice);
   testCase.verifyEqual(string(c.colocation.promice.met_files), ...
      string(prior_case.colocation.promice.met_files));
   testCase.verifyEqual(string(c.colocation.promice.data_files), ...
      string(prior_case.colocation.promice.data_files));
   testCase.verifyTrue(ismember("promice", string(c.forcing_sources)));
   testCase.verifyTrue(ismember("mar3.11", string(c.forcing_sources)));
end

function test_build_observations_false_requires_existing_case(testCase)
   % The fast path must not create a new manifest case without observations.
   [eval_root, input_root] = privateRoots(testCase);

   testCase.verifyError(@() icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", build_observations=false, build_forcing=true, ...
      forcing_sources="mar", ...
      promice_dir=testCase.TestData.promice, ...
      mar_dir=testCase.TestData.mar, ...
      evaluation_data_root=eval_root, ...
       input_data_root=input_root, ...
       icemodel_config_casename="", overwrite=true), ...
      'icemodel:verification:reuseDatasetFamilyCases:missingManifest');

end

function test_colocated_bundle_and_manifest_resolve(testCase)
   % Stage a short single-site bundle and assert every model artifact and
   % the firn manifest entry resolve.

   % A short window inside the staged fixture window keeps the test fast.
   manifest = stageKanm(testCase, "2012-01-01", "2012-01-31");

   % Family manifest shape.
   testCase.verifyEqual(string(manifest.dataset_family), "promice");
   testCase.verifyEqual(numel(manifest.cases), 1);
   testCase.verifyEmpty(manifest.skipped);

   c = manifest.cases(1);
   testCase.verifyEqual(string(c.case_id), "kanm");
   testCase.verifyEqual(string(c.case_type), "firn_observational");
   testCase.verifyEqual(string(c.site_id), "KAN_M");

   % surface_zone is single-sourced from promiceSiteCatalog (KAN_M is in the
   % ablation zone) and must validate against the canonical namelist.
   testCase.verifyEqual(string(c.surface_zone), "ablation");
   testCase.verifyTrue(ismember(string(c.surface_zone), ...
      icemodel.verification.namelists.surfacezone()));

   % permafrost_zone is single-sourced from promiceSiteCatalog: KAN_M is an
   % ice-sheet anchor (not permafrost ground) so it is "none" and validates
   % against the canonical permafrost namelist.
   testCase.verifyEqual(string(c.permafrost_zone), "none");
   testCase.verifyTrue(ismember(string(c.permafrost_zone), ...
      icemodel.verification.namelists.permafrostzone()));

   % eval_target: KAN_M exercises seasonal snow + bare ice.
   testCase.verifyEqual(sort(string(c.eval_target(:))), ...
      sort(["bare_ice"; "seasonal_snow"]));

   % Metadata-only source records (ids, not bundled data). PROMICE native met
   % carries precipitation placeholders, so forcing_ready stays false even
   % though the met file is a selectable forcing source.
   testCase.verifyTrue(ismember("promice", string(c.forcing_sources)));
   testCase.verifyTrue(all(ismember(["mar3.11", "merra2"], ...
      string(c.forcing_sources))));
   testCase.verifyTrue(all(ismember(["promice_obs", "racmo2.3p3"], ...
      string(c.eval_sources))));

   % Site location: WGS84 + EPSG:3413 recorded.
   testCase.verifyEqual(c.site_location.lat_wgs84, 67.067, 'AbsTol', 1e-2);
   testCase.verifyTrue(isfinite(c.site_location.x_epsg3413));
   testCase.verifyTrue(isfinite(c.site_location.y_epsg3413));

   % Colocation metadata: all four model legs present with the right kinds.
   cf = c.colocation;
   testCase.verifyTrue(all(isfield(cf, {'promice', 'mar', 'merra', 'racmo'})));
   testCase.verifyEqual(string(cf.racmo.kind), "point_data_smb_eval");
   testCase.verifyEqual(string(cf.mar.source_id), "mar3.11");
   testCase.verifyEqual(string(cf.merra.source_id), "merra2");
   testCase.verifyEqual(string(cf.racmo.source_id), "racmo2.3p3");
   testCase.verifyTrue(~isempty(cf.promice.met_files));
   testCase.verifyFalse(logical(cf.promice.forcing_ready));
   testCase.verifyTrue(contains(string(cf.promice.forcing_ready_reason), ...
      "ppt"));
   testCase.verifyTrue(isfield(cf.promice, 'forcing_complete_windows'));
   testCase.verifyEmpty(cf.promice.forcing_complete_windows);
   testCase.verifyTrue(~isempty(cf.mar.met_files));
   testCase.verifyTrue(~isempty(cf.merra.met_files));
   testCase.verifyTrue(~isempty(cf.racmo.data_files));
end

function test_staged_files_exist_on_disk(testCase)
   % The individual met / userdata files must be written to disk. The eval is a
   % data-only observations.mat; NO evaluation.mat / reference.mat bundle (the
   % forcing/reference side is never bundled with the eval target).

   stageKanm(testCase, "2012-01-01", "2012-01-31");

   met_dir = fullfile(testCase.TestData.input_root, 'met');
   ud_dir = fullfile(testCase.TestData.input_root, 'userdata');
   eval_dir = fullfile(testCase.TestData.eval_root, 'promice', 'kanm');

   % Forcing/userdata stage into per-product subfolders
   % (met/<source_id>/, userdata/<source_id>/).
   testCase.verifyNotEmpty(dir(fullfile(met_dir, 'promice', 'met_kanm_promice_*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(met_dir, 'mar3.11', ...
      'met_kanm_mar3.11_*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(met_dir, 'merra2', ...
      'met_kanm_merra2_*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(ud_dir, 'promice', 'kanm_promice_*.mat')));
   testCase.verifyNotEmpty(dir(fullfile(ud_dir, 'racmo2.3p3', ...
      'kanm_racmo2.3p3_*.mat')));

   % Eval is the data-only observations.mat (no bundled evaluation/reference copy).
   testCase.verifyTrue(isfile(fullfile(eval_dir, 'observations.mat')));
   testCase.verifyFalse(isfile(fullfile(eval_dir, 'evaluation.mat')));
   testCase.verifyFalse(isfile(fullfile(eval_dir, 'reference.mat')));

   testCase.verifyTrue(isfile(fullfile(testCase.TestData.eval_root, ...
      'promice', 'manifest.json')));
end

function test_observation_refresh_preserves_unrequested_runtime_legs(testCase)
   % A later observation patch leaves every unrequested runtime source untouched.
   prior = stageKanm(testCase, "2012-01-01", "2012-01-31");
   prior_case = prior.cases(string({prior.cases.case_id}) == "kanm");
   prior_promice = jsonencode(prior_case.colocation.promice);

   manifest = icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", build_forcing=false, ...
      promice_dir=testCase.TestData.promice, ...
      startdate="2012-02-01", enddate="2012-02-28", ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      icemodel_config_casename="", overwrite=true);

   c = manifest.cases(string({manifest.cases.case_id}) == "kanm");
   for src = ["mar", "merra", "racmo"]
      testCase.verifyTrue(isfield(c.colocation, src));
      testCase.verifyTrue(logical(c.colocation.(char(src)).staged));
      testCase.verifyEqual(jsonencode(c.colocation.(char(src))), ...
         jsonencode(prior_case.colocation.(char(src))));
   end
   testCase.verifyTrue(ismember("mar3.11", string(c.forcing_sources)));
   testCase.verifyTrue(ismember("merra2", string(c.forcing_sources)));
   testCase.verifyTrue(ismember("racmo2.3p3", string(c.eval_sources)));
   testCase.verifyTrue(ismember("promice", string(c.forcing_sources)));
   testCase.verifyEqual(jsonencode(c.colocation.promice), prior_promice);
   testCase.verifyEqual(string(c.colocation.promice.met_files), ...
      string(prior_case.colocation.promice.met_files));
   testCase.verifyEqual(string(c.colocation.promice.data_files), ...
      string(prior_case.colocation.promice.data_files));
end

function test_colocated_data_reconstitutes_from_individual_files(testCase)
   % The eval target (PROMICE obs) and the RCM reference (RACMO Data) must
   % reconstitute on demand from the individual per-year userdata files the
   % forcing-agnostic manifest declares - no committed forcing bundle needed.

   manifest = stageKanm(testCase, "2012-01-01", "2012-01-31");
   c = manifest.cases(1);
   c.manifest_path = fullfile(testCase.TestData.eval_root, 'promice', ...
      'manifest.json');

   promice = icemodel.verification.helpers.loadColocatedData(c, "promice", ...
      input_data_root=testCase.TestData.input_root);
   testCase.verifyEqual(string(promice.format), "timeseries");
   testCase.verifyTrue(istimetable(promice.data));
   testCase.verifyGreaterThan(height(promice.data), 0);

   racmo = icemodel.verification.helpers.loadColocatedData(c, "racmo", ...
      input_data_root=testCase.TestData.input_root);
   testCase.verifyEqual(string(racmo.format), "timeseries");
   testCase.verifyTrue(istimetable(racmo.data));
   testCase.verifyGreaterThan(height(racmo.data), 0);
end

function test_manifest_json_is_valid(testCase)
   % The written manifest.json must parse back to a struct with the firn
   % case schema fields.

   stageKanm(testCase, "2012-01-01", "2012-01-31");
   manifest_file = fullfile(testCase.TestData.eval_root, ...
      'promice', 'manifest.json');

   decoded = jsondecode(fileread(manifest_file));
   testCase.verifyEqual(string(decoded.dataset_family), "promice");

   needed = icemodel.verification.setup.firnCaseManifestFieldNames();
   have = string(fieldnames(decoded.cases));
   testCase.verifyTrue(all(ismember(needed, have)));
end

function test_data_gated_site_is_skipped_not_fabricated(testCase)
   % A nonexistent station id must be recorded in manifest.skipped with a
   % reason and must NOT create a fabricated case folder. Under merge, any
   % previously staged cases are preserved, so assert the gated site itself is
   % skipped and never fabricated rather than that the manifest is empty.

   manifest = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.importPromiceSites( ...
      case_ids="ZZZ_NOPE", ...
      promice_dir=testCase.TestData.promice, ...
      mar_dir=testCase.TestData.mar, ...
      merra_dir=testCase.TestData.merra, ...
      racmo_dir=testCase.TestData.racmo, ...
      startdate="2012-01-01", enddate="2012-01-31", ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      icemodel_config_casename="", overwrite=true), ...
      'icemodel:verification:importPromiceSites:caseSkipped');

   % No fabricated case for the gated site.
   case_ids = string(arrayfun(@(c) string(c.case_id), manifest.cases));
   testCase.verifyFalse(ismember("zzznope", case_ids));

   % The gated site is recorded in skipped with a reason.
   skip_sites = string(arrayfun(@(s) string(s.site), manifest.skipped));
   testCase.verifyTrue(ismember("zzznope", skip_sites));
   idx = find(skip_sites == "zzznope", 1);
   testCase.verifyTrue(strlength(manifest.skipped(idx).reason) > 0);

   testCase.verifyFalse(isfolder(fullfile(testCase.TestData.eval_root, ...
      'promice', 'zzznope')));
end

function test_skipped_overwrite_preserves_existing_observations(testCase)
   % A skipped same-site overwrite must not clear the previous case directory.
   stageKanm(testCase, "2012-01-01", "2012-01-31");
   obs_file = fullfile(testCase.TestData.eval_root, 'promice', 'kanm', ...
      'observations.mat');
   testCase.verifyTrue(isfile(obs_file));

   icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", promice_dir=testCase.TestData.promice, ...
      startdate="1800-01-01", enddate="1800-01-31", ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      icemodel_config_casename="", overwrite=true), ...
      'icemodel:verification:importPromiceSites:caseSkipped');

   testCase.verifyTrue(isfile(obs_file));
end

function test_staging_second_site_does_not_churn_first(testCase)
   % Staging a SECOND site into a family root that already holds a FIRST site
   % must ADD the second case and leave the first site's case entry + its
   % staged files byte for byte unchanged (the KAN no-churn guarantee). This is
   % the file-level counterpart to test_firn_manifest_merge.

   % Stage KAN_L first (a short window keeps it fast).
   stageSite(testCase, "KAN_L", "2012-01-01", "2012-01-31");

   manifest_file = fullfile(testCase.TestData.eval_root, 'promice', ...
      'manifest.json');
   kanl_manifest_before = fileread(manifest_file);

   % Fingerprint the always-staged KAN_L artifacts: the data-only eval bundle
   % (eval/promice/kanl/observations.mat) and the promice forcing/Data files.
   % writemet/writeuserdata stage into a per-source subfolder (met/<source>/,
   % userdata/<source>/), so glob there. The PROMICE met leg is independently
   % skippable when a station/window lacks a required met channel, so it is
   % folded in only when present rather than required.
   eval_case_dir = fullfile(testCase.TestData.eval_root, 'promice', 'kanl');
   met_dir = fullfile(testCase.TestData.input_root, 'met', 'promice');
   ud_dir = fullfile(testCase.TestData.input_root, 'userdata', 'promice');
   kanl_eval = dir(fullfile(eval_case_dir, 'observations.mat'));
   kanl_met = dir(fullfile(met_dir, 'met_kanl_*'));
   kanl_ud = dir(fullfile(ud_dir, 'kanl_*'));
   testCase.assertNotEmpty(kanl_ud, ...
      'KAN_L staged no promice userdata (eval leg failed)');
   before = fileFingerprints([kanl_eval; kanl_met; kanl_ud], ...
      {eval_case_dir, met_dir, ud_dir});

   % Snapshot the KAN_L case JSON region from the merged manifest.
   m_before = jsondecode(kanl_manifest_before);
   kanl_case_before = encodeCaseById(m_before, "kanl");

   % Stage KAN_M into the SAME roots (merge, default).
   stageSite(testCase, "KAN_M", "2012-01-01", "2012-01-31");

   m_after = jsondecode(fileread(manifest_file));
   ids = string(arrayfun(@(c) string(c.case_id), m_after.cases));
   testCase.verifyTrue(all(ismember(["kanl", "kanm"], ids)));

   % KAN_L case entry is byte-identical after KAN_M was added.
   testCase.verifyEqual(encodeCaseById(m_after, "kanl"), kanl_case_before, ...
      'KAN_L case churned when KAN_M was staged');

   % KAN_L's staged files are byte-identical (size + content hash).
   after = fileFingerprints([dir(fullfile(eval_case_dir, 'observations.mat')); ...
      dir(fullfile(met_dir, 'met_kanl_*')); ...
      dir(fullfile(ud_dir, 'kanl_*'))], {eval_case_dir, met_dir, ud_dir});
   testCase.verifyEqual(after, before, ...
      'KAN_L staged files churned when KAN_M was staged');
end

function test_missing_rcm_legs_still_yield_a_case(testCase)
   % PROMICE met+eval must never be gated by RCM coverage. A site whose
   % MAR/MERRA/RACMO source dirs are absent must still stage a PROMICE case:
   % the RCM legs degrade to skipped legs in colocation, and the site is not
   % wholesale skipped.

   % Resolve PROMICE on its own; this test does not need the RCM sources.
   promice = firstWithData( ...
      string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'promice')), ...
      @(p) ~isempty(dir(fullfile(p, "hour", "*_hour.nc"))));
   testCase.assumeTrue(strlength(promice) > 0, ...
      'PROMICE source not available (no *_hour.nc on disk)');

   % Private roots so nothing committed is touched.
   tmp = tempname;
   mkdir(tmp)
   cleanup = onCleanup(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');
   mkdir(eval_root)
   mkdir(fullfile(input_root, 'met'))
   mkdir(fullfile(input_root, 'userdata'))

   % Point every RCM dir at a path that cannot exist, forcing each builder to
   % throw the "source directory not found" error.
   missing = fullfile(tmp, 'no_such_rcm_source');

   manifest = icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", forcing_sources=["promice", "mar", "merra", "racmo"], ...
      build_forcing=true, ...
      promice_dir=promice, ...
      mar_dir=missing, merra_dir=missing, racmo_dir=missing, ...
      startdate="2012-01-01", enddate="2012-01-31", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      icemodel_config_casename="", overwrite=true);

   % The site STAGED as a case (not skipped) despite every RCM leg failing.
   testCase.verifyEqual(numel(manifest.cases), 1, ...
      'PROMICE site was wrongly skipped when its RCM legs were absent');
   c = manifest.cases(1);
   testCase.verifyEqual(string(c.case_id), "kanm");

   % KAN_M is NOT in skipped - only its RCM legs are.
   if ~isempty(manifest.skipped)
      skip_sites = string(arrayfun(@(s) string(s.site), manifest.skipped));
      testCase.verifyFalse(ismember("KAN_M", skip_sites));
   end

   % PROMICE leg fully stages observations/met, but its native met carries
   % precipitation placeholders and needs filling before a clean run.
   testCase.verifyTrue(ismember("promice", string(c.forcing_sources)));
   testCase.verifyTrue(ismember("promice_obs", string(c.eval_sources)));
   testCase.verifyFalse(logical(c.colocation.promice.forcing_ready));
   testCase.verifyTrue(contains( ...
      string(c.colocation.promice.forcing_ready_reason), "ppt"));
   testCase.verifyTrue(isfield(c.colocation.promice, ...
      'forcing_complete_windows'));
   testCase.verifyEmpty(c.colocation.promice.forcing_complete_windows);

   % Each RCM leg recorded as a skipped leg with a reason, never fabricated.
   cf = c.colocation;
   testCase.verifyFalse(cf.mar.staged);
   testCase.verifyFalse(cf.merra.staged);
   testCase.verifyFalse(cf.racmo.staged);
   testCase.verifyGreaterThan(strlength(string(cf.mar.reason)), 0);

   % RCM legs are absent from the forcing/eval source lists (no fabrication).
   testCase.verifyFalse(ismember("mar3.11", string(c.forcing_sources)));
   testCase.verifyFalse(ismember("merra2", string(c.forcing_sources)));
   testCase.verifyFalse(ismember("racmo2.3p3", string(c.eval_sources)));
end

function test_missing_lwd_promice_site_writes_placeholder_met(testCase)
   % LYN_L has no observed longwave channel in the PROMICE L3 product. The
   % verification importer still needs a native met artifact for runtime
   % swapping, so missing required channels are explicit NaN placeholders
   % instead of whole-site skips.
   lyn_file = fullfile(testCase.TestData.promice, 'hour', 'LYN_L_hour.nc');
   testCase.assumeTrue(isfile(lyn_file), ...
      'LYN_L PROMICE source not available; skipping placeholder-met test.');

   tmp = tempname;
   mkdir(tmp)
   cleanup = onCleanup(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');
   mkdir(eval_root)
   mkdir(fullfile(input_root, 'met'))
   mkdir(fullfile(input_root, 'userdata'))

   manifest = icemodel.verification.setup.importPromiceSites( ...
      case_ids="LYN_L", ...
      promice_dir=testCase.TestData.promice, ...
      build_forcing=true, ...
      startdate="2020-01-01", enddate="2020-01-02", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      icemodel_config_casename="", overwrite=true);

   testCase.verifyEqual(string(manifest.cases.case_id), "lynl");
   testCase.verifyEmpty(manifest.skipped);

   met_files = dir(fullfile(input_root, 'met', 'promice', ...
      'met_lynl_promice_*.mat'));
   testCase.assertNumElements(met_files, 1);

   loaded = load(fullfile(met_files(1).folder, met_files(1).name), 'met');
   testCase.verifyTrue(ismember("lwd", ...
      string(loaded.met.Properties.VariableNames)));
   testCase.verifyTrue(all(isnan(loaded.met.lwd)));
   testCase.verifyTrue(logical( ...
      loaded.met.Properties.UserData.fillwithmissing));
   testCase.verifyWarningFree( ...
      @() icemodel.forcing.helpers.validatemet(loaded.met));
   testCase.verifyTrue(ismember("promice", ...
      string(manifest.cases.forcing_sources)));
   testCase.verifyFalse(logical( ...
      manifest.cases.colocation.promice.forcing_ready));
   testCase.verifyTrue(contains(string( ...
      manifest.cases.colocation.promice.forcing_ready_reason), "lwd"));
   testCase.verifyTrue(isfield(manifest.cases.colocation.promice, ...
      'forcing_complete_windows'));
   testCase.verifyEmpty( ...
      manifest.cases.colocation.promice.forcing_complete_windows);

   clear cleanup
end

function test_one_sample_promice_met_failure_keeps_observations(testCase)
   % A too-short native met window must not discard staged PROMICE targets.
   lyn_file = fullfile(testCase.TestData.promice, 'hour', 'LYN_L_hour.nc');
   testCase.assumeTrue(isfile(lyn_file), ...
      'LYN_L PROMICE source not available; skipping one-sample test.');

   tmp = tempname;
   mkdir(tmp)
   cleanup = onCleanup(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');
   mkdir(eval_root)
   mkdir(fullfile(input_root, 'met'))
   mkdir(fullfile(input_root, 'userdata'))

   manifest = icemodel.verification.setup.importPromiceSites( ...
      case_ids="LYN_L", ...
      promice_dir=testCase.TestData.promice, ...
      build_forcing=true, ...
      startdate="2020-01-01", enddate="2020-01-01", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      icemodel_config_casename="", overwrite=true);

   testCase.verifyEqual(string(manifest.cases.case_id), "lynl");
   testCase.verifyEmpty(manifest.skipped);
   testCase.verifyTrue(isfile(fullfile(eval_root, 'promice', 'lynl', ...
      'observations.mat')));
   testCase.verifyNotEmpty(dir(fullfile(input_root, 'userdata', 'promice', ...
      'lynl_promice_*.mat')));
   testCase.verifyEmpty(dir(fullfile(input_root, 'met', 'promice', ...
      'met_lynl_promice_*.mat')));
   testCase.verifyFalse(ismember("promice", ...
      string(manifest.cases.forcing_sources)));
   testCase.verifyTrue(ismember("promice_obs", ...
      string(manifest.cases.eval_sources)));
   testCase.verifyTrue(manifest.cases.colocation.promice.staged);
   testCase.verifyTrue(manifest.cases.colocation.promice.eval_staged);
   testCase.verifyNotEmpty(manifest.cases.colocation.promice.data_files);
   testCase.verifyEmpty(manifest.cases.colocation.promice.met_files);
   testCase.verifyTrue(isfield(manifest.cases.colocation.promice, ...
      'met_skipped_reason'));
   testCase.verifyTrue(contains(string( ...
      manifest.cases.colocation.promice.met_skipped_reason), ...
      "at least two samples"));

   clear cleanup
end

%% Local helpers
function manifest = stageKanm(testCase, startdate, enddate)
   %STAGEKANM Stage the KAN_M bundle into the private test roots.
   manifest = suppressExpectedWarning( ...
      @() icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", forcing_sources=["promice", "mar", "merra", "racmo"], ...
      build_forcing=true, ...
      promice_dir=testCase.TestData.promice, ...
      mar_dir=testCase.TestData.mar, ...
      merra_dir=testCase.TestData.merra, ...
      racmo_dir=testCase.TestData.racmo, ...
      startdate=startdate, enddate=enddate, ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      icemodel_config_casename="", overwrite=true), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');
end

function [eval_root, input_root] = privateRoots(testCase)
   %PRIVATEROOTS Create isolated eval/input roots for API contract tests.
   tmp = tempname;
   mkdir(tmp)
   testCase.addTeardown(@() rmdir(tmp, 's'));
   eval_root = fullfile(tmp, 'eval');
   input_root = fullfile(tmp, 'input');
   mkdir(eval_root)
   mkdir(fullfile(input_root, 'met'))
   mkdir(fullfile(input_root, 'userdata'))
end

function manifest = stageSite(testCase, site, startdate, enddate)
   %STAGESITE Stage one site's bundle into the private test roots (merge).
   manifest = suppressExpectedWarning( ...
      @() icemodel.verification.setup.importPromiceSites( ...
      case_ids=site, forcing_sources=["promice", "mar", "merra", "racmo"], ...
      build_forcing=true, ...
      promice_dir=testCase.TestData.promice, ...
      mar_dir=testCase.TestData.mar, ...
      merra_dir=testCase.TestData.merra, ...
      racmo_dir=testCase.TestData.racmo, ...
      startdate=startdate, enddate=enddate, ...
      evaluation_data_root=testCase.TestData.eval_root, ...
      input_data_root=testCase.TestData.input_root, ...
      icemodel_config_casename="", overwrite=true), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');
end

function returned = suppressExpectedWarning(fcn, expected_id)
   %SUPPRESSEXPECTEDWARNING Silence incidental warnings covered elsewhere.
   state = warning('off', expected_id);
   cleanup = onCleanup(@() warning(state));
   returned = fcn();
   clear cleanup
end

function fp = fileFingerprints(entries, ~)
   %FILEFINGERPRINTS Map each file to {bytes, checksum} for a churn comparison.
   fp = struct();
   for k = 1:numel(entries)
      e = entries(k);
      f = fullfile(e.folder, e.name);
      key = matlab.lang.makeValidName(e.name);
      fid = fopen(f, 'r');
      raw = fread(fid, Inf, '*uint8');
      fclose(fid);
      fp.(key) = struct('bytes', e.bytes, 'sum', sum(double(raw)));
   end
end

function s = encodeCaseById(manifest, case_id)
   %ENCODECASEBYID Stable JSON encoding of the case with the given id.
   ids = string(arrayfun(@(c) string(c.case_id), manifest.cases));
   s = jsonencode(manifest.cases(ids == case_id), PrettyPrint=true);
end

function path = firstWithData(candidates, hasdata)
   %FIRSTWITHDATA First candidate folder that actually contains data files.
   path = "";
   for c = candidates
      if isfolder(c) && hasdata(c)
         path = c;
         return
      end
   end
end
