function tests = test_retmip_imau_sources
   %TEST_RETMIP_IMAU_SOURCES Verify RetMIP/IMAU metadata and fetch contracts.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install the canonical verification path and allocate an isolated cache.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.cache = tempname;
end

function teardown(testCase)
   % Remove temporary source-cache artifacts.
   if isfolder(testCase.TestData.cache)
      rmdir(testCase.TestData.cache, 's')
   end
   clear testCase.TestData.cleanup
end

function test_prepare_case_root_selects_only_missing_artifacts(testCase)
   % Repeated setup calls may add a missing artifact without replacing a
   % current sibling; overwrite=true remains the explicit replacement signal.
   case_root = fullfile(testCase.TestData.cache, 'artifact-mask');
   icemodel.helpers.ensureDirExists(case_root);
   targets = struct('format', 'legacy_timeseries');
   save(fullfile(case_root, 'observations.mat'), 'targets');
   location = struct('lat_wgs84', 67, 'lon_wgs84', -48);
   prior_case = struct('case_id', 'artifact-mask', ...
      'period', struct('start', '2012-01-01 00:00:00', ...
      'end', '2012-01-31 23:00:00'), 'site_location', location);
   writeText(fullfile(testCase.TestData.cache, 'manifest.json'), ...
      jsonencode(struct('cases', prior_case)));

   selected = icemodel.verification.setup.prepareCaseRoot( ...
      case_root, false, ["observations.mat", "reference.mat"]);

   testCase.verifyEqual(selected, [false, true]);
   covered = struct('period', struct('start', '2012-01-10 00:00:00', ...
      'end', '2012-01-20 23:00:00'), 'site_location', location);
   testCase.verifyFalse(icemodel.verification.setup.prepareCaseRoot( ...
      case_root, false, "observations.mat", covered));
   period_only = rmfield(covered, 'site_location');
   testCase.verifyFalse(icemodel.verification.setup.prepareCaseRoot( ...
      case_root, false, "observations.mat", period_only));

   % A wider or point-incompatible request makes the fixed artifact stale and
   % selects a visible replacement even when overwrite=false.
   wider = covered;
   wider.period.end = '2012-02-29 23:00:00';
   lastwarn("");
   testCase.verifyTrue(icemodel.verification.setup.prepareCaseRoot( ...
      case_root, false, "observations.mat", wider));
   [~, warning_id] = lastwarn;
   testCase.verifyEqual(string(warning_id), ...
      "icemodel:verification:prepareCaseRoot:overwrite");
   moved = covered;
   moved.site_location.lon_wgs84 = -47;
   testCase.verifyTrue(icemodel.verification.setup.prepareCaseRoot( ...
      case_root, false, "observations.mat", moved));

   % Blank bounds can occur in a prior or requested manifest during additive
   % reruns. Either direction is unknowable, so select a rebuild without a
   % zoned-versus-unzoned datetime concatenation error.
   unbounded_prior = prior_case;
   unbounded_prior.period = struct('start', '', 'end', '');
   writeText(fullfile(testCase.TestData.cache, 'manifest.json'), ...
      jsonencode(struct('cases', unbounded_prior)));
   testCase.verifyTrue(icemodel.verification.setup.prepareCaseRoot( ...
      case_root, false, "observations.mat", covered));
   writeText(fullfile(testCase.TestData.cache, 'manifest.json'), ...
      jsonencode(struct('cases', prior_case)));
   unbounded_request = covered;
   unbounded_request.period = struct('start', '', 'end', '');
   testCase.verifyTrue(icemodel.verification.setup.prepareCaseRoot( ...
      case_root, false, "observations.mat", unbounded_request));

   % Without prior manifest evidence an existing fixed artifact is not current
   % for a newly supplied case contract.
   orphan_root = fullfile(testCase.TestData.cache, 'orphan', 'case1');
   icemodel.helpers.ensureDirExists(orphan_root);
   writeText(fullfile(orphan_root, 'observations.mat'), 'orphan');
   testCase.verifyTrue(icemodel.verification.setup.prepareCaseRoot( ...
      orphan_root, false, "observations.mat", covered));

   overwrite_call = @() icemodel.verification.setup.prepareCaseRoot( ...
      case_root, true, ["observations.mat", "reference.mat"]);
   lastwarn("");
   overwrite_selected = overwrite_call();
   [~, warning_id] = lastwarn;
   testCase.verifyEqual(overwrite_selected, [true, true]);
   testCase.verifyEqual(string(warning_id), ...
      "icemodel:verification:prepareCaseRoot:overwrite");
end

function test_prepare_case_root_gates_fixed_artifact_identity(testCase)
   % A covered exact-window observation bundle is reusable only across equal
   % concrete producer identity; legacy-missing facts remain compatible.
   family_root = fullfile(testCase.TestData.cache, 'identity-family');
   case_root = fullfile(family_root, 'case1');
   icemodel.helpers.ensureDirExists(case_root);
   location = struct('lat_wgs84', 67, 'lon_wgs84', -48);
   period = struct('start', '2012-01-01 00:00:00', ...
      'end', '2012-01-31 23:00:00');
   prior_case = struct('case_id', 'case1', 'period', period, ...
      'site_location', location);
   writeText(fullfile(family_root, 'manifest.json'), ...
      jsonencode(struct('cases', prior_case)));

   [fields, matching_values, conflicting_values] = ...
      fixedArtifactIdentityCases();
   metadata = struct();
   for k = 1:numel(fields)
      metadata.(char(fields(k))) = matching_values(k);
   end
   targets = struct('format', 'timeseries', 'data', 1, ...
      'metadata', metadata);
   observation_file = fullfile(case_root, 'observations.mat');
   save(observation_file, 'targets');
   original = fileBytes(observation_file);
   requested = struct('period', period, 'site_location', location, ...
      'artifact_metadata', metadata);

   % Equal identity is an exact no-op and therefore preserves saved bytes.
   testCase.verifyFalse(icemodel.verification.setup.prepareCaseRoot( ...
      case_root, false, "observations.mat", requested));
   testCase.verifyEqual(fileBytes(observation_file), original);

   for k = 1:numel(fields)
      % Hold period, point, and every sibling identity fact constant so each
      % selected rewrite is attributable to exactly one concrete mismatch.
      changed = requested;
      name = char(fields(k));
      changed.artifact_metadata.(name) = conflicting_values(k);
      lastwarn("");
      testCase.verifyTrue(icemodel.verification.setup.prepareCaseRoot( ...
         case_root, false, "observations.mat", changed));
      [~, warning_id] = lastwarn;
      testCase.verifyEqual(string(warning_id), ...
         "icemodel:verification:prepareCaseRoot:overwrite");
      testCase.verifyEqual(fileBytes(observation_file), original);
   end

   % A legacy bundle missing one producer fact remains reusable when every
   % concrete fact it does carry agrees with the request.
   targets.metadata = rmfield(targets.metadata, 'bundle_doi');
   save(observation_file, 'targets');
   legacy = fileBytes(observation_file);
   testCase.verifyFalse(icemodel.verification.setup.prepareCaseRoot( ...
      case_root, false, "observations.mat", requested));
   testCase.verifyEqual(fileBytes(observation_file), legacy);

   % Malformed requested or saved metadata and unreadable fixed artifacts are
   % not legacy-unknown facts; each must select a repair rather than reuse.
   malformed_request = requested;
   malformed_request.artifact_metadata = repmat(metadata, 1, 2);
   testCase.verifyTrue(icemodel.verification.setup.prepareCaseRoot( ...
      case_root, false, "observations.mat", malformed_request));
   targets.metadata = repmat(metadata, 1, 2);
   save(observation_file, 'targets');
   testCase.verifyTrue(icemodel.verification.setup.prepareCaseRoot( ...
      case_root, false, "observations.mat", requested));
   writeText(observation_file, 'corrupt fixed artifact');
   testCase.verifyTrue(icemodel.verification.setup.prepareCaseRoot( ...
      case_root, false, "observations.mat", requested));
end

function test_fetch_retmip_reports_doi_urls(testCase)
   % fetchRetmip must expose the GEUS Dataverse DOI URLs without requiring the
   % large local products to exist.
   [returned, status] = icemodel.verification.setup.fetchRetmip( ...
      cache_dir=testCase.TestData.cache, strict=false, silent=true);

   expected = string(testCase.TestData.cache);
   testCase.verifyEqual(string(returned), expected);
   returned = string({status.product});
   expected = ["forcing", "outputs", "scripts"];
   testCase.verifyTrue(all(ismember(expected, returned)));
   testCase.verifyTrue(any(contains([status.landing_url], ...
      "10.22008/FK2/GZ3CSN")));
   testCase.verifyTrue(any(contains([status.landing_url], ...
      "10.22008/FK2/CVPUJL")));
   testCase.verifyFalse(all([status.present]));
end

function test_fetch_imau_reports_pangaea_urls(testCase)
   % fetchImau must keep the hourly case inventory and daily QA product distinct.
   [~, returned] = icemodel.verification.setup.fetchImau( ...
      cache_dir=testCase.TestData.cache, strict=false, silent=true);

   expected = ["hourly", "daily"];
   testCase.verifyEqual(sort(string({returned.product})), sort(expected));
   testCase.verifyTrue(any(contains([returned.landing_url], ...
      "10.1594/PANGAEA.971647")));
   testCase.verifyTrue(any(contains([returned.landing_url], ...
      "10.1594/PANGAEA.970127")));
end

function test_preview_firn_staging_clamps_long_retmip_case(testCase)
   % RetMIP previews should stage long protocol cases as one-year QA windows.
   cache = makeRetmipProtocolCache(testCase.TestData.cache, "Summit");
   output_root = fullfile(testCase.TestData.cache, 'preview-long');

   returned = icemodel.verification.setup.previewFirnStaging("retmip", ...
      retmip_dir=cache, retmip_case_ids="summit", make_plots=false, ...
      output_root=output_root);

   [t1, t2] = icemodel.verification.setup.periodBounds( ...
      returned.retmip.cases.period);
   testCase.verifyLessThanOrEqual(days(t2 - t1), 366);

   loaded = load(fullfile(output_root, 'eval', 'retmip', ...
      returned.retmip.cases.evaluation_file), 'targets');
   surface = loaded.targets.data.surface;
   testCase.verifyTrue(all(surface.time >= t1 & surface.time <= t2));
end

function test_preview_firn_staging_returns_all_retmip_cases(testCase)
   % Multi-case RetMIP previews should return the final merged manifest.
   cache = makeRetmipProtocolCache(testCase.TestData.cache, "Summit");
   forcing = fullfile(cache, 'forcing');
   outputs = fullfile(cache, 'outputs');
   writeText(fullfile(forcing, 'RetMIP_surface_forcing_KANU.tab'), ...
      sprintf(['time;melt_mmweq;acc_subl_mmweq;Tsurf_K\n' ...
      '"01-May-2012 00:00:00;0;0;260"\n' ...
      '"01-May-2012 03:00:00;1;2;261"\n']));
   touch(fullfile(outputs, 'RetMIP_outputs_KANU.tab'));
   output_root = fullfile(testCase.TestData.cache, 'preview-many');
   missing_rcm = fullfile(testCase.TestData.cache, 'missing-preview-rcm');

   returned = icemodel.verification.setup.previewFirnStaging("retmip", ...
      retmip_dir=cache, retmip_case_ids=["summit", "kanu"], ...
      make_plots=false, output_root=output_root, forcing_sources="mar", ...
      mar_dir=missing_rcm, merra_dir=missing_rcm, racmo_dir=missing_rcm);

   ids = sort(string({returned.retmip.cases.case_id}));
   testCase.verifyEqual(ids, ["kanu", "summit"]);
   cases = returned.retmip.cases;
   for k = 1:numel(cases)
      c = cases(k);
      testCase.verifyTrue(isfield(c.colocation, 'mar'));
      testCase.verifyFalse(isfield(c.colocation, 'merra'));
      testCase.verifyFalse(isfield(c.colocation, 'racmo'));
   end
end

function test_preview_firn_staging_forwards_plot_overwrite(testCase)
   % Preview regeneration should clear obsolete PNGs through the plotter.
   cache = makeRetmipProtocolCache(testCase.TestData.cache, "Summit");
   output_root = fullfile(testCase.TestData.cache, 'preview-overwrite');
   figure_root = fullfile(testCase.TestData.cache, 'preview-figures');
   case_dir = fullfile(figure_root, 'retmip', 'summit');
   obsolete = fullfile(case_dir, 'energy_fluxes.png');
   missing_rcm = fullfile(testCase.TestData.cache, 'missing-preview-rcm');
   icemodel.helpers.ensureDirExists(case_dir);
   fclose(fopen(obsolete, 'w'));

   returned = icemodel.verification.setup.previewFirnStaging("retmip", ...
      retmip_dir=cache, retmip_case_ids="summit", output_root=output_root, ...
      figure_root=figure_root, make_plots=true, overwrite=true, ...
      mar_dir=missing_rcm, merra_dir=missing_rcm, racmo_dir=missing_rcm);
   current = dir(fullfile(case_dir, '*.png'));

   testCase.verifyEqual(string({returned.retmip.cases.case_id}), "summit");
   testCase.verifyFalse(isfile(obsolete));
   testCase.verifyNotEmpty(current);
end

function test_preview_firn_staging_forwards_native_dt_override(testCase)
   % The preview wrapper exposes the same explicit native model-met cadence as
   % its public family importers instead of hiding a staging-only default.
   cache = makeImauSourceCache(testCase.TestData.cache, "S21");
   output_root = fullfile(testCase.TestData.cache, 'preview-native-dt');
   missing_rcm = fullfile(testCase.TestData.cache, 'missing-preview-rcm');

   returned = icemodel.verification.setup.previewFirnStaging("imau", ...
      imau_dir=cache, imau_case_ids="S21", output_root=output_root, ...
      forcing_sources="mar", mar_dir=missing_rcm, ...
      merra_dir=missing_rcm, racmo_dir=missing_rcm, ...
      make_plots=false, dt_out="");

   met_file = string(returned.imau.cases.colocation.imau.met_files(1));
   testCase.verifyTrue(endsWith(met_file, "_1hr.mat"));
end

function test_preview_firn_staging_includes_retmip_native_selector(testCase)
   % RetMIP previews request native runtime files in addition to selected RCMs.
   cache = makeRetmipProtocolCache(testCase.TestData.cache, "Summit");
   gcnet_cache = makeGcnetSurfaceCache(testCase.TestData.cache);
   output_root = fullfile(testCase.TestData.cache, 'preview-retmip-native');
   missing_rcm = fullfile(testCase.TestData.cache, 'missing-preview-rcm');

   returned = icemodel.verification.setup.previewFirnStaging("retmip", ...
      retmip_dir=cache, retmip_case_ids="summit", gcnet_dir=gcnet_cache, ...
      output_root=output_root, forcing_sources="mar", ...
      mar_dir=missing_rcm, merra_dir=missing_rcm, racmo_dir=missing_rcm, ...
      make_plots=false, dt_out="");

   c = returned.retmip.cases;
   testCase.verifyTrue(any(string(c.forcing_sources) == "retmip"));
   testCase.verifyNotEmpty(string(c.colocation.retmip.data_files));
end

function test_import_retmip_half_window_preserves_error_id(testCase)
   % RetMIP window validation should expose a catchable identifier.
   testCase.verifyError(@() icemodel.verification.setup.importRetmip( ...
      dry_run=true, startdate="2012-01-01"), ...
      'icemodel:internal:pairedWindow:invalidWindow');
end

function test_imau_forcing_only_fast_path_requires_manifest(testCase)
   % build_observations=false must not fabricate metadata without staged state.
   root = fullfile(testCase.TestData.cache, 'missing-imau-fast-path');
   testCase.verifyError(@() ...
      icemodel.verification.setup.importImau( ...
      fullfile(root, 'absent-source'), case_ids="S21", output_root=root, ...
      build_observations=false, build_forcing=false), ...
      'icemodel:verification:reuseDatasetFamilyCases:missingManifest');
end

function test_retmip_forcing_only_fast_path_does_not_read_raw_cache(testCase)
   % A forcing-only replay must reuse the staged RetMIP case without reading or
   % creating the raw cache, even when no RCM source is requested.
   root = fullfile(testCase.TestData.cache, 'retmip-forcing-only-root');
   missing_cache = fullfile(root, 'missing-retmip-cache');
   family_root = fullfile(root, 'eval', 'retmip');
   manifest_file = fullfile(family_root, 'manifest.json');
   initial = icemodel.verification.setup.importRetmip( ...
      missing_cache, case_ids="kanu", dry_run=true, ...
      forcing_sources=strings(1, 0), build_forcing=false);
   mkdir(family_root)
   icemodel.verification.setup.writeManifest(manifest_file, initial);

   returned = icemodel.verification.setup.importRetmip( ...
      missing_cache, case_ids="kanu", output_root=root, ...
      forcing_sources=strings(1, 0), build_observations=false, ...
      build_forcing=false);

   testCase.verifyEqual(string({returned.cases.case_id}), "kanu");
   testCase.verifyTrue(isfile(manifest_file));
   testCase.verifyFalse(isfolder(missing_cache));
   testCase.verifyFalse(isfolder(fullfile(root, 'input')));
end

function test_imau_rcm_only_selection_does_not_stage_native_artifacts(testCase)
   % Selecting only an RCM must reuse the case and leave IMAU runtime files absent.
   [mar_dir, ~, ~] = requireRcmFixtures(testCase);
   cache = makeImauSourceCache(testCase.TestData.cache, "S21", ...
      time0=datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'));
   root = fullfile(testCase.TestData.cache, 'imau-rcm-only-root');
   missing_cache = fullfile(root, 'missing-imau-cache');
   icemodel.verification.setup.importImau( ...
      cache, case_ids="S21", output_root=root, ...
      build_forcing=false, overwrite_family=true);

   returned = icemodel.verification.setup.importImau( ...
      missing_cache, case_ids="S21", output_root=root, ...
      forcing_sources="mar", build_observations=false, build_forcing=true, ...
      mar_dir=mar_dir, overwrite=true);

   c = returned.cases;
   testCase.verifyTrue(any(string(c.forcing_sources) == "mar3.11"));
   testCase.verifyFalse(any(string(c.forcing_sources) == "imau"));
   testCase.verifyFalse(isfolder(missing_cache));
   testCase.verifyFalse(isfolder(fullfile(root, 'input', 'met', 'imau')));
   testCase.verifyFalse(isfolder(fullfile(root, 'input', 'userdata', 'imau')));
end

function test_retmip_rcm_only_selection_does_not_stage_native_artifacts(testCase)
   % Selecting only an RCM must reuse the case and leave RetMIP runtime files absent.
   [mar_dir, ~, ~] = requireRcmFixtures(testCase);
   cache = makeRetmipProtocolCache(testCase.TestData.cache, "Dye-2_long");
   root = fullfile(testCase.TestData.cache, 'retmip-rcm-only-selector-root');
   missing_cache = fullfile(root, 'missing-retmip-cache');
   icemodel.verification.setup.importRetmip( ...
      cache, case_ids="dye2", output_root=root, ...
      build_forcing=false, overwrite_family=true);

   returned = icemodel.verification.setup.importRetmip( ...
      missing_cache, case_ids="dye2", output_root=root, ...
      forcing_sources="mar", build_observations=false, build_forcing=true, ...
      mar_dir=mar_dir, overwrite=true);

   c = returned.cases;
   testCase.verifyTrue(any(string(c.forcing_sources) == "mar3.11"));
   testCase.verifyFalse(any(string(c.forcing_sources) == "retmip"));
   testCase.verifyFalse(isfolder(missing_cache));
   testCase.verifyFalse(isfolder(fullfile(root, 'input', 'met', 'retmip')));
   testCase.verifyFalse(isfolder(fullfile(root, 'input', 'userdata', 'retmip')));
end

function test_retmip_rcm_only_replacement_strips_prior_runtime(testCase)
   % Whole-family RCM-only replacement retains the RetMIP protocol contract but
   % removes every omitted native and RCM runtime leg from the manifest.
   [mar_dir, merra_dir, ~] = requireRcmFixtures(testCase);
   cache = makeRetmipProtocolCache(testCase.TestData.cache, "Dye-2_long");
   gcnet_cache = makeGcnetSurfaceCache(testCase.TestData.cache);
   root = fullfile(testCase.TestData.cache, ...
      'retmip-rcm-only-replacement-root');
   initial = icemodel.verification.setup.importRetmip( ...
      cache, case_ids="dye2", output_root=root, ...
      forcing_sources=["retmip", "mar", "merra"], build_forcing=true, ...
      gcnet_dir=gcnet_cache, mar_dir=mar_dir, merra_dir=merra_dir, ...
      overwrite_family=true);
   testCase.verifyTrue(isfield(initial.cases.colocation, 'native_met'));
   testCase.verifyTrue(isfield(initial.cases.colocation, 'merra'));
   missing_cache = fullfile(root, 'missing-retmip-cache');

   replaced = icemodel.verification.setup.importRetmip( ...
      missing_cache, case_ids="dye2", output_root=root, ...
      forcing_sources="mar", build_observations=false, build_forcing=true, ...
      mar_dir=mar_dir, overwrite=true, overwrite_family=true);

   c = replaced.cases;
   testCase.verifyEqual(string(c.forcing_sources), "mar3.11");
   testCase.verifyFalse(isfield(c.colocation, 'native_met'));
   testCase.verifyFalse(isfield(c.colocation, 'merra'));
   testCase.verifyFalse(isfield(c.colocation.retmip, 'native_source'));
   testCase.verifyFalse(isfield( ...
      c.observation_variables, 'native_source'));
   testCase.verifyTrue(ismember("retmip_protocol", string(c.eval_sources)));
   testCase.verifyTrue(ismember("mar3.11", string(c.eval_sources)));
   testCase.verifyFalse(ismember("merra2", string(c.eval_sources)));
   testCase.verifyFalse(isfolder(missing_cache));
end

function test_retmip_strict_rcm_only_import_skips_native_preflight(testCase)
   % Strict RCM-only imports require protocol data but not omitted native caches.
   [mar_dir, ~, ~] = requireRcmFixtures(testCase);
   cache = makeRetmipProtocolCache(testCase.TestData.cache, "Dye-2_long");
   root = fullfile(testCase.TestData.cache, 'retmip-strict-rcm-only-root');
   missing_native = fullfile(root, 'missing-gcnet-cache');

   returned = icemodel.verification.setup.importRetmip( ...
      cache, case_ids="dye2", output_root=root, forcing_sources="mar", ...
      build_forcing=true, skip_missing=false, gcnet_dir=missing_native, ...
      mar_dir=mar_dir, overwrite_family=true);

   c = returned.cases;
   testCase.verifyTrue(logical(c.colocation.mar.staged));
   testCase.verifyFalse(any(string(c.forcing_sources) == "retmip"));
   testCase.verifyFalse(isfolder(missing_native));
   testCase.verifyFalse(isfolder(fullfile(root, 'input', 'met', 'retmip')));
   testCase.verifyFalse(isfolder(fullfile(root, 'input', 'userdata', 'retmip')));
end

function test_retmip_strict_native_only_preflights_all_cases(testCase)
   % Strict native-only refreshes must validate every source before writing any
   % earlier case's runtime artifacts.
   cache = makeRetmipProtocolCache(testCase.TestData.cache, ...
      ["Dye-2_long", "Summit"]);
   gcnet_cache = makeGcnetSurfaceCache(testCase.TestData.cache);
   root = fullfile(testCase.TestData.cache, ...
      'retmip-strict-native-only-preflight-root');
   icemodel.verification.setup.importRetmip( ...
      cache, case_ids=["dye2", "SUM"], output_root=root, ...
      build_forcing=false, overwrite_family=true);
   delete(fullfile(gcnet_cache, 'Summit_surface.nc'));

   f = @() icemodel.verification.setup.importRetmip( ...
      cache, case_ids=["dye2", "SUM"], output_root=root, ...
      forcing_sources="retmip", build_observations=false, ...
      build_forcing=true, skip_missing=false, gcnet_dir=gcnet_cache, ...
      overwrite=true);

   testCase.verifyError(f, ...
      'icemodel:verification:importRetmip:missingGcnetVandecrux');
   testCase.verifyFalse(isfolder(fullfile(root, 'input', 'met', 'retmip')));
   testCase.verifyFalse(isfolder(fullfile(root, 'input', 'userdata', 'retmip')));
end

function test_retmip_strict_native_only_preflights_all_builders(testCase)
   % Strict native-only refreshes must parse every builder before the first
   % requested case can write runtime artifacts.
   cache = makeRetmipProtocolCache(testCase.TestData.cache, ...
      ["Dye-2_long", "Dye-2_16"]);
   gcnet_cache = makeGcnetSurfaceCache(testCase.TestData.cache);
   samimi_cache = makeMalformedSamimiWorkbookCache(testCase.TestData.cache);
   root = fullfile(testCase.TestData.cache, ...
      'retmip-strict-native-only-builder-root');
   icemodel.verification.setup.importRetmip( ...
      cache, case_ids=["dye2", "dye2_16"], output_root=root, ...
      build_forcing=false, overwrite_family=true);

   f = @() icemodel.verification.setup.importRetmip( ...
      cache, case_ids=["dye2", "dye2_16"], output_root=root, ...
      forcing_sources="retmip", build_observations=false, ...
      build_forcing=true, skip_missing=false, gcnet_dir=gcnet_cache, ...
      samimi_dir=samimi_cache, overwrite=true);

   testCase.verifyError(f, ...
      'icemodel:forcing:buildSamimiDye2Data:missingVariables');
   testCase.verifyFalse(isfolder(fullfile(root, 'input', 'met', 'retmip')));
   testCase.verifyFalse(isfolder(fullfile(root, 'input', 'userdata', 'retmip')));
end

function test_retmip_native_only_refresh_updates_native_metadata(testCase)
   % Native-only refreshes must update metadata and cadence with the files.
   [mar_dir, ~, ~] = requireRcmFixtures(testCase);
   cache = makeRetmipProtocolCache(testCase.TestData.cache, "Dye-2_long");
   gcnet_cache = makeGcnetSurfaceCache(testCase.TestData.cache);
   root = fullfile(testCase.TestData.cache, ...
      'retmip-native-only-metadata-root');
   icemodel.verification.setup.importRetmip( ...
      cache, case_ids="dye2", output_root=root, ...
      forcing_sources="mar", build_forcing=true, mar_dir=mar_dir, ...
      overwrite_family=true);

   refreshed = icemodel.verification.setup.importRetmip( ...
      cache, case_ids="dye2", output_root=root, ...
      forcing_sources="retmip", build_observations=false, ...
      build_forcing=true, gcnet_dir=gcnet_cache, overwrite=true, ...
      overwrite_family=true);

   c = refreshed.cases;
   testCase.verifyEqual(string(c.native_timestep), "1hr");
   testCase.verifyEqual( ...
      string(c.observation_variables.native_source.family), ...
      "gcnet_vandecrux");
   testCase.verifyTrue(any(string( ...
      c.observation_variables.native_source.data_variables) == "smb"));
   testCase.verifyFalse(isfield(c.colocation, 'mar'));
   testCase.verifyFalse(any(string(c.forcing_sources) == "mar3.11"));
end

function test_sumup_forcing_only_rejects_blank_source_tokens(testCase)
   % A true forcing build must name at least one runtime artifact.
   root = fullfile(testCase.TestData.cache, 'sumup-blank-source-root');
   family_root = fullfile(root, 'eval', 'sumup');
   manifest_file = fullfile(family_root, 'manifest.json');
   preview = icemodel.verification.setup.importSumup( ...
      "", points=[67.1, -48.1], case_ids="blank_source", ...
      anchors=struct([]), output_root=root, dry_run=true, ...
      forcing_sources=strings(1, 0), build_forcing=false);
   mkdir(family_root)
   icemodel.verification.setup.writeManifest(manifest_file, preview);

   call = @() icemodel.verification.setup.importSumup( ...
      "", case_ids="blank_source", output_root=root, ...
      forcing_sources=["", ""], build_observations=false, ...
      build_forcing=true);

   testCase.verifyError(call, ...
      'icemodel:verification:normalizeForcingSources:emptySelection');
   testCase.verifyTrue(isfile(manifest_file));
end

function test_import_retmip_dry_run_build_forcing_skips_source_probes(testCase)
   % Dry runs should not require native or RCM caches when forcing is selected.
   protocol_cache = makeRetmipProtocolCache( ...
      testCase.TestData.cache, "Summit");
   missing_native = fullfile(testCase.TestData.cache, 'missing-retmip-native');
   returned = icemodel.verification.setup.importRetmip( ...
      protocol_cache, dry_run=true, build_forcing=true, ...
      forcing_sources=["retmip", "mar", "merra", "racmo"], ...
      promice_dir=missing_native, gcnet_dir=missing_native, ...
      samimi_dir=missing_native, imau_dir=missing_native, mar_dir="/no/mar", ...
      merra_dir="/no/merra", racmo_dir="/no/racmo");

   testCase.verifyGreaterThan(numel(returned.cases), 0);
   statuses = string(arrayfun(@(c) ...
      c.colocation.retmip.native_met_status, returned.cases, ...
      'UniformOutput', false));
   testCase.verifyEqual(unique(statuses), "not_checked");
   summit = returned.cases(string({returned.cases.case_id}) == "summit");
   testCase.verifyEqual(string(summit.colocation.retmip.surface_file), "");
   testCase.verifyEmpty(summit.colocation.retmip.profile_files);
   testCase.verifyEmpty(summit.colocation.retmip.model_output_files);
   testCase.verifyFalse(isfolder(missing_native));
end

function test_import_imau_dry_run_build_forcing_skips_source_probes(testCase)
   % Dry runs should not require native or RCM caches when forcing is selected.
   missing_native = fullfile(testCase.TestData.cache, 'missing-imau-native');
   returned = icemodel.verification.setup.importImau(missing_native, ...
      case_ids="S21", dry_run=true, build_forcing=true, ...
      forcing_sources=["imau", "mar"], mar_dir="/no/mar");

   testCase.verifyEqual(string(returned.cases.case_id), "s21");
   testCase.verifyEmpty( ...
      fieldnames(returned.cases.colocation.imau.cache_status));
   testCase.verifyFalse(isfolder(missing_native));
end

function test_import_imau_half_window_preserves_error_id(testCase)
   % IMAU window validation should match the other importers.
   testCase.verifyError(@() icemodel.verification.setup.importImau( ...
      dry_run=true, startdate="2014-01-01"), ...
      'icemodel:internal:pairedWindow:invalidWindow');
end

function test_native_builders_validate_window_before_missing_source_discovery(testCase)
   % Malformed input must win without creating or probing a caller's cache root.
   missing_root = fullfile(testCase.TestData.cache, 'window-precedence');
   calls = { ...
      @() icemodel.forcing.buildImauHourlyData("S21", ...
      source_dir=fullfile(missing_root, 'imau'), startdate="2014-01-01")
      @() icemodel.forcing.buildGcnetVandecruxData("DYE_2", ...
      source_dir=fullfile(missing_root, 'gcnet-data'), startdate="2014-01-01")
      @() icemodel.forcing.buildGcnetVandecruxFirnTemperature("DYE_2", ...
      source_dir=fullfile(missing_root, 'gcnet-firn'), startdate="2014-01-01")
      @() icemodel.forcing.buildSamimiDye2Data( ...
      source_dir=fullfile(missing_root, 'samimi'), startdate="2014-01-01")
      @() icemodel.verification.setup.buildSumupObservations([67, -48], ...
      source_dir=fullfile(missing_root, 'sumup'), startdate="2014-01-01")};

   for k = 1:numel(calls)
      testCase.verifyError(calls{k}, ...
         'icemodel:internal:pairedWindow:invalidWindow');
   end
   testCase.verifyFalse(isfolder(missing_root));
end

function test_set_model_opts_accepts_firn_period_manifest(testCase)
   % Verification manifests use the canonical period window field.
   input_root = fullfile(testCase.TestData.cache, 'period-input');
   mkdir(input_root);
   c = struct('case_id', "s21", 'dataset_family', "imau", ...
      'period', struct('start', '2014-04-12 00:00:00', ...
      'end', '2016-04-13 00:00:00'), ...
      'native_timestep', 'hourly', 'input_data_root', input_root);

   opts = icemodel.test.helpers.setModelOptsForCase(c);

   testCase.verifyEqual(string(opts.forcings), "imau");
   testCase.verifyEqual(string(opts.pathinput), string(input_root));
   testCase.verifyEqual(string(opts.pathuserdata), ...
      string(fullfile(input_root, 'userdata')));
   testCase.verifyTrue(startsWith(string(opts.metfname{1}), ...
      string(fullfile(input_root, 'met'))));
   testCase.verifyEqual(opts.simyears, 2014:2016);
   legacy_override = icemodel.test.helpers.setModelOptsForCase(c, dt=900);
   testCase.verifyEqual(legacy_override.dt, 900);

   c.case_id = "humphrey";
   c.dataset_family = "research_site";
   c.forcing_sources = {'mar'};
   opts = icemodel.test.helpers.setModelOptsForCase(c);

   testCase.verifyEqual(string(opts.forcings), "mar");

   c.forcing_sources = {};
   c.dataset_family = "research_site";
   testCase.verifyError( ...
      @() icemodel.test.helpers.setModelOptsForCase(c), ...
      'icemodel:test:setModelOptsForCase:noRunnableForcing');

   c.dataset_family = "imau";
   c.forcing_sources = {'imau'};
   c.period = struct('start', '', 'end', '');
   testCase.verifyError( ...
      @() icemodel.test.helpers.setModelOptsForCase(c), ...
      'icemodel:test:setModelOptsForCase:unboundedWindow');

   c.dataset_family = "research_site";
   c.forcing_sources = {'mar3.11'};
   c.period = struct('start', '2011-01-01 00:00:00', ...
      'end', '2014-12-31 23:00:00');
   c.colocation = struct('mar', struct('kind', 'point_met', ...
      'staged', true, ...
      'window', struct('start', '2012-01-01 00:00:00', ...
      'end', '2012-12-31 23:00:00')));
   opts = icemodel.test.helpers.setModelOptsForCase(c);
   testCase.verifyEqual(string(opts.forcings), "mar3.11");
   testCase.verifyEqual(opts.simyears, 2012);

   opts = icemodel.test.helpers.setModelOptsForCase(c, ...
      startdate=datetime(2011, 1, 1, TimeZone="UTC"), ...
      enddate=datetime(2014, 12, 31, 23, 0, 0, TimeZone="UTC"));
   testCase.verifyEqual(opts.simyears, 2012);
   testCase.verifyEqual(opts.startdate, ...
      datetime(2012, 1, 1, TimeZone="UTC"));
   testCase.verifyEqual(opts.enddate, ...
      datetime(2012, 12, 31, 23, 0, 0, TimeZone="UTC"));

   c.colocation.mar.window = struct('start', '2009-01-01 00:00:00', ...
      'end', '2019-12-31 23:00:00');
   opts = icemodel.test.helpers.setModelOptsForCase(c);
   testCase.verifyEqual(opts.simyears, 2011:2014);
   testCase.verifyEqual(opts.startdate, ...
      datetime(2011, 1, 1, TimeZone="UTC"));
   testCase.verifyEqual(opts.enddate, ...
      datetime(2014, 12, 31, 23, 0, 0, TimeZone="UTC"));

   c.forcing_sources = {'imau', 'mar3.11'};
   c.period = struct('start', '2011-01-01 00:00:00', ...
      'end', '2011-01-01 23:00:00');
   c.colocation.imau = struct('staged', true, ...
      'met_files', 'imau/met_s21_imau_20110101_20141231_1hr.mat', ...
      'forcing_ready', false, ...
      'forcing_ready_reason', 'required met placeholder channel(s): ppt');
   imau_met_file = fullfile(input_root, 'met', ...
      c.colocation.imau.met_files);
   mkdir(fileparts(imau_met_file));
   met = icemodel.test.fixtures.makeSyntheticMetFile(2011, ...
      nsteps=24, dt_seconds=3600);
   save(imau_met_file, 'met');
   opts = icemodel.test.helpers.setModelOptsForCase(c);
   testCase.verifyEqual(string(opts.forcings), "imau");

   c.forcing_sources = {'imau'};
   opts = icemodel.test.helpers.setModelOptsForCase(c);
   testCase.verifyEqual(string(opts.forcings), "imau");

   c.dataset_family = "retmip";
   c.native_timestep = '15m';
   c.forcing_sources = {'mar3.11'};
   c.period = struct('start', '2011-01-01 00:00:00', ...
      'end', '2011-01-01 23:30:00');
   missing_met_rel = fullfile('mar3.11', ...
      'met_fa_mar3.11_20170101_20171231_1hr.mat');
   covering_met_rel = fullfile('mar3.11', ...
      'met_fa_mar3.11_20100101_20161231_30m.mat');
   outside_met_rel = fullfile('mar3.11', ...
      'met_fa_mar3.11_20150101_20151231_1hr.mat');
   c.colocation.mar.met_files = ...
      {missing_met_rel, covering_met_rel, outside_met_rel};

   % Even when every recorded path exists, select only the member that covers
   % the runtime; a first non-covering/wrong-cadence member cannot control dt.
   first_met_file = fullfile(input_root, 'met', missing_met_rel);
   covering_met_file = fullfile(input_root, 'met', covering_met_rel);
   outside_met_file = fullfile(input_root, 'met', outside_met_rel);
   mkdir(fileparts(covering_met_file));
   met = icemodel.test.fixtures.makeSyntheticMetFile(2011, ...
      nsteps=48, dt_seconds=1800);
   save(covering_met_file, 'met');
   met = icemodel.test.fixtures.makeSyntheticMetFile(2015, ...
      nsteps=24, dt_seconds=3600);
   save(outside_met_file, 'met');
   met = icemodel.test.fixtures.makeSyntheticMetFile(2017, ...
      nsteps=24, dt_seconds=3600);
   save(first_met_file, 'met');
   userdata_rel = fullfile('mar3.11', ...
      'fa_mar3.11_20110101_20141231_30m.mat');
   userdata_file = fullfile(input_root, 'userdata', userdata_rel);
   mkdir(fileparts(userdata_file));
   Data = timetable(datetime(2011, 1, 1, TimeZone="UTC"), 1, ...
      'VariableNames', {'smb'});
   save(userdata_file, 'Data');
   missing_rel = fullfile('mar3.11', ...
      'fa_mar3.11_20150101_20151231_30m.mat');
   c.colocation.mar.data_files = {missing_rel, userdata_rel};
   opts = icemodel.test.helpers.setModelOptsForCase(c);
   testCase.verifyEqual(opts.dt, 1800);
   testCase.verifyEqual(string(opts.metfname), string(covering_met_file));
   equal_override = icemodel.test.helpers.setModelOptsForCase(c, dt=1800);
   testCase.verifyEqual(equal_override.dt, 1800);
   conflict = [];
   try
      icemodel.test.helpers.setModelOptsForCase(c, dt=900);
   catch exception
      conflict = exception;
   end
   testCase.assertNotEmpty(conflict);
   testCase.verifyEqual(string(conflict.identifier), ...
      "icemodel:test:setModelOptsForCase:metCadenceConflict");
   testCase.verifyTrue(all(contains(string(conflict.message), ...
      ["requested dt=900", "1800 seconds", string(covering_met_file)])));
   loaded_met = icemodel.loadmet(opts);
   testCase.verifyEqual(seconds(median(diff(loaded_met.Time))), 1800);
   expected_userdata = string(fullfile(input_root, 'userdata', ...
      [string(missing_rel), string(userdata_rel)]));
   testCase.verifyEqual(string(opts.userdatafname), expected_userdata);

   % A stale explicit list fails during option resolution and must not degrade
   % to cadence-blind legacy discovery of another site/source sibling.
   delete(first_met_file);
   delete(covering_met_file);
   delete(userdata_file);
   caught = [];
   try
      icemodel.test.helpers.setModelOptsForCase(c);
   catch exception
      caught = exception;
   end
   testCase.assertNotEmpty(caught);
   testCase.verifyEqual(string(caught.identifier), ...
      "icemodel:test:setModelOptsForCase:invalidManifestMetFiles");
   testCase.verifyTrue(contains(string(caught.message), ...
      string(fullfile(input_root, 'met', missing_met_rel))));

   c.dataset_family = "imau";
   c.forcing_sources = {'imau'};
   c = rmfield(c, 'colocation');
   c.period = struct('start', '2014-04-12 00:00:00', ...
      'end', '2014-04-13 00:00:00');
   opts = icemodel.test.helpers.setModelOptsForCase(c, ...
      startdate=datetime(2014, 4, 12, TimeZone="UTC"), ...
      enddate=datetime(2014, 4, 13, TimeZone="UTC"));
   testCase.verifyEqual(string(opts.forcings), "imau");
end

function test_set_model_opts_accepts_sorted_contiguous_manifest_met_list(testCase)
   % Adjacent annual files are sorted by saved support and set the exact cadence.
   first = datetime(2011, 12, 31, [22; 23], 0, 0, 'TimeZone', 'UTC');
   second = datetime(2012, 1, 1, [0; 1], 0, 0, 'TimeZone', 'UTC');
   [c, files] = makeManifestMetListCase(testCase, "adjacent", ...
      {second, first}, first(1), second(end));

   opts = icemodel.test.helpers.setModelOptsForCase(c);

   testCase.verifyEqual(string(opts.metfname), string(files([2, 1])));
   testCase.verifyEqual(opts.dt, 3600);
   loaded = icemodel.loadmet(opts);
   testCase.verifyEqual(loaded.Time, [first; second]);
end

function test_set_model_opts_rejects_overlapping_manifest_met_list(testCase)
   % Duplicate boundary rows make an explicit split list ambiguous.
   first = datetime(2012, 1, 1, [0; 1], 0, 0, 'TimeZone', 'UTC');
   second = datetime(2012, 1, 1, [1; 2], 0, 0, 'TimeZone', 'UTC');
   c = makeManifestMetListCase(testCase, "overlap", ...
      {first, second}, first(1), second(end));

   testCase.verifyError( ...
      @() icemodel.test.helpers.setModelOptsForCase(c), ...
      'icemodel:test:setModelOptsForCase:invalidManifestMetFiles');
end

function test_set_model_opts_rejects_gapped_manifest_met_list(testCase)
   % A missing subannual timestep must fail before loadmet concatenation.
   first = datetime(2012, 1, 1, [0; 1], 0, 0, 'TimeZone', 'UTC');
   second = datetime(2012, 1, 1, [3; 4], 0, 0, 'TimeZone', 'UTC');
   c = makeManifestMetListCase(testCase, "gap", ...
      {first, second}, first(1), second(end));

   testCase.verifyError( ...
      @() icemodel.test.helpers.setModelOptsForCase(c), ...
      'icemodel:test:setModelOptsForCase:invalidManifestMetFiles');
end

function test_set_model_opts_rejects_noncovering_manifest_met_list(testCase)
   % A continuous list that ends before the requested runtime is still invalid.
   first = datetime(2012, 1, 1, [0; 1], 0, 0, 'TimeZone', 'UTC');
   second = datetime(2012, 1, 1, [2; 3], 0, 0, 'TimeZone', 'UTC');
   requested_end = datetime(2012, 1, 1, 4, 0, 0, 'TimeZone', 'UTC');
   c = makeManifestMetListCase(testCase, "noncovering", ...
      {first, second}, first(1), requested_end);

   testCase.verifyError( ...
      @() icemodel.test.helpers.setModelOptsForCase(c), ...
      'icemodel:test:setModelOptsForCase:invalidManifestMetFiles');
end

function test_set_model_opts_rejects_mixed_manifest_met_cadence(testCase)
   % Every member of a split runtime must carry the same saved Time cadence.
   first = datetime(2012, 1, 1, [0; 1], 0, 0, 'TimeZone', 'UTC');
   second = datetime(2012, 1, 1, 2, [0; 30], 0, 'TimeZone', 'UTC');
   c = makeManifestMetListCase(testCase, "mixed-cadence", ...
      {first, second}, first(1), second(end));

   testCase.verifyError( ...
      @() icemodel.test.helpers.setModelOptsForCase(c), ...
      'icemodel:test:setModelOptsForCase:invalidManifestMetFiles');
end

function test_retmip_candidate_adapter_returns_protocol_bundle(testCase)
   % RetMIP manifests compare against retmip_protocol_bundle targets, not the
   % SUMup-style subsurface_profile_bundle used by profile observations.
   time = (datetime(2016, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:1))';
   ice1 = struct('Time', time, 'melt', [0; 1], 'Tsfc', [260; 261]);
   ice2 = struct('T', [260 261; 262 263], 'ro_sno', [350 351; 400 401]);
   opts = struct('dz_thermal', 0.1, 'smbmodel', "icemodel", ...
      'sitename', "retmip_test", 'simyears', 2016);
   manifest = struct( ...
      'case_type', "firn_observational", ...
      'dataset_family', "retmip", ...
      'comparison_variables', ["density", "subsurface_temperature", ...
      "melt", "tsfc"], ...
      'eval_sources', "retmip_protocol");

   candidate = icemodel.verification.candidateFromIcemodelOutput( ...
      ice1, ice2, opts, manifest);

   testCase.verifyEqual(string(candidate.format), "retmip_protocol_bundle");
   testCase.verifyTrue(isfield(candidate.data, 'surface'));
   testCase.verifyTrue(isfield(candidate.data.profiles, 'density'));
   testCase.verifyTrue(isfield(candidate.data.profiles, ...
      'subsurface_temperature'));
   testCase.verifyEqual(candidate.data.surface.tsfc, ice1.Tsfc);
end

function test_import_sumup_dry_run_does_not_write_staging_tree(testCase)
   % SUMup dry runs should return a manifest shape without source reads or
   % staged files.
   sumup_dir = fullfile(testCase.TestData.cache, 'missing-sumup-cache');
   site = icemodel.verification.setup.researchSiteCatalog("humphrey");
   root = fullfile(testCase.TestData.cache, 'sumup-dry-run');

   returned = icemodel.verification.setup.importSumup(sumup_dir, ...
      points=[site.lat_wgs84, site.lon_wgs84], case_ids="humphrey", ...
      startdate="2010-07-01", enddate="2012-04-30", ...
      output_root=root, dry_run=true);

   testCase.verifyEqual(string(returned.dataset_family), "sumup");
   testCase.verifyEqual(string(returned.cases.case_id), "humphrey");
   testCase.verifyEqual(string(returned.cases.eval_sources), "sumup_obs");
   testCase.verifyFalse(isfolder(sumup_dir));
   testCase.verifyFalse(isfolder(fullfile(root, 'eval', 'sumup')));
   testCase.verifyFalse(isfolder(fullfile(root, 'input')));
end

function test_profile_bundle_compare_and_plotcase_support_sumup(testCase)
   % Generic compare/plot workflows should accept SUMup profile bundles.
   sumup_dir = requireSumupCache(testCase);
   root = fullfile(testCase.TestData.cache, 'sumup-profile-workflow');
   point = [69.582, -49.295];   % SUMup Humphrey-profile cluster.

   icemodel.verification.setup.importSumup(sumup_dir, ...
      points=point, case_ids="humphrey", ...
      startdate="2010-07-01", enddate="2012-04-30", ...
      output_root=root, overwrite_family=true);

   loaded = load(fullfile(root, 'eval', 'sumup', 'humphrey', ...
      'observations.mat'), 'targets');
   temperature = loaded.targets.data.subsurface_temperature;
   if (istable(temperature) || istimetable(temperature)) ...
         && ismember("subsurface_temperature", ...
         string(temperature.Properties.VariableNames))
      idx = string(temperature.Properties.VariableNames) ...
         == "subsurface_temperature";
      testCase.verifyEqual(string(temperature.Properties.VariableUnits(idx)), ...
         "degC");
   end
   result = icemodel.verification.comparecase("humphrey", ...
      evaluation_data_root=fullfile(root, 'eval'), ...
      candidate=loaded.targets, make_plot=false);
   profile_rows = ismember(string(result.metrics.variable), ...
      ["density", "subsurface_temperature"]);
   testCase.verifyTrue(any(profile_rows));
   testCase.verifyTrue(all(result.metrics.status(profile_rows) == "ok"));

   f = icemodel.verification.plotcase("humphrey", ...
      evaluation_data_root=fullfile(root, 'eval'), visible="off");
   testCase.verifyTrue(isvalid(f));
   close(f)
end

function test_profile_bundle_uses_subsurface_temperature_column(testCase)
   % SUMup temperature axes use the canonical stored variable name.
   eval_root = fullfile(testCase.TestData.cache, 'profile-alias-eval');
   writeTemperatureAliasProfileCase(eval_root);

   loaded = load(fullfile(eval_root, 'sumup', 'tempalias', ...
      'observations.mat'), 'targets');
   candidate = loaded.targets;
   candidate.data.subsurface_temperature.subsurface_temperature = ...
      candidate.data.subsurface_temperature.subsurface_temperature + 1;

   result = icemodel.verification.comparecase("tempalias", ...
      evaluation_data_root=eval_root, candidate=candidate, make_plot=false);

   row = result.metrics(result.metrics.variable == ...
      "subsurface_temperature", :);
   testCase.verifyEqual(row.rmse, 1, 'AbsTol', 1e-12);

   default_result = icemodel.verification.comparecase("tempalias", ...
      evaluation_data_root=eval_root, make_plot=false);
   testCase.verifyEqual(default_result.metrics.status, ...
      "missing_candidate_variable");

   f = icemodel.verification.plotcase("tempalias", ...
      evaluation_data_root=eval_root, source="compare", visible="off");
   testCase.verifyTrue(isvalid(f));
   close(f)
end

function test_profile_comparison_aligns_by_depth_and_family(testCase)
   % Profile comparisons interpolate candidate values onto target depths and use
   % dataset_family to disambiguate shared case ids.
   eval_root = fullfile(testCase.TestData.cache, 'profile-depth-eval');
   writeDepthProfileCase(eval_root, "sumup", "shared", [300; 500]);
   writeDepthProfileCase(eval_root, "promice", "shared", [1300; 1500]);
   loaded = load(fullfile(eval_root, 'sumup', 'shared', ...
      'observations.mat'), 'targets');
   candidate = loaded.targets;
   candidate.data.density = table([0; 1; 2], [300; 400; 500], ...
      'VariableNames', {'depth', 'density'});

   result = icemodel.verification.comparecase("shared", ...
      evaluation_data_root=eval_root, dataset_family="sumup", ...
      candidate=candidate, make_plot=false);

   row = result.metrics(result.metrics.variable == "density", :);
   testCase.verifyEqual(row.rmse, 0, 'AbsTol', 1e-12);
   testCase.verifyEqual(result.aligned.density.depth, [0; 2]);
   testCase.verifyEqual(string(result.manifest.dataset_family), "sumup");
end

function test_profile_comparison_prefers_depth_column_in_timetable(testCase)
   % SUMup profile timetables carry observation time plus a profile-depth
   % column. Depth must be the comparison axis for profile metrics.
   eval_root = fullfile(testCase.TestData.cache, 'profile-timetable-eval');
   writeDepthProfileCase(eval_root, "sumup", "shared", [300; 500]);
   loaded = load(fullfile(eval_root, 'sumup', 'shared', ...
      'observations.mat'), 'targets');
   candidate = loaded.targets;
   Time = datetime(2020, 1, 1:3, TimeZone="UTC")';
   candidate.data.density = timetable(Time, [0; 1; 2], [300; 400; 500], ...
      'VariableNames', {'midpoint', 'density'});

   result = icemodel.verification.comparecase("shared", ...
      evaluation_data_root=eval_root, dataset_family="sumup", ...
      candidate=candidate, make_plot=false);

   row = result.metrics(result.metrics.variable == "density", :);
   testCase.verifyEqual(row.status, "ok");
   testCase.verifyEqual(row.rmse, 0, 'AbsTol', 1e-12);
   testCase.verifyEqual(result.aligned.density.depth, [0; 2]);
end

function test_dated_profiles_pair_exact_dates_without_depth_collapse(testCase)
   % Repeated depths from different dates remain separate; only an exact UTC
   % calendar-date candidate is interpolated onto the observed midpoints.
   eval_root = fullfile(testCase.TestData.cache, 'dated-profile-eval');
   writeDepthProfileCase(eval_root, "sumup", "dated", [300; 500]);
   observation_file = fullfile(eval_root, 'sumup', 'dated', ...
      'observations.mat');
   loaded = load(observation_file, 'targets');
   first_date = datetime(2020, 1, 1, TimeZone="UTC");
   second_date = datetime(2020, 2, 1, TimeZone="UTC");
   % Match the real SUMup schema: one name_key defines the physical profile,
   % while row-level reference/method provenance may vary with depth.
   loaded.targets.data.density = table( ...
      [1; 1; 2; 2], [10; 11; 20; 21], [1; 2; 3; 4], ...
      repmat(70, 4, 1), repmat(-40, 4, 1), ...
      [first_date; first_date; second_date; second_date], ...
      [0.5; 1.5; 0.5; 1.5], [350; 450; 360; 460], ...
      'VariableNames', {'name_key', 'reference_key', 'method_key', ...
      'latitude', 'longitude', 'datetime', 'midpoint', 'density'});
   targets = loaded.targets;
   save(observation_file, 'targets')

   candidate = loaded.targets;
   candidate.data.density = table( ...
      repmat("mar_a", 3, 1), ...
      repmat(first_date + hours(12), 3, 1), [0; 1; 2], [300; 400; 500], ...
      'VariableNames', {'profile_id', 'datetime', 'depth', 'density'});
   result = icemodel.verification.comparecase("dated", ...
      evaluation_data_root=eval_root, dataset_family="sumup", ...
      candidate=candidate, make_plot=false);

   rows = result.metrics(result.metrics.variable == "density", :);
   testCase.verifyEqual(height(rows), 2)
   testCase.verifyEqual(rows.status, ["ok"; "missing_candidate_date"])
   testCase.verifyEqual(rows.rmse(1), 0, AbsTol=1e-12)
   testCase.verifyEqual(height(result.aligned.density), 2)
   testCase.verifyEqual(unique(result.aligned.density.datetime), first_date)

    f = icemodel.verification.plotcase("dated", ...
       evaluation_data_root=eval_root, dataset_family="sumup", ...
       source="compare", candidate=candidate, visible="off");
   testCase.verifyTrue(isvalid(f))
   testCase.verifyGreaterThanOrEqual(numel(findobj(f, 'Type', 'line')), 3)
   close(f)
end

function test_fetch_imau_accepts_flat_cache(testCase)
   % Flat manual caches use the same hourly/daily filename rules as importImau.
   mkdir(testCase.TestData.cache)
   touch(fullfile(testCase.TestData.cache, 'VanTiggelen-etal_2024_S21.tab'));
   touch(fullfile(testCase.TestData.cache, 'GRL_S21_AWS.tab'));

   [~, status] = icemodel.verification.setup.fetchImau( ...
      cache_dir=testCase.TestData.cache, strict=true, silent=true);

   testCase.verifyTrue(all([status.present]));
end

function test_fetch_promice_reports_required_local_files(testCase)
   % fetchPromice must report the local L3 station product and AWS metadata
   % requirements without attempting a download.
   [returned, status] = icemodel.verification.setup.fetchPromice( ...
      cache_dir=testCase.TestData.cache, strict=false, silent=true);

   expected = string(testCase.TestData.cache);
   testCase.verifyEqual(string(returned), expected);
   testCase.verifyEqual(string({status.product}), ["metadata", "hour"]);
   testCase.verifyFalse(all([status.present]));
   metadata = status(string({status.product}) == "metadata");
   testCase.verifyTrue(ismember("AWS_sites_metadata.csv", ...
      metadata.missing_files));
end

function test_fetch_promice_station_filter_passes_full_cache(testCase)
   % A cache with the requested station NetCDF and required metadata satisfies
   % strict mode; station matching is case/underscore tolerant like the reader.
   cache = makePromiceCache(testCase.TestData.cache, ["KAN_M", "KAN_U", ...
      "DY2"]);

   [~, status] = icemodel.verification.setup.fetchPromice( ...
      cache_dir=cache, stations=["kanm", "KAN_U", "dy2"], ...
      strict=true, silent=true);

   testCase.verifyTrue(all([status.present]));
   hour = status(string({status.product}) == "hour");
   testCase.verifyEmpty(hour.missing_stations);
end

function test_fetch_promice_empty_product_selection_is_source_light(testCase)
   % Selecting no products returns the typed empty contract without touching
   % the cache tree or probing the required PROMICE metadata files.
   [returned, status] = icemodel.verification.setup.fetchPromice( ...
      cache_dir=testCase.TestData.cache, products=strings(1, 0));

   testCase.verifyEqual(string(returned), string(testCase.TestData.cache));
   testCase.verifySize(status, [1, 0]);
   testCase.verifyEqual(string(fieldnames(status)), [ ...
      "product"; "landing_url"; "present"; "cache_dir"; ...
      "missing_files"; "missing_stations"]);
   testCase.verifyFalse(isfolder(testCase.TestData.cache));
end

function test_fetch_gcnet_reports_arctic_data_status(testCase)
   % fetchGcnet must expose the three Arctic Data packages and missing
   % RetMIP-relevant station files without loading any NetCDF arrays.
   [returned, status] = icemodel.verification.setup.fetchGcnet( ...
      cache_dir=testCase.TestData.cache, strict=false, silent=true);

   expected = string(testCase.TestData.cache);
   testCase.verifyEqual(string(returned), expected);
   returned = sort(string({status.product}));
   expected = sort(["surface", "firn_temperature", "simulated_firn"]);
   testCase.verifyEqual(returned, expected);
   testCase.verifyTrue(any(contains([status.landing_url], ...
      "10.18739/A2HM52K87")));
   testCase.verifyTrue(any(contains([status.landing_url], ...
      "10.18739/A2833N00P")));
   testCase.verifyTrue(any(contains([status.landing_url], ...
      "10.18739/A2CV4BS43")));
   surface = status(string({status.product}) == "surface");
   testCase.verifyTrue(ismember("DYE_2_surface.nc", ...
      surface.missing_files));
   testCase.verifyTrue(ismember("Summit_surface.nc", ...
      surface.missing_files));
end

function test_fetch_gcnet_empty_products_are_source_light(testCase)
   % The shared typed empty row must not cause a cache mkdir or source scan.
   [returned, status] = icemodel.verification.setup.fetchGcnet( ...
      cache_dir=testCase.TestData.cache, products=strings(1, 0));

   testCase.verifyEqual(string(returned), string(testCase.TestData.cache));
   testCase.verifySize(status, [1, 0]);
   testCase.verifyEqual(string(fieldnames(status)), [ ...
      "product"; "doi"; "landing_url"; "present"; "cache_dir"; ...
      "found_files"; "stations"; "resolved_files"; "missing_files"; ...
      "ambiguous_files"; "provenance_file"; "missing_provenance"]);
   testCase.verifyFalse(isfolder(testCase.TestData.cache));
end

function test_fetch_gcnet_empty_station_selection_uses_defaults(testCase)
   % Empty station selections should not make strict validation vacuously pass.
   mkdir(testCase.TestData.cache)

   testCase.verifyError(@() icemodel.verification.setup.fetchGcnet( ...
      cache_dir=testCase.TestData.cache, stations=strings(1, 0), ...
      products="surface", strict=true, silent=true), ...
      'icemodel:verification:fetchGcnet:missingSources');
end

function test_fetch_gcnet_accepts_package_layout_and_aliases(testCase)
   % The local cache may be flat or unpacked into DOI package folders; station
   % aliases should still validate the RetMIP Dye-2 and Summit files.
   cache = makeGcnetCache(testCase.TestData.cache);

   aliases = ["dye2", "DYE_2", "DY2", "sum", "SUM", "Summit"];
   for alias = aliases
      [~, returned] = icemodel.verification.setup.fetchGcnet( ...
         cache_dir=cache, stations=alias, strict=true, silent=true);

      testCase.verifyTrue(all([returned.present]));
      station_sets = [returned.stations];
      if ismember(alias, ["dye2", "DYE_2", "DY2"])
         testCase.verifyTrue(any(station_sets == "DYE_2"));
      else
         testCase.verifyTrue(any(station_sets == "Summit"));
      end
   end

   [~, returned] = icemodel.verification.setup.fetchGcnet( ...
      cache_dir=cache, stations=["DY2", "SUM"], ...
      strict=true, silent=true);
   simulated = returned(string({returned.product}) == "simulated_firn");
   testCase.verifyTrue(any(contains(simulated.found_files, ...
      "DYE_2_T_ice_bin_1.nc")));
   testCase.verifyEmpty(simulated.missing_files);
end

function test_fetch_gcnet_accepts_flat_cache_without_xml(testCase)
   % Station NetCDF files are sufficient for native staging; XML is provenance.
   cache = fullfile(testCase.TestData.cache, 'gcnet-flat-no-xml');
   mkdir(cache)
   touch(fullfile(cache, 'DYE_2_surface.nc'));

   [~, returned] = icemodel.verification.setup.fetchGcnet( ...
      cache_dir=cache, stations="dye2", products="surface", ...
      strict=true, silent=true);

   testCase.verifyTrue(logical(returned.present));
   testCase.verifyEmpty(returned.missing_files);
   testCase.verifyEqual(returned.missing_provenance, ...
      "Gap_filled_meteorological_data_and_surface_energy.xml");
end

function test_fetch_gcnet_accepts_normalized_station_filenames(testCase)
   % Fetch validation should match the builder's separator/case-insensitive lookup.
   cache = fullfile(testCase.TestData.cache, 'gcnet-normalized');
   mkdir(cache)
   touch(fullfile(cache, 'dye-2 surface.nc'));

   [~, returned] = icemodel.verification.setup.fetchGcnet( ...
      cache_dir=cache, stations="DYE_2", products="surface", ...
      strict=true, silent=true);

   testCase.verifyTrue(logical(returned.present));
   testCase.verifyEmpty(returned.missing_files);
   testCase.verifyTrue(contains(returned.found_files, "dye-2 surface.nc"));
end

function test_fetch_gcnet_rejects_partial_and_ambiguous_matches(testCase)
   % A containing basename is not the required file, and two exact normalized
   % basenames are ambiguous rather than an invitation to select the first one.
   cache = fullfile(testCase.TestData.cache, 'gcnet-partial-ambiguous');
   mkdir(cache)
   touch(fullfile(cache, 'backup_DYE_2_surface.nc'));

   [~, partial] = icemodel.verification.setup.fetchGcnet( ...
      cache_dir=cache, stations="DYE_2", products="surface", ...
      strict=false, silent=true);
   testCase.verifyFalse(logical(partial.present));
   testCase.verifyEmpty(partial.found_files);
   testCase.verifyEqual(partial.missing_files, "DYE_2_surface.nc");
   testCase.verifyEmpty(partial.ambiguous_files);

   % Separator spelling remains normalized, so these two package copies both
   % match the one expected file and must make strict validation fail.
   first = fullfile(cache, 'package-a');
   second = fullfile(cache, 'package-b');
   mkdir(first)
   mkdir(second)
   touch(fullfile(first, 'DYE_2_surface.nc'));
   touch(fullfile(second, 'dye-2 surface.nc'));
   [~, ambiguous] = icemodel.verification.setup.fetchGcnet( ...
      cache_dir=cache, stations="DYE_2", products="surface", ...
      strict=false, silent=true);

   testCase.verifyFalse(logical(ambiguous.present));
   testCase.verifyEmpty(ambiguous.found_files);
   testCase.verifyEqual(ambiguous.missing_files, "DYE_2_surface.nc");
   testCase.verifyNumElements(ambiguous.ambiguous_files, 2);
   testCase.verifyTrue(all(contains(ambiguous.ambiguous_files, "package-")));
   testCase.verifyError(@() icemodel.verification.setup.fetchGcnet( ...
      cache_dir=cache, stations="DYE_2", products="surface", ...
      strict=true, silent=true), ...
      'icemodel:verification:fetchGcnet:missingSources');
end

function test_fetch_gcnet_keeps_multi_product_files_disjoint(testCase)
   % Each requested product receives only its own uniquely resolved files even
   % though the shared product-row helper collects all NetCDF candidates.
   cache = fullfile(testCase.TestData.cache, 'gcnet-multi-product');
   mkdir(cache)
   touch(fullfile(cache, 'DYE_2_surface.nc'));
   touch(fullfile(cache, 'DYE_2_T_firn_obs.nc'));
   for suffix = ["T_ice_bin_1.nc", "rho_bin_1.nc", ...
         "slwc_bin_1.nc", "compaction_bin_1.nc"]
      touch(fullfile(cache, "DYE_2_" + suffix));
   end

   [~, returned] = icemodel.verification.setup.fetchGcnet( ...
      cache_dir=cache, stations="DYE_2", strict=true, silent=true);
   products = string({returned.product});
   surface = returned(products == "surface");
   observed = returned(products == "firn_temperature");
   simulated = returned(products == "simulated_firn");

   testCase.verifyTrue(all([returned.present]));
   testCase.verifyNumElements(surface.found_files, 1);
   testCase.verifyTrue(endsWith(surface.found_files, "DYE_2_surface.nc"));
   testCase.verifyNumElements(observed.found_files, 1);
   testCase.verifyTrue(endsWith( ...
      observed.found_files, "DYE_2_T_firn_obs.nc"));
   testCase.verifyNumElements(simulated.found_files, 4);
   testCase.verifyTrue(all(contains(simulated.found_files, "_bin_1.nc")));
   testCase.verifyNumElements(simulated.resolved_files, 4);
   testCase.verifyEqual([simulated.resolved_files.station], ...
      repmat("DYE_2", 1, 4));
   testCase.verifyEqual([simulated.resolved_files.suffix], ...
      ["_T_ice_bin_", "_rho_bin_", "_slwc_bin_", "_compaction_bin_"]);
end

function test_fetch_gcnet_rejects_ambiguous_numbered_product_class(testCase)
   % Every simulated-firn suffix class resolves one numeric-bin file. Multiple
   % numbered candidates for one class are ambiguous, not additive products.
   cache = fullfile(testCase.TestData.cache, 'gcnet-numbered-ambiguous');
   mkdir(cache)
   for filename = ["DYE_2_T_ice_bin_1.nc", "DYE_2_T_ice_bin_2.nc", ...
         "DYE_2_rho_bin_1.nc", "DYE_2_slwc_bin_1.nc", ...
         "DYE_2_compaction_bin_1.nc"]
      touch(fullfile(cache, filename));
   end

   [~, returned] = icemodel.verification.setup.fetchGcnet( ...
      cache_dir=cache, stations="DYE_2", products="simulated_firn", ...
      strict=false, silent=true);

   testCase.verifyFalse(logical(returned.present));
   testCase.verifyTrue(ismember( ...
      "DYE_2_T_ice_bin_<index>.nc", returned.missing_files));
   testCase.verifyNumElements(returned.ambiguous_files, 2);
   testCase.verifyNumElements(returned.found_files, 3);
   testCase.verifyTrue(all(contains(returned.found_files, ...
      ["rho_bin", "slwc_bin", "compaction_bin"])));
end

function test_gcnet_station_alias_helper_is_shared(testCase)
   % Fetch, inventory, and builders share the same canonical station tokens.
   returned = icemodel.forcing.helpers.gcnetVandecruxStation( ...
      ["dye2", "Dye-2_long", "SUM", "other"]);

   testCase.verifyEqual(returned, ["DYE_2", "DYE_2", "Summit", "other"]);
end

function test_gcnet_inventory_reads_headers_and_aliases(testCase)
   % gcnetInventory should index station/product metadata from headers and
   % XML provenance without reading large NetCDF data arrays.
   cache = makeGcnetHeaderCache(testCase.TestData.cache, "nested");

   [returned, status] = icemodel.verification.setup.gcnetInventory(cache, ...
      stations=["dye2_long", "sum"]);

   testCase.verifyEqual(returned.stations, ["DYE_2", "Summit"]);
   records = returned.records;
   testCase.verifyEqual(numel(records), 12);
   products = unique(string({records.product}), 'stable');
   expected = ["surface", "firn_temperature", "simulated_firn"];
   testCase.verifyEqual(products, expected);

   surface = records(string({records.station}) == "DYE_2" ...
      & string({records.product}) == "surface");
   testCase.verifyTrue(ismember("Ta_2m", string(surface.variables)));
   testCase.verifyEqual(string(surface.units.Ta_2m), "K");
   testCase.verifyEqual(surface.dimensions.time, 3);
   expected = datetime(2000, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
   testCase.verifyEqual(surface.period.start, expected);
   expected = datetime(2000, 1, 1, 2, 0, 0, 'TimeZone', 'UTC');
   testCase.verifyEqual(surface.period.end, expected);
   testCase.verifyEqual(string(surface.cadence), "1hr");
   testCase.verifyTrue(ismember("dye2", surface.aliases));
   testCase.verifyEqual(surface.site_location.elev_m, 2165);
   testCase.verifyEqual(string(surface.provenance.doi), ...
      "10.18739/A2HM52K87");
   surface_status = status(string({status.product}) == "surface");
   testCase.verifyEqual(string(surface.provenance.xml_file), ...
      string(surface_status.provenance_file));

   simulated = records(string({records.station}) == "Summit" ...
      & string({records.kind}) == "simulated_density");
   testCase.verifyTrue(ismember("rho", string(simulated.variables)));
   testCase.verifyEqual(string(simulated.canonical_variables), "density");
   testCase.verifyEqual(simulated.dimensions.level, 2);
   testCase.verifyEqual(string(simulated.depth.variable), "Depth");
   testCase.verifyEqual(string(simulated.depth.units), "m");
   testCase.verifyEqual(simulated.depth.levels, 2);
   testCase.verifyEqual(simulated.depth.sample_min_m, 1, 'AbsTol', 1e-12);
   testCase.verifyEqual(simulated.depth.sample_max_m, 2.2, 'AbsTol', 1e-12);
   testCase.verifyTrue(contains(simulated.depth.sample_policy, ...
      "full depth array not loaded"));
   testCase.verifyEqual(string(simulated.provenance.doi), ...
      "10.18739/A2CV4BS43");
end

function test_gcnet_inventory_uses_fetch_resolved_identity(testCase)
   % A parent directory containing another station token must not override the
   % fetch layer's exact basename classification for either station.
   cache = fullfile(testCase.TestData.cache, 'gcnet-parent-token');
   dye_dir = fullfile(cache, '0-Summit_surface.nc-container');
   summit_dir = fullfile(cache, 'z-package');
   mkdir(dye_dir)
   mkdir(summit_dir)
   dye_file = fullfile(dye_dir, 'DYE_2_surface.nc');
   summit_file = fullfile(summit_dir, 'Summit_surface.nc');
   writeTinyGcnetNetcdf(dye_file, "Ta_2m", "K", false);
   writeTinyGcnetNetcdf(summit_file, "RH_2m", "%", false);

   [inventory, status] = icemodel.verification.setup.gcnetInventory( ...
      cache, stations=["DYE_2", "Summit"], products="surface");

   testCase.verifyTrue(logical(status.present));
   resolved = status.resolved_files;
   testCase.verifyEqual([resolved.station], ["DYE_2", "Summit"]);
   testCase.verifyEqual([resolved.suffix], ["_surface.nc", "_surface.nc"]);
   testCase.verifyEqual([resolved.filename], ...
      [string(dye_file), string(summit_file)]);
   records = inventory.records;
   summit = records(string({records.station}) == "Summit");
   testCase.verifyEqual(string(summit.filename), string(summit_file));
   testCase.verifyTrue(ismember("RH_2m", string(summit.variables)));
   testCase.verifyFalse(ismember("Ta_2m", string(summit.variables)));
end

function test_gcnet_inventory_flat_and_nested_layouts_match(testCase)
   % Complete flat and DOI-style package layouts must expose identical header,
   % kind, alias, and provenance records after removing root-dependent paths.
   layouts = ["flat", "nested"];
   records = cell(size(layouts));
   for k = 1:numel(layouts)
      cache = makeGcnetHeaderCache(testCase.TestData.cache, layouts(k));
      [inventory, status] = icemodel.verification.setup.gcnetInventory( ...
         cache, stations=["dye2_long", "sum"]);
      testCase.verifyTrue(all([status.present]));
      testCase.verifyNumElements(inventory.records, 12);
      records{k} = portableGcnetRecords(inventory.records);
   end

   testCase.verifyEqual(records{1}, records{2});
end

function test_gcnet_inventory_ignores_partial_and_ambiguous_files(testCase)
   % Tolerant inventory must not revive a containing basename or choose between
   % two normalized copies that fetch status has deliberately left unresolved.
   cache = fullfile(testCase.TestData.cache, 'gcnet-inventory-ambiguous');
   first = fullfile(cache, 'package-a');
   second = fullfile(cache, 'package-b');
   mkdir(first)
   mkdir(second)
   touch(fullfile(cache, 'backup_DYE_2_surface.nc'));
   touch(fullfile(first, 'DYE_2_surface.nc'));
   touch(fullfile(second, 'dye-2 surface.nc'));

   [inventory, status] = icemodel.verification.setup.gcnetInventory( ...
      cache, stations="DYE_2", products="surface");

   testCase.verifyFalse(logical(status.present));
   testCase.verifyNumElements(status.ambiguous_files, 2);
   testCase.verifyEmpty(status.found_files);
   testCase.verifyEmpty(inventory.records);
end

function test_gcnet_inventory_keeps_unique_part_of_incomplete_product(testCase)
   % A missing Summit file keeps product status false while the unique DYE_2
   % header remains useful to callers explicitly inspecting a partial cache.
   cache = fullfile(testCase.TestData.cache, 'gcnet-inventory-partial');
   mkdir(cache)
   writeTinyGcnetNetcdf(fullfile(cache, 'DYE_2_surface.nc'), ...
      ["Ta_2m", "RH_2m"], ["K", "%"], false);

   [inventory, status] = icemodel.verification.setup.gcnetInventory( ...
      cache, stations=["DYE_2", "Summit"], products="surface");

   testCase.verifyFalse(logical(status.present));
   testCase.verifyEqual(status.missing_files, "Summit_surface.nc");
   testCase.verifyNumElements(inventory.records, 1);
   testCase.verifyEqual(string(inventory.records.station), "DYE_2");
   testCase.verifyEqual(string(inventory.records.product), "surface");
   testCase.verifyTrue(ismember("Ta_2m", ...
      string(inventory.records.variables)));
end

function test_gcnet_inventory_skips_ambiguous_simulated_class(testCase)
   % One ambiguous temperature class must not suppress the three other unique
   % simulated-firn classes or leak either arbitrary temperature file to records.
   cache = fullfile(testCase.TestData.cache, 'gcnet-inventory-simulated');
   first = fullfile(cache, 'package-a');
   second = fullfile(cache, 'package-b');
   mkdir(first)
   mkdir(second)
   touch(fullfile(first, 'DYE_2_T_ice_bin_1.nc'));
   touch(fullfile(second, 'dye-2 t-ice-bin-1.nc'));
   writeTinyGcnetNetcdf(fullfile(cache, 'DYE_2_rho_bin_1.nc'), ...
      ["rho", "Depth"], ["kg m-3", "m"], true);
   writeTinyGcnetNetcdf(fullfile(cache, 'DYE_2_slwc_bin_1.nc'), ...
      ["slwc", "Depth"], ["m", "m"], true);
   writeTinyGcnetNetcdf(fullfile(cache, 'DYE_2_compaction_bin_1.nc'), ...
      ["compaction", "Depth"], ["m", "m"], true);

   [inventory, status] = icemodel.verification.setup.gcnetInventory( ...
      cache, stations="DYE_2", products="simulated_firn");

   testCase.verifyFalse(logical(status.present));
   testCase.verifyNumElements(status.ambiguous_files, 2);
   testCase.verifyNumElements(inventory.records, 3);
   kinds = sort(string({inventory.records.kind}));
   expected = sort(["simulated_density", "simulated_liquid_water", ...
      "simulated_compaction"]);
   testCase.verifyEqual(kinds, expected);
   testCase.verifyFalse(any(kinds == "simulated_temperature"));
end

function test_import_promice_accepts_positional_source_dir(testCase)
   % importPromiceSites(source_dir, case_ids=...) must parse the source positionally
   % while retaining skip-missing behavior for absent station files.
   cache = makePromiceCache(testCase.TestData.cache, "KAN_U");
   eval_root = fullfile(testCase.TestData.cache, 'eval');
   input_root = fullfile(testCase.TestData.cache, 'input');
   mkdir(eval_root);
   mkdir(input_root);

   returned = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.importPromiceSites(cache, ...
      case_ids="KAN_M", evaluation_data_root=eval_root, ...
      input_data_root=input_root, icemodel_config_casename="", ...
      overwrite_family=true), ...
      'icemodel:verification:importPromiceSites:caseSkipped');

   testCase.verifyEmpty(returned.cases);
   testCase.verifyEqual(string(returned.skipped.site), "kanm");
end

function test_import_promice_keeps_promice_dir_alias(testCase)
   % The existing promice_dir name-value path remains an alias for source_dir.
   cache = makePromiceCache(testCase.TestData.cache, "KAN_U");
   eval_root = fullfile(testCase.TestData.cache, 'eval');
   input_root = fullfile(testCase.TestData.cache, 'input');
   mkdir(eval_root);
   mkdir(input_root);

   returned = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", promice_dir=cache, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      icemodel_config_casename="", overwrite_family=true), ...
      'icemodel:verification:importPromiceSites:caseSkipped');

   testCase.verifyEmpty(returned.cases);
   testCase.verifyEqual(string(returned.skipped.site), "kanm");
end

function test_import_promice_skip_missing_suppresses_fetch_banner(testCase)
   % Optional source misses skip the site without printing retrieval guidance.
   cache = fullfile(testCase.TestData.cache, 'missing-promice');
   output_root = fullfile(testCase.TestData.cache, 'quiet-output');
   f = @() icemodel.verification.setup.importPromiceSites(cache, ...
      case_ids="KAN_M", build_forcing=false, skip_missing=true, ...
      output_root=output_root, overwrite_family=true); %#ok<NASGU>

   output = evalc('manifest = f();');

   testCase.verifyFalse(contains(output, ...
      "PROMICE source cache incomplete"));
   testCase.verifyEmpty(manifest.cases);
end

function test_import_promice_strict_missing_prints_retrieval_guidance(testCase)
   % Fatal source misses keep the actionable PROMICE portal banner visible.
   cache = fullfile(testCase.TestData.cache, 'missing-promice-strict');
   output_root = fullfile(testCase.TestData.cache, 'strict-output');
   f = @() icemodel.verification.setup.importPromiceSites(cache, ...
      case_ids="KAN_M", build_forcing=false, skip_missing=false, ...
      output_root=output_root, overwrite_family=true); %#ok<NASGU>

   output = evalc([ ...
      'testCase.verifyError(f, ' ...
      '''icemodel:verification:fetchPromice:missingSources'');']);

   testCase.verifyTrue(contains(output, ...
      "PROMICE source cache incomplete"));
   testCase.verifyTrue(contains(output, "https://promice.org"));
end

function test_import_promice_dry_run_does_not_write_staging_tree(testCase)
   % PROMICE dry runs should select cases without source reads or staged files.
   promice_dir = fullfile(testCase.TestData.cache, 'missing-promice-cache');
   root = fullfile(testCase.TestData.cache, 'promice-dry-run');

   returned = icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_U", promice_dir=promice_dir, ...
      output_root=root, dry_run=true);

   testCase.verifyEqual(string(returned.dataset_family), "promice");
   testCase.verifyEqual(string(returned.cases.case_id), "kanu");
   testCase.verifyEqual(string(returned.cases.eval_sources), "promice_obs");
   testCase.verifyFalse(isfolder(promice_dir));
   testCase.verifyFalse(isfolder(fullfile(root, 'eval', 'promice')));
   testCase.verifyFalse(isfolder(fullfile(root, 'input')));
end

function test_fetch_retmip_full_cache_passes(testCase)
   % A cache with one file per requested RetMIP product satisfies strict mode.
   for product = ["forcing", "outputs", "scripts"]
      folder = fullfile(testCase.TestData.cache, product);
      mkdir(folder);
      touch(fullfile(folder, product + ".txt"));
   end

   [~, returned] = icemodel.verification.setup.fetchRetmip( ...
      cache_dir=testCase.TestData.cache, strict=true, silent=true);

   testCase.verifyTrue(all([returned.present]));
end

function test_fetch_retmip_accepts_flat_official_outputs(testCase)
   % Official flat NetCDF output filenames should satisfy outputs validation.
   cache = fullfile(testCase.TestData.cache, 'retmip-flat-output-cache');
   mkdir(cache);
   touch(fullfile(cache, 'RetMIP_GEUS_KANU_3hourly_values.nc'));

   [~, returned] = icemodel.verification.setup.fetchRetmip( ...
      cache_dir=cache, products="outputs", strict=true, silent=true);

   testCase.verifyTrue(logical(returned.present));
end

function test_retmip_import_accepts_flat_cache_layout(testCase)
   % Import should accept the flat cache layout that fetchRetmip validates.
   cache = fullfile(testCase.TestData.cache, 'retmip-flat-cache');
   mkdir(cache);
   writeText(fullfile(cache, 'RetMIP_surface_forcing_KANU.tab'), ...
      sprintf(['time;melt_mmweq;acc_subl_mmweq;Tsurf_K\n' ...
      '"01-May-2012 00:00:00;0;0;260"\n' ...
      '"01-May-2012 03:00:00;1;2;261"\n']));
   touch(fullfile(cache, 'RetMIP_outputs_KANU_3hourly_values.nc'));
   touch(fullfile(cache, 'scripts.txt'));
   eval_root = fullfile(testCase.TestData.cache, 'eval-flat-retmip');
   input_root = fullfile(testCase.TestData.cache, 'input-flat-retmip');

   returned = icemodel.verification.setup.importRetmip(cache, ...
      case_ids="kanu", evaluation_data_root=eval_root, ...
      input_data_root=input_root, overwrite_family=true);

   testCase.verifyEqual(string(returned.cases.case_id), "kanu");
   testCase.verifyTrue(isfile(fullfile(eval_root, 'retmip', ...
      'kanu', 'observations.mat')));
end

function test_retmip_case_catalog_uses_protocol_table_and_aliases(testCase)
   % RetMIP catalog rows should use protocol coordinates, windows, and ids.
   returned = icemodel.verification.setup.retmipCaseCatalog( ...
      ["Dye-2_16", "sum"]);

   ids = string({returned.case_id});
   testCase.verifyEqual(ids, ["dye2_2016", "summit"]);
   dye = returned(ids == "dye2_2016");
   testCase.verifyEqual(string(dye.retmip_station_id), "dye2_16");
   testCase.verifyEqual(string(dye.protocol_id), "Dye-2_16");
   testCase.verifyEqual(dye.site_location.lat_wgs84, 66.48001);
   testCase.verifyEqual(dye.site_location.lon_wgs84, -46.27889);
   testCase.verifyEqual(dye.site_location.elev_m, 2165);
   testCase.verifyEqual(string(dye.period.start), "2016-05-02 00:00:00");
   testCase.verifyEqual(string(dye.period.end), "2016-10-28 09:00:00");
   testCase.verifyEqual(string(dye.surface_zone), "percolation");
   summit = returned(ids == "summit");
   testCase.verifyEqual(string(summit.protocol_id), "Summit");
   testCase.verifyEqual(summit.site_location.elev_m, 3254);
   testCase.verifyEqual(string(summit.surface_zone), "accumulation");

   aliases = icemodel.verification.setup.retmipCaseCatalog(["DY2", "SUM"]);
   testCase.verifyEqual(string({aliases.case_id}), ["dye2_long", "summit"]);
   dye_long = aliases(1);
   testCase.verifyEqual(string(dye_long.source_association.family), ...
      "gcnet_vandecrux");
   testCase.verifyEqual(string(dye_long.source_association.source_id), ...
      "DYE_2");
   testCase.verifyEqual(string(dye_long.surface_zone), "percolation");
end

function test_promice_catalog_omits_runtime_selector_aliases(testCase)
   % The public site catalog must not expose obsolete forcing-selector fields.
   catalog = icemodel.verification.setup.promiceSiteCatalog();

   testCase.verifyFalse(isfield(catalog, 'models'));
   testCase.verifyFalse(isfield(catalog, 'rcm_models'));
end

function test_imau_site_catalog_uses_source_accumulation_classes(testCase)
   % PANGAEA identifies S21 as accumulation-zone and S22/S23 as ablation-zone;
   % the manifest inventory must retain those authoritative distinctions.
   returned = icemodel.verification.setup.imauSiteCatalog();

   ids = string({returned.site_id});
   zones = string({returned.surface_zone});
   testCase.verifyEqual(ids, ["S21", "S22", "S23"]);
   testCase.verifyEqual(zones, ["accumulation", "ablation", "ablation"]);
   s21 = returned(ids == "S21");
   testCase.verifyEqual(string(s21.source_association.family), "retmip");
   testCase.verifyEqual(string(s21.source_association.source_id), "fa");
end

function test_retmip_import_dry_run_manifest(testCase)
   % importRetmip dry-run should build the five protocol cases and source links
   % without writing empty observations.
   missing_retmip = fullfile(testCase.TestData.cache, 'missing-retmip-dry');
   missing_gcnet = fullfile(testCase.TestData.cache, 'missing-gcnet-dry');
   returned = icemodel.verification.setup.importRetmip( ...
      missing_retmip, dry_run=true, gcnet_dir=missing_gcnet);

   ids = string({returned.cases.case_id});
   expected = ["kanu", "dye2_long", "dye2_2016", "summit", "fa"];
   testCase.verifyEqual(ids, expected);
   fa = returned.cases(ids == "fa");
   expected = "imau";
   testCase.verifyEqual(string(fa.colocation.source_association.family), ...
      expected);
   expected = "S21";
   testCase.verifyEqual(string(fa.colocation.source_association.source_id), ...
      expected);
   expected = "3hr";
   testCase.verifyEqual(string(fa.native_timestep), expected);
   kanu = returned.cases(ids == "kanu");
   testCase.verifyEmpty(kanu.forcing_sources);
   testCase.verifyEqual(string(kanu.eval_sources), "retmip_protocol");
   testCase.verifyEmpty(fieldnames(kanu.colocation.retmip.cache_status));
   testCase.verifyEqual(kanu.site_location.lat_wgs84, 67.0003);
   testCase.verifyEqual(string(kanu.period.end), "2016-12-31 06:00:00");
   testCase.verifyEqual(string(kanu.comparison_variables(:))', ...
      ["tsfc", "melt", "snowf_subl"]);
   summit = returned.cases(ids == "summit");
   expected = "gcnet_vandecrux";
   testCase.verifyEqual( ...
      string(summit.colocation.source_association.family), expected);
   expected = "Summit";
   testCase.verifyEqual( ...
      string(summit.colocation.source_association.source_id), expected);
   expected = "forcing_disabled";
   testCase.verifyEqual(string(summit.colocation.retmip.native_met_status), ...
      expected);
   testCase.verifyFalse(logical(summit.colocation.native_met.staged));
   testCase.verifyFalse(isfolder(missing_retmip));
   testCase.verifyFalse(isfolder(missing_gcnet));
end

function test_retmip_import_stages_protocol_bundle(testCase)
   % A non-dry RetMIP import stages a real protocol observations.mat bundle
   % from local surface/profile files without creating a normal met timetable.
   cache = fullfile(testCase.TestData.cache, 'retmip-cache');
   forcing = fullfile(cache, 'forcing');
   outputs = fullfile(cache, 'outputs');
   scripts = fullfile(cache, 'scripts');
   mkdir(forcing);
   mkdir(outputs);
   mkdir(scripts);
   writeText(fullfile(forcing, 'RetMIP_surface_forcing_KANU.tab'), ...
      sprintf(['time;melt_mmweq;acc_subl_mmweq;Tsurf_K\n' ...
      '"01-May-2012 00:00:00;0;0;260"\n' ...
      '"01-May-2012 03:00:00;0;0;261"\n']));
   writeText(fullfile(forcing, 'RetMIP_initial_firn_density_KAN-U.tab'), ...
      sprintf(['depth_m;density_kgm3\n' ...
      '"0.1;350"\n']));
   nccreate(fullfile(outputs, 'RetMIP_GEUS_KANU_3hourly_values.nc'), ...
      'temp', 'Dimensions', {'time', 2});
   nccreate(fullfile(outputs, 'RetMIP_GEUS_KANU_3hourly_values.nc'), ...
      'rho', 'Dimensions', {'time', 2});
   touch(fullfile(scripts, 'README.md'));

   eval_root = fullfile(testCase.TestData.cache, 'eval');
   input_root = fullfile(testCase.TestData.cache, 'input-protocol-only');
   missing_promice = fullfile(testCase.TestData.cache, 'missing-promice');
   mkdir(eval_root);
   returned = icemodel.verification.setup.importRetmip(cache, ...
      case_ids="kanu", dry_run=false, evaluation_data_root=eval_root, ...
      input_data_root=input_root, promice_dir=missing_promice, ...
      overwrite_family=true);

   expected = "kanu/observations.mat";
   testCase.verifyEqual(string(returned.cases.evaluation_file), expected);
   observations_file = fullfile(eval_root, 'retmip', expected);
   testCase.verifyTrue(isfile(observations_file));
   testCase.verifyFalse(isfolder(input_root));
   testCase.verifyTrue(logical(returned.cases.colocation.retmip.staged));
   testCase.verifyEqual( ...
      string(returned.cases.colocation.retmip.native_met_status), ...
      "forcing_disabled");
   testCase.verifyEqual( ...
      string(returned.cases.colocation.native_met.reason), ...
      "native forcing disabled because build_forcing=false");
   testCase.verifyTrue(ismember("density", ...
      string(returned.cases.comparison_variables)));
   testCase.verifyFalse(ismember("temp", ...
      string(returned.cases.comparison_variables)));
   testCase.verifyFalse(ismember("rho", ...
      string(returned.cases.comparison_variables)));
   testCase.verifyTrue(ismember("temp", ...
      string(returned.cases.colocation.retmip.model_output_variables)));
   testCase.verifyTrue(ismember("rho", ...
      string(returned.cases.colocation.retmip.model_output_variables)));

   % A repeated protocol import reuses the current atomic observation bundle.
   original_observations = fileBytes(observations_file);
   icemodel.verification.setup.importRetmip(cache, ...
      case_ids="kanu", dry_run=false, evaluation_data_root=eval_root, ...
      input_data_root=input_root, promice_dir=missing_promice, ...
      overwrite=false);
   testCase.verifyEqual(fileBytes(observations_file), original_observations);

   loaded = load(observations_file, 'targets');
   expected_product = string( ...
      returned.cases.observation_variables.protocol_id);
   testCase.verifyEqual(string(loaded.targets.metadata.source_family), ...
      "retmip");
   testCase.verifyEqual(string(loaded.targets.metadata.station), ...
      string(returned.cases.observation_variables.retmip_station_id));
   testCase.verifyEqual(string(loaded.targets.metadata.doi), ...
      "10.22008/FK2/GZ3CSN");
   testCase.verifyEqual(string(loaded.targets.metadata.product), ...
      expected_product);
   testCase.verifyEqual(string(loaded.targets.metadata.schema), ...
      "retmip_protocol_bundle");

   % A concrete fixed-bundle identity conflict forces a same-window rewrite
   % instead of relabeling stale protocol bytes through the fresh manifest.
   loaded.targets.metadata.product = 'stale_protocol';
   targets = loaded.targets;
   save(observations_file, 'targets');
   tampered = fileBytes(observations_file);
   icemodel.verification.setup.importRetmip(cache, ...
      case_ids="kanu", dry_run=false, evaluation_data_root=eval_root, ...
      input_data_root=input_root, promice_dir=missing_promice, ...
      overwrite=false);
   testCase.verifyNotEqual(fileBytes(observations_file), tampered);
   loaded = load(observations_file, 'targets');
   testCase.verifyEqual(string(loaded.targets.metadata.product), ...
      expected_product);

   result = icemodel.verification.comparecase("kanu", ...
      evaluation_data_root=eval_root, candidate=loaded.targets, ...
      make_plot=false);
   testCase.verifyTrue(any(result.metrics.status == "ok"));

   default_result = icemodel.verification.comparecase("kanu", ...
      evaluation_data_root=eval_root, make_plot=false);
   testCase.verifyTrue(any(default_result.metrics.status == ...
      "missing_candidate_variable"));

   candidate = loaded.targets;
   target_surface = loaded.targets.data.surface;
   Time = [target_surface.time(1); target_surface.time(1) + hours(1); ...
      target_surface.time(2)];
   candidate.data.surface = timetable(Time, ...
      [target_surface.tsfc(1); 999; target_surface.tsfc(2)], ...
      'VariableNames', {'tsfc'});
   result = icemodel.verification.comparecase("kanu", ...
      evaluation_data_root=eval_root, candidate=candidate, make_plot=false);
   row = result.metrics(result.metrics.variable == "tsfc", :);
   testCase.verifyEqual(row.rmse, 0, 'AbsTol', 1e-12);

   f = icemodel.verification.plotcase("kanu", ...
      evaluation_data_root=eval_root, source="compare", visible="off");
   testCase.verifyTrue(isvalid(f));
   close(f)

   f = icemodel.verification.plotcase("kanu", ...
      evaluation_data_root=eval_root, variables="tsfc", visible="off");
   lines = findobj(f, 'Type', 'line');
   testCase.verifyNotEmpty(lines);
   testCase.verifyTrue(isdatetime(lines(1).XData));
   close(f)
end

function test_retmip_import_stages_gcnet_native_sources(testCase)
   % Dye-2-long and Summit use Vandecrux GC-Net as native RetMIP met/userdata,
   % while the RetMIP protocol bundle remains the evaluation product.
   cache = makeRetmipProtocolCache(testCase.TestData.cache, ...
      ["Dye-2_long", "Summit"]);
   gcnet_cache = makeGcnetSurfaceCache(testCase.TestData.cache);
   % Mirror the production DYE-2 archive's 01:00 first sample so the manifest
   % must distinguish native support from the 00:00 protocol period.
   dye_file = fullfile(gcnet_cache, 'DYE_2_surface.nc');
   delete(dye_file);
   writeGcnetSurfaceFile(dye_file, ...
      datetime(2000, 1, 1, 1, 0, 0, 'TimeZone', 'UTC'), 0);
   eval_root = fullfile(testCase.TestData.cache, 'eval-gcnet');
   input_root = fullfile(testCase.TestData.cache, 'input-gcnet');
   mkdir(eval_root);
   mkdir(input_root);

   returned = icemodel.verification.setup.importRetmip(cache, ...
      case_ids=["dye2", "SUM"], dry_run=false, gcnet_dir=gcnet_cache, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      forcing_sources="retmip", build_forcing=true, ...
      overwrite_family=true);

   ids = string({returned.cases.case_id});
   testCase.verifyEqual(ids, ["dye2_long", "summit"]);
   for n = 1:numel(returned.cases)
      one = returned.cases(n);
      testCase.verifyEqual(string(one.forcing_sources), "retmip");
      testCase.verifyEqual(string(one.eval_sources), "retmip_protocol");
      testCase.verifyEqual( ...
         string(one.colocation.source_association.family), ...
         "gcnet_vandecrux");
      testCase.verifyEqual( ...
         string(one.colocation.retmip.native_met_status), "staged");
      testCase.verifyFalse(logical(one.colocation.retmip.forcing_ready));
      testCase.verifyTrue(contains( ...
         string(one.colocation.retmip.forcing_ready_reason), "ppt"));
      testCase.verifyTrue(isfield(one.colocation.retmip, ...
         'forcing_complete_windows'));
      testCase.verifyEmpty(one.colocation.retmip.forcing_complete_windows);
      testCase.verifyTrue(logical(one.colocation.native_met.staged));
      testCase.verifyEmpty( ...
         one.colocation.native_met.forcing_complete_windows);
      testCase.verifyFalse(isfield(one.colocation.native_met, 'met_files'));
      testCase.verifyFalse(isfield(one.colocation.native_met, 'data_files'));
      testCase.verifyTrue(contains( ...
         string(one.colocation.retmip.native_source.lwd_policy), ...
         "source-filled"));
      testCase.verifyTrue(contains( ...
         string(one.colocation.retmip.native_source.precip_policy), ...
         "ppt = NaN placeholder"));
      testCase.verifyTrue(any(string( ...
         one.observation_variables.native_source.data_variables) == "smb"));
      testCase.verifyFalse(any(string(one.comparison_variables) == "smb"));

      met_file = fullfile(input_root, 'met', ...
         one.colocation.retmip.met_files(1));
      data_file = fullfile(input_root, 'userdata', ...
         one.colocation.retmip.data_files(1));
      testCase.verifyTrue(isfile(met_file));
      testCase.verifyTrue(isfile(data_file));
      saved_met = load(met_file, 'met');
      testCase.verifyTrue(logical( ...
         saved_met.met.Properties.UserData.fillwithmissing));
      saved_data = load(data_file, 'Data');
      testCase.verifyEqual(string(one.colocation.retmip.window.start), ...
         string(icemodel.verification.setup.formatManifestTime( ...
         saved_data.Data.Time(1))));
      testCase.verifyEqual(string(one.colocation.retmip.window.end), ...
         string(icemodel.verification.setup.formatManifestTime( ...
         saved_data.Data.Time(end))));
   end

   summit = returned.cases(ids == "summit");
   testCase.verifyNotEqual( ...
      string(summit.colocation.retmip.native_met_status), ...
      "pending_gcnet_import");
end

function test_retmip_import_records_missing_gcnet_when_skipping(testCase)
   % Missing Vandecrux files should not fabricate a native met leg, and strict
   % callers can still request a hard failure.
   cache = makeRetmipProtocolCache(testCase.TestData.cache, "Summit");
   missing_gcnet = fullfile(testCase.TestData.cache, 'empty-gcnet');
   mkdir(missing_gcnet);
   eval_root = fullfile(testCase.TestData.cache, 'eval-missing-gcnet');
   input_root = fullfile(testCase.TestData.cache, 'input-missing-gcnet');
   mkdir(eval_root);
   mkdir(input_root);

   returned = icemodel.verification.setup.importRetmip(cache, ...
      case_ids="Summit", dry_run=false, gcnet_dir=missing_gcnet, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      forcing_sources="retmip", build_forcing=true, ...
      skip_missing=true, overwrite_family=true);

   testCase.verifyEmpty(returned.cases.forcing_sources);
   testCase.verifyEqual(string(returned.cases.eval_sources), ...
      "retmip_protocol");
   testCase.verifyEqual( ...
      string(returned.cases.colocation.retmip.native_met_status), ...
      "missing_gcnet_vandecrux");
   testCase.verifyFalse(logical(returned.cases.colocation.native_met.staged));

   f = @() icemodel.verification.setup.importRetmip(cache, ...
      case_ids="Summit", dry_run=false, gcnet_dir=missing_gcnet, ...
      evaluation_data_root=fullfile(testCase.TestData.cache, 'eval-fail'), ...
      input_data_root=fullfile(testCase.TestData.cache, 'input-fail'), ...
      forcing_sources="retmip", build_forcing=true, ...
      skip_missing=false, overwrite_family=true);
   testCase.verifyError(f, ...
      'icemodel:verification:importRetmip:missingGcnetVandecrux');
end

function test_retmip_import_skips_missing_protocol_case(testCase)
   % Missing primary RetMIP protocol files should honor skip_missing=true and
   % still write a manifest with an explicit skipped record.
   cache = fullfile(testCase.TestData.cache, 'empty-retmip-cache');
   mkdir(fullfile(cache, 'forcing'));
   mkdir(fullfile(cache, 'outputs'));
   mkdir(fullfile(cache, 'scripts'));
   root = fullfile(testCase.TestData.cache, 'retmip-missing-root');

   returned = testCase.verifyWarning( ...
      @() icemodel.verification.setup.importRetmip(cache, ...
      case_ids="summit", output_root=root, dry_run=false, ...
      skip_missing=true, overwrite_family=true), ...
      'icemodel:verification:importRetmip:caseSkipped');

   testCase.verifyEmpty(returned.cases);
   testCase.verifyEqual(string(returned.skipped.site), "summit");
   testCase.verifyTrue(contains(string(returned.skipped.reason), ...
      "missing RetMIP surface protocol file"));
   testCase.verifyTrue(isfile(fullfile(root, 'eval', 'retmip', ...
      'manifest.json')));
end

function test_retmip_skipped_case_preserves_existing_case_root(testCase)
   % A skippable missing protocol file must not clear stale case artifacts when
   % overwrite=true is requested for cases that are actually staged.
   cache = fullfile(testCase.TestData.cache, 'partial-retmip-cache');
   mkdir(fullfile(cache, 'forcing'));
   mkdir(fullfile(cache, 'outputs'));
   mkdir(fullfile(cache, 'scripts'));
   root = fullfile(testCase.TestData.cache, 'retmip-stale-root');
   stale_file = fullfile(root, 'eval', 'retmip', 'summit', 'stale.txt');
   mkdir(fileparts(stale_file));
   writeText(stale_file, "keep");

   returned = testCase.verifyWarning( ...
      @() icemodel.verification.setup.importRetmip(cache, ...
      case_ids="summit", output_root=root, dry_run=false, overwrite=true, ...
      skip_missing=true), ...
      'icemodel:verification:importRetmip:caseSkipped');

   testCase.verifyEmpty(returned.cases);
   testCase.verifyTrue(isfile(stale_file));
end

function test_retmip_malformed_surface_preserves_existing_family_root(testCase)
   % A malformed surface table is not a missing-data skip and must fail before
   % overwrite preparation changes any existing family-root byte.
   root = fullfile(testCase.TestData.cache, 'retmip-malformed-surface');
   cache = makeRetmipProtocolCache(fullfile(root, 'source'), "KANU");
   writeText(fullfile(cache, 'forcing', ...
      'RetMIP_surface_forcing_KANU.tab'), ...
      sprintf(['time;melt_mmweq;acc_subl_mmweq;Tsurf_K\n' ...
      '"01-May-2012 00:00:00;0;0;260"\n' ...
      '"01-May-2012 02:00:00;0;0;261"\n']));
   eval_root = fullfile(root, 'eval');
   family_root = fullfile(eval_root, 'retmip');
   stale_file = fullfile(family_root, 'kanu', 'nested', 'sentinel.bin');
   mkdir(fileparts(stale_file));
   writeText(stale_file, "keep-surface");
   before = treeSnapshot(family_root);

   f = @() icemodel.verification.setup.importRetmip(cache, ...
      case_ids="kanu", evaluation_data_root=eval_root, ...
      input_data_root=fullfile(root, 'input'), overwrite=true, ...
      skip_missing=true);

   testCase.verifyError(f, ...
      'icemodel:verification:readRetmipProtocolTable:badCadence');
   testCase.verifyEqual(treeSnapshot(family_root), before);
end

function test_retmip_malformed_later_profile_preserves_existing_root(testCase)
   % Every discovered profile must parse in memory before case-root preparation;
   % a valid first profile cannot hide a malformed later profile.
   root = fullfile(testCase.TestData.cache, 'retmip-malformed-profile');
   cache = makeRetmipProtocolCache(fullfile(root, 'source'), "KANU");
   writeText(fullfile(cache, 'forcing', ...
      'RetMIP_initial_firn_density_KAN-U.tab'), ...
      sprintf('depth_m;density_kgm3\n"0.1;350"\n'));
   writeText(fullfile(cache, 'forcing', ...
      'RetMIP_initial_firn_temperature_KAN-U.tab'), ...
      sprintf('height_m;label\n"0.1;bad"\n'));
   eval_root = fullfile(root, 'eval');
   family_root = fullfile(eval_root, 'retmip');
   stale_file = fullfile(family_root, 'kanu', 'sentinel.bin');
   mkdir(fileparts(stale_file));
   writeText(stale_file, "keep-profile");
   before = treeSnapshot(family_root);

   f = @() icemodel.verification.setup.importRetmip(cache, ...
      case_ids="kanu", evaluation_data_root=eval_root, ...
      input_data_root=fullfile(root, 'input'), overwrite=true, ...
      skip_missing=true);

   testCase.verifyError(f, ...
      'icemodel:verification:readRetmipProfileTable:badProfile');
   testCase.verifyEqual(treeSnapshot(family_root), before);
end

function test_retmip_empty_clamped_surface_preserves_case_root(testCase)
   % Metadata overlap is not enough: a clamp with no protocol rows is skippable,
   % and overwrite preparation must leave the existing case root untouched.
   root = fullfile(testCase.TestData.cache, 'retmip-empty-clamp');
   cache = makeRetmipProtocolCache(fullfile(root, 'source'), "KANU");
   family_root = fullfile(root, 'eval', 'retmip');
   stale_file = fullfile(family_root, 'kanu', 'sentinel.bin');
   mkdir(fileparts(stale_file));
   writeText(stale_file, "keep-empty-window");
   before = treeSnapshot(fullfile(family_root, 'kanu'));

   returned = testCase.verifyWarning( ...
      @() icemodel.verification.setup.importRetmip(cache, ...
      case_ids="kanu", evaluation_data_root=fullfile(root, 'eval'), ...
      input_data_root=fullfile(root, 'input'), ...
      startdate="2012-05-01", enddate="2012-05-01", ...
      overwrite=true, skip_missing=true), ...
      'icemodel:verification:importRetmip:caseSkipped');

   testCase.verifyEmpty(returned.cases);
   testCase.verifyEqual(string(returned.skipped.site), "kanu");
   testCase.verifyTrue(contains(string(returned.skipped.reason), ...
      "no rows in requested window"));
   testCase.verifyEqual(treeSnapshot(fullfile(family_root, 'kanu')), before);
end

function test_retmip_strict_later_parse_failure_precedes_all_mutation(testCase)
   % Strict multi-case staging must validate a later case before writing an
   % earlier valid case or replacing the existing manifest.
   root = fullfile(testCase.TestData.cache, 'retmip-strict-preflight');
   cache = makeRetmipProtocolCache(fullfile(root, 'source'), ...
      ["KANU", "Summit"]);
   writeText(fullfile(cache, 'forcing', ...
      'RetMIP_surface_forcing_Summit.tab'), ...
      sprintf(['time;melt_mmweq;acc_subl_mmweq;Tsurf_K\n' ...
      '"01-Jan-2001 00:00:00;0;0;260"\n' ...
      '"01-Jan-2001 02:00:00;1;2;261"\n']));
   family_root = fullfile(root, 'eval', 'retmip');
   for case_id = ["kanu", "summit"]
      stale_file = fullfile(family_root, case_id, 'sentinel.bin');
      mkdir(fileparts(stale_file));
      writeText(stale_file, "keep-" + case_id);
   end
   writeText(fullfile(family_root, 'manifest.json'), "keep-manifest");
   before = treeSnapshot(family_root);

   f = @() icemodel.verification.setup.importRetmip(cache, ...
      case_ids=["kanu", "summit"], ...
      evaluation_data_root=fullfile(root, 'eval'), ...
      input_data_root=fullfile(root, 'input'), overwrite=true, ...
      overwrite_family=true, skip_missing=false);

   testCase.verifyError(f, ...
      'icemodel:verification:readRetmipProtocolTable:badCadence');
   testCase.verifyEqual(treeSnapshot(family_root), before);
end

function test_retmip_strict_later_empty_window_precedes_all_mutation(testCase)
   % Strict preflight applies the requested clamp to every case before an earlier
   % valid case can mutate its root or the family manifest.
   root = fullfile(testCase.TestData.cache, 'retmip-strict-empty-window');
   cache = makeRetmipProtocolCache(fullfile(root, 'source'), ...
      ["KANU", "Summit"]);
   writeText(fullfile(cache, 'forcing', ...
      'RetMIP_surface_forcing_KANU.tab'), ...
      sprintf(['time;melt_mmweq;acc_subl_mmweq;Tsurf_K\n' ...
      '"01-May-2012 00:00:00;0;0;260"\n' ...
      '"01-May-2012 03:00:00;1;2;261"\n']));
   family_root = fullfile(root, 'eval', 'retmip');
   for case_id = ["kanu", "summit"]
      stale_file = fullfile(family_root, case_id, 'sentinel.bin');
      mkdir(fileparts(stale_file));
      writeText(stale_file, "keep-" + case_id);
   end
   writeText(fullfile(family_root, 'manifest.json'), "keep-manifest");
   before = treeSnapshot(family_root);

   f = @() icemodel.verification.setup.importRetmip(cache, ...
      case_ids=["kanu", "summit"], ...
      evaluation_data_root=fullfile(root, 'eval'), ...
      input_data_root=fullfile(root, 'input'), ...
      startdate="2012-05-01", enddate="2012-05-01", ...
      overwrite=true, overwrite_family=true, skip_missing=false);

   testCase.verifyError(f, ...
      'icemodel:verification:importRetmip:emptyProtocolWindow');
   testCase.verifyEqual(treeSnapshot(family_root), before);
end

function test_retmip_missing_promice_retains_protocol_provenance(testCase)
   % Missing native PROMICE data is a recorded forcing gap, not grounds to drop
   % the RetMIP protocol case or invent staged met/userdata files.
   root = fullfile(testCase.TestData.cache, 'retmip-missing-promice');
   cache = makeRetmipProtocolCache(fullfile(root, 'source'), "KANU");
   missing_promice = fullfile(root, 'empty-promice');
   returned = icemodel.verification.setup.importRetmip(cache, ...
      case_ids="kanu", evaluation_data_root=fullfile(root, 'eval'), ...
      input_data_root=fullfile(root, 'input'), promice_dir=missing_promice, ...
      forcing_sources="retmip", build_forcing=true, ...
      skip_missing=true, overwrite_family=true);

   c = returned.cases;
   testCase.verifyEqual(string(c.case_id), "kanu");
   testCase.verifyEmpty(c.forcing_sources);
   testCase.verifyEqual(string(c.eval_sources), "retmip_protocol");
   testCase.verifyEqual(string(c.colocation.retmip.native_met_status), ...
      "missing_promice");
   retmip_reason = string(c.colocation.retmip.native_met_skipped_reason);
   native_reason = string(c.colocation.native_met.reason);
   testCase.verifyTrue(contains(retmip_reason, "KAN_U") ...
      && contains(retmip_reason, "hour"));
   testCase.verifyEqual(native_reason, retmip_reason);
   testCase.verifyEqual(string(c.colocation.native_met.source_family), ...
      "promice");
   testCase.verifyEqual(string(c.colocation.native_met.source_id), "KAN_U");
   testCase.verifyEqual(string(c.colocation.retmip.native_source.family), ...
      "promice");
   testCase.verifyEqual(string(c.colocation.retmip.native_source.source_id), ...
      "KAN_U");
   testCase.verifyFalse(isfield(c.colocation.retmip, 'met_files'));
   testCase.verifyFalse(isfield(c.colocation.retmip, 'data_files'));
   input_listing = dir(fullfile(root, 'input', '**', '*'));
   testCase.verifyFalse(any(~[input_listing.isdir]));

   % Strict handling raises the established identifier before touching either
   % the existing case root or manifest.
   strict_family_root = fullfile(root, 'strict-eval', 'retmip');
   strict_stale = fullfile(strict_family_root, 'kanu', 'sentinel.bin');
   mkdir(fileparts(strict_stale));
   writeText(strict_stale, "keep-strict-native");
   icemodel.verification.setup.writeManifest( ...
      fullfile(strict_family_root, 'manifest.json'), returned);
   before = treeSnapshot(strict_family_root);
   f = @() icemodel.verification.setup.importRetmip(cache, ...
      case_ids="kanu", evaluation_data_root=fullfile(root, 'strict-eval'), ...
      input_data_root=fullfile(root, 'strict-input'), ...
      promice_dir=missing_promice, forcing_sources="retmip", ...
      build_forcing=true, skip_missing=false, overwrite=true);
   testCase.verifyError(f, ...
      'icemodel:verification:importRetmip:missingNativeSource');
   testCase.verifyEqual(treeSnapshot(strict_family_root), before);
end

function test_samimi_dye2_builder_preserves_half_hourly_native_data(testCase)
   % The source adapter converts units but leaves cadence to artifact writers.
   cache = makeSamimiWorkbookCache(testCase.TestData.cache);

   [met, metadata, Data] = icemodel.forcing.buildSamimiDye2Met( ...
      source_dir=cache, startdate="2016-05-02 00:00:00", ...
      enddate="2016-05-02 01:30:00");

   testCase.verifyEqual(height(Data), 4);
   testCase.verifyEqual(minutes(diff(Data.Time)), repmat(30, 3, 1));
   testCase.verifyEqual(Data.tair(1), 260.15, 'AbsTol', 1e-12);
   testCase.verifyTrue(metadata.fillwithmissing);
   testCase.verifyEqual(Data.psfc(1), 78000, 'AbsTol', 1e-12);
   testCase.verifyEqual(Data.swd(1), 0, 'AbsTol', 1e-12);
   testCase.verifyEqual(Data.swu(1), 0, 'AbsTol', 1e-12);
   testCase.verifyTrue(isnan(Data.albedo(1)));
   testCase.verifyTrue(all(isnan(Data.albedo(2:4))));
   testCase.verifyEqual(Data.snow_depth, [0.95; 0.96; 0.97; 0.98], ...
      'AbsTol', 1e-12);
   testCase.verifyFalse(ismember("tsfc", ...
      string(Data.Properties.VariableNames)));
   testCase.verifyEqual(metadata.shortwave_negative_counts.swd, 1);
   testCase.verifyEqual(metadata.shortwave_negative_counts.swu, 1);
   testCase.verifyEqual(string(metadata.channel_map.snow_depth), "dsnow");
   testCase.verifyEqual(metadata.albedo_qc_counts.low_solar_elevation, 3);
   testCase.verifyEqual(metadata.albedo_qc_counts.total, 3);
   testCase.verifyTrue(all(isnan(Data.rainf)));
   testCase.verifyTrue(all(isnan(Data.snowf)));
   testCase.verifyTrue(all(isnan(met.ppt)));
   testCase.verifyTrue(all(ismember( ...
      icemodel.forcing.helpers.metvariables(), ...
      string(met.Properties.VariableNames))));
   testCase.verifyWarningFree(@() icemodel.forcing.helpers.validatemet(met));
   testCase.verifyEqual(string(metadata.source_family), "samimi");
   testCase.verifyEqual(size(metadata.source_variables, 2), 1);
   testCase.verifyEqual(string(met.Properties.UserData.source_family), ...
      "samimi");
   testCase.verifyTrue(isfield(met.Properties.UserData, 'met_variables'));
   testCase.verifyTrue(all(isfinite([metadata.site_location.x_epsg3413, ...
      metadata.site_location.y_epsg3413])));
   expected_policy = "rainf = NaN placeholder; snowf = NaN placeholder; " + ...
      "ppt = NaN placeholder via data2met because no precipitation channel exists";
   testCase.verifyEqual(string(metadata.precip_policy), expected_policy);
end

function test_samimi_dye2_builder_finds_nested_workbook(testCase)
   % DOI/package cache folders should be accepted just like flat manual caches.
   cache = makeSamimiWorkbookCache(testCase.TestData.cache);
   nested = fullfile(cache, 'doi-package');
   mkdir(nested)
   movefile(fullfile(cache, 'Dye2_AWS_Summer2016.xlsx'), ...
      fullfile(nested, 'Samimi_Dye_AWS_2016_fixture.xlsx'));

   Data = icemodel.forcing.buildSamimiDye2Data(source_dir=cache, ...
      startdate="2016-05-02 00:00:00", enddate="2016-05-02 01:30:00");

   testCase.verifyEqual(height(Data), 4);
   testCase.verifyEqual(Data.tair(1), 260.15, 'AbsTol', 1e-12);
end

function test_samimi_dye2_builder_requires_wind_direction(testCase)
   % The workbook validator should catch an omitted wind direction column
   % before MATLAB raises a generic table-variable error.
   cache = makeSamimiWorkbookWithoutWindDirection(testCase.TestData.cache);

   f = @() icemodel.forcing.buildSamimiDye2Data(source_dir=cache);

   testCase.verifyError(f, ...
      'icemodel:forcing:buildSamimiDye2Data:missingVariables');
end

function test_samimi_dye2_builder_supports_strict_mode(testCase)
   % fillwithmissing=false is exposed and leaves validation strict.
   cache = makeSamimiWorkbookCache(testCase.TestData.cache);

   [met, metadata] = icemodel.forcing.buildSamimiDye2Met( ...
      source_dir=cache, startdate="2016-05-02 00:00:00", ...
      enddate="2016-05-02 01:30:00", fillwithmissing=false);

   testCase.verifyFalse(metadata.fillwithmissing);
   testCase.verifyWarningFree(@() icemodel.forcing.helpers.validatemet(met));
   testCase.verifyFalse(met.Properties.UserData.fillwithmissing);
end

function test_samimi_dye2_builder_preserves_native_wind_direction(testCase)
   % Circular resampling belongs to the writer; the builder keeps source rows.
   cache = makeSamimiWorkbookCache(testCase.TestData.cache);
   filename = fullfile(cache, 'Dye2_AWS_Summer2016.xlsx');
   T = readtable(filename, 'VariableNamingRule', 'preserve');
   T.winddir = [350; 10; 90; 90];
   writetable(T, filename, 'Sheet', 'WS_TDR_Samira_Greenland_2016_Ha');

   Data = icemodel.forcing.buildSamimiDye2Data(source_dir=cache, ...
      startdate="2016-05-02 00:00:00", ...
      enddate="2016-05-02 01:30:00");

   testCase.verifyEqual(Data.wdir, [350; 10; 90; 90]);
end

function test_retmip_import_stages_samimi_dye2_native_source(testCase)
   % RetMIP Dye-2 2016 should stage Samimi native met/userdata using the
   % retmip runtime source label while preserving Samimi provenance.
   retmip_cache = makeRetmipProtocolCache(testCase.TestData.cache, ...
      "Dye-2_16");
   samimi_cache = makeSamimiWorkbookCache(testCase.TestData.cache);
   movefile(fullfile(samimi_cache, 'Dye2_AWS_Summer2016.xlsx'), ...
      fullfile(samimi_cache, 'Samimi_Dye_AWS_2016_fixture.xlsx'));
   eval_root = fullfile(testCase.TestData.cache, 'eval-samimi');
   input_root = fullfile(testCase.TestData.cache, 'input-samimi');
   mkdir(eval_root);
   mkdir(input_root);

   returned = icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids="dye2_16", dry_run=false, samimi_dir=samimi_cache, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      forcing_sources="retmip", build_forcing=true, ...
      overwrite_family=true);

   c = returned.cases;
   testCase.verifyEqual(string(c.case_id), "dye2_2016");
   testCase.verifyEqual(string(c.forcing_sources), "retmip");
   testCase.verifyEqual(string(c.eval_sources), "retmip_protocol");
   testCase.verifyEqual(string(c.colocation.source_association.family), ...
      "samimi");
   testCase.verifyEqual(string(c.colocation.retmip.native_met_status), ...
      "staged");
   testCase.verifyEqual(string(c.colocation.native_met.source_family), ...
      "samimi");
   testCase.verifyFalse(logical(c.colocation.retmip.forcing_ready));
   testCase.verifyTrue(contains( ...
      string(c.colocation.retmip.forcing_ready_reason), "ppt"));
   testCase.verifyEqual(string(c.native_timestep), "00:30:00");
   testCase.verifyTrue(any(string( ...
      c.observation_variables.native_source.data_variables) == "snow_depth"));
   testCase.verifyFalse(any(string(c.comparison_variables) == "snow_depth"));

   met_file = fullfile(input_root, 'met', c.colocation.retmip.met_files(1));
   data_file = fullfile(input_root, 'userdata', ...
      c.colocation.retmip.data_files(1));
   testCase.verifyTrue(isfile(met_file));
   testCase.verifyTrue(isfile(data_file));
   saved_met = load(met_file, 'met');
   saved_data = load(data_file, 'Data');
   testCase.verifyEqual(seconds(median(diff(saved_met.met.Time))), 900);
   testCase.verifyEqual(seconds(median(diff(saved_data.Data.Time))), 3600);
end

function test_retmip_samimi_native_override_writes_and_resolves_30m(testCase)
   % dt_out="" must preserve the proven Samimi cadence through runtime lookup.
   retmip_cache = makeRetmipProtocolCache(testCase.TestData.cache, ...
      "Dye-2_16");
   writeText(fullfile(retmip_cache, 'forcing', ...
      'RetMIP_surface_forcing_Dye-2_16.tab'), ...
      sprintf(['time;melt_mmweq;acc_subl_mmweq;Tsurf_K\n' ...
      '"02-May-2016 00:00:00;0;0;260"\n' ...
      '"02-May-2016 03:00:00;1;2;261"\n']));
   samimi_cache = makeSamimiWorkbookCache(testCase.TestData.cache);
   eval_root = fullfile(testCase.TestData.cache, 'eval-samimi-native');
   input_root = fullfile(testCase.TestData.cache, 'input-samimi-native');

   returned = icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids="dye2_16", samimi_dir=samimi_cache, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      forcing_sources="retmip", build_forcing=true, dt_out="", ...
      startdate="2016-05-02 00:00:00", ...
      enddate="2016-05-02 01:30:00", ...
      overwrite_family=true);

   c = returned.cases;
   met_file = fullfile(input_root, 'met', c.colocation.retmip.met_files(1));
   saved = load(met_file, 'met');
   testCase.verifyTrue(endsWith(string(met_file), "_30m.mat"));
   testCase.verifyEqual(seconds(median(diff(saved.met.Time))), 1800);

   % Exercise the same manifest-to-configureRun path used by verification runs.
   c.dataset_family = returned.dataset_family;
   c.input_data_root = string(input_root);
   opts = icemodel.test.helpers.setModelOptsForCase(c);
   testCase.verifyEqual(opts.dt, 1800);
   testCase.verifyTrue(isfile(opts.metfname{1}));
   testCase.verifyEqual(string(opts.metfname{1}), string(met_file));
   loaded = icemodel.loadmet(opts);
   testCase.verifyEqual(seconds(median(diff(loaded.Time))), 1800);
end

function test_retmip_protocol_refresh_preserves_native_runtime_leg(testCase)
   % An ordinary protocol-only refresh must not erase an earlier native build.
   retmip_cache = makeRetmipProtocolCache(testCase.TestData.cache, ...
      "Dye-2_16");
   alternate_root = fullfile(testCase.TestData.cache, ...
      'alternate-retmip-protocol');
   alternate_retmip_cache = makeRetmipProtocolCache( ...
      alternate_root, "Dye-2_16");
   samimi_cache = makeSamimiWorkbookCache(testCase.TestData.cache);
   eval_root = fullfile(testCase.TestData.cache, 'eval-retmip-additive');
   input_root = fullfile(testCase.TestData.cache, 'input-retmip-additive');

   first = icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids="dye2_16", samimi_dir=samimi_cache, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      forcing_sources="retmip", build_forcing=true, ...
      overwrite_family=true);
   prior = first.cases;
   prior_met = string(prior.colocation.retmip.met_files);
   prior_data = string(prior.colocation.retmip.data_files);
   met_file = fullfile(input_root, 'met', prior_met(1));
   data_file = fullfile(input_root, 'userdata', prior_data(1));
   met_bytes = fileBytes(met_file);
   data_bytes = fileBytes(data_file);

   % The second call refreshes the protocol observation contract only and does
   % not receive the native Samimi cache, so preservation cannot depend on a
   % hidden native-source probe or rewrite.
   second = icemodel.verification.setup.importRetmip(alternate_retmip_cache, ...
      case_ids="dye2_16", evaluation_data_root=eval_root, ...
      input_data_root=input_root, forcing_sources=strings(1, 0), ...
      build_forcing=false, overwrite=true);
   current = second.cases;

   testCase.verifyEqual(string(current.colocation.retmip.met_files), prior_met);
   testCase.verifyEqual(string(current.colocation.retmip.data_files), prior_data);
   testCase.verifyEqual(string(current.forcing_sources), ...
      string(prior.forcing_sources));
   testCase.verifyEqual(jsonencode(current.colocation.retmip.native_source), ...
      jsonencode(prior.colocation.retmip.native_source));
   testCase.verifyEqual(jsonencode(current.colocation.retmip.window), ...
      jsonencode(prior.colocation.retmip.window));
   testCase.verifyEqual(string(current.colocation.retmip.source_dir), ...
      string(alternate_retmip_cache));
   testCase.verifyTrue(contains(string(current.colocation.retmip.surface_file), ...
      string(alternate_retmip_cache)));
   testCase.verifyEqual(jsonencode(current.colocation.native_met), ...
      jsonencode(prior.colocation.native_met));
   testCase.verifyEqual(string(current.native_timestep), ...
      string(prior.native_timestep));
   testCase.verifyEqual(fileBytes(met_file), met_bytes);
   testCase.verifyEqual(fileBytes(data_file), data_bytes);
end

function test_retmip_protocol_refresh_rejects_changed_native_identity(testCase)
   % A protocol-only refresh cannot attach files from another native producer.
   retmip_cache = makeRetmipProtocolCache(testCase.TestData.cache, ...
      "Dye-2_16");
   samimi_cache = makeSamimiWorkbookCache(testCase.TestData.cache);
   eval_root = fullfile(testCase.TestData.cache, 'eval-retmip-identity');
   input_root = fullfile(testCase.TestData.cache, 'input-retmip-identity');

   first = icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids="dye2_16", samimi_dir=samimi_cache, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      forcing_sources="retmip", build_forcing=true, ...
      overwrite_family=true);
   prior = first.cases;
   manifest_file = fullfile(eval_root, 'retmip', 'manifest.json');
   prior_manifest = ...
      icemodel.verification.helpers.readFamilyManifest(manifest_file);
   fresh_identity = prior.colocation.retmip.native_source;
   met_file = fullfile(input_root, 'met', ...
      prior.colocation.retmip.met_files(1));
   data_file = fullfile(input_root, 'userdata', ...
      prior.colocation.retmip.data_files(1));
   met_bytes = fileBytes(met_file);
   data_bytes = fileBytes(data_file);

   % Exercise every stable native identity field independently. Missing legacy
   % facts remain compatible, but a known conflict must require a fresh build.
   fields = ["family", "source_id", "relationship"];
   replacements = ["promice", "KAN_U", "evaluation_source"];
   for k = 1:numel(fields)
      mutated = prior_manifest;
      name = char(fields(k));
      mutated.cases.colocation.retmip.native_source.(name) = ...
         char(replacements(k));
      mutated.cases.colocation.retmip.source_dir = 'stale-protocol-cache';
      mutated.cases.colocation.retmip.surface_file = 'stale-surface.tab';
      writeText(manifest_file, jsonencode(mutated));

      refreshed = icemodel.verification.setup.importRetmip(retmip_cache, ...
         case_ids="dye2_16", evaluation_data_root=eval_root, ...
         input_data_root=input_root, forcing_sources=strings(1, 0), ...
         build_forcing=false, overwrite=true);
      current = refreshed.cases;

      testCase.verifyFalse(isfield(current.colocation.retmip, 'met_files'));
      testCase.verifyFalse(isfield(current.colocation.retmip, 'data_files'));
      testCase.verifyFalse(isfield( ...
         current.observation_variables, 'native_source'));
      testCase.verifyEmpty(current.forcing_sources);
      testCase.verifyEqual(string(current.native_timestep), "3hr");
      testCase.verifyEqual( ...
         string(current.colocation.retmip.native_met_status), ...
         "identity_changed_requires_rebuild");
      testCase.verifyTrue(contains(string( ...
         current.colocation.retmip.native_met_skipped_reason), ...
         "build_forcing=true"));
      testCase.verifyFalse(logical(current.colocation.native_met.staged));
      testCase.verifyEqual(string(current.colocation.retmip.source_dir), ...
         string(retmip_cache));
      testCase.verifyTrue(contains( ...
         string(current.colocation.retmip.surface_file), ...
         string(retmip_cache)));
      testCase.verifyFalse(isfield( ...
         current.colocation.retmip.native_source, 'source_file'));
      for identity_field = fields
         identity_name = char(identity_field);
         testCase.verifyEqual(string( ...
            current.colocation.retmip.native_source.(identity_name)), ...
            string(fresh_identity.(identity_name)));
      end
      testCase.verifyEqual(fileBytes(met_file), met_bytes);
      testCase.verifyEqual(fileBytes(data_file), data_bytes);
   end
end

function test_retmip_import_rethrows_malformed_native_source(testCase)
   % Malformed native files are regressions, not missing-source skips.
   retmip_cache = makeRetmipProtocolCache(testCase.TestData.cache, ...
      "Dye-2_16");
   samimi_cache = makeMalformedSamimiWorkbookCache(testCase.TestData.cache);
   eval_root = fullfile(testCase.TestData.cache, 'eval-samimi-malformed');
   input_root = fullfile(testCase.TestData.cache, 'input-samimi-malformed');
   mkdir(eval_root);
   mkdir(input_root);

   f = @() icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids="dye2_16", dry_run=false, samimi_dir=samimi_cache, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      forcing_sources="retmip", build_forcing=true, ...
      skip_missing=true, overwrite_family=true);

   testCase.verifyError(f, ...
      'icemodel:forcing:buildSamimiDye2Data:missingVariables');
end

function test_retmip_import_stages_promice_kanu_native_source(testCase)
   % RetMIP KAN_U should stage PROMICE-derived native met/userdata with
   % placeholder precipitation under the retmip runtime source label.
   promice_dir = requirePromiceVerificationCache(testCase, "KAN_U");
   retmip_cache = makeRetmipProtocolCache(testCase.TestData.cache, "KANU");
   eval_root = fullfile(testCase.TestData.cache, 'eval-promice-kanu');
   input_root = fullfile(testCase.TestData.cache, 'input-promice-kanu');
   mkdir(eval_root);
   mkdir(input_root);

   returned = icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids="kanu", dry_run=false, promice_dir=promice_dir, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      forcing_sources="retmip", build_forcing=true, ...
      overwrite_family=true);

   c = returned.cases;
   testCase.verifyEqual(string(c.case_id), "kanu");
   testCase.verifyEqual(string(c.forcing_sources), "retmip");
   testCase.verifyEqual(string(c.eval_sources), "retmip_protocol");
   testCase.verifyEqual(string(c.colocation.source_association.family), ...
      "promice");
   testCase.verifyEqual(string(c.colocation.retmip.native_met_status), ...
      "staged");
   testCase.verifyEqual(string(c.colocation.native_met.source_family), ...
      "promice");
   testCase.verifyFalse(logical(c.colocation.retmip.forcing_ready));
   testCase.verifyTrue(contains( ...
      string(c.colocation.retmip.forcing_ready_reason), "ppt"));

   met_file = fullfile(input_root, 'met', c.colocation.retmip.met_files(1));
   data_file = fullfile(input_root, 'userdata', ...
      c.colocation.retmip.data_files(1));
   testCase.verifyTrue(isfile(met_file));
   testCase.verifyTrue(isfile(data_file));
   loaded = load(met_file);
   testCase.verifyTrue(all(ismember(icemodel.forcing.helpers.metvariables(), ...
      string(loaded.met.Properties.VariableNames))));
   testCase.verifyWarningFree( ...
      @() icemodel.forcing.helpers.validatemet(loaded.met));
   testCase.verifyTrue(contains( ...
      string(c.colocation.retmip.native_source.precip_policy), ...
      "ppt = NaN placeholder"));
end

function test_imau_import_dry_run_manifest(testCase)
   % importImau dry-run should expose S21/S22/S23 and the S21->RetMIP FA link.
   missing_imau = fullfile(testCase.TestData.cache, 'missing-imau-dry');
   returned = icemodel.verification.setup.importImau( ...
      missing_imau, dry_run=true);
   testCase.verifyEqual(string(returned.source_version), ...
      "hourly-s21-s22-s23+daily-qa");

   ids = string({returned.cases.site_id});
   expected = ["S21", "S22", "S23"];
   testCase.verifyEqual(ids, expected);
   s21 = returned.cases(ids == "S21");
   expected = "retmip";
   testCase.verifyEqual(string(s21.colocation.source_association.family), ...
      expected);
   expected = "fa";
   testCase.verifyEqual(string(s21.colocation.source_association.source_id), ...
      expected);
   expected = "1hr";
   testCase.verifyEqual(string(s21.native_timestep), expected);
   testCase.verifyEmpty(s21.forcing_sources);
   testCase.verifyEqual(string(s21.eval_sources), "imau_obs");
   testCase.verifyEmpty(fieldnames(s21.colocation.imau.cache_status));
   testCase.verifyEqual(unique(string( ...
      {s21.colocation.daily_qa.records.reason})), ...
      "daily QA cache not read for this call");
   testCase.verifyFalse(isfolder(missing_imau));
end

function test_imau_import_without_forcing_writes_no_input_artifacts(testCase)
   % Observation staging may read native data but build_forcing=false owns no input.
   cache = makeImauSourceCache(testCase.TestData.cache, "S21");
   eval_root = fullfile(testCase.TestData.cache, 'eval-imau-observations');
   input_root = fullfile(testCase.TestData.cache, 'input-imau-disabled');

   returned = icemodel.verification.setup.importImau(cache, ...
      case_ids="S21", evaluation_data_root=eval_root, ...
      input_data_root=input_root, build_forcing=false, ...
      overwrite_family=true);

   testCase.verifyTrue(isfile(fullfile(eval_root, 'imau', ...
      returned.cases.evaluation_file)));
   testCase.verifyEmpty(returned.cases.forcing_sources);
   testCase.verifyEmpty(returned.cases.colocation.imau.met_files);
   testCase.verifyEmpty(returned.cases.colocation.imau.data_files);
   testCase.verifyEqual( ...
      string(returned.cases.colocation.imau.forcing_ready_reason), ...
      "native forcing disabled because build_forcing=false");
   testCase.verifyFalse(isfolder(input_root));
end

function test_imau_observation_refresh_preserves_native_artifacts(testCase)
   % Fresh observation provenance wins while compatible native artifacts remain.
   cache = makeImauSourceCache(testCase.TestData.cache, "S21");
   alternate_root = fullfile(testCase.TestData.cache, 'alternate-imau-source');
   mkdir(alternate_root);
   alternate_cache = makeImauSourceCache(alternate_root, "S21", ...
      daily_sites=strings(1, 0));
   appendImauHourlyFixtureRow(alternate_cache, "S21", ...
      datetime(2014, 4, 12, 2, 0, 0, 'TimeZone', 'UTC'));
   eval_root = fullfile(testCase.TestData.cache, 'eval-imau-refresh');
   input_root = fullfile(testCase.TestData.cache, 'input-imau-refresh');

   icemodel.verification.setup.importImau(cache, ...
      case_ids="S21", evaluation_data_root=eval_root, ...
      input_data_root=input_root, forcing_sources="imau", ...
      build_forcing=true, overwrite_family=true);

   % Persist an artifact-identity sentinel to prove the selective merge retains
   % forcing provenance without retaining stale observation-source fields.
   manifest_file = fullfile(eval_root, 'imau', 'manifest.json');
   prior_manifest = ...
      icemodel.verification.helpers.readFamilyManifest(manifest_file);
   prior_manifest.cases.colocation.imau.artifact_metadata = ...
      struct('source_id', 'S21', 'artifact_cadence_seconds', 900);
   prior_manifest.cases.colocation.imau.evaluation_file = ...
      'stale/observations.mat';
   writeText(manifest_file, jsonencode(prior_manifest));
   prior = prior_manifest.cases;
   prior_met = string(prior.colocation.imau.met_files);
   prior_data = string(prior.colocation.imau.data_files);
   prior_sources = string(prior.forcing_sources);
   met_file = fullfile(input_root, 'met', prior_met(1));
   data_file = fullfile(input_root, 'userdata', prior_data(1));
   met_bytes = fileBytes(met_file);
   data_bytes = fileBytes(data_file);

   refreshed = icemodel.verification.setup.importImau(alternate_cache, ...
      case_ids="S21", evaluation_data_root=eval_root, ...
      input_data_root=input_root, forcing_sources=strings(1, 0), ...
      startdate="2014-04-12 00:00:00", ...
      enddate="2014-04-12 02:00:00", ...
      build_forcing=false, overwrite=true);
   current = refreshed.cases;

   testCase.verifyEqual(string(current.colocation.imau.met_files), prior_met);
   testCase.verifyEqual(string(current.colocation.imau.data_files), prior_data);
   testCase.verifyEqual(string(current.forcing_sources), prior_sources);
   testCase.verifyEqual(current.colocation.imau.forcing_ready, ...
      prior.colocation.imau.forcing_ready);
   testCase.verifyEqual(string(current.colocation.imau.forcing_ready_reason), ...
      string(prior.colocation.imau.forcing_ready_reason));
   testCase.verifyEqual(current.colocation.imau.forcing_complete_windows, ...
      prior.colocation.imau.forcing_complete_windows);
   testCase.verifyEqual(current.colocation.imau.artifact_metadata, ...
      prior.colocation.imau.artifact_metadata);
   testCase.verifyEqual(string(current.colocation.imau.source_dir), ...
      string(alternate_cache));
   testCase.verifyNotEqual(string(current.colocation.imau.source_dir), ...
      string(prior.colocation.imau.source_dir));
   testCase.verifyTrue(contains(string(current.colocation.imau.source_file), ...
      string(alternate_cache)));
   testCase.verifyEqual(string(current.colocation.imau.evaluation_file), ...
      string(current.evaluation_file));
   testCase.verifyNotEqual(string(current.colocation.imau.evaluation_file), ...
      string(prior.colocation.imau.evaluation_file));
   testCase.verifyEqual(string(current.period.end), ...
      "2014-04-12 02:00:00");
   testCase.verifyEqual(string(current.colocation.imau.window.start), ...
      "2014-04-12 00:00:00");
   testCase.verifyEqual(string(current.colocation.imau.window.end), ...
      "2014-04-12 02:00:00");
   testCase.verifyFalse(logical(current.colocation.daily_qa.validated));
   testCase.verifyFalse(logical( ...
      current.observation_variables.daily_qa_validated));
   testCase.verifyEqual(fileBytes(met_file), met_bytes);
   testCase.verifyEqual(fileBytes(data_file), data_bytes);
end

function test_imau_observation_refresh_rejects_native_identity_conflicts(testCase)
   % Conflicting producer metadata unlinks prior native refs without deleting files.
   cache = makeImauSourceCache(testCase.TestData.cache, "S21");
   eval_root = fullfile(testCase.TestData.cache, 'eval-imau-identity');
   input_root = fullfile(testCase.TestData.cache, 'input-imau-identity');
   first = icemodel.verification.setup.importImau(cache, ...
      case_ids="S21", evaluation_data_root=eval_root, ...
      input_data_root=input_root, forcing_sources="imau", ...
      build_forcing=true, overwrite_family=true);
   prior = first.cases;
   manifest_file = fullfile(eval_root, 'imau', 'manifest.json');
   prior_manifest = ...
      icemodel.verification.helpers.readFamilyManifest(manifest_file);
   met_file = fullfile(input_root, 'met', ...
      prior.colocation.imau.met_files(1));
   data_file = fullfile(input_root, 'userdata', ...
      prior.colocation.imau.data_files(1));
   met_bytes = fileBytes(met_file);
   data_bytes = fileBytes(data_file);

   % One variant conflicts on producer station; the other keeps S21 but moves
   % the finite artifact point. Matching case/DOIs cannot authorize either ref.
   for variant = ["source", "point"]
      mutated = prior_manifest;
      metadata = struct('source_id', 'S21', ...
         'artifact_cadence_seconds', 900);
      if variant == "source"
         metadata.source_id = 'S22';
      else
         metadata.lat_wgs84 = prior.site_location.lat_wgs84 + 1;
         metadata.lon_wgs84 = prior.site_location.lon_wgs84;
      end
      mutated.cases.colocation.imau.artifact_metadata = metadata;
      writeText(manifest_file, jsonencode(mutated));

      refreshed = icemodel.verification.setup.importImau(cache, ...
         case_ids="S21", evaluation_data_root=eval_root, ...
         input_data_root=input_root, forcing_sources=strings(1, 0), ...
         build_forcing=false, overwrite=true);
      current = refreshed.cases;

      testCase.verifyEqual(string(current.colocation.imau.kind), ...
         "hourly_aws_eval");
      testCase.verifyEmpty(current.colocation.imau.met_files);
      testCase.verifyEmpty(current.colocation.imau.data_files);
      testCase.verifyFalse(isfield( ...
         current.colocation.imau, 'artifact_metadata'));
      testCase.verifyFalse(logical(current.colocation.imau.forcing_ready));
      testCase.verifyTrue(contains(string( ...
         current.colocation.imau.forcing_ready_reason), ...
         "build_forcing=true"));
      testCase.verifyEmpty(current.forcing_sources);
      testCase.verifyEqual(string(current.colocation.imau.source_dir), ...
         string(cache));
      testCase.verifyEqual(fileBytes(met_file), met_bytes);
      testCase.verifyEqual(fileBytes(data_file), data_bytes);
   end
end

function test_imau_hourly_builder_maps_native_met(testCase)
   % The IMAU hourly adapter should map corrected AWS channels, attach location
   % metadata, and leave unsupported precipitation channels missing.
   cache = makeImauSourceCache(testCase.TestData.cache, "S21", ...
      time0=datetime(2014, 4, 12, 14, 0, 0, 'TimeZone', 'UTC'));

   [met, metadata, Data] = icemodel.forcing.buildImauHourlyMet( ...
      "S21", source_dir=cache);

   testCase.verifyEqual(height(Data), 2);
   testCase.verifyEqual(hours(diff(Data.Time)), 1);
   testCase.verifyEqual(Data.tair(1), 273.15);
   testCase.verifyEqual(Data.psfc(1), 90000);
   testCase.verifyTrue(all(isnan(Data.rainf)));
   testCase.verifyTrue(all(isnan(Data.snowf)));
   testCase.verifyTrue(all(isnan(met.ppt)));
   testCase.verifyEqual( ...
      Data.Properties.CustomProperties.Lat, 66.181304);
   testCase.verifyTrue(all(ismember( ...
      icemodel.forcing.helpers.metvariables(), ...
      string(met.Properties.VariableNames))));
   testCase.verifyWarningFree(@() icemodel.forcing.helpers.validatemet(met));
   testCase.verifyEqual(string(metadata.source_family), "imau");
   testCase.verifyEqual(size(metadata.source_variables, 2), 1);
   testCase.verifyEqual(string(met.Properties.UserData.station), "S21");
   testCase.verifyTrue(isfield(met.Properties.UserData, 'met_variables'));
   testCase.verifyTrue(all(isfinite([metadata.site_location.x_epsg3413, ...
      metadata.site_location.y_epsg3413])));
   testCase.verifyTrue(contains(metadata.precip_policy, "rainf = NaN"));
   testCase.verifyTrue(contains(metadata.precip_policy, "placeholder"));
   testCase.verifyEqual(Data.albedo, [0.7; NaN]);
   testCase.verifyEqual(Data.swd, [100; 110]);
   testCase.verifyEqual(Data.swu, [50; 55]);
   testCase.verifyEqual(Data.swn(1), 50);
   testCase.verifyEqual(Data.swn(2), 55);
   testCase.verifyEqual(Data.netr(2), 5);
   testCase.verifyTrue(contains(metadata.albedo_policy, ...
      "0.2 lower-bound code"));
   testCase.verifyEqual(metadata.albedo_qc_counts.low_light, 0);
   testCase.verifyEqual(metadata.albedo_qc_counts.low_solar_elevation, 0);
   testCase.verifyEqual(metadata.albedo_qc_counts.source_floor, 1);
   testCase.verifyEqual(metadata.albedo_qc_counts.total, 1);
   testCase.verifyEqual(metadata.albedo_qc_counts.derived_balance, 0);
   testCase.verifyEqual( ...
      metadata.albedo_transient_qc.episode_day_count, 0);
   testCase.verifyEmpty(metadata.albedo_transient_qc.episode_dates);
   testCase.verifyTrue(contains(metadata.shortwave_balance_policy, ...
      "dailyAlbedoAnomalyFlags"));
end

function test_imau_hourly_builder_masks_collapsed_floor_coded_balance(testCase)
   % A bright source-floor sample is excluded from derived balances only when
   % its raw reflected-shortwave ratio independently confirms the collapse.
   cache = makeImauSourceCache(testCase.TestData.cache, "S21");
   filename = icemodel.forcing.helpers.locateImauHourlyFile(cache, "S21");
   lines = readlines(filename);
   header = find(startsWith(lines, "Date/Time"), 1, 'first');
   fields = split(lines(header + 2), sprintf('\t')).';
   fields(15) = "11";
   lines(header + 2) = strjoin(fields, sprintf('\t'));
   writelines(lines, filename);

   [~, metadata, Data] = icemodel.forcing.buildImauHourlyMet( ...
      "S21", source_dir=cache);

   testCase.verifyEqual(Data.swu(2), 11);
   testCase.verifyTrue(isnan(Data.albedo(2)));
   testCase.verifyTrue(isnan(Data.swn(2)));
   testCase.verifyTrue(isnan(Data.netr(2)));
   testCase.verifyTrue(contains(metadata.shortwave_balance_policy, ...
      "raw swu/swd <= 0.2"));
   testCase.verifyEqual(metadata.albedo_qc_counts.derived_balance, 1);
end

function test_imau_hourly_builder_ignores_surface_temperature_fill(testCase)
   % A source fill must not make metchecks treat the remaining Kelvin samples as
   % Celsius and clamp every valid surface temperature to zero.
   cache = makeImauSourceCache(testCase.TestData.cache, "S22");
   filename = icemodel.forcing.helpers.locateImauHourlyFile(cache, "S22");
   lines = readlines(filename);
   header = find(startsWith(lines, "Date/Time"), 1, 'first');
   fields = split(lines(header + 2), sprintf('\t')).';
   fields(26) = "-1273.05";
   lines(header + 2) = strjoin(fields, sprintf('\t'));
   writelines(lines, filename);

   [met, ~, Data] = icemodel.forcing.buildImauHourlyMet( ...
      "S22", source_dir=cache);

   testCase.verifyEqual(Data.tsfc(1), 268.15);
   testCase.verifyTrue(isnan(Data.tsfc(2)));
   testCase.verifyEqual(met.tsfc, Data.tsfc);
end

function test_imau_hourly_builder_enforces_paired_utc_window(testCase)
   % The direct builder rejects every malformed pair and preserves zoned instants.
   cache = makeImauSourceCache(testCase.TestData.cache, "S21");
   error_id = 'icemodel:internal:pairedWindow:invalidWindow';
   testCase.verifyError(@() icemodel.forcing.buildImauHourlyData( ...
      "S21", source_dir=cache, startdate="2014-04-12"), error_id);
   testCase.verifyError(@() icemodel.forcing.buildImauHourlyData( ...
      "S21", source_dir=cache, enddate="2014-04-12"), error_id);
   testCase.verifyError(@() icemodel.forcing.buildImauHourlyData( ...
      "S21", source_dir=cache, startdate="2014-04-13", ...
      enddate="2014-04-12"), error_id);

   full_data = icemodel.forcing.buildImauHourlyData("S21", source_dir=cache);
   expected = full_data.Time(1:min(2, height(full_data)));
   window_start = expected(1);
   window_end = expected(end);
   window_start.TimeZone = 'America/New_York';
   window_end.TimeZone = 'America/New_York';
   Data = icemodel.forcing.buildImauHourlyData("S21", source_dir=cache, ...
      startdate=window_start, enddate=window_end);

   testCase.verifyEqual(Data.Time.TimeZone, 'UTC');
   testCase.verifyEqual(Data.Time, expected);
end

function test_imau_hourly_builder_regularizes_source_gaps(testCase)
   % Real IMAU hourly files can contain outage gaps. The builder should insert
   % missing rows on the hourly axis instead of interpolating or failing the met
   % contract.
   cache = makeImauSourceCache(testCase.TestData.cache, "S21");
   filename = icemodel.forcing.helpers.locateImauHourlyFile(cache, "S21");
   lines = readlines(filename);
   header = find(startsWith(lines, "Date/Time"), 1, 'first');
   fields = split(lines(header + 2), sprintf('\t')).';
   fields(1) = "2014-04-12T02:00:00";
   lines(header + 2) = strjoin(fields, sprintf('\t'));
   writelines(lines, filename);

   [met, ~, Data] = icemodel.forcing.buildImauHourlyMet( ...
      "S21", source_dir=cache);

   testCase.verifyEqual(height(Data), 3);
   testCase.verifyEqual(hours(diff(Data.Time)), ones(2, 1));
   testCase.verifyTrue(isnan(Data.tair(2)));
   testCase.verifyWarningFree(@() icemodel.forcing.helpers.validatemet(met));
end

function test_imau_import_stages_native_hourly_sources(testCase)
   % Non-dry IMAU import should stage only the three hourly first-pass cases,
   % with daily 19-site content recorded as QA provenance.
   cache = makeImauSourceCache(testCase.TestData.cache, ...
      ["S21", "S22", "S23"]);
   eval_root = fullfile(testCase.TestData.cache, 'eval-imau');
   input_root = fullfile(testCase.TestData.cache, 'input-imau');
   mkdir(eval_root);
   mkdir(input_root);

   returned = icemodel.verification.setup.importImau(cache, ...
      case_ids=["S21", "S22", "S23"], ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      forcing_sources="imau", build_forcing=true, ...
      overwrite_family=true);

   ids = string({returned.cases.case_id});
   testCase.verifyEqual(ids, ["s21", "s22", "s23"]);
   s21 = returned.cases(ids == "s21");
   testCase.verifyEqual(string(s21.forcing_sources), "imau");
   testCase.verifyEqual(string(s21.eval_sources), "imau_obs");
   testCase.verifyFalse(logical(s21.colocation.imau.forcing_ready));
   testCase.verifyTrue(contains( ...
      string(s21.colocation.imau.forcing_ready_reason), "ppt"));
   testCase.verifyTrue(isfield(s21.colocation.imau, ...
      'forcing_complete_windows'));
   testCase.verifyEmpty(s21.colocation.imau.forcing_complete_windows);
   testCase.verifyEqual( ...
      string(s21.colocation.source_association.family), "retmip");
   testCase.verifyEqual( ...
      string(s21.colocation.source_association.source_id), "fa");
   testCase.verifyTrue(logical(s21.colocation.daily_qa.validated));
   testCase.verifyTrue(logical(s21.observation_variables.daily_qa_validated));
   testCase.verifyTrue(any(string(s21.comparison_variables) == ...
      "surface_height"));

   met_file = fullfile(input_root, 'met', s21.colocation.imau.met_files(1));
   data_file = fullfile(input_root, 'userdata', ...
      s21.colocation.imau.data_files(1));
   testCase.verifyTrue(isfile(fullfile(eval_root, 'imau', ...
      s21.evaluation_file)));
   testCase.verifyTrue(isfile(met_file));
   testCase.verifyTrue(isfile(data_file));
   loaded = load(met_file);
   testCase.verifyWarningFree( ...
      @() icemodel.forcing.helpers.validatemet(loaded.met));
   testCase.verifyTrue(logical( ...
      loaded.met.Properties.UserData.fillwithmissing));
   testCase.verifyTrue(endsWith(string(met_file), "_15m.mat"));
   testCase.verifyEqual(seconds(median(diff(loaded.met.Time))), 900);
   native = load(data_file, 'Data');
   testCase.verifyEqual(seconds(median(diff(native.Data.Time))), 3600);

   % A subset repeat preserves its own native artifacts and an unrelated site.
   observation_file = fullfile(eval_root, 'imau', s21.evaluation_file);
   observation = load(observation_file, 'targets');
   testCase.verifyEqual(string(observation.targets.metadata.source_family), ...
      "imau");
   testCase.verifyEqual(string(observation.targets.metadata.station), "S21");
   testCase.verifyEqual(string(observation.targets.metadata.doi), ...
      "10.1594/PANGAEA.969585");
   testCase.verifyEqual(string(observation.targets.metadata.bundle_doi), ...
      "10.1594/PANGAEA.971647");
   s22 = returned.cases(ids == "s22");
   unrelated_file = fullfile(eval_root, 'imau', s22.evaluation_file);
   artifacts = [string(observation_file); string(met_file); ...
      string(data_file); string(unrelated_file)];
   original = cellfun(@fileBytes, cellstr(artifacts), 'UniformOutput', false);
   icemodel.verification.setup.importImau(cache, case_ids="S21", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      forcing_sources="imau", build_forcing=true, overwrite=false);
   repeated = cellfun(@fileBytes, cellstr(artifacts), 'UniformOutput', false);
   testCase.verifyEqual(repeated, original);

   % A same-window producer mismatch rewrites only the fixed observation bundle;
   % native input files and unrelated cases remain governed by their own gates.
   observation.targets.metadata.source_family = 'samimi';
   targets = observation.targets;
   save(observation_file, 'targets');
   tampered = fileBytes(observation_file);
   icemodel.verification.setup.importImau(cache, case_ids="S21", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      forcing_sources="imau", build_forcing=true, overwrite=false);
   testCase.verifyNotEqual(fileBytes(observation_file), tampered);
   repaired = load(observation_file, 'targets');
   testCase.verifyEqual(string(repaired.targets.metadata.source_family), ...
      "imau");

   % The same public importer forwards an explicit native model-met cadence.
   native_eval = fullfile(testCase.TestData.cache, 'eval-imau-native');
   native_input = fullfile(testCase.TestData.cache, 'input-imau-native');
   native_manifest = icemodel.verification.setup.importImau(cache, ...
      case_ids="S21", evaluation_data_root=native_eval, ...
      input_data_root=native_input, forcing_sources="imau", ...
      build_forcing=true, dt_out="", overwrite_family=true);
   native_case = native_manifest.cases;
   native_met_file = fullfile(native_input, 'met', ...
      native_case.colocation.imau.met_files(1));
   native_met = load(native_met_file, 'met');
   testCase.verifyTrue(endsWith(string(native_met_file), "_1hr.mat"));
   testCase.verifyEqual(seconds(median(diff(native_met.met.Time))), 3600);
end

function test_imau_import_reports_missing_daily_qa(testCase)
   % Missing daily QA files should be manifest provenance, not extra cases or a
   % reason to lose the hourly native source when skip_missing=true.
   cache = makeImauSourceCache(testCase.TestData.cache, "S21", ...
      daily_sites=strings(1, 0));
   eval_root = fullfile(testCase.TestData.cache, 'eval-imau-missing-daily');
   input_root = fullfile(testCase.TestData.cache, 'input-imau-missing-daily');
   mkdir(eval_root);
   mkdir(input_root);

   returned = icemodel.verification.setup.importImau(cache, ...
      case_ids="S21", dry_run=false, evaluation_data_root=eval_root, ...
      input_data_root=input_root, overwrite_family=true);

   testCase.verifyEqual(numel(returned.cases), 1);
   c = returned.cases;
   testCase.verifyEqual(string(c.case_id), "s21");
   testCase.verifyTrue(logical(c.colocation.imau.staged));
   testCase.verifyFalse(logical(c.colocation.daily_qa.validated));
   testCase.verifyEqual( ...
      string(c.colocation.daily_qa.missing_sites), "S21");
end

function test_imau_strict_native_only_refresh_does_not_require_daily_qa(testCase)
   % Native-only refreshes require hourly input but preserve the prior daily QA
   % record instead of reopening the observation-only daily cache.
   [mar_dir, ~, ~] = requireRcmFixtures(testCase);
   cache = makeImauSourceCache(testCase.TestData.cache, "S21", ...
      daily_sites=strings(1, 0));
   root = fullfile(testCase.TestData.cache, ...
      'imau-strict-native-only-refresh-root');
   initial = icemodel.verification.setup.importImau( ...
      cache, case_ids="S21", output_root=root, ...
      forcing_sources="mar", build_forcing=true, mar_dir=mar_dir, ...
      overwrite_family=true);

   refreshed = icemodel.verification.setup.importImau( ...
      cache, case_ids="S21", output_root=root, forcing_sources="imau", ...
      build_observations=false, build_forcing=true, skip_missing=false, ...
      overwrite=true, overwrite_family=true);

   c = refreshed.cases;
   testCase.verifyEqual(jsonencode(c.colocation.daily_qa), ...
      jsonencode(initial.cases.colocation.daily_qa));
   testCase.verifyEqual(string(c.forcing_sources), "imau");
   testCase.verifyNotEmpty(c.colocation.imau.met_files);
   testCase.verifyNotEmpty(c.colocation.imau.data_files);
   testCase.verifyFalse(isfield(c.colocation, 'mar'));
   testCase.verifyFalse(any(string(c.forcing_sources) == "mar3.11"));
end

function test_imau_import_strict_non_dry_fails_missing_daily_qa(testCase)
   % Strict non-dry IMAU staging should fail before writing partial site data.
   cache = makeImauSourceCache(testCase.TestData.cache, "S21", ...
      daily_sites="S22");
   eval_root = fullfile(testCase.TestData.cache, 'eval-imau-strict-daily');
   input_root = fullfile(testCase.TestData.cache, 'input-imau-strict-daily');
   mkdir(eval_root);
   mkdir(input_root);

   f = @() icemodel.verification.setup.importImau(cache, ...
      case_ids="S21", dry_run=false, skip_missing=false, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      overwrite_family=true);

   testCase.verifyError(f, ...
      'icemodel:verification:importImau:missingDailyQa');
end

function test_imau_import_strict_preflights_hourly_files(testCase)
   % Strict non-dry IMAU staging should fail before writing partial cases.
   cache = makeImauSourceCache(testCase.TestData.cache, "S21", ...
      daily_sites=["S21", "S22"]);
   eval_root = fullfile(testCase.TestData.cache, 'eval-imau-strict-hourly');
   input_root = fullfile(testCase.TestData.cache, 'input-imau-strict-hourly');
   mkdir(eval_root);
   mkdir(input_root);

   f = @() icemodel.verification.setup.importImau(cache, ...
      case_ids=["S21", "S22"], dry_run=false, skip_missing=false, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      overwrite_family=true);

   testCase.verifyError(f, ...
      'icemodel:verification:importImau:missingHourlySource');
   testCase.verifyFalse(isfolder(fullfile(eval_root, 'imau', 's21')));
end

function test_imau_hourly_lookup_ignores_flat_daily_qa(testCase)
   % A flat cache may contain hourly and daily PANGAEA files side by side; the
   % hourly locator must not hand the daily QA file to the hourly parser.
   cache = fullfile(testCase.TestData.cache, 'imau-flat-cache');
   mkdir(cache);
   writeImauHourlyFixture(fullfile(cache, ...
      'VanTiggelen-etal_2024_S21.tab'), "S21", ...
      datetime(2014, 4, 12, 0, 0, 0, 'TimeZone', 'UTC'));
   writeImauDailyFixture(fullfile(cache, 'GRL_S21_AWS.tab'), "S21", ...
      datetime(2014, 4, 12, 0, 0, 0, 'TimeZone', 'UTC'));

   filename = icemodel.forcing.helpers.locateImauHourlyFile(cache, "S21");

   testCase.verifyTrue(contains(filename, "VanTiggelen"));
end

function test_imau_import_reports_malformed_daily_qa(testCase)
   % A daily QA artifact with no valid table should be reported in the manifest
   % while the hourly S21 staging path remains intact.
   cache = makeImauSourceCache(testCase.TestData.cache, "S21");
   writeText(fullfile(cache, 'daily', 'GRL_S21_AWS.tab'), ...
      sprintf(['/* DATA DESCRIPTION:\n' ...
      'Citation:\tbad daily QA fixture\n*/\nnot a table\n']));
   eval_root = fullfile(testCase.TestData.cache, 'eval-imau-bad-daily');
   input_root = fullfile(testCase.TestData.cache, 'input-imau-bad-daily');
   mkdir(eval_root);
   mkdir(input_root);

   returned = icemodel.verification.setup.importImau(cache, ...
      case_ids="S21", dry_run=false, evaluation_data_root=eval_root, ...
      input_data_root=input_root, overwrite_family=true);

   c = returned.cases;
   testCase.verifyTrue(logical(c.colocation.imau.staged));
   testCase.verifyFalse(logical(c.colocation.daily_qa.validated));
   testCase.verifyEqual( ...
      string(c.colocation.daily_qa.invalid_sites), "S21");
end

function test_imau_import_stages_rcm_fixtures(testCase)
   % The shared dataset-family RCM helper should attach MAR/MERRA/RACMO legs to
   % an IMAU case when a fast one-point fixture window is available.
   [mar_dir, merra_dir, racmo_dir] = requireRcmFixtures(testCase);
   cache = makeImauSourceCache(testCase.TestData.cache, "S21", ...
      time0=datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'));
   root = fullfile(testCase.TestData.cache, 'imau-rcm-root');

   returned = icemodel.verification.setup.importImau(cache, ...
      case_ids="S21", output_root=root, dry_run=false, ...
      build_forcing=true, ...
      forcing_sources=["imau", "mar", "merra", "racmo"], ...
      mar_dir=mar_dir, merra_dir=merra_dir, ...
      racmo_dir=racmo_dir, overwrite_family=true);

   c = returned.cases;
   testCase.verifyEqual(string(c.case_id), "s21");
   testCase.verifyTrue(any(string(c.eval_sources) == "imau_obs"));
   testCase.verifyTrue(any(string(c.eval_sources) == "racmo2.3p3"));
   testCase.verifyTrue(any(string(c.forcing_sources) == "imau"));
   testCase.verifyTrue(any(string(c.forcing_sources) == "mar3.11"));
   testCase.verifyTrue(any(string(c.forcing_sources) == "merra2"));
   testCase.verifyFalse(any(string(c.forcing_sources) == "racmo2.3p3"));
   prior_merra = jsonencode(c.colocation.merra);
   prior_racmo = jsonencode(c.colocation.racmo);

   % A MAR-only refresh is additive: omitted existing MERRA/RACMO stay put.
   absent_cache = fullfile(root, 'absent-imau-cache');
   refreshed = icemodel.verification.setup.importImau( ...
      absent_cache, ...
      case_ids="S21", output_root=root, dry_run=false, ...
      forcing_sources="mar", build_forcing=true, ...
      build_observations=false, mar_dir=mar_dir, ...
      merra_dir=merra_dir, racmo_dir=racmo_dir, overwrite=true);
   c = refreshed.cases;
   testCase.verifyTrue(logical(c.colocation.merra.staged));
   testCase.verifyTrue(logical(c.colocation.racmo.staged));
   testCase.verifyEqual(jsonencode(c.colocation.merra), prior_merra);
   testCase.verifyEqual(jsonencode(c.colocation.racmo), prior_racmo);
   testCase.verifyTrue(any(string(c.forcing_sources) == "imau"));
   testCase.verifyTrue(any(string(c.eval_sources) == "imau_obs"));
   testCase.verifyTrue(any(string(c.forcing_sources) == "merra2"));
   testCase.verifyTrue(any(string(c.eval_sources) == "racmo2.3p3"));
   testCase.verifyFalse(isfolder(absent_cache));
end

function test_retmip_import_stages_imau_fa_native_source(testCase)
   % RetMIP FA should reuse IMAU S21 native met/Data under the retmip runtime
   % label while remaining a distinct case from the standalone IMAU s21 import.
   retmip_cache = makeRetmipProtocolCache(testCase.TestData.cache, "FA");
   imau_cache = makeImauSourceCache(testCase.TestData.cache, "S21");
   retmip_eval = fullfile(testCase.TestData.cache, 'eval-retmip-fa');
   retmip_input = fullfile(testCase.TestData.cache, 'input-retmip-fa');
   imau_eval = fullfile(testCase.TestData.cache, 'eval-imau-fa');
   imau_input = fullfile(testCase.TestData.cache, 'input-imau-fa');
   mkdir(retmip_eval);
   mkdir(retmip_input);
   mkdir(imau_eval);
   mkdir(imau_input);

   returned = icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids="fa", dry_run=false, imau_dir=imau_cache, ...
      evaluation_data_root=retmip_eval, input_data_root=retmip_input, ...
      forcing_sources="retmip", build_forcing=true, ...
      overwrite_family=true);
   imau = icemodel.verification.setup.importImau(imau_cache, ...
      case_ids="S21", dry_run=false, evaluation_data_root=imau_eval, ...
      input_data_root=imau_input, forcing_sources="imau", ...
      build_forcing=true, overwrite_family=true);

   c = returned.cases;
   testCase.verifyEqual(string(c.case_id), "fa");
   testCase.verifyEqual(string(c.forcing_sources), "retmip");
   testCase.verifyEqual(string(c.eval_sources), "retmip_protocol");
   testCase.verifyEqual(string(c.colocation.source_association.family), ...
      "imau");
   testCase.verifyEqual(string(c.colocation.source_association.source_id), ...
      "S21");
   testCase.verifyEqual(string(c.colocation.retmip.native_met_status), ...
      "staged");
   testCase.verifyEqual(string(c.colocation.native_met.source_family), ...
      "imau");
   testCase.verifyFalse(logical(c.colocation.retmip.forcing_ready));
   testCase.verifyTrue(contains( ...
      string(c.colocation.retmip.forcing_ready_reason), "ppt"));
   testCase.verifyTrue(any(string( ...
      c.observation_variables.native_source.data_variables) == ...
      "surface_height"));
   testCase.verifyFalse(any(string(c.comparison_variables) == ...
      "surface_height"));

   retmip_met = string(c.colocation.retmip.met_files(1));
   imau_met = string(imau.cases.colocation.imau.met_files(1));
   testCase.verifyTrue(startsWith(retmip_met, "retmip"));
   testCase.verifyTrue(startsWith(imau_met, "imau"));
   testCase.verifyTrue(endsWith(retmip_met, "_15m.mat"));
   retmip_met_file = fullfile(retmip_input, 'met', retmip_met);
   retmip_data_file = fullfile(retmip_input, 'userdata', ...
      c.colocation.retmip.data_files(1));
   saved_met = load(retmip_met_file, 'met');
   saved_data = load(retmip_data_file, 'Data');
   testCase.verifyEqual(seconds(median(diff(saved_met.met.Time))), 900);
   testCase.verifyEqual(seconds(median(diff(saved_data.Data.Time))), 3600);
   testCase.verifyEqual(string(imau.cases.case_id), "s21");
   testCase.verifyNotEqual(string(c.case_id), string(imau.cases.case_id));

   % RetMIP's protocol bundle and direct native writer are both idempotent.
   observations_file = fullfile(retmip_eval, 'retmip', c.evaluation_file);
   artifacts = [string(observations_file); string(retmip_met_file); ...
      string(retmip_data_file)];
   original = cellfun(@fileBytes, cellstr(artifacts), 'UniformOutput', false);
   icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids="fa", dry_run=false, imau_dir=imau_cache, ...
      evaluation_data_root=retmip_eval, input_data_root=retmip_input, ...
      forcing_sources="retmip", build_forcing=true, overwrite=false);
   repeated = cellfun(@fileBytes, cellstr(artifacts), 'UniformOutput', false);
   testCase.verifyEqual(repeated, original);
end

function test_retmip_imau_availability_ignores_flat_daily_qa(testCase)
   % Real RetMIP staging must not mistake a flat daily QA file for S21 hourly
   % data when another station makes the hourly product cache present.
   retmip_cache = makeRetmipProtocolCache(testCase.TestData.cache, "FA");
   imau_cache = fullfile(testCase.TestData.cache, 'imau-flat-daily-only');
   output_root = fullfile(testCase.TestData.cache, ...
      'retmip-imau-flat-daily-root');
   mkdir(imau_cache);
   writeImauHourlyFixture(fullfile(imau_cache, ...
      'VanTiggelen-etal_2024_S22.tab'), "S22", ...
      datetime(2014, 4, 12, 0, 0, 0, 'TimeZone', 'UTC'));
   writeImauDailyFixture(fullfile(imau_cache, 'GRL_S21_AWS.tab'), "S21", ...
      datetime(2014, 4, 12, 0, 0, 0, 'TimeZone', 'UTC'));

   returned = icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids="fa", dry_run=false, imau_dir=imau_cache, ...
      output_root=output_root, ...
      forcing_sources="retmip", build_forcing=true, ...
      overwrite_family=true);

   c = returned.cases;
   testCase.verifyEqual(string(c.colocation.retmip.native_met_status), ...
      "missing_imau");
   testCase.verifyTrue(contains( ...
      string(c.colocation.retmip.native_met_skipped_reason), ...
      "missing IMAU hourly source product"));
end

function test_retmip_imau_availability_accepts_flat_hourly(testCase)
   % Real RetMIP staging accepts a flat IMAU hourly file for RetMIP FA.
   retmip_cache = makeRetmipProtocolCache(testCase.TestData.cache, "FA");
   imau_cache = fullfile(testCase.TestData.cache, 'imau-flat-hourly');
   output_root = fullfile(testCase.TestData.cache, ...
      'retmip-imau-flat-hourly-root');
   mkdir(imau_cache);
   writeImauHourlyFixture(fullfile(imau_cache, ...
      'VanTiggelen-etal_2024_S21.tab'), "S21", ...
      datetime(2014, 4, 12, 0, 0, 0, 'TimeZone', 'UTC'));

   returned = icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids="fa", dry_run=false, imau_dir=imau_cache, ...
      output_root=output_root, ...
      forcing_sources="retmip", build_forcing=true, ...
      overwrite_family=true);

   testCase.verifyEqual( ...
      string(returned.cases.colocation.retmip.native_met_status), ...
      "staged");
end

function test_retmip_import_stages_rcm_fixtures(testCase)
   % RetMIP convergence should write native/protocol artifacts first, then add
   % RCM legs through the shared staging helper and refresh source lists.
   [mar_dir, merra_dir, racmo_dir] = requireRcmFixtures(testCase);
   retmip_cache = makeRetmipProtocolCache(testCase.TestData.cache, ...
      ["Dye-2_long", "FA"]);
   gcnet_cache = makeGcnetSurfaceCache(testCase.TestData.cache);
   imau_cache = makeImauSourceCache(testCase.TestData.cache, "S21");
   root = fullfile(testCase.TestData.cache, 'retmip-rcm-root');

   returned = icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids=["dye2", "fa"], output_root=root, dry_run=false, ...
      build_forcing=true, ...
      forcing_sources=["retmip", "mar", "merra", "racmo"], ...
      gcnet_dir=gcnet_cache, imau_dir=imau_cache, ...
      mar_dir=mar_dir, merra_dir=merra_dir, racmo_dir=racmo_dir, ...
      overwrite_family=true);

   ids = string({returned.cases.case_id});
   testCase.verifyEqual(ids, ["dye2_long", "fa"]);
   dye = returned.cases(ids == "dye2_long");
   testCase.verifyTrue(any(string(dye.forcing_sources) == "retmip"));
   testCase.verifyTrue(any(string(dye.forcing_sources) == "mar3.11"));
   testCase.verifyTrue(any(string(dye.forcing_sources) == "merra2"));
   testCase.verifyFalse(any(string(dye.forcing_sources) == "racmo2.3p3"));
   testCase.verifyTrue(any(string(dye.eval_sources) == "retmip_protocol"));
   testCase.verifyTrue(any(string(dye.eval_sources) == "racmo2.3p3"));
   testCase.verifyTrue(logical(dye.colocation.mar.staged));
   testCase.verifyTrue(logical(dye.colocation.merra.staged));
   testCase.verifyTrue(logical(dye.colocation.racmo.staged));
   testCase.verifyFalse(logical(dye.colocation.retmip.forcing_ready));
   prior_merra = jsonencode(dye.colocation.merra);
   prior_racmo = jsonencode(dye.colocation.racmo);

   fa = returned.cases(ids == "fa");
   testCase.verifyTrue(any(string(fa.forcing_sources) == "retmip"));
   testCase.verifyTrue(any(string(fa.eval_sources) == "retmip_protocol"));
   testCase.verifyTrue(isfield(fa.colocation, 'mar'));
   testCase.verifyTrue(isfield(fa.colocation, 'merra'));
   testCase.verifyTrue(isfield(fa.colocation, 'racmo'));
   testCase.verifyFalse(logical(fa.colocation.mar.staged));
   testCase.verifyFalse(logical(fa.colocation.merra.staged));
   testCase.verifyFalse(logical(fa.colocation.racmo.staged));
   testCase.verifyNotEmpty(string(fa.colocation.racmo.reason));

   % Refreshing only MAR for the same case must preserve omitted cached legs.
   absent_cache = fullfile(root, 'absent-retmip-cache');
   refreshed = icemodel.verification.setup.importRetmip( ...
      absent_cache, ...
      case_ids="dye2", output_root=root, dry_run=false, ...
      forcing_sources="mar", build_observations=false, ...
      build_forcing=true, gcnet_dir=gcnet_cache, imau_dir=imau_cache, ...
      mar_dir=mar_dir, merra_dir=merra_dir, racmo_dir=racmo_dir, ...
      overwrite=true);
   dye = refreshed.cases(string({refreshed.cases.case_id}) == "dye2_long");
   testCase.verifyTrue(logical(dye.colocation.merra.staged));
   testCase.verifyTrue(logical(dye.colocation.racmo.staged));
   testCase.verifyEqual(jsonencode(dye.colocation.merra), prior_merra);
   testCase.verifyEqual(jsonencode(dye.colocation.racmo), prior_racmo);
   testCase.verifyTrue(any(string(dye.forcing_sources) == "retmip"));
   testCase.verifyTrue(any(string(dye.eval_sources) == "retmip_protocol"));
   testCase.verifyTrue(any(string(dye.forcing_sources) == "merra2"));
   testCase.verifyTrue(any(string(dye.eval_sources) == "racmo2.3p3"));
   testCase.verifyFalse(isfolder(absent_cache));
end

function test_retmip_partial_rcm_refresh_only_updates_requested_case(testCase)
   % A build_forcing refresh for one RetMIP case must not stage RCM legs for
   % preserved manifest cases outside the current request.
   [mar_dir, merra_dir, racmo_dir] = requireRcmFixtures(testCase);
   retmip_cache = makeRetmipProtocolCache(testCase.TestData.cache, ...
      ["Dye-2_long", "FA"]);
   gcnet_cache = makeGcnetSurfaceCache(testCase.TestData.cache);
   imau_cache = makeImauSourceCache(testCase.TestData.cache, "S21");
   root = fullfile(testCase.TestData.cache, 'retmip-partial-rcm-root');

   icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids=["dye2", "fa"], output_root=root, dry_run=false, ...
      gcnet_dir=gcnet_cache, imau_dir=imau_cache, overwrite_family=true);
   returned = icemodel.verification.setup.importRetmip(retmip_cache, ...
      case_ids="fa", output_root=root, dry_run=false, build_forcing=true, ...
      forcing_sources=["mar", "merra", "racmo"], ...
      gcnet_dir=gcnet_cache, imau_dir=imau_cache, ...
      mar_dir=mar_dir, merra_dir=merra_dir, racmo_dir=racmo_dir, ...
      overwrite=true);

   ids = string({returned.cases.case_id});
   dye = returned.cases(ids == "dye2_long");
   fa = returned.cases(ids == "fa");
   testCase.verifyFalse(isfield(dye.colocation, 'mar'));
   testCase.verifyFalse(isfield(dye.colocation, 'merra'));
   testCase.verifyFalse(isfield(dye.colocation, 'racmo'));
   testCase.verifyTrue(isfield(fa.colocation, 'mar'));
   testCase.verifyTrue(isfield(fa.colocation, 'merra'));
   testCase.verifyTrue(isfield(fa.colocation, 'racmo'));
end

function test_sumup_colocation_records_mixed_anchor_metadata(testCase)
   % sumupColocation should return nearest family/source metadata from explicit
   % mixed anchors instead of assuming every anchor is PROMICE.
   anchors = [ ...
      struct('site', "FA", 'family', "retmip", 'source_id', "fa", ...
         'x_epsg3413', 0, 'y_epsg3413', 0)
      struct('site', "S21", 'family', "imau", 'source_id', "S21", ...
         'x_epsg3413', 50000, 'y_epsg3413', 0)];

   [tf, returned, distance_km] = ...
      icemodel.verification.helpers.sumupColocation(1000, 0, ...
      anchors=anchors, threshold_km=7.5);

   testCase.verifyTrue(tf);
   expected = 1;
   testCase.verifyEqual(distance_km, expected);
   expected = "retmip";
   testCase.verifyEqual(string(returned.family), expected);
   expected = "fa";
   testCase.verifyEqual(string(returned.source_id), expected);

   [tf2, returned2, distance2_km] = ...
      icemodel.verification.setup.anchorColocation(1000, 0, ...
      anchors=anchors, threshold_km=7.5);

   % The SUMup compatibility helper and the setup-level generic helper should
   % report the same nearest mixed anchor.
   testCase.verifyEqual(tf2, tf);
   testCase.verifyEqual(distance2_km, distance_km);
   testCase.verifyEqual(string(returned2.family), "retmip");
end

function test_import_sumup_inherits_mixed_anchor_regime_fields(testCase)
   % SUMup cases colocated with non-PROMICE anchors should inherit the staged
   % anchor's regime fields instead of falling back to unknown.
   proj = icemodel.forcing.helpers.psnProjection();
   point = [66.2, -39.0];
   [x3413, y3413] = projfwd(proj, point(1), point(2));
   anchors = struct('site', "FA", 'family', "retmip", ...
      'source_id', "fa", 'case_id', "fa", ...
      'lat_wgs84', point(1), 'lon_wgs84', point(2), ...
      'x_epsg3413', x3413, 'y_epsg3413', y3413, ...
      'surface_zone', "percolation", 'eval_target', "firn", ...
      'permafrost_zone', "none");

   returned = icemodel.verification.setup.importSumup( ...
      fullfile(testCase.TestData.cache, 'missing-sumup'), ...
      points=point, anchors=anchors, dry_run=true);

   testCase.verifyEqual(string(returned.cases.surface_zone), "percolation");
   testCase.verifyEqual(string(returned.cases.eval_target), "firn");
   testCase.verifyEqual(string(returned.cases.permafrost_zone), "none");
end

function test_sumup_colocation_keeps_nearest_when_outside_threshold(testCase)
   % The nearest anchor is still useful metadata even when the point is outside
   % the colocation threshold.
   anchors = struct('site', "HUMPHREY", 'family', "research_site", ...
      'source_id', "humphrey", 'x_epsg3413', 0, 'y_epsg3413', 0);

   [tf, returned, distance_km] = ...
      icemodel.verification.helpers.sumupColocation(20000, 0, ...
      anchors=anchors, threshold_km=7.5);

   testCase.verifyFalse(tf);
   expected = 20;
   testCase.verifyEqual(distance_km, expected);
   expected = "HUMPHREY";
   testCase.verifyEqual(string(returned.site), expected);
   expected = "research_site";
   testCase.verifyEqual(string(returned.family), expected);
end

function test_research_site_catalog_defines_humphrey_catchall(testCase)
   % Humphrey belongs to the generic research_site family, not a bespoke anchor
   % category or a SUMup subfamily.
   returned = icemodel.verification.setup.researchSiteCatalog("humphrey");

   expected = "research_site";
   testCase.verifyEqual(string(returned.family), expected);
   expected = "humphrey";
   testCase.verifyEqual(string(returned.source_id), expected);
   testCase.verifyEqual(string(returned.period.start), ...
      "2007-01-01 00:00:00");
   testCase.verifyEqual(string(returned.period.end), ...
      "2008-12-31 23:00:00");
   testCase.verifyEqual(string(returned.surface_zone), "percolation");
   testCase.verifyEqual(string(returned.permafrost_zone), "none");
   testCase.verifyGreaterThan(returned.lat_wgs84, 69);
   testCase.verifyLessThan(returned.lon_wgs84, -48);
end

function test_mixed_anchor_catalog_reads_staged_manifests(testCase)
   % mixedAnchorCatalog should tolerate partial eval trees and expose source
   % availability from the staged manifest rows it does find.
   eval_root = fullfile(testCase.TestData.cache, 'catalog-eval');
   writeTinyManifest(eval_root, "promice", "cp1", "CP1", 69.0, -49.0, ...
      1, 2, {'promice'}, {'promice'}, {'promice_obs'});
   writeTinyManifest(eval_root, "retmip", "fa", "FA", 66.2, -39.0, ...
      5, 6, {}, {}, {'retmip_protocol'});
   writeTinyManifest(eval_root, "imau", "s21", "S21", 66.2, -39.1, ...
      7, 8, {}, {}, {'imau_obs'});
   writeTinyManifest(eval_root, "research_site", "humphrey", "humphrey", ...
      69.725714, -48.190512, 3, 4, {}, {}, {'sumup_obs'});

   anchors = icemodel.verification.setup.mixedAnchorCatalog( ...
      evaluation_data_root=eval_root);

   testCase.verifyEqual(numel(anchors), 4);
   families = sort(string({anchors.family}));
   testCase.verifyEqual(families, ["imau", "promice", "research_site", ...
      "retmip"]);
   cp1 = anchors(string({anchors.case_id}) == "cp1");
   testCase.verifyEqual(string(cp1.site), "CP1");
   testCase.verifyEqual(string(cp1.met_sources), "promice");
   testCase.verifyEqual(string(cp1.userdata_sources), "promice");
   fa = anchors(string({anchors.case_id}) == "fa");
   testCase.verifyEqual(string(fa.family), "retmip");
   testCase.verifyEqual(string(fa.eval_sources), "retmip_protocol");
   testCase.verifyEqual(string(fa.surface_zone), "percolation");
   testCase.verifyEqual(string(fa.permafrost_zone), "none");
end

function test_mixed_anchor_catalog_derives_projected_coords(testCase)
   % Finite lat/lon should keep anchors usable when manifests omit EPSG:3413.
   eval_root = fullfile(testCase.TestData.cache, 'catalog-derived-eval');
   writeTinyManifest(eval_root, "retmip", "summit", "Summit", ...
      72.57972, -38.50454, NaN, NaN, {}, {}, {'retmip_protocol'});

   anchors = icemodel.verification.setup.mixedAnchorCatalog( ...
      evaluation_data_root=eval_root);
   testCase.verifyEqual(numel(anchors), 1);
   summit = anchors(1);
   proj = icemodel.forcing.helpers.psnProjection();
   [expected_x, expected_y] = projfwd(proj, 72.57972, -38.50454);

   testCase.verifyTrue(isfinite(summit.x_epsg3413));
   testCase.verifyTrue(isfinite(summit.y_epsg3413));
   testCase.verifyEqual(summit.x_epsg3413, expected_x, 'AbsTol', 1e-6);
   testCase.verifyEqual(summit.y_epsg3413, expected_y, 'AbsTol', 1e-6);
end

function test_import_research_sites_dry_run_is_source_free(testCase)
   % A clean-machine preview should need no SUMup, PROMICE, RCM, or stage roots.
   root = fullfile(testCase.TestData.cache, 'research-dry-source-free');
   missing_sumup = fullfile(root, 'sumup');
   eval_root = fullfile(root, 'eval');
   input_root = fullfile(root, 'input');
   mar_root = fullfile(root, 'mar');
   merra_root = fullfile(root, 'merra');
   racmo_root = fullfile(root, 'racmo');
   returned = icemodel.verification.setup.importResearchSites(missing_sumup, ...
       family="research_site", observation_source="sumup", ...
       case_ids="humphrey", dry_run=true, build_forcing=true, ...
       skip_missing=false, evaluation_data_root=eval_root, ...
       input_data_root=input_root, mar_dir=mar_root, ...
       merra_dir=merra_root, racmo_dir=racmo_root);

   testCase.verifyEqual(string(returned.dataset_family), "research_site");
   testCase.verifyEqual(string(returned.cases.case_id), "humphrey");
   testCase.verifyEqual(string(fieldnames(returned.cases)), ...
      icemodel.verification.setup.firnCaseManifestFieldNames());
   testCase.verifyEqual(string(returned.cases.period.start), ...
      "2007-01-01 00:00:00");
   testCase.verifyEqual(string(returned.cases.period.end), ...
      "2008-12-31 23:00:00");
   testCase.verifyEqual(string(returned.cases.surface_zone), "percolation");
   testCase.verifyEqual(string(returned.cases.permafrost_zone), "none");
   testCase.verifyEmpty(returned.cases.forcing_sources);
   testCase.verifyEqual(string(returned.cases.eval_sources), "sumup_obs");
   testCase.verifyTrue(any(string(returned.cases.comparison_variables) == "smb"));
   testCase.verifyTrue(isfield(returned.cases.observation_variables, ...
      'density'));
   testCase.verifyTrue(isfield(returned.cases.colocation, ...
      'nearest_noncolocated_promice'));
   nearest = returned.cases.colocation.nearest_noncolocated_promice;
   testCase.verifyEqual(string(nearest.nearest_anchor), "");
   testCase.verifyEqual(string(nearest.nearest_source_id), "");
   testCase.verifyTrue(isnan(nearest.distance_km));
   testCase.verifyEqual(nearest.threshold_km, 0);
   testCase.verifyFalse(logical(nearest.is_colocated));
   testCase.verifyFalse(logical(returned.cases.colocation.research_site.staged));
   testCase.verifyFalse(logical(returned.cases.colocation.sumup.staged));
   testCase.verifyEqual(string(returned.cases.evaluation_file), "");
   testCase.verifyEqual(string(returned.cases.colocation.sumup.obs_file), "");
   supplied_roots = [root, missing_sumup, eval_root, input_root, mar_root, ...
      merra_root, racmo_root];
   testCase.verifyFalse(any(isfolder(supplied_roots)));
end

function test_import_research_sites_invalid_window_precedes_root_resolution(testCase)
   % Invalid paired input wins before a real import can create output/cache roots.
   output_root = fullfile(testCase.TestData.cache, ...
      'research-invalid-window-output');
   testCase.verifyError(@() ...
      icemodel.verification.setup.importResearchSites("", ...
      output_root=output_root, startdate="2012-01-01"), ...
      'icemodel:internal:pairedWindow:invalidWindow');
   testCase.verifyFalse(isfolder(output_root));
end

function test_import_research_sites_dry_run_skips_rcm_probe(testCase)
   % Dry runs should not require RCM archives even when build_forcing is true.
   eval_root = fullfile(testCase.TestData.cache, 'research-dry-rcm-eval');
   writePromiceAnchorManifest(eval_root);

   returned = icemodel.verification.setup.importResearchSites("", ...
      family="research_site", observation_source="sumup", ...
      case_ids="humphrey", dry_run=true, build_forcing=true, ...
      evaluation_data_root=eval_root, ...
      input_data_root=fullfile(testCase.TestData.cache, 'research-dry-rcm-input'), ...
      mar_dir="/no/mar", merra_dir="/no/merra", racmo_dir="/no/racmo");

   testCase.verifyEqual(string(returned.dataset_family), "research_site");
   testCase.verifyEqual(string(returned.cases.case_id), "humphrey");
end

function test_import_research_sites_strict_preflights_sumup_cache(testCase)
   % Non-dry research-site import should fail before partial staging when the
   % backing SUMup cache is absent and skip_missing is disabled.
   source_dir = fullfile(testCase.TestData.cache, 'empty-sumup-cache');
   eval_root = fullfile(testCase.TestData.cache, 'research-strict-eval');
   input_root = fullfile(testCase.TestData.cache, 'research-strict-input');

   f = @() icemodel.verification.setup.importResearchSites(source_dir, ...
      case_ids="humphrey", evaluation_data_root=eval_root, ...
      input_data_root=input_root, skip_missing=false); %#ok<NASGU>
   evalc("testCase.verifyError(f, " + ...
      """icemodel:verification:fetchSumup:missingSources"");");

   testCase.verifyFalse(isfolder(fullfile(eval_root, 'research_site')));
end

function test_import_research_sites_stages_humphrey_observations(testCase)
   % Non-dry Humphrey import should stage SUMup-derived observations and keep
   % native research_site met absent.
   sumup_dir = requireSumupCache(testCase);
   eval_root = fullfile(testCase.TestData.cache, 'research-obs-eval');
   input_root = fullfile(testCase.TestData.cache, 'research-obs-input');
   writePromiceAnchorManifest(eval_root);

   returned = icemodel.verification.setup.importResearchSites(sumup_dir, ...
      case_ids="humphrey", startdate="2010-07-01", ...
      enddate="2012-04-30", evaluation_data_root=eval_root, ...
      input_data_root=input_root, overwrite_family=true);

   c = returned.cases;
   testCase.verifyEqual(string(c.case_id), "humphrey");
   testCase.verifyTrue(isfile(fullfile(eval_root, 'research_site', ...
      'humphrey', 'observations.mat')));
   testCase.verifyEmpty(c.forcing_sources);
   testCase.verifyEqual(string(c.eval_sources), "sumup_obs");
   testCase.verifyTrue(any(string(c.comparison_variables) == "smb"));
   testCase.verifyTrue(logical(c.colocation.sumup.staged));
   testCase.verifyFalse(logical(c.colocation.research_site.staged));
   nearest = c.colocation.nearest_noncolocated_promice;
   testCase.verifyGreaterThan(nearest.distance_km, nearest.threshold_km);
   testCase.verifyTrue(any(string(nearest.nearest_anchor) == ...
      ["CP1", "SWC", "JAR"]));

   % Research-site observation refreshes are additive for an existing case.
   observation_file = fullfile(eval_root, 'research_site', ...
      'humphrey', 'observations.mat');
   original = fileBytes(observation_file);
   icemodel.verification.setup.importResearchSites(sumup_dir, ...
      case_ids="humphrey", startdate="2010-07-01", ...
      enddate="2012-04-30", evaluation_data_root=eval_root, ...
      input_data_root=input_root, overwrite=false);
   testCase.verifyEqual(fileBytes(observation_file), original);
end

function test_import_research_sites_skips_empty_sumup_window(testCase)
   % Existing SUMup files with no selected records should produce a skipped case,
   % not a staged but empty observations.mat target.
   sumup_dir = requireSumupCache(testCase);
   eval_root = fullfile(testCase.TestData.cache, 'research-empty-eval');
   input_root = fullfile(testCase.TestData.cache, 'research-empty-input');
   writePromiceAnchorManifest(eval_root);

   returned = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.importResearchSites(sumup_dir, ...
      case_ids="humphrey", startdate="1800-01-01", enddate="1800-01-02", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      overwrite_family=true, skip_missing=true), ...
      'icemodel:verification:importResearchSites:caseSkipped');

   testCase.verifyEmpty(returned.cases);
   testCase.verifyEqual(string(returned.skipped.site), "humphrey");
   testCase.verifyTrue(contains(string(returned.skipped.reason), ...
      "no SUMup observations found"));
   testCase.verifyFalse(isfile(fullfile(eval_root, 'research_site', ...
      'humphrey', 'observations.mat')));
   testCase.verifyFalse(isfolder(fullfile(eval_root, 'research_site', ...
      'humphrey')));
end

function test_import_sumup_skips_empty_observation_window(testCase)
   % SUMup imports should skip all-empty observation bundles, not stage a case.
   sumup_dir = requireSumupCache(testCase);
   site = icemodel.verification.setup.researchSiteCatalog("humphrey");
   eval_root = fullfile(testCase.TestData.cache, 'sumup-empty-eval');
   input_root = fullfile(testCase.TestData.cache, 'sumup-empty-input');

   returned = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.importSumup(sumup_dir, ...
      points=[site.lat_wgs84, site.lon_wgs84], case_ids="humphrey", ...
      startdate="1800-01-01", enddate="1800-01-02", ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      overwrite_family=true, skip_missing=true), ...
      'icemodel:verification:importSumup:caseSkipped');

   testCase.verifyEmpty(returned.cases);
   testCase.verifyEqual(string(returned.skipped.site), "humphrey");
   testCase.verifyTrue(contains(string(returned.skipped.reason), ...
      "no SUMup observations found"));
   testCase.verifyFalse(isfile(fullfile(eval_root, 'sumup', ...
      'humphrey', 'observations.mat')));
end

function test_import_research_sites_stages_rcm_fixtures(testCase)
   % RCM convenience staging should add MAR/MERRA met+Data and RACMO Data while
   % preserving the SUMup observation target.
   sumup_dir = requireSumupCache(testCase);
   [mar_dir, merra_dir, racmo_dir] = requireRcmFixtures(testCase);
   root = fullfile(testCase.TestData.cache, 'research-rcm-root');
   writePromiceAnchorManifest(fullfile(root, 'eval'));

   returned = icemodel.verification.setup.importResearchSites(sumup_dir, ...
      case_ids="humphrey", output_root=root, startdate="2010-07-01", ...
      enddate="2012-04-30", build_forcing=true, ...
      mar_dir=mar_dir, merra_dir=merra_dir, racmo_dir=racmo_dir, ...
      overwrite_family=true);

   c = returned.cases;
   testCase.verifyEqual(string(c.case_id), "humphrey");
   testCase.verifyTrue(any(string(c.eval_sources) == "sumup_obs"));
   testCase.verifyTrue(any(string(c.eval_sources) == "racmo2.3p3"));
   testCase.verifyTrue(any(string(c.forcing_sources) == "mar3.11"));
   testCase.verifyTrue(any(string(c.forcing_sources) == "merra2"));
   testCase.verifyFalse(any(string(c.forcing_sources) == "racmo2.3p3"));

   for src = ["mar", "merra"]
      leg = c.colocation.(char(src));
      testCase.verifyTrue(logical(leg.staged));
      testCase.verifyTrue(isfield(leg, 'met_files') && ~isempty(leg.met_files));
      testCase.verifyTrue(isfield(leg, 'data_files') && ~isempty(leg.data_files));
      met_file = fullfile(root, 'input', 'met', leg.met_files(1));
      data_file = fullfile(root, 'input', 'userdata', leg.data_files(1));
      saved_met = load(met_file, 'met');
      saved_data = load(data_file, 'Data');
      testCase.verifyTrue(endsWith(string(met_file), "_15m.mat"));
      testCase.verifyEqual(seconds(median(diff(saved_met.met.Time))), 900);
      testCase.verifyEqual(seconds(median(diff(saved_data.Data.Time))), 3600);
   end
   testCase.verifyTrue(logical(c.colocation.racmo.staged));
   testCase.verifyTrue(isfield(c.colocation.racmo, 'data_files') ...
      && ~isempty(c.colocation.racmo.data_files));
   testCase.verifyFalse(isfield(c.colocation.racmo, 'met_files'));
   prior_merra = jsonencode(c.colocation.merra);
   prior_racmo = jsonencode(c.colocation.racmo);

   % A MAR-only update must not erase omitted existing MERRA/RACMO legs.
   absent_cache = fullfile(root, 'absent-sumup-cache');
   refreshed = icemodel.verification.setup.importResearchSites( ...
      absent_cache, ...
      case_ids="humphrey", output_root=root, startdate="2010-07-01", ...
      enddate="2012-04-30", forcing_sources=["", "mar"], ...
      build_forcing=true, build_observations=false, ...
      mar_dir=mar_dir, merra_dir=merra_dir, racmo_dir=racmo_dir, ...
      overwrite=true);
   c = refreshed.cases;
   testCase.verifyTrue(logical(c.colocation.merra.staged));
   testCase.verifyTrue(logical(c.colocation.racmo.staged));
   testCase.verifyEqual(jsonencode(c.colocation.merra), prior_merra);
   testCase.verifyEqual(jsonencode(c.colocation.racmo), prior_racmo);
   testCase.verifyTrue(any(string(c.eval_sources) == "sumup_obs"));
   testCase.verifyTrue(any(string(c.forcing_sources) == "merra2"));
   testCase.verifyTrue(any(string(c.eval_sources) == "racmo2.3p3"));
   testCase.verifyFalse(isfolder(absent_cache));
end

function test_import_sumup_humphrey_matches_research_site_observations(testCase)
   % importSumup remains the dataset importer, but at the Humphrey point it
   % should stage materially the same SUMup observation data as research_site.
   sumup_dir = requireSumupCache(testCase);
   site = icemodel.verification.setup.researchSiteCatalog("humphrey");
   point = [site.lat_wgs84, site.lon_wgs84];
   root = fullfile(testCase.TestData.cache, 'humphrey-equivalence');
   research_eval = fullfile(root, 'research-eval');
   research_input = fullfile(root, 'research-input');
   sumup_eval = fullfile(root, 'sumup-eval');
   sumup_input = fullfile(root, 'sumup-input');
   writePromiceAnchorManifest(research_eval);

   research = icemodel.verification.setup.importResearchSites(sumup_dir, ...
      case_ids="humphrey", startdate="2010-07-01", ...
      enddate="2012-04-30", evaluation_data_root=research_eval, ...
      input_data_root=research_input, overwrite_family=true);
   sumup = icemodel.verification.setup.importSumup(sumup_dir, ...
      points=point, anchors=researchSiteAnchor(site), case_ids="humphrey", ...
      startdate="2010-07-01", enddate="2012-04-30", ...
      evaluation_data_root=sumup_eval, input_data_root=sumup_input, ...
      overwrite_family=true);

   rc = research.cases;
   sc = sumup.cases;
   testCase.verifyEqual(string(rc.case_id), "humphrey");
   testCase.verifyEqual(string(sc.case_id), "humphrey");
   testCase.verifyEqual(string(sc.colocation.anchor.nearest_family), ...
      "research_site");
   testCase.verifyTrue(logical(sc.colocation.anchor.is_colocated));
   testCase.verifyEqual(string(sc.comparison_variables), "smb");
   testCase.verifyEqual(string(sc.observation_variables.present_groups), ...
      "smb");

   r = load(fullfile(research_eval, 'research_site', rc.evaluation_file));
   sumup_file = fullfile(sumup_eval, 'sumup', sc.evaluation_file);
   s = load(sumup_file);
   testCase.verifyEqual(s.targets.data.smb, r.targets.data.smb);

   % SUMup's shared family orchestration reuses an existing observation bundle.
   original = fileBytes(sumup_file);
   icemodel.verification.setup.importSumup(sumup_dir, ...
      points=point, anchors=researchSiteAnchor(site), case_ids="humphrey", ...
      startdate="2010-07-01", enddate="2012-04-30", ...
      evaluation_data_root=sumup_eval, input_data_root=sumup_input, ...
      overwrite=false);
   testCase.verifyEqual(fileBytes(sumup_file), original);
end

function test_import_sumup_rcm_only_reuses_observations_without_source_cache( ...
      testCase)
   % RCM-only SUMup refreshes reuse the case without probing the omitted cache.
   sumup_dir = requireSumupCache(testCase);
   [mar_dir, merra_dir, ~] = requireRcmFixtures(testCase);
   site = icemodel.verification.setup.researchSiteCatalog("humphrey");
   point = [site.lat_wgs84, site.lon_wgs84];
   root = fullfile(testCase.TestData.cache, 'sumup-rcm-only-root');
   icemodel.verification.setup.importSumup(sumup_dir, ...
      points=point, case_ids="humphrey", output_root=root, ...
      startdate="2010-07-01", enddate="2012-04-30", ...
      forcing_sources=["mar", "merra"], build_forcing=true, ...
      mar_dir=mar_dir, merra_dir=merra_dir, overwrite_family=true);
   absent_cache = fullfile(root, 'absent-sumup-cache');

   refreshed = icemodel.verification.setup.importSumup(absent_cache, ...
      points=point, case_ids="humphrey", output_root=root, ...
      startdate="2010-07-01", enddate="2012-04-30", ...
      forcing_sources="mar", build_observations=false, ...
      build_forcing=true, mar_dir=mar_dir, overwrite=true, ...
      overwrite_family=true);

   c = refreshed.cases;
   testCase.verifyTrue(logical(c.colocation.mar.staged));
   testCase.verifyTrue(any(string(c.forcing_sources) == "mar3.11"));
   testCase.verifyTrue(any(string(c.eval_sources) == "sumup_obs"));
   testCase.verifyFalse(isfield(c.colocation, 'merra'));
   testCase.verifyFalse(any(string(c.forcing_sources) == "merra2"));
   testCase.verifyEqual(string(c.site_id), string(site.site_id));
   testCase.verifyFalse(isfolder(absent_cache));
end

function test_import_sumup_skipped_point_preserves_stale_manifest_case(testCase)
   % An ordinary failed SUMup point restage records the skip but preserves the
   % last valid case; overwrite_family is the explicit destructive boundary.
   cache = fullfile(testCase.TestData.cache, 'empty-sumup-cache');
   eval_root = fullfile(testCase.TestData.cache, 'sumup-stale-eval');
   input_root = fullfile(testCase.TestData.cache, 'sumup-stale-input');
   family_root = fullfile(eval_root, 'sumup');
   mkdir(cache);
   mkdir(family_root);
   mkdir(input_root);

   stale = sumupManifestEntry("sumup01");
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "sumup", "", "", "stale-test", string(datetime('today')), stale, ...
      struct('site', {}, 'reason', {}));
   icemodel.verification.setup.writeManifest( ...
      fullfile(family_root, 'manifest.json'), manifest);

   returned = testCase.verifyWarning( ...
      @() icemodel.verification.setup.importSumup(cache, ...
      points=[69, -48], evaluation_data_root=eval_root, ...
      input_data_root=input_root, skip_missing=true), ...
      'icemodel:verification:importSumup:caseSkipped');

   expected = jsondecode(jsonencode(stale));
   testCase.verifyEqual(jsonencode(returned.cases), jsonencode(expected));
   testCase.verifyEqual(string(returned.source_version), "sumup2025");
   testCase.verifyEqual(string(returned.skipped.site), "sumup01");
   testCase.verifyEqual(numel(returned.skipped), 1);
end

function test_import_sumup_default_points_use_requested_eval_root(testCase)
   % Default mixed anchors should come from the caller-selected eval root.
   cache = fullfile(testCase.TestData.cache, 'empty-sumup-default-cache');
   eval_root = fullfile(testCase.TestData.cache, 'sumup-default-eval');
   input_root = fullfile(testCase.TestData.cache, 'sumup-default-input');
   mkdir(cache);
   writePromiceAnchorManifest(eval_root);

   returned = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.importSumup(cache, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      skip_missing=true), ...
      'icemodel:verification:importSumup:caseSkipped');

   testCase.verifyEmpty(returned.cases);
   testCase.verifyGreaterThan(numel(returned.skipped), 0);
end

function test_import_sumup_default_points_use_non_promice_anchors(testCase)
   % Default SUMup staging should use any staged firn anchor family, not only
   % PROMICE.
   cache = fullfile(testCase.TestData.cache, 'empty-sumup-retmip-cache');
   eval_root = fullfile(testCase.TestData.cache, 'sumup-retmip-eval');
   input_root = fullfile(testCase.TestData.cache, 'sumup-retmip-input');
   mkdir(cache);
   writeTinyManifest(eval_root, "retmip", "dye2_2016", "Dye-2", ...
      66.48, -46.28, 1, 2, {}, {}, {'retmip_protocol'});

   returned = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.importSumup(cache, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      skip_missing=true), ...
      'icemodel:verification:importSumup:caseSkipped');

   testCase.verifyEmpty(returned.cases);
   testCase.verifyGreaterThan(numel(returned.skipped), 0);
end

function test_import_sumup_default_points_deduplicate_anchor_locations(testCase)
   % Multiple source-family anchors at one location should request one SUMup case.
   cache = fullfile(testCase.TestData.cache, 'empty-sumup-duplicate-cache');
   eval_root = fullfile(testCase.TestData.cache, 'sumup-duplicate-eval');
   input_root = fullfile(testCase.TestData.cache, 'sumup-duplicate-input');
   mkdir(cache);
   writeTinyManifestEntries(eval_root, "retmip", [ ...
      tinyManifestEntry("dye2_long", "Dye-2_long", 66.48001, -46.27889, ...
      1, 2, {}, {}, {'retmip_protocol'})
      tinyManifestEntry("dye2_2016", "Dye-2_16", 66.48001, -46.27889, ...
      1, 2, {}, {}, {'retmip_protocol'})]);

   returned = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.importSumup(cache, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      skip_missing=true), ...
      'icemodel:verification:importSumup:caseSkipped');

   testCase.verifyEmpty(returned.cases);
   testCase.verifyEqual(numel(returned.skipped), 1);
end

function test_import_sumup_default_points_deduplicate_anchor_case_ids(testCase)
   % Default mixed-anchor staging should not request duplicate case ids when
   % multiple families describe one logical site with slightly different
   % coordinates.
   cache = fullfile(testCase.TestData.cache, 'sumup-duplicate-id-cache');
   eval_root = fullfile(testCase.TestData.cache, 'sumup-duplicate-id-eval');
   input_root = fullfile(testCase.TestData.cache, 'sumup-duplicate-id-input');
   writeTinyManifestEntries(eval_root, "promice", ...
      tinyManifestEntry("kanu", "KAN_U", 67.0003, -47.0243, ...
      NaN, NaN, {'promice'}, {'promice'}, {'promice_obs'}));
   writeTinyManifestEntries(eval_root, "retmip", ...
      tinyManifestEntry("kanu", "KAN_U", 67.0020, -47.0260, ...
      NaN, NaN, {'retmip'}, {'retmip'}, {'retmip_protocol'}));

   returned = icemodel.verification.setup.importSumup(cache, ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      dry_run=true);

   testCase.verifyEqual(numel(returned.cases), 1);
   testCase.verifyEqual(string(returned.cases.case_id), "kanu");
end

function test_import_sumup_default_catalog_uses_2025_canonical_map(testCase)
   % The fully implicit release catalog drops only the audited redundant cases,
   % external Humphrey ownership, and anchors without nearby SUMup records while
   % retaining both Dye-2 targets.
   eval_root = fullfile(testCase.TestData.cache, 'sumup-canonical-eval');
   sites = ["SER_B", "MIT", "TAS_U", "TAS_L", "THU_L", "THU_L2", ...
       "THU_U", "ZAC_A", "ZAC_L", "ZAC_U", "HUMPHREY", "DY2", ...
       "Dye-2_long", "KANT", "S23", "LYN_L", "LYN_T", "NUK_B"];
   writeSumupCatalogManifest(eval_root, sites);

   returned = icemodel.verification.setup.importSumup("", ...
      evaluation_data_root=eval_root, ...
      input_data_root=fullfile(testCase.TestData.cache, 'sumup-canonical-input'), ...
      dry_run=true);

   actual = sort(string({returned.cases.case_id}));
   expected = sort(["mit", "tasl", "thuu", "zacu", "dy2", ...
      "dye-2long", "kant"]);
   testCase.verifyEqual(actual, expected);
end

function test_import_sumup_default_catalog_retains_loser_without_keeper(testCase)
   % Static-map losers are still valid targets when their audited keeper is absent.
   eval_root = fullfile(testCase.TestData.cache, 'sumup-no-keeper-eval');
   sites = ["SER_B", "TAS_U", "THU_L", "THU_L2", "ZAC_A", "ZAC_L", ...
      "DY2", "Dye-2_long"];
   writeSumupCatalogManifest(eval_root, sites);

   returned = icemodel.verification.setup.importSumup("", ...
      evaluation_data_root=eval_root, ...
      input_data_root=fullfile(testCase.TestData.cache, 'sumup-no-keeper-input'), ...
      dry_run=true);

   actual = sort(string({returned.cases.case_id}));
   expected = sort(lower(erase(sites, "_")));
   testCase.verifyEqual(actual, expected);
end

function test_import_sumup_explicit_requests_bypass_2025_canonical_map(testCase)
   % Explicit points/case ids and explicit anchor catalogs retain even mapped
   % losers, no-data anchors, and the SUMup Humphrey alias; only the fully
   % implicit path is pruned.
   anchors = makeSumupCatalogAnchors( ...
      ["SER_B", "MIT", "HUMPHREY", "S23", "LYN_L", "LYN_T", "NUK_B"]);
   points = [[anchors.lat_wgs84].', [anchors.lon_wgs84].'];
   eval_root = fullfile(testCase.TestData.cache, 'sumup-explicit-eval');

   explicit_points = icemodel.verification.setup.importSumup("", ...
      points=points, case_ids=["serb", "mit", "humphrey", "s23", ...
      "lynl", "lynt", "nukb"], ...
      evaluation_data_root=eval_root, dry_run=true);
   testCase.verifyEqual(sort(string({explicit_points.cases.case_id})), ...
      ["humphrey", "lynl", "lynt", "mit", "nukb", "s23", "serb"]);

   explicit_anchors = icemodel.verification.setup.importSumup("", ...
      anchors=anchors, evaluation_data_root=eval_root, dry_run=true);
   testCase.verifyEqual(sort(string({explicit_anchors.cases.case_id})), ...
      ["humphrey", "lynl", "lynt", "mit", "nukb", "s23", "serb"]);
end

function touch(filename)
   %TOUCH Create a tiny file for source-cache presence tests.
   fid = fopen(filename, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, 'placeholder\n');
   clear cleaner
end

function writeText(filename, text)
   %WRITETEXT Write one tiny source-cache fixture.
   fid = fopen(filename, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', text);
   clear cleaner
end

function cache = makeRetmipProtocolCache(root, tokens)
   %MAKERETMIPPROTOCOLCACHE Create tiny RetMIP forcing files for importer tests.
   cache = fullfile(root, 'retmip-protocol-cache');
   forcing = fullfile(cache, 'forcing');
   outputs = fullfile(cache, 'outputs');
   scripts = fullfile(cache, 'scripts');
   mkdir(forcing);
   mkdir(outputs);
   mkdir(scripts);
   touch(fullfile(scripts, 'README.md'));

   for token = reshape(tokens, 1, [])
      writeText(fullfile(forcing, ...
         "RetMIP_surface_forcing_" + token + ".tab"), ...
         sprintf(['time;melt_mmweq;acc_subl_mmweq;Tsurf_K\n' ...
         '"01-Jan-2001 00:00:00;0;0;260"\n' ...
         '"01-Jan-2001 03:00:00;1;2;261"\n']));
      touch(fullfile(outputs, "RetMIP_outputs_" + token + ".tab"));
   end
end

function cache = makeGcnetSurfaceCache(root)
   %MAKEGCNETSURFACECACHE Create tiny Vandecrux surface NetCDF fixtures.
   cache = fullfile(root, 'gcnet-surface-cache');
   mkdir(cache);
   touch(fullfile(cache, ...
      "Gap_filled_meteorological_data_and_surface_energy.xml"));
   writeGcnetSurfaceFile(fullfile(cache, 'DYE_2_surface.nc'), ...
      datetime(2000, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'), 0);
   writeGcnetSurfaceFile(fullfile(cache, 'Summit_surface.nc'), ...
      datetime(2001, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'), 10);
end

function cache = makeSamimiWorkbookCache(root)
   %MAKESAMIMIWORKBOOKCACHE Create a tiny Samimi Dye-2 AWS workbook fixture.
   cache = fullfile(root, 'samimi-cache');
   mkdir(cache);
   Time = datetime(2016, 5, 2, 0, 0, 0) + minutes([0; 30; 60; 90]);
   T = table(Time, [-13; -12; -11; -10], [80; 82; 84; 86], ...
      [780; 781; 782; 783], [-1; 200; 300; 400], ...
      [-2; 100; 120; 160], [250; 251; 252; 253], ...
      [240; 241; 242; 243], [5; 6; 7; 8], [180; 181; 182; 183], ...
      [95; 96; 97; 98], ...
      'VariableNames', {'Var1', 'airT', 'RH', 'P', 'SW in', ...
      'SW out', 'LW in', 'LW out', 'wind', 'winddir', 'dsnow'});
   writetable(T, fullfile(cache, 'Dye2_AWS_Summer2016.xlsx'), ...
      'Sheet', 'WS_TDR_Samira_Greenland_2016_Ha');
end

function cache = makeMalformedSamimiWorkbookCache(root)
   %MAKEMALFORMEDSAMIMIWORKBOOKCACHE Create a Samimi workbook missing columns.
   cache = fullfile(root, 'samimi-malformed-cache');
   mkdir(cache);
   T = table(datetime(2016, 5, 2, 0, 0, 0), ...
      'VariableNames', {'Var1'});
   writetable(T, fullfile(cache, 'Dye2_AWS_Summer2016.xlsx'), ...
      'Sheet', 'WS_TDR_Samira_Greenland_2016_Ha');
end

function cache = makeSamimiWorkbookWithoutWindDirection(root)
   %MAKESAMIMIWORKBOOKWITHOUTWINDDIRECTION Omit only the winddir column.
   cache = fullfile(root, 'samimi-missing-winddir-cache');
   mkdir(cache);
   Time = datetime(2016, 5, 2, 0, 0, 0);
   T = table(Time, -13, 80, 780, 100, 60, 250, 240, 5, 95, ...
      'VariableNames', {'Var1', 'airT', 'RH', 'P', 'SW in', ...
      'SW out', 'LW in', 'LW out', 'wind', 'dsnow'});
   writetable(T, fullfile(cache, 'Dye2_AWS_Summer2016.xlsx'), ...
      'Sheet', 'WS_TDR_Samira_Greenland_2016_Ha');
end

function [c, files] = makeManifestMetListCase( ...
      testCase, label, time_axes, requested_start, requested_end)
   %MAKEMANIFESTMETLISTCASE Write an explicit saved-Time met-list fixture.
   input_root = fullfile(testCase.TestData.cache, "manifest-met-" + label);
   met_root = fullfile(input_root, 'met', 'mar3.11');
   mkdir(met_root);
   files = strings(1, numel(time_axes));
   relative = strings(1, numel(time_axes));
   for n = 1:numel(time_axes)
      times = reshape(time_axes{n}, [], 1);
      cadence = seconds(times(2) - times(1));
      suffix = icemodel.forcing.helpers.metTimestepSuffix(cadence);
      relative(n) = fullfile('mar3.11', sprintf( ...
         'met_fa_mar3.11_%s_part%d_%s.mat', label, n, suffix));
      files(n) = fullfile(input_root, 'met', relative(n));
      met = icemodel.test.fixtures.makeSyntheticMetFile(year(times(1)), ...
         nsteps=numel(times), dt_seconds=cadence);
      met.Time = times;
      save(files(n), 'met');
   end

   period = struct('start', char(requested_start, 'yyyy-MM-dd HH:mm:ss'), ...
      'end', char(requested_end, 'yyyy-MM-dd HH:mm:ss'));
   leg = struct('staged', true, 'met_files', {cellstr(relative)});
   c = struct('case_id', "fa", 'dataset_family', "retmip", ...
      'period', period, 'native_timestep', '15m', ...
      'input_data_root', input_root, 'forcing_sources', {{'mar3.11'}}, ...
      'colocation', struct('mar', leg));
end

function cache = makeImauSourceCache(root, site_ids, kwargs)
   %MAKEIMAUSOURCECACHE Create tiny IMAU hourly/daily PANGAEA fixtures.
   arguments
      root (1, 1) string
      site_ids (1, :) string
      kwargs.daily_sites (1, :) string = site_ids
      kwargs.time0 (1, 1) datetime = ...
         datetime(2014, 4, 12, 0, 0, 0, 'TimeZone', 'UTC')
   end

   cache = fullfile(root, 'imau-cache');
   hourly = fullfile(cache, 'hourly');
   daily = fullfile(cache, 'daily');
   mkdir(hourly);
   mkdir(daily);
   for site = reshape(site_ids, 1, [])
      writeImauHourlyFixture(fullfile(hourly, ...
         "VanTiggelen-etal_2024_" + site + ".tab"), site, kwargs.time0);
   end
   for site = reshape(kwargs.daily_sites, 1, [])
      writeImauDailyFixture(fullfile(daily, "GRL_" + site + "_AWS.tab"), ...
         site, kwargs.time0);
   end
end

function writeImauHourlyFixture(filename, site_id, time0)
   %WRITEIMAUHOURLYFIXTURE Write a two-row PANGAEA-shaped hourly table.
   station_doi = imauHourlyDoi(site_id);
   t1 = char(time0, 'yyyy-MM-dd''T''HH:mm:ss');
   t2 = char(time0 + hours(1), 'yyyy-MM-dd''T''HH:mm:ss');
   metadata = sprintf(['/* DATA DESCRIPTION:\n' ...
      'Citation:\tVan Tiggelen et al. (2024): Hourly data at %s ' ...
      '[dataset]. PANGAEA, https://doi.org/%s,\n' ...
      '\tIn: bundled publication]. PANGAEA, ' ...
      'https://doi.org/10.1594/PANGAEA.971647\n' ...
      'Event(s):\tGRL_%s (%s) * LATITUDE: 66.360000 * ' ...
      'LONGITUDE: -39.310000 * DATE/TIME START: %s * DATE/TIME END: %s ' ...
      '* ELEVATION: 1615.0 m\n' ...
      '*/\n'], site_id, station_doi, site_id, site_id, t1, t2);

   rows = [ ...
      [t1, "2", "10", "-1", "0", "1.1", "1.2", "85", "80", ...
      "4", "5", "180", "900", "100", "50", "200", "250", ...
      "3.4", "3.3", "0.1", "0.2", "3.5", "66.181304", ...
      "-39.042994", "0.7", "-5"]
      [t2, "2", "10", "-2", "1", "1.0", "1.1", "84", "81", ...
      "5", "6", "181", "901", "110", "55", "201", "251", ...
      "3.5", "3.4", "0.2", "0.3", "3.6", "66.181304", ...
      "-39.042994", "0.2", "-4"]];

   text = metadata + strjoin(imauHourlyHeaders(), sprintf('\t')) + newline;
   for k = 1:size(rows, 1)
      text = text + strjoin(rows(k, :), sprintf('\t')) + newline;
   end
   writeText(filename, text);
end

function appendImauHourlyFixtureRow(cache, site_id, time)
   %APPENDIMAUHOURLYFIXTUREROW Extend one hourly fixture without changing DOI.
   filename = fullfile(cache, 'hourly', ...
      "VanTiggelen-etal_2024_" + site_id + ".tab");
   values = [char(time, 'yyyy-MM-dd''T''HH:mm:ss'), "2", "10", "-3", ...
      "2", "0.9", "1.0", "83", "82", "6", "7", "182", "902", ...
      "120", "60", "202", "252", "3.6", "3.5", "0.3", "0.4", ...
      "3.7", "66.181304", "-39.042994", "0.5", "-3"];
   text = fileread(filename);
   writeText(filename, string(text) + strjoin(values, sprintf('\t')) + newline);
end

function headers = imauHourlyHeaders()
   %IMAUHOURLYHEADERS Return the PANGAEA hourly header columns used in tests.
   headers = ["Date/Time", "Height [m] (height 1)", ...
      "Height [m] (height 2)", "TTT [degC]", ...
      "TTT [degC] (corrected at 2m height)", ...
      "Humidity spec [g/kg]", ...
      "Humidity spec [g/kg] (corrected at 2m height)", ...
      "RH [%]", "RH [%] (corrected at 2m height)", "ff [m/s]", ...
      "FF10 [m/s] (corrected at 10m height)", "dd [deg]", ...
      "PPPP [hPa]", "SWD [W/m**2]", "SWU [W/m**2]", ...
      "LWD [W/m**2]", "LWU [W/m**2]", "Surf elev [m] raw", ...
      "Surf elev [m] filtered", "Surf elev change [m] sonic", ...
      "Surf elev change [m] sonic and ablation", "Height [m] boom", ...
      "Latitude", "Longitude", "Alb frac", "Surf temp [degC]"];
end

function writeImauDailyFixture(filename, site_id, time0)
   %WRITEIMAUDAILYFIXTURE Write a minimal daily QA PANGAEA table.
   day = char(time0, 'yyyy-MM-dd');
   station_doi = imauDailyDoi(site_id);
   metadata = sprintf(['/* DATA DESCRIPTION:\n' ...
      'Citation:\tVan Tiggelen et al. (2024): Daily SEB at %s ' ...
      '[dataset]. PANGAEA, https://doi.org/%s,\n' ...
      '\tIn: bundled publication]. PANGAEA, ' ...
      'https://doi.org/10.1594/PANGAEA.970127\n' ...
      'Event(s):\tGRL_%s (%s) * LATITUDE: 66.360000 * ' ...
      'LONGITUDE: -39.310000 * DATE/TIME START: %sT00:00:00 * ' ...
      'DATE/TIME END: %sT00:00:00 * ELEVATION: 1615.0 m\n' ...
      '*/\n'], site_id, station_doi, site_id, site_id, day, day);
   header = strjoin(["Height [m] (height 1)", ...
      "Height [m] (height 2)", "Date/Time", "TTT day m [degC]"], ...
      sprintf('\t'));
   row = strjoin(["2", "10", string(day), "-20"], sprintf('\t'));
   writeText(filename, metadata + header + newline + row + newline);
end

function doi = imauHourlyDoi(site_id)
   %IMAUHOURLYDOI Return the per-station hourly PANGAEA DOI.
   switch string(site_id)
      case "S21"
         doi = "10.1594/PANGAEA.969585";
      case "S22"
         doi = "10.1594/PANGAEA.969629";
      case "S23"
         doi = "10.1594/PANGAEA.969631";
   end
end

function doi = imauDailyDoi(site_id)
   %IMAUDAILYDOI Return the per-station daily QA PANGAEA DOI.
   switch string(site_id)
      case "S21"
         doi = "10.1594/PANGAEA.970148";
      case "S22"
         doi = "10.1594/PANGAEA.970149";
      case "S23"
         doi = "10.1594/PANGAEA.970150";
   end
end

function writeGcnetSurfaceFile(filename, t0, offset)
   %WRITEGCNETSURFACEFILE Write the source channels used by the RetMIP import.
   base = datetime(1900, 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
   time = days(t0 + hours(0:2) - base);
   writeGcnetSurfaceVar(filename, "time", time, ...
      "days since 1900-1-1 0:0:0", "Time (UTC)");
   writeGcnetSurfaceVar(filename, "SRin", [100 200 300] + offset, ...
      "Wm-2", "Incoming shortwave radiation");
   writeGcnetSurfaceVar(filename, "SRout", [70 120 150] + offset, ...
      "Wm-2", "Outgoing shortwave radiation");
   writeGcnetSurfaceVar(filename, "LRin", [250 251 252] + offset, ...
      "Wm-2", "Downward longwave radiation flux");
   writeGcnetSurfaceVar(filename, "LRout", [240 241 242] + offset, ...
      "Wm-2", "Upward longwave radiation flux");
   writeGcnetSurfaceVar(filename, "SHF", [1 2 3], "Wm-2", ...
      "Sensible heat flux");
   writeGcnetSurfaceVar(filename, "LHF", [4 5 6], "Wm-2", ...
      "Latent heat flux");
   writeGcnetSurfaceVar(filename, "H_surf_obs", [10 10.1 10.2], "m", ...
      "Observed surface height");
   writeGcnetSurfaceVar(filename, "SMB", [0.01 0.02 0.03], "m_weq", ...
      "Surface mass balance");
   writeGcnetSurfaceVar(filename, "melt", [0.1 0.2 0.3], "m_weq", "Melt");
   writeGcnetSurfaceVar(filename, "sublimation", [-0.01 -0.02 -0.03], ...
      "m_weq", "Net sublimation");
   writeGcnetSurfaceVar(filename, "snowfall", [0 0.001 0.002], "m_weq", ...
      "Snowfall");
   writeGcnetSurfaceVar(filename, "Tsurf", [259 260 261], "K", ...
      "Surface temperature");
   writeGcnetSurfaceVar(filename, "Ta_2m", [260 261 262], "K", ...
      "2m air temperature");
   writeGcnetSurfaceVar(filename, "RH_2m", [70 71 72], "%", ...
      "2m relative humidity with regards to water");
   writeGcnetSurfaceVar(filename, "WS_10m", [5 6 7], "ms-1", ...
      "10 m wind speed");
   writeGcnetSurfaceVar(filename, "Pres", [800 801 802], "hPa", ...
      "air pressure at the surface");
end

function writeGcnetSurfaceVar(filename, name, values, units, long_name)
   %WRITEGCNETSURFACEVAR Write one one-dimensional NetCDF fixture variable.
   nccreate(filename, name, 'Dimensions', {'time', numel(values)});
   ncwrite(filename, name, values);
   ncwriteatt(filename, name, 'units', units);
   ncwriteatt(filename, name, 'long_name', long_name);
end

function writeTinyManifest(eval_root, family, case_id, site_id, lat, lon, x, y, ...
      forcing_sources, data_sources, eval_sources)
   %WRITETINYMANIFEST Write a minimal firn family manifest for catalog tests.
   entry = tinyManifestEntry(case_id, site_id, lat, lon, x, y, ...
      forcing_sources, data_sources, eval_sources);
   if ismember(string(family), ["retmip", "imau", "research_site"])
      % Mixed-anchor fixtures should carry the same curated regime fields that
      % real RetMIP/IMAU/research_site manifests provide.
      entry.surface_zone = "percolation";
      entry.permafrost_zone = "none";
   end
   writeTinyManifestEntries(eval_root, family, entry);
end

function writeTemperatureAliasProfileCase(eval_root)
   %WRITETEMPERATUREALIASPROFILECASE Write a profile bundle alias fixture.
   family_root = fullfile(eval_root, 'sumup');
   case_root = fullfile(family_root, 'tempalias');
   mkdir(case_root);

   Time = datetime(2010, 1, 1:3, 'TimeZone', 'UTC')';
   reference_key = [101; 102; 103];
   observed = [-10; -9; -8];
   subsurface_temperature = timetable(Time, reference_key, observed, ...
      'VariableNames', {'reference_key', 'subsurface_temperature'});
   subsurface_temperature.Properties.DimensionNames{1} = 'Time';
   targets = struct('format', 'subsurface_profile_bundle', ...
      'data', struct('format', 'subsurface_profile_bundle', ...
      'subsurface_temperature', subsurface_temperature), ...
      'metadata', struct('source', 'test'));
   save(fullfile(case_root, 'observations.mat'), 'targets');

   values = { ...
      'tempalias'
      'firn_observational'
      'tempalias'
      'tempalias'
      'unknown'
      {'firn'}
      'unknown'
      struct('lat_wgs84', 69, 'lon_wgs84', -48, ...
      'x_epsg3413', 0, 'y_epsg3413', 0, 'elev_m', 1000)
      struct('start', '2010-01-01 00:00:00', ...
      'end', '2010-01-03 00:00:00')
      char(fullfile('tempalias', 'observations.mat'))
      {}
      {'sumup_obs'}
      {'subsurface_temperature'}
      struct('subsurface_temperature', 'test observation')
      struct('sumup', struct('kind', 'firn_profile_obs', ...
      'staged', true, 'obs_file', char(fullfile('tempalias', ...
      'observations.mat'))))
      'irregular'
      'test profile alias case'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
   writeTinyManifestEntries(eval_root, "sumup", entry);
end

function writeDepthProfileCase(eval_root, family, case_id, density_values)
   %WRITEDEPTHPROFILECASE Write a tiny profile bundle with a density axis.
   family_root = fullfile(eval_root, family);
   case_root = fullfile(family_root, case_id);
   mkdir(case_root);

   density = table([0; 2], density_values(:), ...
      'VariableNames', {'depth', 'density'});
   targets = struct('format', 'subsurface_profile_bundle', ...
      'data', struct('format', 'subsurface_profile_bundle', ...
      'density', density), ...
      'metadata', struct('source', char(family)));
   save(fullfile(case_root, 'observations.mat'), 'targets');

   eval_source = string(family) + "_obs";
   colocation = struct();
   colocation.(char(family)) = struct('kind', 'firn_profile_obs', ...
      'staged', true, 'obs_file', char(fullfile(case_id, ...
      'observations.mat')));
   values = { ...
      char(case_id)
      'firn_observational'
      char(case_id)
      char(case_id)
      'unknown'
      {'firn'}
      'unknown'
      struct('lat_wgs84', 69, 'lon_wgs84', -48, ...
      'x_epsg3413', 0, 'y_epsg3413', 0, 'elev_m', 1000)
      struct('start', '', 'end', '')
      char(fullfile(case_id, 'observations.mat'))
      {}
      cellstr(eval_source)
      {'density'}
      struct('density', 'test density profile')
      colocation
      'irregular'
      'test depth profile case'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
   writeTinyManifestEntries(eval_root, family, entry);
end

function writeTinyManifestEntries(eval_root, family, entries)
   %WRITETINYMANIFESTENTRIES Write a minimal firn manifest with entries.
   family_root = fullfile(eval_root, family);
   if ~isfolder(family_root)
      mkdir(family_root);
   end
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      family, "", "", "test", string(datetime('today')), entries, ...
      struct('site', {}, 'reason', {}));
   icemodel.verification.setup.writeManifest( ...
      fullfile(family_root, "manifest.json"), manifest);
end

function entry = tinyManifestEntry(case_id, site_id, lat, lon, x, y, ...
      forcing_sources, data_sources, eval_sources)
   %TINYMANIFESTENTRY Build one minimal firn manifest entry for catalog tests.
   colocation = struct();
   for k = 1:numel(data_sources)
      src = char(data_sources{k});
      colocation.(src) = struct('kind', src, 'staged', true, ...
         'met_files', {{sprintf('%s-met.mat', src)}}, ...
         'data_files', {{sprintf('%s-data.mat', src)}});
   end
   values = { ...
      char(case_id)
      'firn_observational'
      char(site_id)
      char(site_id)
      'unknown'
      {'firn'}
      'unknown'
      struct('lat_wgs84', lat, 'lon_wgs84', lon, ...
         'x_epsg3413', x, 'y_epsg3413', y, 'elev_m', 1000)
      struct('start', '', 'end', '')
      ''
      forcing_sources
      eval_sources
      {}
      struct()
      colocation
      'irregular'
      'test manifest'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
end

function writeSumupCatalogManifest(eval_root, sites)
   %WRITESUMUPCATALOGMANIFEST Seed an implicit mixed-anchor catalog fixture.
   anchors = makeSumupCatalogAnchors(sites);
   anchor_families = repmat("promice", 1, numel(anchors));
   anchor_families(string({anchors.site}) == "HUMPHREY") = "research_site";
   for family = ["promice", "research_site"]
      selected = find(anchor_families == family);
      if isempty(selected)
         continue
      end
      entries = repmat(tinyManifestEntry("", "", NaN, NaN, NaN, NaN, ...
         {}, {}, {}), 1, numel(selected));
      for n = 1:numel(selected)
         k = selected(n);
         entries(n) = tinyManifestEntry( ...
            anchors(k).case_id, anchors(k).site, ...
            anchors(k).lat_wgs84, anchors(k).lon_wgs84, ...
            anchors(k).x_epsg3413, anchors(k).y_epsg3413, ...
            {}, {}, {'sumup_obs'});
      end
      writeTinyManifestEntries(eval_root, family, entries);
   end
end

function anchors = makeSumupCatalogAnchors(sites)
   %MAKESUMUPCATALOGANCHORS Build distinct valid Greenland catalog anchors.
   proj = icemodel.forcing.helpers.psnProjection();
   template = struct('site', "", 'family', "promice", 'source_id', "", ...
      'case_id', "", 'lat_wgs84', NaN, 'lon_wgs84', NaN, ...
      'x_epsg3413', NaN, 'y_epsg3413', NaN, ...
      'surface_zone', "unknown", 'eval_target', "firn", ...
      'permafrost_zone', "unknown");
   anchors = repmat(template, 1, numel(sites));
   for k = 1:numel(sites)
      lat = 64.5 + 0.35 * k;
      lon = -49 + 0.2 * k;
      [x, y] = projfwd(proj, lat, lon);
      case_id = lower(erase(sites(k), "_"));
      anchors(k) = struct('site', sites(k), 'family', "promice", ...
         'source_id', case_id, 'case_id', case_id, ...
         'lat_wgs84', lat, 'lon_wgs84', lon, ...
         'x_epsg3413', x, 'y_epsg3413', y, ...
         'surface_zone', "unknown", 'eval_target', "firn", ...
         'permafrost_zone', "unknown");
   end
end

function writePromiceAnchorManifest(eval_root)
   %WRITEPROMICEANCHORMANIFEST Write CP1/SWC/JAR anchors near Humphrey.
   site = icemodel.verification.setup.researchSiteCatalog("humphrey");
   proj = icemodel.forcing.helpers.psnProjection();
   [x0, y0] = projfwd(proj, site.lat_wgs84, site.lon_wgs84);
   ids = ["CP1", "SWC", "JAR"];
   dx = [20000, 30000, 40000];
   entries = repmat(emptyPromiceAnchorEntry(), 1, numel(ids));
   for k = 1:numel(ids)
      entries(k) = promiceAnchorEntry(ids(k), x0 + dx(k), y0);
   end

   family_root = fullfile(eval_root, "promice");
   mkdir(family_root);
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "promice", "", "", "test-promice-anchors", ...
      string(datetime('today')), entries, ...
      struct('site', {}, 'reason', {}));
   icemodel.verification.setup.writeManifest( ...
      fullfile(family_root, "manifest.json"), manifest);
end

function anchor = researchSiteAnchor(site)
   %RESEARCHSITEANCHOR Build one mixed-anchor row for a research_site target.
   proj = icemodel.forcing.helpers.psnProjection();
   [x3413, y3413] = projfwd(proj, site.lat_wgs84, site.lon_wgs84);
   anchor = struct('site', string(site.site_id), ...
      'family', "research_site", ...
      'source_id', string(site.source_id), ...
      'case_id', string(site.case_id), ...
      'lat_wgs84', site.lat_wgs84, ...
      'lon_wgs84', site.lon_wgs84, ...
      'x_epsg3413', x3413, ...
      'y_epsg3413', y3413, ...
      'elev_m', site.elev_m);
end

function entry = emptyPromiceAnchorEntry()
   %EMPTYPROMICEANCHORENTRY Prototype PROMICE anchor manifest entry.
   entry = promiceAnchorEntry("ZZZ", NaN, NaN);
end

function entry = promiceAnchorEntry(site_id, x3413, y3413)
   %PROMICEANCHORENTRY Build one staged PROMICE anchor row for tests.
   colocation = struct('promice', struct('kind', 'station_met_and_eval', ...
      'staged', true, 'met_files', {{'promice-met.mat'}}, ...
      'data_files', {{'promice-data.mat'}}));
   values = { ...
      char(lower(site_id))
      'firn_observational'
      char(site_id)
      char(site_id)
      'unknown'
      {'firn'}
      'unknown'
      struct('lat_wgs84', 69, 'lon_wgs84', -48, ...
         'x_epsg3413', x3413, 'y_epsg3413', y3413, 'elev_m', 1500)
      struct('start', '', 'end', '')
      char(fullfile(lower(site_id), 'observations.mat'))
      {'promice'}
      {'promice_obs'}
      {}
      struct()
      colocation
      'daily'
      'test PROMICE anchor'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
end

function entry = sumupManifestEntry(case_id)
   %SUMUPMANIFESTENTRY Build one stale SUMup manifest row for merge tests.
   values = { ...
      char(case_id)
      'firn_observational'
      char(case_id)
      char(case_id)
      ''
      {'firn'}
      ''
      struct('lat_wgs84', 69, 'lon_wgs84', -48, ...
         'x_epsg3413', 0, 'y_epsg3413', 0, 'elev_m', NaN)
      struct('start', '', 'end', '')
      char(fullfile(case_id, 'observations.mat'))
      {}
      {'sumup_obs'}
      {'density'}
      struct()
      struct('sumup', struct('kind', 'firn_profile_obs', ...
         'staged', true, 'obs_file', char(fullfile(case_id, ...
         'observations.mat'))))
      'irregular'
      'stale test row'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
end

function source_dir = requireSumupCache(testCase)
   %REQUIRESUMUPCACHE Return the local SUMup cache or skip the test cleanly.
   source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'sumup'));
   testCase.assumeTrue(isfile(fullfile(source_dir, ...
      'SUMup_2025_SMB_greenland.nc')), ...
      'SUMup verification cache absent; skipping research_site import test.');
end

function source_dir = requirePromiceVerificationCache(testCase, station)
   %REQUIREPROMICEVERIFICATIONCACHE Return the local PROMICE cache or skip.
   source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'promice'));
   required = ["AWS_sites_metadata.csv", "AWS_stations_metadata.csv", ...
      "AWS_variables.csv"];
   has_metadata = all(isfile(fullfile(source_dir, required)));
   has_station = isfile(fullfile(source_dir, 'hour', station + "_hour.nc"));
   testCase.assumeTrue(has_metadata && has_station, ...
      sprintf('PROMICE %s verification cache absent; skipping.', station));
end

function [mar_dir, merra_dir, racmo_dir] = requireRcmFixtures(testCase)
   %REQUIRERCMFIXTURES Return fast RCM fixtures or skip the test cleanly.
   root = string(fullfile(icemodel.internal.fullpath('data'), ...
      'test', 'forcing'));
   mar_dir = firstWithData(fullfile(root, 'mar'), ...
      @(p) ~isempty(dir(fullfile(p, "MARv3.11*.nc"))));
   merra_dir = firstWithData(fullfile(root, 'merra2'), ...
      @(p) ~isempty(dir(fullfile(p, "slv", "*_Nx.*.nc4*"))));
   racmo_dir = firstWithData(fullfile(root, 'racmo'), ...
      @(p) ~isempty(dir(fullfile(p, "*.RACMO*.nc"))));
   testCase.assumeTrue(strlength(mar_dir) > 0 ...
      && strlength(merra_dir) > 0 && strlength(racmo_dir) > 0, ...
      'RCM fast fixtures absent; skipping research_site RCM test.');
end

function p = firstWithData(candidates, hasData)
   %FIRSTWITHDATA First candidate dir that exists and satisfies hasData.
   p = "";
   for c = reshape(string(candidates), 1, [])
      if isfolder(c) && hasData(c)
         p = c;
         return
      end
   end
end

function cache = makePromiceCache(root, station)
   %MAKEPROMICECACHE Create the smallest PROMICE-shaped cache for unit tests.
   cache = fullfile(root, 'promice-cache');
   hour = fullfile(cache, 'hour');
   mkdir(hour);
   stations = reshape(string(station), 1, []);
   rows = join(stations + ",ice sheet", newline);
   writeText(fullfile(cache, "AWS_sites_metadata.csv"), ...
      sprintf('site_id,location_type\n%s\n', rows));
   touch(fullfile(cache, "AWS_stations_metadata.csv"));
   touch(fullfile(cache, "AWS_variables.csv"));
   for one = stations
      touch(fullfile(hour, one + "_hour.nc"));
   end
end

function cache = makeGcnetCache(root)
   %MAKEGCNETCACHE Create a tiny package-subfolder Vandecrux cache fixture.
   cache = fullfile(root, 'gcnet-cache');
   surface = fullfile(cache, 'gap-filled-surface');
   firn_temperature = fullfile(cache, 'firn-temperature');
   simulated = fullfile(cache, 'simulated-firn');
   mkdir(surface);
   mkdir(firn_temperature);
   mkdir(simulated);

   % Package metadata files prove each Arctic Data DOI package is represented.
   touch(fullfile(surface, ...
      "Gap_filled_meteorological_data_and_surface_energy.xml"));
   touch(fullfile(firn_temperature, ...
      "Firn_temperatures_and_measurement_depths_at_nine.xml"));
   touch(fullfile(simulated, ...
      "Simulated_firn_density_temperature_liquid_water.xml"));

   % Station files are deliberately split across product folders to exercise
   % recursive validation rather than flat-cache-only matching.
   for station = ["DYE_2", "Summit"]
      touch(fullfile(surface, station + "_surface.nc"));
      touch(fullfile(firn_temperature, station + "_T_firn_obs.nc"));
      touch(fullfile(simulated, station + "_T_ice_bin_1.nc"));
      touch(fullfile(simulated, station + "_rho_bin_1.nc"));
      touch(fullfile(simulated, station + "_slwc_bin_1.nc"));
      touch(fullfile(simulated, station + "_compaction_bin_1.nc"));
   end
end

function cache = makeGcnetHeaderCache(root, layout)
   %MAKEGCNETHEADERCACHE Create flat or nested valid inventory test headers.
   cache = fullfile(root, "gcnet-header-cache-" + layout);
   if layout == "flat"
      mkdir(cache)
      surface = cache;
      firn_temperature = cache;
      simulated = cache;
   else
      surface = fullfile(cache, 'surface');
      firn_temperature = fullfile(cache, 'firn-temperature');
      simulated = fullfile(cache, 'simulated-firn');
      mkdir(surface);
      mkdir(firn_temperature);
      mkdir(simulated);
   end

   writeGcnetXml(fullfile(surface, ...
      "Gap_filled_meteorological_data_and_surface_energy.xml"), ...
      "10.18739/A2HM52K87", "Gap-filled meteorological data");
   writeGcnetXml(fullfile(firn_temperature, ...
      "Firn_temperatures_and_measurement_depths_at_nine.xml"), ...
      "10.18739/A2833N00P", "Firn temperatures");
   writeGcnetXml(fullfile(simulated, ...
      "Simulated_firn_density_temperature_liquid_water.xml"), ...
      "10.18739/A2CV4BS43", "Simulated firn state");

   for station = ["DYE_2", "Summit"]
      writeTinyGcnetNetcdf(fullfile(surface, station + "_surface.nc"), ...
         ["Ta_2m", "RH_2m"], ["K", "%"], false);
      writeTinyGcnetNetcdf(fullfile(firn_temperature, ...
         station + "_T_firn_obs.nc"), ["T_firn", "Depth"], ...
         ["degC", "m"], true);
      writeTinyGcnetNetcdf(fullfile(simulated, ...
         station + "_T_ice_bin_1.nc"), ["T_ice", "Depth"], ...
         ["K", "m"], true);
      writeTinyGcnetNetcdf(fullfile(simulated, ...
         station + "_rho_bin_1.nc"), ["rho", "Depth"], ...
         ["kg m-3", "m"], true);
      writeTinyGcnetNetcdf(fullfile(simulated, ...
         station + "_slwc_bin_1.nc"), ["slwc", "Depth"], ...
         ["m", "m"], true);
      writeTinyGcnetNetcdf(fullfile(simulated, ...
         station + "_compaction_bin_1.nc"), ["compaction", "Depth"], ...
         ["m", "m"], true);
   end
end

function records = portableGcnetRecords(records)
   %PORTABLEGCNETRECORDS Remove layout-dependent roots before parity checks.
   for k = 1:numel(records)
      [~, name, extension] = fileparts(records(k).filename);
      records(k).filename = string(name) + string(extension);
      [~, name, extension] = fileparts(records(k).provenance.xml_file);
      records(k).provenance.xml_file = string(name) + string(extension);
   end
end

function writeGcnetXml(filename, doi, title)
   %WRITEGCNETXML Write the XML snippets used by provenance parsing tests.
   text = sprintf(['<?xml version="1.0" encoding="UTF-8"?>\n' ...
      '<eml:eml packageId="doi:%s">\n' ...
      '<dataset><title>%s</title><abstract><para>Vandecrux et al. 2020 ' ...
      'test citation.</para></abstract></dataset>\n' ...
      '</eml:eml>\n'], doi, title);
   writeText(filename, text);
end

function writeTinyGcnetNetcdf(filename, variables, units, has_level)
   %WRITETINYGCNETNETCDF Write tiny valid coordinate and data-variable headers.
   nccreate(filename, 'time', 'Dimensions', {'time', 3});
   ncwrite(filename, 'time', [36524 36524 + 1 / 24 36524 + 2 / 24]);
   ncwriteatt(filename, 'time', 'units', 'days since 1900-1-1 0:0:0');

   if has_level
      nccreate(filename, 'level', 'Dimensions', {'level', 2});
      ncwrite(filename, 'level', [1 2]);
      dims = {'time', 3, 'level', 2};
   else
      dims = {'time', 3};
   end

   for k = 1:numel(variables)
      nccreate(filename, variables(k), 'Dimensions', dims);
      if has_level && variables(k) == "Depth"
         data = [1.0 2.0; 1.1 2.1; 1.2 2.2];
      elseif has_level
         data = reshape(1:6, 3, 2);
      else
         data = 1:3;
      end
      ncwrite(filename, variables(k), data);
      ncwriteatt(filename, variables(k), 'units', units(k));
   end
end

function bytes = fileBytes(filename)
   %FILEBYTES Read one staged binary artifact for no-churn assertions.
   fid = fopen(filename, 'r');
   cleanup = onCleanup(@() fclose(fid));
   bytes = fread(fid, Inf, '*uint8');
end

function snapshot = treeSnapshot(root)
   %TREESNAPSHOT Capture relative paths and exact bytes below one staging root.
   listing = dir(fullfile(root, '**', '*'));
   listing = listing(~ismember(string({listing.name}), [".", ".."]));
   fullnames = string(fullfile({listing.folder}, {listing.name}));
   fullnames = fullnames(:);
   paths = extractAfter(fullnames, strlength(string(root)) + 1);
   [paths, order] = sort(paths);
   listing = listing(order);
   fullnames = fullnames(order);

   % Directory entries participate in equality; file entries additionally carry
   % their exact binary payload so content-only changes are visible.
   bytes = cell(numel(listing), 1);
   for k = 1:numel(listing)
      if listing(k).isdir
         bytes{k} = zeros(0, 1, 'uint8');
      else
         bytes{k} = fileBytes(fullnames(k));
      end
   end
   snapshot = struct('paths', paths, ...
      'isdir', reshape([listing.isdir], [], 1), 'bytes', {bytes});
end

function [fields, matching, conflicting] = fixedArtifactIdentityCases()
   %FIXEDARTIFACTIDENTITYCASES Return documented scalar reuse facts.
   fields = ["source_family", "station", "doi", "bundle_doi", ...
      "product", "schema", "method"];
   matching = ["imau", "S21", "10.1594/PANGAEA.969585", ...
      "10.1594/PANGAEA.971647", "hourly_aws", ...
   "verification_timeseries", "nearest"];
   conflicting = ["samimi", "S22", "10.1594/PANGAEA.969629", ...
      "10.1594/PANGAEA.970127", "protocol_bundle", ...
      "retmip_protocol_bundle", "natural"];
end
