function tests = test_dataset_family_import_helpers
   %TEST_DATASET_FAMILY_IMPORT_HELPERS Verify neutral dataset-family staging glue.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Allocate a manifest root for helper-only tests.
   testCase.TestData.tmp = tempname;
   mkdir(testCase.TestData.tmp)
end

function teardown(testCase)
   % Remove temporary manifests.
   if isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
end

function test_selectSiteCatalogEntries_returns_full_catalog_for_empty_ids(testCase)
   % Empty ids select the complete source catalog without reordering it.
   sites = tinySiteCatalog();

   returned = icemodel.verification.setup.selectSiteCatalogEntries( ...
      sites, strings(1, 0), "test:catalog:unknown", "test site");

   testCase.verifyEqual(string({returned.site_id}), ["A", "B", "C"]);
end

function test_selectSiteCatalogEntries_preserves_catalog_order(testCase)
   % Requested ids filter the catalog while retaining its canonical order.
   sites = tinySiteCatalog();

   returned = icemodel.verification.setup.selectSiteCatalogEntries( ...
      sites, ["C", "A"], "test:catalog:unknown", "test site");

   testCase.verifyEqual(string({returned.site_id}), ["A", "C"]);
end

function test_selectSiteCatalogEntries_rejects_unknown_ids(testCase)
   % The shared selector preserves each caller's catchable error identifier.
   f = @() icemodel.verification.setup.selectSiteCatalogEntries( ...
      tinySiteCatalog(), ["B", "missing"], ...
      "test:catalog:unknown", "test site");

   testCase.verifyError(f, 'test:catalog:unknown');
end

function test_loadPriorDatasetFamilyCases_reads_additive_refresh_cases(testCase)
   % An ordinary observation refresh reads the prior normalized case entries.
   manifest_file = writePriorFamilyManifest(testCase);

   returned = ...
      icemodel.verification.setup.loadPriorDatasetFamilyCases(manifest_file);

   testCase.verifyEqual(numel(returned), 1);
   testCase.verifyEqual(string(returned.case_id), "case1");
end

function test_loadPriorDatasetFamilyCases_distinguishes_replacement_paths(testCase)
   % Forcing-only replacement reads its requested case; observation replacement
   % starts without prior state.
   manifest_file = writePriorFamilyManifest(testCase);
   malformed_file = fullfile(testCase.TestData.tmp, 'malformed-prior.json');
   touch(malformed_file);

   forcing_rebuild = ...
      icemodel.verification.setup.loadPriorDatasetFamilyCases(manifest_file);
   forcing_replacement = ...
      icemodel.verification.setup.loadPriorDatasetFamilyCases( ...
      manifest_file, overwrite_family=true, build_observations=false);
   observation_replacement = ...
      icemodel.verification.setup.loadPriorDatasetFamilyCases( ...
      malformed_file, overwrite_family=true, build_observations=true);

   testCase.verifyEqual(string(forcing_rebuild.case_id), "case1");
   testCase.verifyEqual(string(forcing_replacement.case_id), "case1");
   testCase.verifyEmpty(observation_replacement);
end

function test_loadPriorDatasetFamilyCases_tolerates_missing_manifest(testCase)
   % A new family root has no prior cases and does not fail additive staging.
   missing_file = fullfile(testCase.TestData.tmp, 'missing-prior.json');

   returned = ...
      icemodel.verification.setup.loadPriorDatasetFamilyCases(missing_file);

   testCase.verifyEmpty(returned);
end

function test_priorCaseById_returns_exact_match_or_empty(testCase)
   % Shared importer lookup preserves one exact prior entry and misses cleanly.
   cases = struct('case_id', {'first', 'second'}, 'value', {1, 2});

   returned = icemodel.verification.setup.priorCaseById(cases, "second");
   missing = icemodel.verification.setup.priorCaseById(cases, "missing");
   legacy = icemodel.verification.setup.priorCaseById( ...
      rmfield(cases, 'case_id'), "second");

   testCase.verifyEqual(returned, cases(2));
   testCase.verifyEmpty(missing);
   testCase.verifyEmpty(legacy);
end

function test_datasetFamilyStagingPaths_returns_canonical_layout(testCase)
   % Every importer must derive the same family and runtime output paths.
   eval_root = fullfile(testCase.TestData.tmp, 'eval');
   input_root = fullfile(testCase.TestData.tmp, 'input');

   [family_root, manifest_file, met_outdir, userdata_outdir] = ...
      icemodel.verification.setup.datasetFamilyStagingPaths( ...
      eval_root, input_root, "family");

   testCase.verifyEqual(family_root, string(fullfile(eval_root, 'family')));
   testCase.verifyEqual(manifest_file, ...
      string(fullfile(eval_root, 'family', 'manifest.json')));
   testCase.verifyEqual(met_outdir, string(fullfile(input_root, 'met')));
   testCase.verifyEqual(userdata_outdir, ...
      string(fullfile(input_root, 'userdata')));
end

function test_prepareReplacementCaseEntry_clears_only_runtime_state(testCase)
   % Whole-family forcing replacement retains observations while removing every
   % prior native/RCM runtime leg.
   rcm_leg = struct('staged', true, 'met_files', {{'old-met.mat'}}, ...
      'data_files', {{'old-data.mat'}});

   sumup_entry = struct('colocation', struct( ...
      'sumup', struct('staged', true), 'mar', rcm_leg), ...
      'forcing_sources', {{'mar3.11'}}, 'eval_sources', {{'sumup_obs'}});
   sumup_entry = ...
      icemodel.verification.setup.prepareReplacementCaseEntry( ...
      sumup_entry, "sumup");
   testCase.verifyFalse(isfield(sumup_entry.colocation, 'mar'));
   testCase.verifyTrue(isfield(sumup_entry.colocation, 'sumup'));
   testCase.verifyEmpty(sumup_entry.forcing_sources);
   testCase.verifyEqual(string(sumup_entry.eval_sources), "sumup_obs");

   promice_leg = struct('kind', 'station_met_and_eval', 'staged', true, ...
      'eval_staged', true, 'met_files', {{'old-met.mat'}}, ...
      'data_files', {{'old-data.mat'}}, 'forcing_ready', true, ...
      'forcing_ready_reason', 'ready', ...
      'forcing_complete_windows', struct('start_time', '2012-01-01'), ...
      'met_skipped_reason', 'old failure');
   promice_entry = struct('colocation', struct( ...
      'promice', promice_leg, 'racmo', rcm_leg), ...
      'forcing_sources', {{'promice', 'racmo2.3p2'}}, ...
      'eval_sources', {{'promice_obs'}});
   promice_entry = ...
      icemodel.verification.setup.prepareReplacementCaseEntry( ...
      promice_entry, "promice");
   testCase.verifyFalse(logical(promice_entry.colocation.promice.staged));
   testCase.verifyTrue(logical(promice_entry.colocation.promice.eval_staged));
   testCase.verifyEmpty(promice_entry.colocation.promice.met_files);
   testCase.verifyEmpty(promice_entry.colocation.promice.data_files);
   testCase.verifyFalse(isfield( ...
      promice_entry.colocation.promice, 'forcing_ready'));
   testCase.verifyFalse(isfield(promice_entry.colocation, 'racmo'));
   testCase.verifyEmpty(promice_entry.forcing_sources);
   testCase.verifyEqual(string(promice_entry.eval_sources), "promice_obs");

   imau_leg = struct('kind', 'hourly_aws_met_and_eval', 'staged', true, ...
      'met_files', {{'old-met.mat'}}, 'data_files', {{'old-data.mat'}}, ...
      'forcing_ready', true, 'forcing_ready_reason', 'ready', ...
      'forcing_complete_windows', struct('start_time', '2012-01-01'), ...
      'artifact_metadata', struct('source_id', 'S21'));
   imau_entry = struct('colocation', struct('imau', imau_leg, 'mar', rcm_leg), ...
      'forcing_sources', {{'imau', 'mar3.11'}}, ...
      'eval_sources', {{'imau_obs'}});
   imau_entry = ...
      icemodel.verification.setup.prepareReplacementCaseEntry( ...
      imau_entry, "imau");
   testCase.verifyFalse(logical(imau_entry.colocation.imau.staged));
   testCase.verifyTrue(logical(imau_entry.colocation.imau.eval_staged));
   testCase.verifyEmpty(imau_entry.colocation.imau.met_files);
   testCase.verifyEmpty(imau_entry.colocation.imau.data_files);
   testCase.verifyFalse(isfield( ...
      imau_entry.colocation.imau, 'artifact_metadata'));
   testCase.verifyFalse(isfield(imau_entry.colocation, 'mar'));
   testCase.verifyEmpty(imau_entry.forcing_sources);
   testCase.verifyEqual(string(imau_entry.eval_sources), "imau_obs");

   retmip_leg = struct('kind', 'protocol_userdata_and_native_met', ...
      'staged', true, 'protocol_timestep', '3hr', ...
      'native_source', struct('family', 'gcnet_vandecrux'), ...
      'native_met_status', 'staged', 'met_files', {{'old-met.mat'}}, ...
      'data_files', {{'old-data.mat'}}, 'window', struct(), ...
      'forcing_ready', true, 'forcing_ready_reason', 'ready');
   retmip_entry = struct('colocation', struct('retmip', retmip_leg, ...
      'native_met', struct('staged', true), 'merra', rcm_leg), ...
      'forcing_sources', {{'retmip', 'merra2'}}, ...
      'eval_sources', {{'retmip_protocol'}}, ...
      'observation_variables', struct('protocol_id', 'dye2', ...
      'native_source', struct('family', 'gcnet_vandecrux')), ...
      'native_timestep', '1hr');
   retmip_entry = ...
      icemodel.verification.setup.prepareReplacementCaseEntry( ...
      retmip_entry, "retmip");
   testCase.verifyFalse(isfield(retmip_entry.colocation, 'native_met'));
   testCase.verifyFalse(isfield(retmip_entry.colocation, 'merra'));
   testCase.verifyFalse(isfield( ...
      retmip_entry.colocation.retmip, 'native_source'));
   testCase.verifyFalse(isfield( ...
      retmip_entry.observation_variables, 'native_source'));
   testCase.verifyEqual(string(retmip_entry.native_timestep), "3hr");
   testCase.verifyEmpty(retmip_entry.forcing_sources);
   testCase.verifyEqual(string(retmip_entry.eval_sources), "retmip_protocol");

   plain = struct('case_id', 'plain');
   testCase.verifyEqual( ...
      icemodel.verification.setup.prepareReplacementCaseEntry( ...
      plain, "plain"), plain);
end

function test_rcmSourceSelection_accepts_optional_native_source(testCase)
   % Native-plus-RCM importers use the validator with their family source.
   validator = @icemodel.verification.validators.mustBeRcmSourceSelection;

   testCase.verifyWarningFree(@() validator(["mar", "racmo"]));
   testCase.verifyWarningFree(@() validator(["promice", "merra"], "promice"));
   testCase.verifyWarningFree(@() validator(["imau", "merra"], "imau"));
   testCase.verifyWarningFree(@() validator(["retmip", "racmo"], "retmip"));
   testCase.verifyWarningFree(@() validator("", "promice"));
   testCase.verifyError(@() validator("promice"), ...
      'icemodel:validators:mustBeRcmSourceSelection');
   testCase.verifyError(@() validator("unknown", "promice"), ...
      'icemodel:validators:mustBeRcmSourceSelection');
end

function test_normalizeForcingSources_enforces_ordered_build_selection(testCase)
   % Every non-atomic importer shares blank filtering, stable order, and gating.
   normalize = @icemodel.verification.setup.normalizeForcingSources;

   testCase.verifyEqual(normalize(["", "mar", "mar", "racmo"], true), ...
      ["mar", "racmo"]);
   testCase.verifyEmpty(normalize(["", ""], false));
   testCase.verifyError(@() normalize(["", ""], true), ...
      'icemodel:verification:normalizeForcingSources:emptySelection');
end

function test_stateCaseEntry_preserves_atomic_state(testCase)
   % Atomic family state has no colocation graph and returns its entry verbatim.
   expected = struct('case_id', 'atomic', 'notes', 'complete bundle');
   state = struct('case_id', "atomic", 'entry', expected);

   returned = icemodel.verification.setup.stateCaseEntry(state);

   testCase.verifyEqual(returned, expected);
end

function test_laughCaseValidator_accepts_vectors_and_reports_unknown(testCase)
   % The canonical plural selector validates every requested atomic case.
   validator = @icemodel.verification.validators.mustBeLaughTestCase;
   valid = icemodel.verification.namelists.caseid("laugh_tests");

   testCase.verifyWarningFree(@() validator(valid));
   testCase.verifyError(@() validator([valid, "unknown"]), ...
      'icemodel:verification:validators:mustBeLaughTestCase');
end

function test_importers_use_case_ids_and_canonical_manifest_fields(testCase)
   % Every family accepts case_ids and emits one of the two persisted schemas.
   output_root = fullfile(testCase.TestData.tmp, 'api-contract');
   atomic = { ...
      icemodel.verification.setup.importEsmSnowmip( ...
      case_ids="cdp", output_root=output_root, dry_run=true), ...
      icemodel.verification.setup.importLaughTests( ...
      case_ids="colbeck1976", output_root=output_root, dry_run=true)};
   firn = { ...
      icemodel.verification.setup.importPromiceSites( ...
      case_ids="KAN_M", output_root=output_root, dry_run=true, ...
      forcing_sources="promice", build_forcing=false), ...
      icemodel.verification.setup.importImau( ...
      case_ids="S21", output_root=output_root, dry_run=true, ...
      forcing_sources=strings(1, 0)), ...
      icemodel.verification.setup.importRetmip( ...
      case_ids="kanu", output_root=output_root, dry_run=true, ...
      forcing_sources=strings(1, 0)), ...
      icemodel.verification.setup.importResearchSites( ...
      case_ids="humphrey", output_root=output_root, dry_run=true, ...
      forcing_sources=strings(1, 0)), ...
      icemodel.verification.setup.importSumup( ...
      points=[67.1, -48.1], case_ids="api_contract", ...
      anchors=struct([]), output_root=output_root, dry_run=true, ...
      forcing_sources=strings(1, 0))};

   family_fields = string( ...
      icemodel.verification.setup.familyManifestFieldNames());
   atomic_fields = ...
      icemodel.verification.setup.caseManifestFieldNames();
   firn_fields = ...
      icemodel.verification.setup.firnCaseManifestFieldNames();
   for manifest = [atomic, firn]
      testCase.verifyEqual(string(fieldnames(manifest{1})), family_fields);
   end
   for manifest = atomic
      testCase.verifyEqual(string(fieldnames(manifest{1}.cases)), atomic_fields);
   end
   for manifest = firn
      testCase.verifyEqual(string(fieldnames(manifest{1}.cases)), firn_fields);
   end
   testCase.verifyFalse(isfolder(output_root));
end

function test_pairedWindow_normalizes_optional_utc_bounds(testCase)
   % The shared helper supports blank, string, and datetime paired inputs.
   [t1, t2, enabled] = ...
      icemodel.internal.pairedWindow("", "");
   testCase.verifyFalse(enabled);
   testCase.verifyTrue(isnat(t1) && isnat(t2));
   testCase.verifyEqual(t1.TimeZone, 'UTC');

   [t1, t2, enabled] = icemodel.internal.pairedWindow( ...
      "2012-01-01", "2012-01-02");
   testCase.verifyTrue(enabled);
   testCase.verifyEqual(t1, datetime(2012, 1, 1, 'TimeZone', 'UTC'));
   testCase.verifyEqual(t2, datetime(2012, 1, 2, 'TimeZone', 'UTC'));

   nat = NaT('TimeZone', 'UTC');
   [~, ~, enabled] = icemodel.internal.pairedWindow(nat, nat);
   testCase.verifyFalse(enabled);

   % Zoned inputs retain their instant while sharing one UTC representation.
   ny_start = datetime(2012, 1, 1, 0, 0, 0, ...
      'TimeZone', 'America/New_York');
   ny_end = ny_start + hours(1);
   [t1, t2, enabled] = icemodel.internal.pairedWindow(ny_start, ny_end);
   testCase.verifyTrue(enabled);
   testCase.verifyEqual(t1.TimeZone, 'UTC');
   testCase.verifyEqual(t1, datetime(2012, 1, 1, 5, 0, 0, ...
      'TimeZone', 'UTC'));
   testCase.verifyEqual(t2, t1 + hours(1));
end

function test_pairedWindow_rejects_half_reversed_and_vector_bounds(testCase)
   % Invalid public windows share one stable error before downstream reads.
   error_id = 'icemodel:internal:pairedWindow:invalidWindow';
   testCase.verifyError(@() ...
      icemodel.internal.pairedWindow("2012-01-01", ""), ...
      error_id);
   testCase.verifyError(@() ...
      icemodel.internal.pairedWindow("", "2012-01-02"), ...
      error_id);
   testCase.verifyError(@() ...
      icemodel.internal.pairedWindow( ...
      "2012-01-02", "2012-01-01"), error_id);
   testCase.verifyError(@() ...
      icemodel.internal.pairedWindow( ...
      ["2012-01-01", "2012-01-02"], "2012-01-03"), error_id);
   testCase.verifyError(@() ...
      icemodel.internal.pairedWindow("not-a-date", "still-not-a-date"), ...
      error_id);
end

function test_preserveRcmLegs_keeps_only_compatible_failed_refresh(testCase)
   % The shared preservation helper covers both dataset and manifest staging.
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   prior = struct('mar', existingLeg("mar"));
   period = struct('start', '2012-01-01', 'end', '2012-01-02');
   skipped = struct('mar', struct('staged', false, 'reason', 'no source'));

   preserved = icemodel.verification.setup.preserveRcmLegs( ...
      prior, skipped, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root);
   testCase.verifyFalse(isfield(preserved, 'mar'));

   degraded = struct('mar', rmfield(prior.mar, 'data_files'));
   preserved = icemodel.verification.setup.preserveRcmLegs( ...
      prior, degraded, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root);
   testCase.verifyFalse(isfield(preserved, 'mar'));

   conflict = skipped;
   conflict.mar.reason = 'artifact metadata does not match requested point';
   retained = icemodel.verification.setup.preserveRcmLegs( ...
      prior, conflict, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root);
   testCase.verifyTrue(isfield(retained, 'mar'));
   testCase.verifyTrue(logical( ...
      retained.mar.replace_prior_artifacts));

   explicit = skipped;
   explicit.mar.replace_prior_artifacts = true;
   retained = icemodel.verification.setup.preserveRcmLegs( ...
      prior, explicit, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root);
   testCase.verifyTrue(isfield(retained, 'mar'));
   testCase.verifyTrue(logical( ...
      retained.mar.replace_prior_artifacts));

   retained = icemodel.verification.setup.preserveRcmLegs( ...
      prior, skipped, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root, overwrite=true);
   testCase.verifyFalse(isfield(retained, 'mar'));

   retained = icemodel.verification.setup.preserveRcmLegs( ...
      prior, skipped, "mar", struct('start', '', 'end', ''), ...
      met_outdir=met_root, userdata_outdir=data_root);
   testCase.verifyTrue(isfield(retained, 'mar'));
   testCase.verifyFalse(isfield( ...
      retained.mar, 'replace_prior_artifacts'));

   retained = icemodel.verification.setup.preserveRcmLegs( ...
      prior, skipped, "mar", period, met_outdir="missing", ...
      userdata_outdir="missing");
   testCase.verifyTrue(isfield(retained, 'mar'));
   testCase.verifyTrue(logical( ...
      retained.mar.replace_prior_artifacts));

   prior.mar.sample_method = 'natural';
   retained = icemodel.verification.setup.preserveRcmLegs( ...
      prior, skipped, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root, method="nearest");
   testCase.verifyTrue(isfield(retained, 'mar'));
   testCase.verifyTrue(logical( ...
      retained.mar.replace_prior_artifacts));

   successful = struct('mar', prior.mar);
   retained = icemodel.verification.setup.preserveRcmLegs( ...
      prior, successful, ["mar", "merra"], period, ...
      met_outdir=met_root, userdata_outdir=data_root);
   testCase.verifyTrue(isfield(retained, 'mar'));
   testCase.verifyTrue(logical( ...
      retained.mar.replace_prior_artifacts));
end

function test_preserveRcmLegs_keeps_exact_successful_replay(testCase)
   % An overwrite=false replay that resolves the same staged leg is a no-op.
   % Reuse-only notes must not churn the manifest, while changed references or
   % coverage remain material patches for the additive merge.
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   prior = struct('mar', existingLeg("mar"));
   period = struct('start', '2012-01-01', 'end', '2012-01-02');
   replay = prior;
   replay.mar.note = ...
      'Existing staged met/Data files fully cover requested window.';

   preserved = icemodel.verification.setup.preserveRcmLegs( ...
      prior, replay, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root, method="nearest", point=[67, -48]);
   testCase.verifyFalse(isfield(preserved, 'mar'));

   % Newly computed readiness must survive an otherwise exact cache replay;
   % once persisted, JSON char/string and empty-array differences remain a no-op.
   upgraded = replay;
   upgraded.mar.forcing_ready = false;
   upgraded.mar.forcing_ready_reason = ...
      'required met placeholder/gap channel(s): albedo';
   upgraded.mar.forcing_complete_windows = struct( ...
      'start_time', "2012-01-01T01:00:00Z", ...
      'end_time', "2012-01-02T00:00:00Z", 'sample_count', 24);
   retained = icemodel.verification.setup.preserveRcmLegs( ...
      prior, upgraded, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root, method="nearest", point=[67, -48]);
   testCase.verifyTrue(isfield(retained, 'mar'));
   testCase.verifyFalse(logical(retained.mar.forcing_ready));

   decoded = jsondecode(jsonencode(upgraded));
   preserved = icemodel.verification.setup.preserveRcmLegs( ...
      decoded, upgraded, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root, method="nearest", point=[67, -48]);
   testCase.verifyFalse(isfield(preserved, 'mar'));

   changed_reference = replay;
   changed_reference.mar.met_files = ...
      "mar3.11/met_case1_mar3.11_20120101_20131231_1hr.mat";
   retained = icemodel.verification.setup.preserveRcmLegs( ...
      prior, changed_reference, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root, method="nearest", point=[67, -48]);
   testCase.verifyTrue(isfield(retained, 'mar'));
   testCase.verifyEqual(retained.mar.met_files, ...
      changed_reference.mar.met_files);

   changed_window = replay;
   changed_window.mar.window.end = '2012-01-03';
   retained = icemodel.verification.setup.preserveRcmLegs( ...
      prior, changed_window, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root, method="nearest", point=[67, -48]);
   testCase.verifyTrue(isfield(retained, 'mar'));
   testCase.verifyEqual(retained.mar.window, changed_window.mar.window);

   changed_product = replay;
   changed_product.mar.source_id = 'mar-other';
   retained = icemodel.verification.setup.preserveRcmLegs( ...
      prior, changed_product, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root, method="nearest", point=[67, -48]);
   testCase.verifyTrue(isfield(retained, 'mar'));
end

function test_preserveRcmLegs_validates_referenced_artifact_metadata(testCase)
   % Concrete saved method/point conflicts invalidate a manifest fallback.
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   prior = struct('mar', existingLeg("mar"));
   period = struct('start', '2012-01-01', 'end', '2012-01-02');
   skipped = struct('mar', struct('staged', false, 'reason', 'no source'));

   matching = icemodel.verification.setup.preserveRcmLegs( ...
      prior, skipped, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root, method="nearest", point=[67, -48]);
   testCase.verifyFalse(isfield(matching, 'mar'));

   wrong_point = icemodel.verification.setup.preserveRcmLegs( ...
      prior, skipped, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root, method="nearest", point=[68, -48]);
   testCase.verifyTrue(isfield(wrong_point, 'mar'));
   testCase.verifyTrue(logical( ...
      wrong_point.mar.replace_prior_artifacts));

   % Clear the manifest-level method so the saved artifact method itself owns
   % this conflict check.
   prior.mar.sample_method = '';
   wrong_method = icemodel.verification.setup.preserveRcmLegs( ...
      prior, skipped, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root, method="natural", point=[67, -48]);
   testCase.verifyTrue(isfield(wrong_method, 'mar'));
   testCase.verifyTrue(logical( ...
      wrong_method.mar.replace_prior_artifacts));

   % A partial coordinate pair is malformed concrete metadata, not an
   % unstamped legacy artifact.
   met_file = fullfile(met_root, prior.mar.met_files);
   saved = load(met_file, 'met', 'artifact_metadata', '-mat');
   saved.artifact_metadata = rmfield(saved.artifact_metadata, 'lon_wgs84');
   met = saved.met;
   artifact_metadata = saved.artifact_metadata;
   save(met_file, 'met', 'artifact_metadata')
   prior.mar.sample_method = 'nearest';
   partial_point = icemodel.verification.setup.preserveRcmLegs( ...
      prior, skipped, "mar", period, met_outdir=met_root, ...
      userdata_outdir=data_root, method="nearest", point=[67, -48]);
   testCase.verifyTrue(isfield(partial_point, 'mar'));
   testCase.verifyTrue(logical( ...
      partial_point.mar.replace_prior_artifacts));
end

function test_readRcmArtifactMetadata_handles_saved_and_legacy_files(testCase)
   % The lightweight reader distinguishes one valid struct from legacy forms.
   missing = fullfile(testCase.TestData.tmp, 'missing-artifact.mat');
   testCase.verifyEmpty(fieldnames( ...
      icemodel.verification.setup.readRcmArtifactMetadata(missing)));

   plain_file = fullfile(testCase.TestData.tmp, 'plain-artifact.mat');
   value = 1;
   save(plain_file, 'value')
   testCase.verifyEmpty(fieldnames( ...
      icemodel.verification.setup.readRcmArtifactMetadata(plain_file)));

   malformed_file = fullfile(testCase.TestData.tmp, ...
      'malformed-artifact.mat');
   artifact_metadata = "bad";
   save(malformed_file, 'artifact_metadata')
   testCase.verifyEmpty(fieldnames( ...
      icemodel.verification.setup.readRcmArtifactMetadata(malformed_file)));

   nonscalar_file = fullfile(testCase.TestData.tmp, ...
      'nonscalar-artifact.mat');
   artifact_metadata = repmat(struct('sample_method', 'nearest'), 1, 2);
   save(nonscalar_file, 'artifact_metadata')
   testCase.verifyEmpty(fieldnames( ...
      icemodel.verification.setup.readRcmArtifactMetadata(nonscalar_file)));

   valid_file = fullfile(testCase.TestData.tmp, 'valid-artifact.mat');
   artifact_metadata = struct('sample_method', 'nearest', ...
      'lat_wgs84', 67, 'lon_wgs84', -48);
   save(valid_file, 'artifact_metadata')
   returned = ...
      icemodel.verification.setup.readRcmArtifactMetadata(valid_file);
   testCase.verifyEqual(returned, artifact_metadata);
end

function test_rcmArtifactOutputDirs_resolves_defaults_and_overrides(testCase)
   % Shared staging/preservation roots match writer defaults and pass overrides.
   [default_met, default_data] = ...
      icemodel.verification.setup.rcmArtifactOutputDirs("", "");
   testCase.verifyEqual(default_met, ...
      string(fullfile(icemodel.getpath('input'), 'met')));
   testCase.verifyEqual(default_data, string(icemodel.getpath('userdata')));

   [explicit_met, explicit_data] = ...
      icemodel.verification.setup.rcmArtifactOutputDirs("met-root", ...
      "data-root");
   testCase.verifyEqual(explicit_met, "met-root");
   testCase.verifyEqual(explicit_data, "data-root");
end

function test_reuseDatasetFamilyCases_preserves_entry_and_resolves_leg(testCase)
   % The forcing-only state preserves its staged entry and narrows only the leg.
   manifest_file = fullfile(testCase.TestData.tmp, 'reuse-manifest.json');
   state = helperState();
   state.colocation.mar = existingLeg("mar");
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);
   prior = jsondecode(fileread(manifest_file));
   coverage.mar = struct('years', 2012, 'year_min', 2012, ...
      'year_max', 2012, 'reason', "");

   [reused, alive, skipped] = ...
      icemodel.verification.setup.reuseDatasetFamilyCases( ...
      manifest_file, "case1", reusePrototype(true), ...
      forcing_sources="mar", coverage=coverage, ...
      startdate="2012-01-01", enddate="2012-01-01 23:00:00");

   testCase.verifyTrue(alive);
   testCase.verifyEmpty(skipped);
   testCase.verifyTrue(reused.reuse_entry);
   testCase.verifyEqual(jsonencode(reused.entry), ...
      jsonencode(prior.cases));
   testCase.verifyEqual(jsonencode(reused.colocation), ...
      jsonencode(prior.cases.colocation));
   testCase.verifyEqual(reused.point, [67, -48]);
   testCase.verifyEqual(string(reused.site_id), "CASE_1");
   testCase.verifyEqual(string(reused.site_name), "Case 1");
   testCase.verifyEqual(reused.site_location, prior.cases.site_location);
   testCase.verifyEqual(string(reused.evaluation_file_rel), ...
      string(prior.cases.evaluation_file));
   testCase.verifyTrue(logical(reused.leg.mar.staged));
   testCase.verifyEqual(reused.leg.mar.start, ...
      icemodel.verification.setup.ensureUtc('2012-01-01'));
   testCase.verifyEqual(reused.leg.mar.end, ...
      icemodel.verification.setup.ensureUtc('2012-01-01 23:00:00'));
end

function test_reuseDatasetFamilyCases_accepts_independent_forcing_window(testCase)
   % SUMup forcing years may be independent without relabeling observations.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'reuse-independent-window-manifest.json');
   state = helperState();
   state.period = struct('start', '', 'end', '');
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);
   coverage.mar = struct('years', 2013, 'year_min', 2013, ...
      'year_max', 2013, 'reason', "");
   t1 = datetime(2013, 1, 1, 'TimeZone', 'UTC');
   t2 = datetime(2013, 12, 31, 23, 0, 0, 'TimeZone', 'UTC');

   reused = icemodel.verification.setup.reuseDatasetFamilyCases( ...
      manifest_file, "case1", reusePrototype(false), ...
      forcing_sources="mar", coverage=coverage, ...
      forcing_startdate=t1, forcing_enddate=t2);

   testCase.verifyEqual(string(reused.period.start), "");
   testCase.verifyEqual(string(reused.period.end), "");
   testCase.verifyTrue(logical(reused.leg.mar.staged));
   testCase.verifyEqual(reused.leg.mar.start, t1);
   testCase.verifyEqual(reused.leg.mar.end, t2);
end

function test_reuseDatasetFamilyCases_keeps_unbounded_period_without_forcing(testCase)
   % A metadata-only fast path retains a prior unbounded period verbatim.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'reuse-unbounded-manifest.json');
   state = helperState();
   state.period = struct('start', '', 'end', '');
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   nat = NaT('TimeZone', 'UTC');
   reused = icemodel.verification.setup.reuseDatasetFamilyCases( ...
      manifest_file, "case1", reusePrototype(false), ...
      startdate="NaT", enddate="NaT", ...
      forcing_startdate=nat, forcing_enddate=nat);

   testCase.verifyEqual(string(reused.period.start), "");
   testCase.verifyEqual(string(reused.period.end), "");
   testCase.verifyEmpty(fieldnames(reused.leg));
end

function test_reuseDatasetFamilyCases_rejects_missing_state(testCase)
   % Fast attachment requires both the family manifest and every requested case.
   missing = fullfile(testCase.TestData.tmp, 'missing-manifest.json');
   testCase.verifyError(@() ...
      icemodel.verification.setup.reuseDatasetFamilyCases( ...
      missing, "case1", reusePrototype(false)), ...
      'icemodel:verification:reuseDatasetFamilyCases:missingManifest');

   manifest_file = fullfile(testCase.TestData.tmp, 'reuse-missing-case.json');
   state = helperState();
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);
   testCase.verifyError(@() ...
      icemodel.verification.setup.reuseDatasetFamilyCases( ...
      manifest_file, "other", reusePrototype(false)), ...
      'icemodel:verification:reuseDatasetFamilyCases:missingCase');

   empty_file = fullfile(testCase.TestData.tmp, 'reuse-empty-family.json');
   empty_manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "helper", "", "test", "v1", string(datetime('today')), struct([]), ...
      struct('site', {}, 'reason', {}));
   icemodel.verification.setup.writeManifest(empty_file, empty_manifest);
   testCase.verifyError(@() ...
      icemodel.verification.setup.reuseDatasetFamilyCases( ...
      empty_file, "case1", reusePrototype(false)), ...
      'icemodel:verification:reuseDatasetFamilyCases:missingCase');

   missing_cases_file = fullfile(testCase.TestData.tmp, ...
      'reuse-no-cases-field.json');
   icemodel.verification.setup.writeManifest( ...
      missing_cases_file, rmfield(empty_manifest, 'cases'));
   testCase.verifyError(@() ...
      icemodel.verification.setup.reuseDatasetFamilyCases( ...
      missing_cases_file, "case1", reusePrototype(false)), ...
      'icemodel:verification:reuseDatasetFamilyCases:missingCase');
end

function test_reuseDatasetFamilyCases_rejects_invalid_windows(testCase)
   % Half, reversed, unbounded, and out-of-period windows fail before mutation.
   missing_manifest = fullfile(testCase.TestData.tmp, ...
      'reuse-invalid-before-missing.json');
   testCase.verifyError(@() ...
      icemodel.verification.setup.reuseDatasetFamilyCases( ...
      missing_manifest, "case1", reusePrototype(false), ...
      startdate="2012-01-01"), ...
      'icemodel:internal:pairedWindow:invalidWindow');
   testCase.verifyFalse(isfile(missing_manifest));

   manifest_file = fullfile(testCase.TestData.tmp, 'reuse-window-errors.json');
   state = helperState();
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);
   testCase.verifyError(@() ...
      icemodel.verification.setup.reuseDatasetFamilyCases( ...
      manifest_file, "case1", reusePrototype(false), ...
      startdate="2012-01-01"), ...
      'icemodel:internal:pairedWindow:invalidWindow');
   testCase.verifyError(@() ...
      icemodel.verification.setup.reuseDatasetFamilyCases( ...
      manifest_file, "case1", reusePrototype(false), ...
      startdate="2012-01-02", enddate="2012-01-01"), ...
      'icemodel:internal:pairedWindow:invalidWindow');
   testCase.verifyError(@() ...
      icemodel.verification.setup.reuseDatasetFamilyCases( ...
      manifest_file, "case1", reusePrototype(false), ...
      startdate="2011-12-31", enddate="2012-01-01"), ...
      'icemodel:verification:reuseDatasetFamilyCases:windowConflict');

   unbounded_file = fullfile(testCase.TestData.tmp, ...
      'reuse-window-unbounded.json');
   state.period = struct('start', '', 'end', '');
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=unbounded_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);
   testCase.verifyError(@() ...
      icemodel.verification.setup.reuseDatasetFamilyCases( ...
      unbounded_file, "case1", reusePrototype(false), ...
      startdate="2012-01-01", enddate="2012-01-02"), ...
      'icemodel:verification:reuseDatasetFamilyCases:unboundedPeriod');
end

function test_buildDatasetFamilyManifest_refreshes_source_lists(testCase)
   % The neutral manifest builder should use the family callback each time, so
   % source lists refresh when the same case gains an RCM colocation leg.
   manifest_file = fullfile(testCase.TestData.tmp, 'manifest.json');
   state = helperState();
   alive = true;

   manifest = icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", skipped=struct('site', "missing", ...
      'reason', "no source"), source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   testCase.verifyEqual(string(manifest.cases.eval_sources), "sumup_obs");
   testCase.verifyEmpty(manifest.cases.forcing_sources);
   testCase.verifyEqual(string(manifest.skipped.site), "missing");

   state.colocation.mar = struct('kind', 'point_met', 'staged', true, ...
      'source', 'mar', 'source_id', 'mar3.11', ...
      'met_files', "mar3.11/met_case1_mar3.11_2012_1hr.mat");
   manifest = icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", skipped=struct('site', {}, 'reason', {}), ...
      source_url="test", source_version="v2", entry_callback=@helperEntry);

   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.forcing_sources)));
   testCase.verifyEqual(string(manifest.skipped.site), "missing");
end

function test_stageDatasetRcmForcing_persists_after_each_source(testCase)
   % Skipped legs avoid archive reads but still exercise one-source-at-a-time
   % delegation and the per-source persist callback.
   state = helperState();
   state.colocation.merra = existingLeg("merra");
   alive = true;
   calls = 0;
   snapshots = cell(1, 2);

   returned = icemodel.verification.setup.stageDatasetRcmForcing( ...
      state, alive, dataset_family="helper", ...
      forcing_sources=["mar", "racmo"], ...
      leg_callback=@skippedLeg, persist_callback=@persistSnapshot);

   testCase.verifyEqual(calls, 2);
   testCase.verifyEqual(numel(snapshots), 2);
   testCase.verifyTrue(isfield(returned.colocation, 'mar'));
   testCase.verifyFalse(logical(returned.colocation.mar.staged));
   testCase.verifyTrue(isfield(returned.colocation, 'racmo'));
   testCase.verifyFalse(logical(returned.colocation.racmo.staged));
   testCase.verifyTrue(isfield(returned.colocation, 'merra'));
   testCase.verifyTrue(logical(returned.colocation.merra.staged));

   function persistSnapshot(st)
      %PERSISTSNAPSHOT Capture each kill-safe per-source state.
      calls = calls + 1;
      snapshots{calls} = st;
   end
end

function test_stageDatasetRcmForcing_runs_after_source_callback_before_persist( ...
      testCase)
   % Optional family products attach to the completed source leg before the
   % kill-safe manifest snapshot observes that state.
   state = helperState();
   callback_sources = strings(0, 1);
   persisted = cell(1, 2);
   persist_count = 0;

   returned = icemodel.verification.setup.stageDatasetRcmForcing( ...
      state, true, dataset_family="helper", ...
      forcing_sources=["mar", "racmo"], leg_callback=@skippedLeg, ...
      after_source_callback=@attachOptional, persist_callback=@capturePersist);

   testCase.verifyEqual(callback_sources, ["mar"; "racmo"])
   testCase.verifyEqual(string(returned.colocation.racmo.optional_product), ...
      "racmo")
   testCase.verifyEqual(string(persisted{1}.colocation.mar.optional_product), ...
      "mar")
   testCase.verifyEqual(string(persisted{2}.colocation.racmo.optional_product), ...
      "racmo")

   function st = attachOptional(st, live, source)
      %ATTACHOPTIONAL Stamp a source-local marker on each live state record.
      callback_sources(end + 1, 1) = source;
      for idx = reshape(live, 1, [])
         st(idx).colocation.(char(source)).optional_product = char(source);
      end
   end

   function capturePersist(st)
      %CAPTUREPERSIST Prove callback changes precede per-source persistence.
      persist_count = persist_count + 1;
      persisted{persist_count} = st;
   end
end

function test_stageDatasetRcmForcing_deduplicates_sources_stably(testCase)
   % Duplicate selectors must not repeat source staging or persistence callbacks.
   state = helperState();
   calls = 0;
   returned = icemodel.verification.setup.stageDatasetRcmForcing( ...
      state, true, dataset_family="helper", ...
      forcing_sources=["mar", "mar", "racmo", "mar"], ...
      leg_callback=@skippedLeg, persist_callback=@countPersist);

   testCase.verifyEqual(calls, 2);
   testCase.verifyTrue(isfield(returned.colocation, 'mar'));
   testCase.verifyTrue(isfield(returned.colocation, 'racmo'));

   function countPersist(~)
      %COUNTPERSIST Count one callback for each unique source in stable order.
      calls = calls + 1;
   end
end

function test_shared_rcm_method_validation_is_side_effect_free(testCase)
   % Invalid interpolation choices fail at both shared boundaries before writes.
   state = helperState();
   met_root = fullfile(testCase.TestData.tmp, 'invalid-method-met');
   testCase.verifyError(@() ...
      icemodel.verification.setup.stageDatasetRcmForcing( ...
      state, true, forcing_sources="mar", leg_callback=@skippedLeg, ...
      met_outdir=met_root, method="cubic"), ...
      'MATLAB:validators:mustBeMember');
   testCase.verifyFalse(isfolder(met_root));

   manifest_root = fullfile(testCase.TestData.tmp, 'invalid-method-runner');
   manifest_file = fullfile(manifest_root, 'manifest.json');
   testCase.verifyError(@() ...
      icemodel.verification.setup.runDatasetFamilyImport( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", entry_callback=@helperEntry, ...
      build_forcing=true, forcing_sources="mar", ...
      leg_callback=@skippedLeg, method="cubic"), ...
      'MATLAB:validators:mustBeMember');
   testCase.verifyFalse(isfolder(manifest_root));

   source_root = fullfile(testCase.TestData.tmp, 'invalid-source-runner');
   testCase.verifyError(@() ...
      icemodel.verification.setup.runDatasetFamilyImport( ...
      state, true, dataset_family="helper", ...
      manifest_file=fullfile(source_root, 'manifest.json'), ...
      requested_ids="case1", entry_callback=@helperEntry, ...
      forcing_sources="not_an_rcm"), ...
      'icemodel:validators:mustBeRcmSourceSelection');
   testCase.verifyFalse(isfolder(source_root));
end

function test_stageDatasetRcmForcing_noop_does_not_require_callback(testCase)
   % Empty live/source selections return unchanged before callback validation.
   state = helperState();

   no_live = icemodel.verification.setup.stageDatasetRcmForcing( ...
      state, false, forcing_sources="mar");
   no_sources = icemodel.verification.setup.stageDatasetRcmForcing( ...
      state, true, forcing_sources=strings(1, 0));

   testCase.verifyEqual(no_live, state);
   testCase.verifyEqual(no_sources, state);
end

function test_stageDatasetRcmForcing_reattaches_cached_requested_source(testCase)
   % A blank optional storage alias preserves the ordinary cache identity.
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   state = helperState();
   state.storage_alias = "";
   alive = true;

   returned = icemodel.verification.setup.stageDatasetRcmForcing( ...
      state, alive, dataset_family="helper", ...
      forcing_sources=["", "mar", ""], ...
      leg_callback=@skippedLeg, met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyTrue(logical(returned.colocation.mar.staged));
   testCase.verifyTrue(isfield(returned.colocation.mar, 'met_files'));
   testCase.verifyEqual(string(returned.colocation.mar.window.start), ...
      "2012-01-01 00:00:00");
   testCase.verifyEqual(string(returned.colocation.mar.window.end), ...
      "2012-01-02 00:00:00");
end

function test_stageDatasetRcmForcing_uses_optional_storage_alias(testCase)
   % A collision-safe cache identity changes only RCM file discovery/naming.
   met_root = fullfile(testCase.TestData.tmp, 'storage-alias-met');
   data_root = fullfile(testCase.TestData.tmp, 'storage-alias-userdata');
   writeTaggedExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_retmip_kanu_mar3.11_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'retmip_kanu_mar3.11_20120101_20121231.mat'), ...
      "nearest", [67, -48]);
   state = helperState();
   state.storage_alias = "retmip_kanu";

   returned = icemodel.verification.setup.stageDatasetRcmForcing( ...
      state, true, dataset_family="retmip", forcing_sources="mar", ...
      leg_callback=@skippedLeg, met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyEqual(string(returned.case_id), "case1");
   testCase.verifyTrue(contains( ...
      string(returned.colocation.mar.data_files), "retmip_kanu_mar3.11"));
   testCase.verifyFalse(contains( ...
      string(returned.colocation.mar.data_files), "/case1_mar3.11"));
end

function test_rcmStorageAlias_scopes_only_proved_collisions(testCase)
   % Exact family/case mappings leave every ordinary cache identity unchanged.
   testCase.verifyEqual( ...
      icemodel.verification.setup.rcmStorageAlias("retmip", "KANU"), ...
      "retmip_kanu");
   testCase.verifyEqual( ...
      icemodel.verification.setup.rcmStorageAlias( ...
      "sumup", "firn aquifer (fa)"), "fa");
   testCase.verifyEqual( ...
      icemodel.verification.setup.rcmStorageAlias("sumup", "kanu"), "kanu");
end

function test_stageDatasetRcmForcing_uses_entry_period_for_cache_discovery(testCase)
   % RetMIP/research_site states can carry discovery period under entry.
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   state = rmfield(helperState(), 'period');
   state.entry = struct('period', struct('start', '2012-01-01 00:00:00', ...
      'end', '2012-01-02 00:00:00'));
   alive = true;

   returned = icemodel.verification.setup.stageDatasetRcmForcing( ...
      state, alive, dataset_family="helper", forcing_sources="mar", ...
      leg_callback=@skippedLeg, met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyTrue(logical(returned.colocation.mar.staged));
   testCase.verifyTrue(isfield(returned.colocation.mar, 'data_files'));
end

function test_stageDatasetRcmForcing_ranks_staged_leg_over_full_case(testCase)
   % Staged family legs retain the full case only as their discovery window.
   met_root = fullfile(testCase.TestData.tmp, 'dataset-discovery-met');
   data_root = fullfile(testCase.TestData.tmp, 'dataset-discovery-data');
   point = [67, -48];
   writeTaggedExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20100101_20121231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20100101_20121231.mat'), "nearest", point);
   writeTaggedExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20141231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20141231.mat'), "nearest", point);
   state = helperState();
   state.period = struct('start', '2010-01-01', 'end', '2012-12-31');

   returned = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.stageDatasetRcmForcing( ...
      state, true, dataset_family="helper", forcing_sources="mar", ...
      leg_callback=@fixedClippedLeg, met_outdir=met_root, ...
      userdata_outdir=data_root, dt_out=""), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   testCase.verifyTrue(contains( ...
      string(returned.colocation.mar.met_files), "20100101_20121231"));
   testCase.verifyTrue(contains( ...
      string(returned.colocation.mar.data_files), "20100101_20121231"));
   testCase.verifyEqual(string(returned.colocation.mar.window.start), ...
      "2012-01-01 00:00:00");
end

function test_runDatasetFamilyImport_preserves_skips_and_final_sources(testCase)
   % The runner writes the native manifest, stages skipped RCM legs through the
   % shared path, and returns the final manifest with refreshed source lists.
   manifest_file = fullfile(testCase.TestData.tmp, 'runner-manifest.json');
   state = helperState();
   alive = true;
   skipped = struct('site', "case2", 'reason', "missing observation");

   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids=["case1", "case2"], skipped=skipped, ...
      source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, build_forcing=true, ...
      forcing_sources=["mar", "merra"], leg_callback=@skippedLeg);

   testCase.verifyEqual(string(manifest.skipped.site), "case2");
   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyTrue(isfield(manifest.cases.colocation, 'merra'));
   testCase.verifyFalse(logical(manifest.cases.colocation.mar.staged));
   testCase.verifyEmpty(manifest.cases.forcing_sources);
end

function test_listcases_handles_skip_only_manifest(testCase)
   % A skip-only import writes cases=[] plus skipped records; readers should treat
   % it as an empty case family instead of indexing into the empty array.
   eval_root = fullfile(testCase.TestData.tmp, 'eval-root');
   manifest_file = fullfile(eval_root, 'promice', 'manifest.json');
   mkdir(fileparts(manifest_file))
   state = helperState();
   alive = false;
   skipped = struct('site', "case1", 'reason', "missing source");

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="promice", manifest_file=manifest_file, ...
      requested_ids="case1", skipped=skipped, source_url="test", ...
      source_version="v1", entry_callback=@helperEntry, overwrite_family=true);

   cases = icemodel.verification.listcases( ...
      evaluation_data_root=eval_root, dataset_family="promice");

   testCase.verifyEmpty(cases);
end

function test_stageRcmForcing_struct_manifest_preserves_skipped(testCase)
   % RCM-only manifest rewrites must not erase missing-observation records.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'struct-manifest-stage-rcm.json');
   state = helperState();
   alive = true;
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids=["case1", "case2"], ...
      skipped=struct('site', "case2", 'reason', "missing observation"), ...
      source_url="test", source_version="v1", entry_callback=@helperEntry, ...
      overwrite_family=true);
   manifest = jsondecode(fileread(manifest_file));
   missing_source = fullfile(testCase.TestData.tmp, 'missing-rcm-source');

   manifest = icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=manifest, manifest_file=manifest_file, ...
      forcing_sources="mar", ...
      met_outdir=fullfile(testCase.TestData.tmp, 'met'), ...
      userdata_outdir=fullfile(testCase.TestData.tmp, 'userdata'), ...
      mar_dir=missing_source, overwrite_family=true);

   testCase.verifyEqual(string(manifest.skipped.site), "case2");
   testCase.verifyEqual(string(manifest.skipped.reason), ...
      "missing observation");
end

function test_stageRcmForcing_manifest_does_not_sanitize_omitted_rcm(testCase)
   % A manifest patch leaves an unrequested leg untouched, even after point edit.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'manifest-preserve-omitted-rcm.json');
   met_root = fullfile(testCase.TestData.tmp, 'omitted-compatible-met');
   data_root = fullfile(testCase.TestData.tmp, 'omitted-compatible-userdata');
   writeTaggedExistingFiles(fullfile(met_root, 'merra2', ...
      'met_case1_merra2_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'merra2', ...
      'case1_merra2_20120101_20121231.mat'), "nearest", [67, -48]);
   state = helperState();
   state.colocation.merra = existingLeg("merra");
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", entry_callback=@helperEntry, ...
      overwrite_family=true);
   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.site_location.lat_wgs84 = 68.0;

   returned = icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=manifest, manifest_file=manifest_file, ...
      forcing_sources=strings(1, 0), met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyTrue(isfield(returned.cases.colocation, 'merra'));
   testCase.verifyTrue(logical(returned.cases.colocation.merra.staged));
   testCase.verifyTrue(ismember("merra2", ...
      string(returned.cases.forcing_sources)));
end

function test_resolveStagingRoots_pairs_explicit_eval_root(testCase)
   % An explicit eval root implies the sibling input root unless overridden.
   eval_root = fullfile(testCase.TestData.tmp, 'custom', 'eval');

   [returned_eval, returned_input] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      evaluation_data_root=eval_root);

   testCase.verifyEqual(string(returned_eval), string(eval_root));
   testCase.verifyEqual(string(returned_input), ...
      string(fullfile(testCase.TestData.tmp, 'custom', 'input')));
end

function test_resolveStagingRoots_explicit_input_overrides_pair(testCase)
   % A caller-provided input root remains authoritative.
   eval_root = fullfile(testCase.TestData.tmp, 'custom', 'eval');
   input_root = fullfile(testCase.TestData.tmp, 'other-input');

   [~, returned_input] = icemodel.verification.setup.resolveStagingRoots( ...
      evaluation_data_root=eval_root, input_data_root=input_root);

   testCase.verifyEqual(string(returned_input), string(input_root));
end

function test_runDatasetFamilyImport_build_forcing_preserves_omitted_rcm(testCase)
   % A partial RCM refresh preserves same-location unrequested RCM sources.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-preserve-manifest.json');
   met_root = fullfile(testCase.TestData.tmp, 'preserve-omitted-met');
   data_root = fullfile(testCase.TestData.tmp, 'preserve-omitted-userdata');
   writeTaggedExistingFiles(fullfile(met_root, 'merra2', ...
      'met_case1_merra2_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'merra2', ...
      'case1_merra2_20120101_20121231.mat'), "nearest", [67, -48]);
   writeTaggedExistingData(fullfile(data_root, 'racmo2.3p3', ...
      'case1_racmo2.3p3_20120101_20121231.mat'), "nearest", [67, -48]);
   state = helperState();
   alive = true;
   state.colocation.mar = existingLeg("mar");
   state.colocation.merra = existingLeg("merra");
   state.colocation.racmo = existingLeg("racmo");

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);
   prior = jsondecode(fileread(manifest_file));
   prior_merra = jsonencode(prior.cases.colocation.merra);
   prior_racmo = jsonencode(prior.cases.colocation.racmo);

   state = helperState();
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, build_forcing=true, ...
      forcing_sources="mar", leg_callback=@skippedLeg, ...
      met_outdir=met_root, userdata_outdir=data_root);

   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyFalse(logical(manifest.cases.colocation.mar.staged));
   testCase.verifyFalse(isfield( ...
      manifest.cases.colocation.mar, 'replace_prior_artifacts'));
   testCase.verifyTrue(isfield(manifest.cases.colocation, 'merra'));
   testCase.verifyTrue(isfield(manifest.cases.colocation, 'racmo'));
   testCase.verifyEqual(jsonencode(manifest.cases.colocation.merra), ...
      prior_merra);
   testCase.verifyEqual(jsonencode(manifest.cases.colocation.racmo), ...
      prior_racmo);
   testCase.verifyTrue(any(string(manifest.cases.forcing_sources) ...
      == "merra2"));
   testCase.verifyTrue(any(string(manifest.cases.eval_sources) ...
      == "racmo2.3p3"));
end

function test_runDatasetFamilyImport_overwrite_preserves_prior_when_raw_unavailable(testCase)
   % overwrite=true keeps a compatible prior leg when no raw leg can rebuild it.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-overwrite-requested-rcm-manifest.json');
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   state = helperState();
   alive = true;
   state.colocation.mar = existingLeg("mar");

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   state = helperState();
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, build_forcing=true, ...
      forcing_sources="mar", leg_callback=@skippedLeg, overwrite=true, ...
      met_outdir=met_root, userdata_outdir=data_root);

   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyTrue(logical(manifest.cases.colocation.mar.staged));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.forcing_sources)));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.eval_sources)));
end

function test_runDatasetFamilyImport_persists_reused_rcm_readiness(testCase)
   % A canonical exact-cache replay must upgrade a legacy leg with readiness,
   % persist it through the family merge, then converge byte-for-byte.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-readiness-upgrade-manifest.json');
   met_root = fullfile(testCase.TestData.tmp, 'readiness-upgrade-met');
   data_root = fullfile(testCase.TestData.tmp, 'readiness-upgrade-userdata');
   writeTaggedExistingFiles(fullfile(met_root, 'merra2', ...
      'met_case1_merra2_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'merra2', ...
      'case1_merra2_20120101_20121231.mat'), ...
      "nearest", [67, -48], true);
   state = helperState();
   state.colocation.merra = existingLeg("merra");

   % Seed the exact pre-fix manifest shape: staged files without readiness.
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);
   state = helperState();
   run = @() icemodel.verification.setup.runDatasetFamilyImport( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, build_forcing=true, ...
      forcing_sources="merra", leg_callback=@exactLeg, ...
      met_outdir=met_root, userdata_outdir=data_root, dt_out="");
   manifest = icemodel.test.helpers.captureExpectedWarning(testCase, run, ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   leg = manifest.cases.colocation.merra;
   testCase.verifyTrue(isfield(leg, 'forcing_ready'));
   testCase.verifyFalse(logical(leg.forcing_ready));
   testCase.verifyTrue(contains(string(leg.forcing_ready_reason), "albedo"));
   testCase.verifyEqual(double(leg.forcing_complete_windows.sample_count), 1);
   persisted = jsondecode(fileread(manifest_file));
   testCase.verifyTrue(isfield( ...
      persisted.cases.colocation.merra, 'forcing_ready'));
   before = fileBytes(manifest_file);

   % The same public replay now observes semantically identical diagnostics
   % across JSON types and must leave the durable manifest byte-identical.
   manifest = icemodel.test.helpers.captureExpectedWarning(testCase, run, ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');
   testCase.verifyFalse(logical( ...
      manifest.cases.colocation.merra.forcing_ready));
   testCase.verifyEqual(fileBytes(manifest_file), before);

   function L = exactLeg(~, ~)
      %EXACTLEG Match the legacy leg window so readiness is the only patch.
      L = struct('staged', true, 'years', 2012, ...
         'start', datetime(2012, 1, 1, 'TimeZone', 'UTC'), ...
         'end', datetime(2012, 1, 2, 'TimeZone', 'UTC'), 'reason', "");
   end
end

function test_runDatasetFamilyImport_location_patch_preserves_omitted_rcm(testCase)
   % Editing a case point does not sanitize RCM sources omitted from this patch.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-location-stale-manifest.json');
   state = helperState();
   alive = true;
   state.colocation.mar = existingLeg("mar");
   state.colocation.merra = existingLeg("merra");

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);
   prior = jsondecode(fileread(manifest_file));
   prior_mar = jsonencode(prior.cases.colocation.mar);
   prior_merra = jsonencode(prior.cases.colocation.merra);

   state = helperState();
   state.point = [68.0, -49.0];
   state.site_location.lat_wgs84 = 68.0;
   state.site_location.lon_wgs84 = -49.0;
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry);

   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyTrue(isfield(manifest.cases.colocation, 'merra'));
   testCase.verifyEqual(jsonencode(manifest.cases.colocation.mar), prior_mar);
   testCase.verifyEqual(jsonencode(manifest.cases.colocation.merra), ...
      prior_merra);
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.forcing_sources)));
end

function test_runDatasetFamilyImport_preserves_omitted_skipped_rcm(testCase)
   % An omitted source leg is additive state even when its prior build skipped.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-drop-skipped-manifest.json');
   state = helperState();
   alive = true;
   state.colocation.mar = struct('kind', 'point_met', ...
      'staged', false, 'reason', 'old MAR failure');

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   state = helperState();
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry);

   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyFalse(logical(manifest.cases.colocation.mar.staged));
   testCase.verifyFalse(isfield( ...
      manifest.cases.colocation.mar, 'replace_prior_artifacts'));
   testCase.verifyEqual(string(manifest.cases.colocation.mar.reason), ...
      "old MAR failure");
end

function test_runDatasetFamilyImport_partial_refresh_preserves_omitted_rcm(testCase)
   % A requested MAR refresh cannot sanitize an omitted existing MERRA leg.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-omitted-stale-manifest.json');
   state = helperState();
   alive = true;
   state.colocation.mar = existingLeg("mar");
   state.colocation.merra = existingLeg("merra");

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);
   prior = jsondecode(fileread(manifest_file));
   prior_merra = jsonencode(prior.cases.colocation.merra);

   state = helperState();
   state.period = struct('start', '2013-01-01', 'end', '2013-01-02');
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, build_forcing=true, ...
      forcing_sources="mar", leg_callback=@skippedLeg);

   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyFalse(logical(manifest.cases.colocation.mar.staged));
   testCase.verifyTrue(isfield(manifest.cases.colocation, 'merra'));
   testCase.verifyEqual(jsonencode(manifest.cases.colocation.merra), ...
      prior_merra);
   testCase.verifyTrue(ismember("merra2", ...
      string(manifest.cases.forcing_sources)));
end

function test_runDatasetFamilyImport_period_patch_preserves_omitted_rcm(testCase)
   % A new entry period does not delete a source omitted from the current patch.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-entry-period-manifest.json');
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   state = helperState();
   alive = true;
   state.colocation.mar = existingLeg("mar");

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);
   prior = jsondecode(fileread(manifest_file));
   prior_mar = jsonencode(prior.cases.colocation.mar);

   state = helperState();
   state.entry = struct('period', struct('start', '2013-01-01', ...
      'end', '2013-01-02'));
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyEqual(jsonencode(manifest.cases.colocation.mar), prior_mar);
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.forcing_sources)));
end

function test_runDatasetFamilyImport_observation_only_preserves_rcm(testCase)
   % Observation-only refreshes keep same-location RCM legs in the manifest.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-observation-only-manifest.json');
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   state = helperState();
   alive = true;
   state.colocation.mar = existingLeg("mar");

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   state = helperState();
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, build_forcing=true, ...
      forcing_sources="mar", leg_callback=@skippedLeg, ...
      met_outdir=met_root, userdata_outdir=data_root);

   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyTrue(logical(manifest.cases.colocation.mar.staged));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.forcing_sources)));
end

function test_runDatasetFamilyImport_observation_only_does_not_attach_cache(testCase)
   % build_forcing=false must not discover or attach unrequested cache files.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-observation-cache-manifest.json');
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   state = helperState();
   alive = true;

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   state = helperState();
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyFalse(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyEmpty(manifest.cases.forcing_sources);
end

function test_runDatasetFamilyImport_overwrite_family_keeps_cache_reuse(testCase)
   % overwrite_family rewrites the manifest, not existing expensive RCM files.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-overwrite-family-cache-manifest.json');
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   state = helperState();
   alive = true;

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   state = helperState();
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, build_forcing=true, ...
      forcing_sources="mar", leg_callback=@skippedLeg, ...
      met_outdir=met_root, userdata_outdir=data_root, ...
      overwrite_family=true);

   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyTrue(logical(manifest.cases.colocation.mar.staged));
   testCase.verifyTrue(isfield(manifest.cases.colocation.mar, 'data_files'));
   testCase.verifyTrue(isfield(manifest.cases.colocation.mar, 'met_files'));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.forcing_sources)));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.eval_sources)));
end

function test_runDatasetFamilyImport_overwrite_preserves_cache_when_raw_unavailable(testCase)
   % A compatible discovered cache survives when overwrite has no raw leg.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-overwrite-cache-manifest.json');
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   state = helperState();
   alive = true;

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   state = helperState();
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, build_forcing=true, ...
      forcing_sources="mar", leg_callback=@skippedLeg, ...
      met_outdir=met_root, userdata_outdir=data_root, overwrite=true);

   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyTrue(logical(manifest.cases.colocation.mar.staged));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.forcing_sources)));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.eval_sources)));
end

function test_runDatasetFamilyImport_observation_only_preserves_clipped_rcm(testCase)
   % Observation-only refreshes preserve source-clipped artifacts for the site.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-observation-stale-manifest.json');
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   state = helperState();
   alive = true;
   state.colocation.mar = existingLeg("mar");

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   state = helperState();
   state.period = struct('start', '2011-01-01', 'end', '2013-01-02');
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyTrue(logical(manifest.cases.colocation.mar.staged));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.forcing_sources)));
end

function test_runDatasetFamilyImport_blank_period_does_not_rediscover_rcm(testCase)
   % A blank refreshed period cannot define a safe cached-artifact overlap.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-observation-blank-manifest.json');
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   state = helperState();
   alive = true;

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   state = helperState();
   state.period = struct('start', '', 'end', '');
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyFalse(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyFalse(ismember("mar3.11", ...
      string(manifest.cases.forcing_sources)));
end

function test_runDatasetFamilyImport_preserves_unrequested_clipped_rcm(testCase)
   % A metadata-only patch preserves an omitted RCM leg without widening it.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-observation-clipped-manifest.json');
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   state = helperState();
   alive = true;
   state.colocation.mar = existingLeg("mar");

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   state = helperState();
   state.period = struct('start', '2010-01-01', 'end', '2020-01-01');
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyTrue(logical(manifest.cases.colocation.mar.staged));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.forcing_sources)));
   testCase.verifyEqual(string(manifest.cases.colocation.mar.window.start), ...
      "2012-01-01");
   testCase.verifyEqual(string(manifest.cases.colocation.mar.window.end), ...
      "2012-01-02");
end

function test_runDatasetFamilyImport_preserves_prior_unstamped_rcm(testCase)
   % Prior manifest refs can carry compatible unstamped expensive artifacts.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-unstamped-prior-manifest.json');
   met_root = fullfile(testCase.TestData.tmp, 'unstamped-met');
   data_root = fullfile(testCase.TestData.tmp, 'unstamped-userdata');
   writeBareExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'));
   state = helperState();
   alive = true;
   state.colocation.mar = existingLeg("mar");

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   state = helperState();
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyTrue(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyTrue(logical(manifest.cases.colocation.mar.staged));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.forcing_sources)));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.eval_sources)));
end

function test_runDatasetFamilyImport_preserves_prior_data_on_partial_cache(testCase)
   % Met-only cache rediscovery must not drop prior manifest Data refs.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-partial-cache-prior-manifest.json');
   met_root = fullfile(testCase.TestData.tmp, 'runner-partial-cache-met');
   data_root = fullfile(testCase.TestData.tmp, 'runner-partial-cache-userdata');
   writeTaggedExistingMet(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat'), ...
      "nearest", [67, -48]);
   writeBareExistingData(fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'));
   state = helperState();
   alive = true;
   state.colocation.mar = existingLeg("mar");

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   state = helperState();
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyTrue(isfield(manifest.cases.colocation.mar, 'met_files'));
   testCase.verifyTrue(isfield(manifest.cases.colocation.mar, 'data_files'));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(manifest.cases.eval_sources)));
end

function test_runDatasetFamilyImport_does_not_attach_cache_without_request(testCase)
   % A fresh native manifest ignores unrelated compatible input artifacts.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-new-manifest-cache.json');
   [met_root, data_root] = seedExistingRcmArtifacts(testCase);
   state = helperState();
   alive = true;

   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyFalse(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyFalse(ismember("mar3.11", ...
      string(manifest.cases.forcing_sources)));
   testCase.verifyFalse(ismember("mar3.11", ...
      string(manifest.cases.eval_sources)));
end

function test_runDatasetFamilyImport_does_not_derive_met_when_disabled(testCase)
   % build_forcing=false leaves cached Data untouched and writes no met file.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-data-only-cache.json');
   met_root = fullfile(testCase.TestData.tmp, 'data-only-met');
   data_root = fullfile(testCase.TestData.tmp, 'data-only-userdata');
   writeTaggedMetContractData(fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'), "nearest", [67, -48]);
   state = helperState();
   alive = true;

   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, met_outdir=met_root, ...
      userdata_outdir=data_root);

   testCase.verifyFalse(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyFalse(isfolder(met_root));
   testCase.verifyEmpty(manifest.cases.forcing_sources);
end

function test_runDatasetFamilyImport_overwrite_discards_existing_rcm(testCase)
   % A full family overwrite starts from the new native state and does not
   % carry stale RCM legs from the previous manifest.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'runner-overwrite-manifest.json');
   state = helperState();
   alive = true;
   state.colocation.mar = struct('kind', 'point_met', 'staged', true, ...
      'source', 'mar', 'source_id', 'mar3.11', ...
      'sample_method', 'nearest', ...
      'met_files', "mar3.11/met_case1_mar3.11_2012_1hr.mat");

   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   state = helperState();
   manifest = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v2", ...
      entry_callback=@helperEntry, build_forcing=true, ...
      forcing_sources=strings(1, 0), overwrite_family=true);

   testCase.verifyFalse(isfield(manifest.cases.colocation, 'mar'));
   testCase.verifyEmpty(manifest.cases.forcing_sources);
end

function test_stampArtifactMetadata_preserves_source_units(testCase)
   % Source-specific observational metadata should survive recursive canonical
   % fill, which labels blank fields but does not overwrite curated fields.
   smb = table(0.2, "SUMup", 'VariableNames', {'smb', 'source'});
   smb.Properties.VariableUnits = {'m w.e.', ''};
   smb.Properties.VariableDescriptions = {'interval water equivalent', ''};
   smb = addprop(smb, "StandardNames", "table");
   smb.Properties.CustomProperties.StandardNames = ["", ""];

   artifact = struct('targets', struct('smb', smb));
   artifact = icemodel.verification.setup.stampArtifactMetadata(artifact);
   returned = artifact.targets.smb;

   testCase.verifyEqual(string(returned.Properties.VariableUnits( ...
      string(returned.Properties.VariableNames) == "smb")), "m w.e.");
   testCase.verifyEqual(string(returned.Properties.VariableDescriptions( ...
      string(returned.Properties.VariableNames) == "smb")), ...
      "interval water equivalent");
   testCase.verifyEqual(returned.Properties.CustomProperties.StandardNames( ...
      string(returned.Properties.VariableNames) == "smb"), "");
   testCase.verifyEqual(string(returned.Properties.VariableUnits( ...
      string(returned.Properties.VariableNames) == "source")), "");
end

function test_stampArtifactMetadata_labels_interval_smb_as_accumulated(testCase)
   % Blank interval-SMB metadata should not inherit the model rate-unit meaning.
   t1 = datetime(2012, 1, 1, 'TimeZone', 'UTC');
   t2 = datetime(2012, 1, 31, 'TimeZone', 'UTC');
   smb = table(t1, t2, 0.2, ...
      'VariableNames', {'start_date', 'end_date', 'smb'});

   artifact = struct('targets', struct('smb', smb));
   artifact = icemodel.verification.setup.stampArtifactMetadata(artifact);
   returned = artifact.targets.smb;
   names = string(returned.Properties.VariableNames);

   testCase.verifyEqual(string(returned.Properties.VariableUnits( ...
      names == "smb")), "m w.e.");
   testCase.verifyTrue(contains(string(returned.Properties.VariableDescriptions( ...
      names == "smb")), "accumulated"));
   testCase.verifyEqual(returned.Properties.CustomProperties.StandardNames( ...
      names == "smb"), "");
end

function test_repairRcmArtifactMetadata_reports_ambiguous_alias(testCase)
   % Alias fallback must not guess when current manifests disagree on location.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-ambiguous-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-ambiguous-input');
   writeTinyRepairManifest(eval_root, "promice", "kanu", 67.0, -47.0, "");
   writeTinyRepairManifest(eval_root, "retmip", "kanu", 67.1, -47.0, "");
   artifact = fullfile(input_root, 'userdata', 'mar3.11', ...
      'kanu_mar3.11_2012.mat');
   ensureParent(artifact);
   touch(artifact);

   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root);
   statuses = string({report.records.status});

   testCase.verifyTrue(any(statuses == "ambiguous"));
end

function test_repairRcmArtifactMetadata_accepts_yearly_referenced_artifact(testCase)
   % Current demo fixtures can be yearly userdata files rather than windows.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-yearly-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-yearly-input');
   rel = "mar3.11/case1_mar3.11_2012.mat";
   writeTinyRepairManifest(eval_root, "promice", "case1", 67.0, -48.0, rel);
   writeTaggedExistingData(fullfile(input_root, 'userdata', rel), ...
      "nearest", [67, -48]);

   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root);
   statuses = string({report.records.status});

   testCase.verifyTrue(any(statuses == "would_repair"));
end

function test_repairRcmArtifactMetadata_matches_tmp_symlink_scoped_keys(testCase)
   % A /tmp-spelled family scope must match canonical /private/tmp dir results
   % exactly without falling back to another family's artifact alias.
   [~, token] = fileparts(tempname);
   root_name = "icemodel-repair-path-" + string(token);
   alias_root = fullfile("/tmp", root_name);
   canonical_root = fullfile("/private/tmp", root_name);
   mkdir(alias_root)
   root_cleanup = onCleanup(@() rmdir(alias_root, 's'));
   testCase.assumeTrue(isfolder(canonical_root));

   % Keep the dry-run proof independent of any workstation MAR source setting.
   prior_mar_dir = getenv('ICEMODEL_MAR_DIR');
   env_cleanup = onCleanup(@() setenv('ICEMODEL_MAR_DIR', prior_mar_dir));
   setenv('ICEMODEL_MAR_DIR', '');

   % The selected manifest spelling contains a missing component plus '..'; the
   % actual file and the out-of-scope RetMIP file both live in the canonical tree.
   alias_eval = fullfile(alias_root, 'eval');
   alias_input = fullfile(alias_root, 'input');
   canonical_eval = fullfile(canonical_root, 'eval');
   canonical_input = fullfile(canonical_root, 'input');
   actual_rel = "mar3.11/case1_mar3.11_2012.mat";
   manifest_rel = "mar3.11/not-created/../case1_mar3.11_2012.mat";
   outside_rel = "mar3.11/case2_mar3.11_2012.mat";
   writeTinyRepairManifest(alias_eval, "promice", "case1", ...
      67.0, -48.0, manifest_rel);
   writeTinyRepairManifest(alias_eval, "retmip", "case2", ...
      68.0, -48.0, outside_rel);
   selected_file = fullfile(alias_input, 'userdata', actual_rel);
   outside_file = fullfile(alias_input, 'userdata', outside_rel);
   writeBareExistingData(selected_file);
   writeBareExistingData(outside_file);
   selected_before = fileBytes(selected_file);
   outside_before = fileBytes(outside_file);

   % Alias and canonical roots must classify the identical exact PROMICE record
   % while remaining source-light and leaving the lexical component nonexistent.
   alias_report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      alias_input, eval_root=alias_eval, dataset_family="promice");
   canonical_report = ...
      icemodel.verification.setup.repairRcmArtifactMetadata( ...
      canonical_input, eval_root=canonical_eval, dataset_family="promice");
   testCase.verifyEqual(numel(alias_report.records), 1);
   testCase.verifyEqual(numel(canonical_report.records), 1);
   testCase.verifyEqual(string({alias_report.records.filename}), ...
      string({canonical_report.records.filename}));
   testCase.verifyEqual(string({alias_report.records.status}), ...
      string({canonical_report.records.status}));
   testCase.verifyEqual(string(alias_report.records.actions), ...
      string(canonical_report.records.actions));
   testCase.verifyEqual(string(alias_report.records.status), "would_repair");
   testCase.verifyEqual(string(alias_report.records.filename), ...
      string(fullfile(canonical_input, 'userdata', actual_rel)));

   % Family scoping and dry-run source-light behavior are part of the regression,
   % not incidental consequences of matching one path spelling.
   for report = {alias_report, canonical_report}
      testCase.verifyEqual(report{1}.mar_source_reads, 0);
      testCase.verifyEqual(report{1}.merra_source_reads, 0);
      testCase.verifyEqual(report{1}.modis_source_reads, 0);
   end
   testCase.verifyEqual(fileBytes(selected_file), selected_before);
   testCase.verifyEqual(fileBytes(outside_file), outside_before);
   testCase.verifyFalse(isfolder(fullfile(alias_input, 'userdata', ...
      'mar3.11', 'not-created')));
   clear env_cleanup root_cleanup
end

function test_repairRcmArtifactMetadata_reports_exact_file_conflict(testCase)
   % Exact manifest refs must not be overwritten when case locations disagree.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-conflict-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-conflict-input');
   rel = "mar3.11/shared_mar3.11_2012.mat";
   writeTinyRepairManifest(eval_root, "promice", "case1", 67.0, -48.0, rel);
   writeTinyRepairManifest(eval_root, "retmip", "case2", 68.0, -48.0, rel);
   writeTaggedExistingData(fullfile(input_root, 'userdata', rel), ...
      "nearest", [67, -48]);

   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root);
   statuses = string({report.records.status});

   testCase.verifyTrue(any(statuses == "ambiguous"));
end

function test_repairRcmArtifactMetadata_uses_manifest_sample_method(testCase)
   % A manifest exact ref decides the missing sample_method during repair.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-method-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-method-input');
   rel = "mar3.11/case1_mar3.11_2012.mat";
   data_file = fullfile(input_root, 'userdata', rel);
   writeTinyRepairManifest(eval_root, "promice", "case1", 67.0, -48.0, ...
      rel, "natural");
   writeBareExistingData(data_file);

   icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dry_run=false);
   loaded = load(data_file, 'artifact_metadata');

   testCase.verifyEqual(string(loaded.artifact_metadata.sample_method), ...
      "natural");
end

function test_repairRcmArtifactMetadata_accepts_15m_met_artifact(testCase)
   % Fifteen-minute RCM met files are current-cache artifacts too.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-15m-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-15m-input');
   rel = "mar3.11/met_case1_mar3.11_20120101_20121231_15m.mat";
   met_file = fullfile(input_root, 'met', rel);
   writeTinyRepairManifest(eval_root, "promice", "case1", 67.0, -48.0, ...
      rel, "nearest");
   writeTaggedExistingMet(met_file, "nearest", [67, -48]);

   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root);
   statuses = string({report.records.status});

   testCase.verifyTrue(any(statuses == "would_repair"));
end

function test_repairRcmArtifactMetadata_applies_mar_daily_idempotently(testCase)
   % Native daily sector-1 RU/SMB repair applies only to hourly source Data.
   % The paired 15-minute met payload must remain numerically byte-for-byte
   % equivalent while generic metadata is synchronized, and pass two reads no RU/SMB.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-mar-daily-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-mar-daily-input');
   mar_dir = fullfile(testCase.TestData.tmp, 'repair-mar-daily-source');
   met_rel = "mar3.11/met_case1_mar3.11_20120101_20120102_15m.mat";
   data_rel = "mar3.11/case1_mar3.11_20120101_20120102.mat";
   met_file = fullfile(input_root, 'met', met_rel);
   data_file = fullfile(input_root, 'userdata', data_rel);
   writeTinyRepairPairManifest(eval_root, "case1", 67.0, -48.0, ...
      met_rel, data_rel);
   writeLegacyMarPair(met_file, data_file);
   met_before = makeLegacyMarMet15m(met_file);
   writeTinyMarDailySource(mar_dir, 2012);

   dry = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      mar_dir=string(mar_dir));
   testCase.verifyTrue(all(string({dry.records.status}) == "would_repair"));
   is_met = contains(string({dry.records.filename}), filesep + "met" + filesep);
   testCase.verifyEqual(nnz(is_met), 1);
   testCase.verifyFalse(any(string(dry.records(is_met).actions) == ...
      "apply_mar_daily_constrained_qc"));
   testCase.verifyTrue(any(string(dry.records(~is_met).actions) == ...
      "apply_mar_daily_constrained_qc"));
   testCase.verifyEqual(dry.records(is_met).mar_qc_status, "not_applicable");
   % RU and SMB are independent native variables, hence two source reads for
   % the one Data payload; the met payload performs no additional native read.
   testCase.verifyEqual(dry.mar_source_reads, 2);
   testCase.verifyEqual(dry.mar_cached_series, 1);

   written = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      mar_dir=string(mar_dir), dry_run=false);
   testCase.verifyTrue(all(string({written.records.status}) == "repaired"));
   loaded_met = load(met_file, 'met', 'artifact_metadata');
   loaded_data = load(data_file, 'Data', 'artifact_metadata');
   expected_runoff = [ones(24, 1) * 0.001; ones(24, 1) * 0.002];
   expected_smb = -expected_runoff;
   testCase.verifyEqual(loaded_met.met.Time, met_before.Time);
   testCase.verifyEqual(loaded_met.met.Properties.VariableNames, ...
      met_before.Properties.VariableNames);
   testCase.verifyEqual(loaded_met.met.Variables, met_before.Variables, ...
      'AbsTol', 0);
   testCase.verifyEqual(loaded_data.Data.runoff, expected_runoff, ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(loaded_data.Data.smb, expected_smb, ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(loaded_data.Data.melt, (1:48)');
   testCase.verifyEqual(loaded_met.met.Properties.UserData, ...
      loaded_met.artifact_metadata);
   testCase.verifyEqual(loaded_data.Data.Properties.UserData, ...
      loaded_data.artifact_metadata);
   testCase.verifyEqual(string( ...
      loaded_data.artifact_metadata.mar_qc_method), ...
      "daily_constrained_hourly");
   testCase.verifyEqual(loaded_data.artifact_metadata.mar_qc_sector, 1);

   second = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      mar_dir=string(mar_dir));
   testCase.verifyTrue(all(string({second.records.status}) == "unchanged"));
   testCase.verifyEqual(second.mar_source_reads, 0);
end

function test_repairRcmArtifactMetadata_aligns_overlong_mar_ledgers_source_light( ...
      testCase)
   % A clipped hourly Data/15-minute met pair with full-calendar provenance is
   % repaired metadata-only, preserves all payload values, and is idempotent.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-mar-align-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-mar-align-input');
   met_rel = "mar3.11/met_case1_mar3.11_20120110_20120112_15m.mat";
   data_rel = "mar3.11/case1_mar3.11_20120110_20120112.mat";
   met_file = fullfile(input_root, 'met', met_rel);
   data_file = fullfile(input_root, 'userdata', data_rel);
   writeTinyRepairPairManifest(eval_root, "case1", 67.0, -48.0, ...
      met_rel, data_rel);
   [data_before, met_before, sentinel] = ...
      writeOverlongMarLedgerPair(met_file, data_file);

   dry = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice");
   testCase.verifyTrue(all(string({dry.records.status}) == "would_repair"));
   testCase.verifyTrue(all(arrayfun(@(record) ...
      any(string(record.actions) == "align_mar_daily_metadata"), ...
      dry.records)));
   testCase.verifyFalse(any(arrayfun(@(record) ...
      any(string(record.actions) == "apply_mar_daily_constrained_qc"), ...
      dry.records)));
   testCase.verifyEqual(dry.mar_source_reads, 0);

   written = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      dry_run=false);
   testCase.verifyTrue(all(string({written.records.status}) == "repaired"));
   data = load(data_file, 'Data', 'artifact_metadata', 'sentinel_payload');
   derived = load(met_file, 'met', 'artifact_metadata', 'sentinel_payload');

   % Numeric values, axes, variable classes/names, and unrelated MAT payloads
   % are invariant; only timetable/top-level metadata is synchronized.
   testCase.verifyEqual(data.Data.Time, data_before.Time);
   testCase.verifyEqual(data.Data.Variables, data_before.Variables, 'AbsTol', 0);
   testCase.verifyEqual(data.Data.Properties.VariableNames, ...
      data_before.Properties.VariableNames);
   testCase.verifyEqual(class(data.Data), class(data_before));
   testCase.verifyEqual(derived.met.Time, met_before.Time);
   testCase.verifyEqual(derived.met.Variables, met_before.Variables, 'AbsTol', 0);
   testCase.verifyEqual(derived.met.Properties.VariableNames, ...
      met_before.Properties.VariableNames);
   testCase.verifyEqual(class(derived.met), class(met_before));
   testCase.verifyEqual(data.sentinel_payload, sentinel);
   testCase.verifyEqual(derived.sentinel_payload, sentinel);
   testCase.verifyEqual(data.Data.Properties.UserData, data.artifact_metadata);
   testCase.verifyEqual(derived.met.Properties.UserData, ...
      derived.artifact_metadata);

   % Jan 10/12 are partial, Jan 11 is complete. Both artifacts carry the same
   % clipped native ledger while cumulative and non-day provenance survives.
   for metadata = {data.artifact_metadata, derived.artifact_metadata}
      md = metadata{1};
      testCase.verifyEqual(md.mar_qc_runoff_day_status, uint8([3; 1; 3]));
      testCase.verifyEqual(md.mar_qc_smb_day_status, uint8([3; 1; 3]));
      testCase.verifyEqual(md.mar_diagnostic_melt_day_status, ...
         uint8([3; 1; 3]));
      testCase.verifyTrue(all(isnan( ...
         md.mar_qc_runoff_daily_reference_mwe([1 3]))));
      testCase.verifyTrue(all(isnan( ...
         md.mar_diagnostic_melt_residual_mwe_day([1 3]))));
      testCase.verifyEqual(md.mar_qc_complete_utc_day_count, 1);
      testCase.verifyEqual(md.mar_qc_partial_utc_day_count, 2);
      testCase.verifyEqual(md.mar_qc_replaced_runoff_count, 17);
      testCase.verifyEqual(md.mar_qc_replaced_smb_count, 23);
      testCase.verifyEqual(md.source_files, "full-mar-2012.nc");
      testCase.verifyEqual(md.sentinel_policy, "preserve");
   end

   second = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice");
   testCase.verifyTrue(all(string({second.records.status}) == "unchanged"));
   testCase.verifyTrue(all(arrayfun(@(record) isempty(record.actions), ...
      second.records)));
   testCase.verifyEqual(second.mar_source_reads, 0);
end

function test_repairRcmArtifactMetadata_rejects_ambiguous_mar_ledger_axis( ...
      testCase)
   % Repair must not infer a calendar mapping when saved per-day vectors disagree.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-mar-ambiguous-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-mar-ambiguous-input');
   met_rel = "mar3.11/met_case1_mar3.11_20120110_20120112_15m.mat";
   data_rel = "mar3.11/case1_mar3.11_20120110_20120112.mat";
   met_file = fullfile(input_root, 'met', met_rel);
   data_file = fullfile(input_root, 'userdata', data_rel);
   writeTinyRepairPairManifest(eval_root, "case1", 67.0, -48.0, ...
      met_rel, data_rel);
   writeOverlongMarLedgerPair(met_file, data_file);

   % Corrupt both metadata copies consistently so parity is not the reason for
   % failure; only the unresolvable source-day lengths should block repair.
   for item = {met_file, 'met'; data_file, 'Data'}'
      loaded = load(item{1});
      artifact_metadata = loaded.artifact_metadata;
      artifact_metadata.mar_qc_runoff_daily_reference_mwe(end) = [];
      payload = loaded.(item{2});
      payload.Properties.UserData = artifact_metadata;
      loaded.(item{2}) = payload;
      save(item{1}, '-struct', 'loaded')
   end

   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice");
   statuses = string({report.records.status});
   reasons = string({report.records.reason});
   testCase.verifyTrue(all(statuses == "restage_required"), ...
      strjoin(statuses + ":" + reasons, newline));
   testCase.verifyTrue(all(contains(reasons, "ledger lengths disagree")));
   testCase.verifyEqual(report.mar_source_reads, 0);
end

function test_repairRcmArtifactMetadata_signed_mar_rz_source_light( ...
      testCase)
   % Legacy roundoff-only RZ metadata is replaced on both artifact cadences
   % without source reads or any numeric/time mutation, then remains idempotent.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-mar-rz-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-mar-rz-input');
   met_rel = "mar3.11/met_case1_mar3.11_20120101_20120102_15m.mat";
   data_rel = "mar3.11/case1_mar3.11_20120101_20120102.mat";
   met_file = fullfile(input_root, 'met', met_rel);
   data_file = fullfile(input_root, 'userdata', data_rel);
   writeTinyRepairPairManifest(eval_root, "case1", 67.0, -48.0, ...
      met_rel, data_rel);
   [data_before, met_before] = writeLegacySignedMarPair(met_file, data_file);

   dry = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice");
   dry_status = string({dry.records.status});
   dry_reason = string({dry.records.reason});
   testCase.verifyEqual(dry_status, ["would_repair", "would_repair"], ...
      strjoin(dry_status + ":" + dry_reason, newline));
   testCase.verifyTrue(all(arrayfun(@(record) any(string(record.actions) ...
      == "stamp_mar_refreeze_signed_metadata"), dry.records)));
   testCase.verifyEqual(dry.mar_source_reads, 0);

   written = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      dry_run=false);
   written_status = string({written.records.status});
   written_reason = string({written.records.reason});
   testCase.verifyEqual(written_status, ["repaired", "repaired"], ...
      strjoin(written_status + ":" + written_reason, newline));
   loaded_data = load(data_file, 'Data', 'artifact_metadata');
   loaded_met = load(met_file, 'met', 'artifact_metadata');
   testCase.verifyEqual(loaded_data.Data.Properties.RowTimes, ...
      data_before.Properties.RowTimes);
   testCase.verifyEqual(loaded_data.Data.Variables, data_before.Variables, ...
      'AbsTol', 0);
   testCase.verifyEqual(loaded_met.met.Properties.RowTimes, ...
      met_before.Properties.RowTimes);
   testCase.verifyEqual(loaded_met.met.Variables, met_before.Variables, ...
      'AbsTol', 0);
   for item = {loaded_data, 'Data'; loaded_met, 'met'}'
      metadata = item{1}.artifact_metadata;
      payload = item{1}.(item{2});
      testCase.verifyEqual(payload.Properties.UserData, metadata);
      testCase.verifyEqual(string( ...
         metadata.mar_diagnostic_refreeze_deposition_sign), ...
         "native_signed_combined_term_preserved_not_pure_refreeze");
      testCase.verifyEqual( ...
         metadata.mar_diagnostic_refreeze_negative_day_count, 1);
      testCase.verifyEqual( ...
         metadata.mar_diagnostic_refreeze_negative_minimum_mwe_h, -2e-4);
      testCase.verifyEqual( ...
         metadata.mar_diagnostic_refreeze_material_negative_threshold_mwe_h, ...
         1e-8);
      testCase.verifyEqual( ...
         metadata.mar_diagnostic_refreeze_material_negative_day_count, 1);
      testCase.verifyEqual( ...
         metadata.mar_diagnostic_refreeze_material_negative_minimum_mwe_h, ...
         -2e-4);
      testCase.verifyFalse(isfield(metadata, ...
         'mar_diagnostic_refreeze_negative_tolerance_mwe_h'));
   end

   second = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice");
   testCase.verifyEqual(string({second.records.status}), ...
      ["unchanged", "unchanged"]);
   testCase.verifyTrue(all(arrayfun(@(record) isempty(record.actions), ...
      second.records)));
   testCase.verifyEqual(second.mar_source_reads, 0);
end

function test_repairRcmArtifactMetadata_requires_restage_for_constant_daily_mar( ...
      testCase)
   % A formerly canonical constant-daily artifact cannot recover discarded
   % hourly structure from RU/SMB alone and must not be relabelled as hybrid.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-mar-legacy-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-mar-legacy-input');
   mar_dir = fullfile(testCase.TestData.tmp, 'repair-mar-legacy-source');
   met_rel = "mar3.11/met_case1_mar3.11_20120101_20120102_1hr.mat";
   data_rel = "mar3.11/case1_mar3.11_20120101_20120102.mat";
   met_file = fullfile(input_root, 'met', met_rel);
   data_file = fullfile(input_root, 'userdata', data_rel);
   writeTinyRepairPairManifest(eval_root, "case1", 67.0, -48.0, ...
      met_rel, data_rel);
   writeLegacyMarPair(met_file, data_file);
   writeTinyMarDailySource(mar_dir, 2012);

   loaded = load(met_file);
   loaded.artifact_metadata.mar_qc_method = 'native_daily';
   loaded.met.Properties.UserData = loaded.artifact_metadata;
   met = loaded.met;
   artifact_metadata = loaded.artifact_metadata;
   save(met_file, 'met', 'artifact_metadata')
   loaded = load(data_file);
   loaded.artifact_metadata.mar_qc_method = 'native_daily';
   loaded.Data.Properties.UserData = loaded.artifact_metadata;
   Data = loaded.Data;
   artifact_metadata = loaded.artifact_metadata;
   save(data_file, 'Data', 'artifact_metadata')

   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      mar_dir=string(mar_dir));

   statuses = string({report.records.status});
   reasons = string({report.records.reason});
   is_met = contains(string({report.records.filename}), filesep + "met" + filesep);
   diagnostic = strjoin(statuses + ":" + reasons, newline);
   testCase.verifyEqual(nnz(is_met), 1);
   testCase.verifyEqual(statuses(~is_met), "restage_required", diagnostic);
   testCase.verifyTrue(contains(reasons(~is_met), "requires full restage"), ...
      diagnostic);
   testCase.verifyNotEqual(statuses(is_met), "restage_required", diagnostic);
   testCase.verifyEqual(report.mar_source_reads, 0);
end

function test_repairRcmArtifactMetadata_canonicalizes_racmo_idempotently(testCase)
   % RACMO precip is renamed once while every unrelated value and time is kept.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-racmo-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-racmo-input');
   rel = "racmo2.3p3/case1_racmo2.3p3_2012.mat";
   data_file = fullfile(input_root, 'userdata', rel);
   writeTinyRepairManifest(eval_root, "promice", "case1", 67.0, -48.0, rel);
   writeLegacyRacmoData(data_file);

   dry = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice");
   testCase.verifyEqual(string(dry.records.status), "would_repair");
   testCase.verifyTrue(ismember("rename_racmo_precip_to_ppt", ...
      string(dry.records.actions)));

   icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      dry_run=false);
   loaded = load(data_file, 'Data', 'auxiliary');
   testCase.verifyTrue(ismember("ppt", ...
      string(loaded.Data.Properties.VariableNames)));
   testCase.verifyFalse(ismember("precip", ...
      string(loaded.Data.Properties.VariableNames)));
   testCase.verifyEqual(loaded.Data.ppt, [1; 2] * 1e-6);
   testCase.verifyEqual(loaded.Data.runoff, [3; 4]);
   testCase.verifyEqual(loaded.auxiliary, ...
      struct('label', "preserve", 'values', [1, 2, 3]));

   second = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice");
   testCase.verifyEqual(string(second.records.status), "unchanged");
end

function test_repairRcmArtifactMetadata_flips_marked_merra_once(testCase)
   % An explicit legacy positive-upward marker is inverted and replaced by the
   % canonical marker; the next pass must not invert the values again.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-merra-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-merra-input');
   rel = "merra2/case1_merra2_2012.mat";
   data_file = fullfile(input_root, 'userdata', rel);
   writeTinyRepairManifest(eval_root, "promice", "case1", 67.0, -48.0, rel);
   writeLegacyMerraData(data_file);

   first = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      dry_run=false);
   loaded = load(data_file, 'Data', 'artifact_metadata');
   testCase.verifyEqual(string(first.records.merra_flux_orientation), ...
      "positive_upward");
   testCase.verifyEqual(loaded.Data.shf, [-10; 5]);
   testCase.verifyEqual(loaded.Data.lhf, [-20; 4]);
   testCase.verifyEqual( ...
      string(loaded.artifact_metadata.merra_flux_sign_convention), ...
      "positive_toward_surface");

   second = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      dry_run=false);
   reloaded = load(data_file, 'Data');
   testCase.verifyEqual(string(second.records.status), "unchanged");
   testCase.verifyEqual(reloaded.Data.shf, loaded.Data.shf);
   testCase.verifyEqual(reloaded.Data.lhf, loaded.Data.lhf);
end

function test_repairRcmArtifactMetadata_caches_modis_and_is_idempotent(testCase)
   % Paired met/Data artifacts at one point share one GEUS read; a current
   % coverage marker makes the next pass source-light and non-mutating.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-modis-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-modis-input');
   modis_dir = fullfile(testCase.TestData.tmp, ...
      'scratch_2825_parent', 'repair-modis-source');
   met_rel = "mar3.11/met_case1_mar3.11_20120101_20120102_1hr.mat";
   data_rel = "mar3.11/case1_mar3.11_20120101_20120102.mat";
   met_file = fullfile(input_root, 'met', met_rel);
   data_file = fullfile(input_root, 'userdata', data_rel);
   writeTinyRepairPairManifest(eval_root, "case1", 67.0, -48.0, ...
      met_rel, data_rel);
   writeBareModisPair(met_file, data_file);
   writeTinyGeusModis(modis_dir, 2012, [0.5, 0.6]);

   first = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      modis_dir=string(modis_dir), dry_run=false);
   met_loaded = load(met_file, 'met');
   data_loaded = load(data_file, 'Data');
   testCase.verifyEqual(numel(first.records), 2);
   testCase.verifyEqual(first.modis_source_reads, 1);
   testCase.verifyGreaterThan(nnz(isfinite(met_loaded.met.modis)), 0);
   testCase.verifyGreaterThan(nnz(isfinite(data_loaded.Data.modis)), 0);

   second = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      modis_dir=string(modis_dir));
   testCase.verifyTrue(all(string({second.records.status}) == "unchanged"));
   testCase.verifyEqual(second.modis_source_reads, 0);
end

function test_repairRcmArtifactMetadata_batches_modis_points_by_year(testCase)
   % Distinct requested sites in one source year come from one bounding GEUS
   % read, not one full source-grid read per artifact or location.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-modis-batch-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-modis-batch-input');
   modis_dir = fullfile(testCase.TestData.tmp, 'repair-modis-batch-source');
   rel1 = "mar3.11/case1_mar3.11_2012.mat";
   rel2 = "mar3.11/case2_mar3.11_2012.mat";
   file1 = fullfile(input_root, 'userdata', rel1);
   file2 = fullfile(input_root, 'userdata', rel2);
   writeTinyTwoPointRepairManifest(eval_root, rel1, rel2);
   writeBareExistingData(file1);
   writeBareExistingData(file2);
   writeTinyGeusModis(modis_dir, 2012, [0.5, 0.6]);

   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      modis_dir=string(modis_dir), dry_run=false);
   loaded1 = load(file1, 'Data');
   loaded2 = load(file2, 'Data');

   testCase.verifyEqual(report.modis_source_reads, 1);
   testCase.verifyEqual(report.modis_cached_series, 2);
   testCase.verifyGreaterThan(nnz(isfinite(loaded1.Data.modis)), 0);
   testCase.verifyGreaterThan(nnz(isfinite(loaded2.Data.modis)), 0);
end

function test_repairRcmArtifactMetadata_distinguishes_no_modis_coverage(testCase)
   % An all-NaN temporary MODIS column outside the source years is removed and
   % classified as absent coverage, not as a failed external read.
   eval_root = fullfile(testCase.TestData.tmp, 'repair-modis-gap-eval');
   input_root = fullfile(testCase.TestData.tmp, 'repair-modis-gap-input');
   modis_dir = fullfile(testCase.TestData.tmp, 'repair-modis-gap-source');
   rel = "mar3.11/case1_mar3.11_1999.mat";
   data_file = fullfile(input_root, 'userdata', rel);
   writeTinyRepairManifest(eval_root, "promice", "case1", 67.0, -48.0, rel);
   writeAllNanModisData(data_file, 1999);
   writeTinyGeusModis(modis_dir, 2012, [0.5, 0.6]);

   report = icemodel.verification.setup.repairRcmArtifactMetadata( ...
      input_root, eval_root=eval_root, dataset_family="promice", ...
      modis_dir=string(modis_dir), dry_run=false);
   loaded = load(data_file, 'Data');

   testCase.verifyEqual(string(report.records.modis_status), ...
      "no_source_coverage");
   testCase.verifyFalse(ismember("modis", ...
      string(loaded.Data.Properties.VariableNames)));
   testCase.verifyEqual(report.modis_source_reads, 0);
end

function test_stageRcmForcing_marks_existing_sample_method(testCase)
   % Reused window files keep their saved method when the sampled point matches.
   met_root = fullfile(testCase.TestData.tmp, 'met');
   data_root = fullfile(testCase.TestData.tmp, 'userdata');
   mkdir(fullfile(met_root, 'mar3.11'));
   mkdir(fullfile(data_root, 'mar3.11'));
   writeTaggedExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'), "nearest", [67, -48]);

   L = struct('staged', true, 'years', 2012, ...
      'start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'end', datetime(2012, 6, 2, 'TimeZone', 'UTC'), ...
      'reason', '');
   legspec = struct('alias', "case1", 'mar', L);

   colocation = testCase.verifyWarning( ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      [67, -48], legspec=legspec, forcing_sources=["", "mar", ""], ...
      met_outdir=met_root, userdata_outdir=data_root, method="nearest"), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   testCase.verifyEqual(string(colocation{1}.mar.sample_method), "nearest");
   testCase.verifyFalse(isfield(colocation{1}.mar, 'requested_sample_method'));
end

function test_stageRcmForcing_reuses_existing_when_raw_leg_is_skipped(testCase)
   % Cached artifacts remain usable even if a fresh raw archive probe is absent.
   met_root = fullfile(testCase.TestData.tmp, 'met-skipped-raw');
   data_root = fullfile(testCase.TestData.tmp, 'userdata-skipped-raw');
   mkdir(fullfile(met_root, 'mar3.11'));
   mkdir(fullfile(data_root, 'mar3.11'));
   writeTaggedExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'), "nearest", [67, -48]);

   L = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
      'reason', 'MAR absent', ...
      'discovery_start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'discovery_end', datetime(2012, 6, 2, 'TimeZone', 'UTC'));
   legspec = struct('alias', "case1", 'mar', L);

   colocation = testCase.verifyWarning( ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      [67, -48], legspec=legspec, forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, method="nearest"), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   testCase.verifyTrue(logical(colocation{1}.mar.staged));
   testCase.verifyTrue(isfield(colocation{1}.mar, 'met_files'));
   testCase.verifyEqual(string(colocation{1}.mar.window.start), ...
      "2012-06-01 00:00:00");
   testCase.verifyEqual(string(colocation{1}.mar.window.end), ...
      "2012-06-02 00:00:00");
end

function test_stageRcmForcing_rejects_unproven_racmo_mask_cache(testCase)
   % A legacy coastal RACMO artifact without mask-selection provenance must not
   % bypass the current selector merely because its point/method metadata match.
   data_root = fullfile(testCase.TestData.tmp, 'userdata-racmo-mask');
   data_file = fullfile(data_root, 'racmo2.3p3', ...
      'case1_racmo2.3p3_20120101_20121231.mat');
   writeTaggedExistingData(data_file, "nearest", [67, -48]);
   saved = load(data_file, 'Data', 'artifact_metadata');
   saved.artifact_metadata = rmfield(saved.artifact_metadata, ...
      {'racmo_ice_mask_applied', 'racmo_point_max_distance_m'});
   saved.Data.Properties.UserData = saved.artifact_metadata;
   Data = saved.Data;
   artifact_metadata = saved.artifact_metadata;
   save(data_file, 'Data', 'artifact_metadata')

   L = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
      'reason', 'RACMO absent', ...
      'discovery_start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'discovery_end', datetime(2012, 6, 2, 'TimeZone', 'UTC'));
   legspec = struct('alias', "case1", 'racmo', L);
   colocation = icemodel.verification.setup.stageRcmForcing( ...
      [67, -48], legspec=legspec, forcing_sources="racmo", ...
      userdata_outdir=data_root, method="nearest", ...
      racmo_dir=fullfile(testCase.TestData.tmp, 'missing-racmo'));

   testCase.verifyFalse(logical(colocation{1}.racmo.staged));
   testCase.verifyFalse(isfield(colocation{1}.racmo, 'data_files'));
   testCase.verifySubstring(string(colocation{1}.racmo.reason), ...
      "artifact metadata does not match");
end

function test_stageRcmForcing_preserves_met_only_when_build_fails(testCase)
   % MAR/MERRA raw failure keeps compatible met-only forcing available.
   met_root = fullfile(testCase.TestData.tmp, 'met-only-buildable');
   data_root = fullfile(testCase.TestData.tmp, 'userdata-only-buildable');
   writeTaggedExistingMet(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_15m.mat'), ...
      "nearest", [67, -48]);

   L = struct('staged', true, 'years', 2012, ...
      'start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'end', datetime(2012, 6, 2, 'TimeZone', 'UTC'), ...
      'reason', '');
   legspec = struct('alias', "case1", 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      [67, -48], legspec=legspec, forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, method="nearest", ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'));

   testCase.verifyTrue(logical(colocation{1}.mar.staged));
   testCase.verifyTrue(isfield(colocation{1}.mar, 'met_files'));
   testCase.verifyFalse(isfield(colocation{1}.mar, 'data_files'));
   testCase.verifyTrue(contains(string(colocation{1}.mar.note), ...
      "Raw refresh failed"));
end

function test_stageRcmForcing_selects_compatible_cached_met_data_pair(testCase)
   % Pair selection should not drop Data because a longer disjoint cache exists.
   met_root = fullfile(testCase.TestData.tmp, 'met-pair');
   data_root = fullfile(testCase.TestData.tmp, 'userdata-pair');
   writeTaggedExistingData(fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'), "nearest", [67, -48]);
   writeTaggedExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20130601_20131231_15m.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20130601_20131231.mat'), "nearest", [67, -48]);

   L = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
      'reason', 'MAR absent', ...
      'discovery_start', datetime(2012, 1, 1, 'TimeZone', 'UTC'), ...
      'discovery_end', datetime(2013, 12, 31, 23, 0, 0, 'TimeZone', 'UTC'));
   legspec = struct('alias', "case1", 'mar', L);

   colocation = testCase.verifyWarning( ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      [67, -48], legspec=legspec, forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, method="nearest"), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   testCase.verifyTrue(isfield(colocation{1}.mar, 'data_files'));
   testCase.verifyTrue(contains(string(colocation{1}.mar.data_files), ...
      "20130601_20131231"));
   testCase.verifyEqual(string(colocation{1}.mar.window.start), ...
      "2013-06-01 00:00:00");
end

function test_stageRcmForcing_prefers_required_cover_before_discovery_rank(testCase)
   % A broad partial cache cannot outrank an exact required-window cache.
   met_root = fullfile(testCase.TestData.tmp, 'met-required-rank');
   data_root = fullfile(testCase.TestData.tmp, 'userdata-required-rank');
   point = [67, -48];
   writeTaggedExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20100101_20120630_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20100101_20120630.mat'), "nearest", point);
   writeTaggedExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'), "nearest", point);
   L = fixedClippedLeg(struct(), "mar");
   L.discovery_start = datetime(2010, 1, 1, 'TimeZone', 'UTC');
   L.discovery_end = datetime(2014, 12, 31, 23, 0, 0, ...
      'TimeZone', 'UTC');
   legspec = struct('alias', "case1", 'mar', L);

   colocation = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      point, legspec=legspec, forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, dt_out="", ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar')), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   testCase.verifyTrue(contains( ...
      string(colocation{1}.mar.met_files), "20120101_20121231"));
   testCase.verifyTrue(contains( ...
      string(colocation{1}.mar.data_files), "20120101_20121231"));
   testCase.verifyFalse(contains( ...
      string(colocation{1}.mar.note), "Raw refresh failed"));
end

function test_stageRcmForcing_rejects_existing_point_mismatch(testCase)
   % Existing window files sampled at another point must not be reused by alias.
   met_root = fullfile(testCase.TestData.tmp, 'met-point');
   data_root = fullfile(testCase.TestData.tmp, 'userdata-point');
   met_file = fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat');
   data_file = fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat');
   writeTaggedExistingFiles(met_file, data_file, "nearest", [67, -48]);

   L = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
      'reason', 'MAR absent', ...
      'discovery_start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'discovery_end', datetime(2012, 6, 2, 'TimeZone', 'UTC'));
   legspec = struct('alias', "case1", 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      [68, -48], legspec=legspec, forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, method="nearest", ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'));

   testCase.verifyFalse(logical(colocation{1}.mar.staged));
   testCase.verifyFalse(isfield(colocation{1}.mar, 'met_files'));
   testCase.verifyTrue(contains(string(colocation{1}.mar.reason), ...
      "artifact metadata does not match"));
end

function test_stageRcmForcing_unstamped_cache_does_not_suppress_build(testCase)
   % Missing cache metadata should not by itself skip a requested raw build.
   met_root = fullfile(testCase.TestData.tmp, 'met-unstamped-buildable');
   data_root = fullfile(testCase.TestData.tmp, 'userdata-unstamped-buildable');
   writeBareExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'));

   L = struct('staged', true, 'years', 2012, ...
      'start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'end', datetime(2012, 6, 2, 'TimeZone', 'UTC'), ...
      'reason', '');
   legspec = struct('alias', "case1", 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      [67, -48], legspec=legspec, forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, method="nearest", ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'));

   testCase.verifyFalse(logical(colocation{1}.mar.staged));
   testCase.verifyFalse(contains(string(colocation{1}.mar.reason), ...
      "artifact metadata is missing"));
   testCase.verifyTrue(contains(string(colocation{1}.mar.reason), ...
      "source directory not found"));
end

function test_stageRcmForcing_conflicts_mismatched_cache_before_rebuild(testCase)
   % A same-name cache for another point must not be overwritten by rebuild.
   met_root = fullfile(testCase.TestData.tmp, 'met-point-buildable');
   data_root = fullfile(testCase.TestData.tmp, 'userdata-point-buildable');
   met_file = fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat');
   data_file = fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat');
   writeTaggedExistingFiles(met_file, data_file, "nearest", [67, -48]);
   before = dir(met_file);

   L = struct('staged', true, 'years', 2012, ...
      'start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'end', datetime(2012, 6, 2, 'TimeZone', 'UTC'), ...
      'reason', '');
   legspec = struct('alias', "case1", 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      [68, -48], legspec=legspec, forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, method="nearest", ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'));
   after = dir(met_file);

   testCase.verifyFalse(logical(colocation{1}.mar.staged));
   testCase.verifyTrue(contains(string(colocation{1}.mar.reason), ...
      "artifact metadata does not match"));
   testCase.verifyEqual(after.datenum, before.datenum);
end

function test_stageRcmForcing_manifest_conflict_replaces_prior_leg(testCase)
   % Manifest refresh must not preserve old legs after metadata conflicts.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'manifest-conflict-replace.json');
   met_root = fullfile(testCase.TestData.tmp, 'met-manifest-conflict');
   data_root = fullfile(testCase.TestData.tmp, 'userdata-manifest-conflict');
   writeTaggedExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'), "nearest", [67, -48]);
   state = helperState();
   state.colocation.mar = existingLeg("mar");
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.site_location.lat_wgs84 = 68.0;
   returned = icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=manifest, manifest_file=manifest_file, ...
      forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'), ...
      overwrite=true);

   testCase.verifyFalse(logical(returned.cases.colocation.mar.staged));
   testCase.verifyTrue(contains(string(returned.cases.colocation.mar.reason), ...
      "artifact metadata does not match"));
end

function test_stageRcmForcing_manifest_preserves_prior_unstamped_after_raw_failure(testCase)
   % A referenced unstamped prior leg survives an unavailable raw refresh.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'manifest-missing-metadata-preserve.json');
   met_root = fullfile(testCase.TestData.tmp, ...
      'met-manifest-missing-metadata');
   data_root = fullfile(testCase.TestData.tmp, ...
      'userdata-manifest-missing-metadata');
   writeBareExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'));
   state = helperState();
   state.colocation.mar = existingLeg("mar");
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   % Directory discovery cannot newly attach the unstamped files, but the prior
   % manifest references remain the compatibility evidence after raw failure.
   manifest = jsondecode(fileread(manifest_file));
   returned = icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=manifest, manifest_file=manifest_file, ...
      forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'), ...
      overwrite=true);

   testCase.verifyTrue(logical(returned.cases.colocation.mar.staged));
   testCase.verifyEqual(string(returned.cases.colocation.mar.met_files), ...
      "mar3.11/met_case1_mar3.11_20120101_20121231_1hr.mat");
   testCase.verifyEqual(string(returned.cases.colocation.mar.data_files), ...
      "mar3.11/case1_mar3.11_20120101_20121231.mat");
end

function test_stageRcmForcing_manifest_ignores_unreferenced_unstamped_collision(testCase)
   % An unknown directory cache cannot displace a different compatible prior.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'manifest-unreferenced-unstamped.json');
   met_root = fullfile(testCase.TestData.tmp, ...
      'met-unreferenced-unstamped');
   data_root = fullfile(testCase.TestData.tmp, ...
      'userdata-unreferenced-unstamped');

   % The legacy-path files are explicitly referenced by the prior manifest.
   prior = existingLeg("mar");
   prior.met_files = ...
      "legacy/met_case1_mar3.11_20120101_20121231_1hr.mat";
   prior.data_files = ...
      "legacy/case1_mar3.11_20120101_20121231.mat";
   writeBareExistingFiles(fullfile(met_root, prior.met_files), ...
      fullfile(data_root, prior.data_files));

   % Same-alias files in the current storage folder are unreferenced and lack
   % metadata, so discovery must not attach them when the raw source is absent.
   writeBareExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'));
   state = helperState();
   state.colocation.mar = prior;
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   manifest = jsondecode(fileread(manifest_file));
   returned = icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=manifest, manifest_file=manifest_file, ...
      forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'));

   testCase.verifyTrue(logical(returned.cases.colocation.mar.staged));
   testCase.verifyEqual(string(returned.cases.colocation.mar.met_files), ...
      prior.met_files);
   testCase.verifyEqual(string(returned.cases.colocation.mar.data_files), ...
      prior.data_files);
end

function test_stageRcmForcing_manifest_rejects_legacy_prior_point_mismatch(testCase)
   % Referenced files outside discovery still honor concrete saved coordinates.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'manifest-legacy-point-mismatch.json');
   met_root = fullfile(testCase.TestData.tmp, ...
      'met-legacy-point-mismatch');
   data_root = fullfile(testCase.TestData.tmp, ...
      'userdata-legacy-point-mismatch');
   prior = existingLeg("mar");
   prior.met_files = ...
      "legacy/met_case1_mar3.11_20120101_20121231_1hr.mat";
   prior.data_files = ...
      "legacy/case1_mar3.11_20120101_20121231.mat";
   writeTaggedExistingFiles(fullfile(met_root, prior.met_files), ...
      fullfile(data_root, prior.data_files), "nearest", [66, -47]);
   state = helperState();
   state.colocation.mar = prior;
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   manifest = jsondecode(fileread(manifest_file));
   returned = icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=manifest, manifest_file=manifest_file, ...
      forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'));

   testCase.verifyFalse(logical(returned.cases.colocation.mar.staged));
   testCase.verifyFalse(isfield(returned.cases.colocation.mar, 'met_files'));
   testCase.verifyFalse(isfield( ...
      returned.cases.colocation.mar, 'replace_prior_artifacts'));
   testCase.verifyGreaterThan(strlength( ...
      string(returned.cases.colocation.mar.reason)), 0);
end

function test_stageRcmForcing_manifest_preserves_prior_data_on_partial_cache(testCase)
   % A met-only rediscovery must not drop prior manifest-backed Data refs.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'manifest-partial-cache.json');
   met_root = fullfile(testCase.TestData.tmp, 'met-partial-cache');
   data_root = fullfile(testCase.TestData.tmp, 'userdata-partial-cache');
   writeTaggedExistingMet(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat'), ...
      "nearest", [67, -48]);
   writeBareExistingData(fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'));
   state = helperState();
   state.colocation.mar = existingLeg("mar");
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   manifest = jsondecode(fileread(manifest_file));
   returned = icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=manifest, manifest_file=manifest_file, ...
      forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'));

   testCase.verifyTrue(isfield(returned.cases.colocation.mar, 'data_files'));
   testCase.verifyTrue(isfield(returned.cases.colocation.mar, 'met_files'));
   testCase.verifyTrue(ismember("mar3.11", ...
      string(returned.cases.eval_sources)));
end

function test_stageRcmForcing_manifest_drops_prior_leg_with_missing_files(testCase)
   % Prior staged legs are preserved only while their artifacts still resolve.
   manifest_file = fullfile(testCase.TestData.tmp, ...
      'manifest-prior-missing-files.json');
   state = helperState();
   state.colocation.mar = existingLeg("mar");
   icemodel.verification.setup.buildDatasetFamilyManifest( ...
      state, true, dataset_family="helper", manifest_file=manifest_file, ...
      requested_ids="case1", source_url="test", source_version="v1", ...
      entry_callback=@helperEntry, overwrite_family=true);

   manifest = jsondecode(fileread(manifest_file));
   returned = icemodel.verification.setup.stageRcmForcing( ...
      obs_manifest=manifest, manifest_file=manifest_file, ...
      forcing_sources="mar", ...
      met_outdir=fullfile(testCase.TestData.tmp, 'missing-met'), ...
      userdata_outdir=fullfile(testCase.TestData.tmp, 'missing-userdata'), ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'));

   testCase.verifyFalse(logical(returned.cases.colocation.mar.staged));
   testCase.verifyFalse(ismember("mar3.11", ...
      string(returned.cases.forcing_sources)));
   testCase.verifyFalse(isfield( ...
      returned.cases.colocation.mar, 'replace_prior_artifacts'));
end

function test_stageRcmForcing_bad_data_only_cache_does_not_skip_rebuild(testCase)
   % A corrupt Data-only cache should not suppress raw rebuild eligibility.
   data_root = fullfile(testCase.TestData.tmp, 'userdata-bad-data-only');
   data_file = fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat');
   writeBadTaggedExistingData(data_file, "nearest", [67, -48]);

   L = struct('staged', true, 'years', 2012, ...
      'start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'end', datetime(2012, 6, 2, 'TimeZone', 'UTC'), ...
      'reason', '');
   legspec = struct('alias', "case1", 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      [67, -48], legspec=legspec, forcing_sources="mar", ...
      userdata_outdir=data_root, method="nearest", ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'));

   testCase.verifyFalse(logical(colocation{1}.mar.staged));
   testCase.verifyFalse(isfield(colocation{1}.mar, 'data_files'));
   testCase.verifyFalse(contains(string(colocation{1}.mar.reason), "Data"));
end

function test_stageRcmForcing_rejects_existing_method_mismatch(testCase)
   % A natural-sampling request must not reuse a nearest-sampled staged artifact.
   met_root = fullfile(testCase.TestData.tmp, 'met-mismatch');
   data_root = fullfile(testCase.TestData.tmp, 'userdata-mismatch');
   met_file = fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat');
   data_file = fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat');
   writeTaggedExistingFiles(met_file, data_file, "nearest");

   L = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
      'reason', 'MAR absent', ...
      'discovery_start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'discovery_end', datetime(2012, 6, 2, 'TimeZone', 'UTC'));
   legspec = struct('alias', "case1", 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      [67, -48], legspec=legspec, forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, method="natural", ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'));

   testCase.verifyFalse(logical(colocation{1}.mar.staged));
   testCase.verifyTrue(contains(string(colocation{1}.mar.reason), ...
      "artifact metadata does not match"));
end

function test_stageRcmForcing_conflicts_stale_cache_when_buildable(testCase)
   % A stale same-alias cache must not be overwritten by a rebuild attempt.
   met_root = fullfile(testCase.TestData.tmp, 'met-stale-buildable');
   data_root = fullfile(testCase.TestData.tmp, 'userdata-stale-buildable');
   writeTaggedExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'), "nearest", [68, -48]);

   L = struct('staged', true, 'years', 2012, ...
      'start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'end', datetime(2012, 6, 2, 'TimeZone', 'UTC'), ...
      'reason', '');
   legspec = struct('alias', "case1", 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      [67, -48], legspec=legspec, forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, method="nearest", ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'));

   testCase.verifyFalse(logical(colocation{1}.mar.staged));
   testCase.verifyTrue(contains(string(colocation{1}.mar.reason), ...
      "artifact metadata does not match"));
end

function test_stageRcmForcing_derives_missing_met_from_existing_data(testCase)
   % A reusable MAR/MERRA Data file is enough to recreate the met companion.
   met_root = fullfile(testCase.TestData.tmp, 'met-derived');
   data_root = fullfile(testCase.TestData.tmp, 'userdata-derived');
   data_file = fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat');
   writeTaggedExistingData(data_file, "nearest", [67, -48]);

   L = struct('staged', true, 'years', 2012, ...
      'start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'end', datetime(2012, 6, 2, 'TimeZone', 'UTC'), ...
      'reason', '');
   legspec = struct('alias', "case1", 'mar', L);

   colocation = icemodel.verification.setup.stageRcmForcing( ...
      [67, -48], legspec=legspec, forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, method="nearest", ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mar'));

   met_file = fullfile(met_root, colocation{1}.mar.met_files(1));
   testCase.verifyTrue(logical(colocation{1}.mar.staged));
   testCase.verifyTrue(isfile(met_file));
   testCase.verifyTrue(contains(string(colocation{1}.mar.data_files), ...
      "case1_mar3.11_20120101_20121231.mat"));
end

function test_stageRcmForcing_ignores_incompatible_met_when_data_is_reusable(testCase)
   % A compatible Data side may rebuild met without linking or replacing a
   % different-point/method met artifact discovered in the same window.
   met_root = fullfile(testCase.TestData.tmp, 'mixed-compatible-data');
   data_root = fullfile(testCase.TestData.tmp, 'mixed-compatible-data-userdata');
   data_file = fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat');
   incompatible_met = fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20100101_20131231_15m.mat');
   writeTaggedExistingData(data_file, "nearest", [67, -48]);
   writeTaggedExistingMet(incompatible_met, "natural", [68, -47]);
   incompatible_before = fileBytes(incompatible_met);

   L = struct('staged', true, 'years', 2012, ...
      'start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'end', datetime(2012, 6, 2, 'TimeZone', 'UTC'), 'reason', '');
   legspec = struct('alias', "case1", 'mar', L);
   colocation = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      [67, -48], legspec=legspec, forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, method="nearest", ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mixed-mar')), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   leg = colocation{1}.mar;
   returned_met = fullfile(met_root, leg.met_files(1));
   testCase.verifyTrue(logical(leg.staged));
   testCase.verifyTrue(contains(string(leg.data_files), ...
      "case1_mar3.11_20120101_20121231.mat"));
   testCase.verifyNotEqual(string(returned_met), string(incompatible_met));
   metadata = icemodel.verification.setup.readRcmArtifactMetadata(returned_met);
   testCase.verifyEqual(string(metadata.sample_method), "nearest");
   testCase.verifyEqual([metadata.lat_wgs84, metadata.lon_wgs84], [67, -48]);
   testCase.verifyEqual(fileBytes(incompatible_met), incompatible_before);
end

function test_stageRcmForcing_ignores_incompatible_data_in_met_only_fallback(testCase)
   % A compatible met side remains usable when raw refresh fails; an incompatible
   % Data sibling must remain byte-identical and absent from the returned leg.
   met_root = fullfile(testCase.TestData.tmp, 'mixed-compatible-met');
   data_root = fullfile(testCase.TestData.tmp, 'mixed-compatible-met-userdata');
   met_file = fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_15m.mat');
   incompatible_data = fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20100101_20131231.mat');
   writeTaggedExistingMet(met_file, "nearest", [67, -48]);
   writeTaggedExistingData(incompatible_data, "natural", [68, -47]);
   incompatible_before = fileBytes(incompatible_data);

   L = struct('staged', true, 'years', 2012, ...
      'start', datetime(2012, 6, 1, 'TimeZone', 'UTC'), ...
      'end', datetime(2012, 6, 2, 'TimeZone', 'UTC'), 'reason', '');
   legspec = struct('alias', "case1", 'mar', L);
   colocation = icemodel.test.helpers.captureExpectedWarning(testCase, ...
      @() icemodel.verification.setup.stageRcmForcing( ...
      [67, -48], legspec=legspec, forcing_sources="mar", ...
      met_outdir=met_root, userdata_outdir=data_root, method="nearest", ...
      mar_dir=fullfile(testCase.TestData.tmp, 'missing-mixed-mar')), ...
      'icemodel:verification:stageRcmForcing:existingWindowFile');

   leg = colocation{1}.mar;
   returned_met = fullfile(met_root, leg.met_files(1));
   testCase.verifyTrue(logical(leg.staged));
   testCase.verifyFalse(isfield(leg, 'data_files'));
   testCase.verifyEqual(string(returned_met), string(met_file));
   metadata = icemodel.verification.setup.readRcmArtifactMetadata(returned_met);
   testCase.verifyEqual(string(metadata.sample_method), "nearest");
   testCase.verifyEqual([metadata.lat_wgs84, metadata.lon_wgs84], [67, -48]);
   testCase.verifyEqual(fileBytes(incompatible_data), incompatible_before);
end

function test_stageDatasetFamilyCases_rethrows_unexpected_errors(testCase)
   % skip_missing is for absent data/windows, not malformed staging regressions.
   f = @() icemodel.verification.setup.stageDatasetFamilyCases( ...
      1, helperState(), @(~, ~) malformedStage(), ...
      warning_id="test:skip", skip_missing=true);

   testCase.verifyError(f, 'test:malformedStage');
end

function test_stageDatasetFamilyCases_skips_label_callback_missing_data(testCase)
   % Label lookup can touch source catalogs and should honor skip_missing.
   [~, alive, skipped] = icemodel.verification.setup.stageDatasetFamilyCases( ...
      "lynl", helperState(), @(~, ~) stageShouldNotRun(), ...
      warning_id="test:skip", skip_missing=true, ...
      label_callback=@missingPromiceLabel);

   testCase.verifyFalse(alive);
   testCase.verifyEqual(string(skipped.site), "lynl");
   testCase.verifyTrue(contains(string(skipped.reason), "missing station"));
end

function test_stageDatasetFamilyCases_skips_empty_protocol_window(testCase)
   % Empty protocol windows are per-case data gaps, not family aborts.
   [~, alive, skipped] = icemodel.verification.setup.stageDatasetFamilyCases( ...
      "site1", helperState(), @(~, ~) emptyProtocolWindow(), ...
      warning_id="test:skip", skip_missing=true);

   testCase.verifyFalse(alive);
   testCase.verifyEqual(string(skipped.site), "site1");
   testCase.verifyTrue(contains(string(skipped.reason), ...
      "no rows in requested window"));
end

function test_stageDatasetFamilyCases_skips_promice_source_not_found(testCase)
   % Missing PROMICE source roots are absent-data gaps under skip_missing.
   [~, alive, skipped] = icemodel.verification.setup.stageDatasetFamilyCases( ...
      "kanm", helperState(), @(~, ~) missingPromiceSource(), ...
      warning_id="test:skip", skip_missing=true);

   testCase.verifyFalse(alive);
   testCase.verifyEqual(string(skipped.site), "kanm");
   testCase.verifyTrue(contains(string(skipped.reason), ...
      "source missing"));
end

%% Local fixture helpers
function state = reusePrototype(include_location)
   %REUSEPROTOTYPE Minimal state schema required by the guarded reuse helper.
   state = struct('case_id', "", 'site_id', "", 'site_name', "", ...
      'point', [NaN NaN], ...
      'evaluation_file_rel', "", ...
      'entry', struct(), 'period', struct('start', '', 'end', ''), ...
      'colocation', struct(), 'leg', struct(), 'reuse_entry', false, ...
      'dry_run', false);
   if include_location
      state.site_location = struct();
   end
end

function state = helperState()
   %HELPERSTATE Minimal state record consumed by the shared helper callbacks.
   state = struct( ...
      'case_id', "case1", ...
      'site_id', "CASE_1", ...
      'site_name', "Case 1", ...
      'point', [67.0, -48.0], ...
      'site_location', struct('lat_wgs84', 67.0, 'lon_wgs84', -48.0, ...
         'x_epsg3413', 1, 'y_epsg3413', 2, 'elev_m', 1000), ...
      'period', struct('start', '2012-01-01', 'end', '2012-01-02'), ...
      'evaluation_file_rel', 'case1/observations.mat', ...
      'colocation', struct('sumup', ...
         struct('kind', 'firn_profile_obs', 'staged', true)), ...
      'comparison_variables', "density", ...
      'observation_variables', struct('density', 'profile density'));
end

function L = existingLeg(source)
   %EXISTINGLEG Return a staged RCM manifest leg for preservation tests.
   product = icemodel.verification.namelists.rcmProductIds(source);
   L = struct('kind', 'point_met', 'staged', true, ...
      'source', char(source), 'source_id', char(product), ...
      'sample_method', 'nearest', ...
      'window', struct('start', '2012-01-01', 'end', '2012-01-02'), ...
      'data_files', product + "/case1_" + product + "_20120101_20121231.mat");
   if source ~= "racmo"
      L.met_files = product + "/met_case1_" + product ...
         + "_20120101_20121231_1hr.mat";
   end
end

function [met_root, data_root] = seedExistingRcmArtifacts(testCase)
   %SEEDEXISTINGRCMARTIFACTS Write tiny current-contract files for preservation.
   met_root = fullfile(testCase.TestData.tmp, 'preserved-met');
   data_root = fullfile(testCase.TestData.tmp, 'preserved-userdata');
   writeTaggedExistingFiles(fullfile(met_root, 'mar3.11', ...
      'met_case1_mar3.11_20120101_20121231_1hr.mat'), ...
      fullfile(data_root, 'mar3.11', ...
      'case1_mar3.11_20120101_20121231.mat'), "nearest", [67, -48]);
end

function writeTaggedExistingFiles( ...
      met_file, data_file, method, point, missing_albedo)
   %WRITETAGGEDEXISTINGFILES Create tiny staged RCM artifacts with method metadata.
   if nargin < 4
      point = [67, -48];
   end
   if nargin < 5
      % Ordinary cache fixtures are strictly ready unless a test opts out.
      missing_albedo = false;
   end
   Data = timetable((datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:1))');
   artifact_metadata = currentRacmoCacheMetadata(method, point);
   Data.tair = [260; 261];
   Data.Properties.UserData = artifact_metadata;

   required = icemodel.forcing.helpers.metvariables();
   met = array2timetable(ones(height(Data), numel(required)), ...
      RowTimes=Data.Time, VariableNames=cellstr(required));
   met = icemodel.forcing.helpers.stampMetadata(met);
   if missing_albedo
      % Model a source-faithful MERRA darkness gap on one required channel.
      met.albedo(1) = NaN;
   end
   met.Properties.UserData = artifact_metadata;

   ensureParent(met_file);
   ensureParent(data_file);
   save(met_file, 'met', 'artifact_metadata')
   save(data_file, 'Data', 'artifact_metadata')
end

function writeTaggedExistingMet(met_file, method, point)
   %WRITETAGGEDEXISTINGMET Create a tiny staged RCM met artifact.
   time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:1))';
   artifact_metadata = currentRacmoCacheMetadata(method, point);

   % Cache reuse now validates and diagnoses the exact selected met bytes, so
   % this fixture must model a structurally valid saved forcing artifact.
   required = icemodel.forcing.helpers.metvariables();
   met = array2timetable(ones(numel(time), numel(required)), ...
      RowTimes=time, VariableNames=cellstr(required));
   met = icemodel.forcing.helpers.stampMetadata(met);
   met.Properties.UserData = artifact_metadata;

   ensureParent(met_file);
   save(met_file, 'met', 'artifact_metadata')
end

function metadata = currentRacmoCacheMetadata(method, point)
   %CURRENTRACMOCACHEMETADATA Minimal proof carried by mask-aware cache fixtures.
   metadata = struct('sample_method', char(method), ...
      'lat_wgs84', point(1), 'lon_wgs84', point(2), ...
      'racmo_ice_mask_applied', true, ...
      'racmo_point_max_distance_m', 15e3);
end

function writeTaggedExistingData(data_file, method, point)
   %WRITETAGGEDEXISTINGDATA Create a tiny staged RCM Data artifact.
   Data = timetable((datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:1))');
   artifact_metadata = currentRacmoCacheMetadata(method, point);
   Data.tair = [260; 261];
   Data.ppt = [0; 0];
   Data.Properties.UserData = artifact_metadata;

   ensureParent(data_file);
   save(data_file, 'Data', 'artifact_metadata')
end

function writeTaggedMetContractData(data_file, method, point)
   %WRITETAGGEDMETCONTRACTDATA Create Data convertible to a met artifact.
   time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:1))';
   vars = icemodel.forcing.helpers.metvariables();
   values = ones(numel(time), numel(vars));
   Data = array2timetable(values, 'RowTimes', time, ...
      'VariableNames', cellstr(vars));
   Data.tair(:) = 260;
   Data.ppt(:) = 0;
   artifact_metadata = struct('sample_method', char(method), ...
      'lat_wgs84', point(1), 'lon_wgs84', point(2));
   Data.Properties.UserData = artifact_metadata;

   ensureParent(data_file);
   save(data_file, 'Data', 'artifact_metadata')
end

function writeBadTaggedExistingData(data_file, method, point)
   %WRITEBADTAGGEDEXISTINGDATA Save matching metadata with no Data payload.
   artifact_metadata = struct('sample_method', char(method), ...
      'lat_wgs84', point(1), 'lon_wgs84', point(2));
   notData = "bad cache";

   ensureParent(data_file);
   save(data_file, 'notData', 'artifact_metadata')
end

function writeTinyRepairManifest(eval_root, family, case_id, lat, lon, ...
      data_file, method)
   %WRITETINYREPAIRMANIFEST Write a current-schema manifest for repair tests.
   if nargin < 7
      method = "nearest";
   end
   family_root = fullfile(eval_root, family);
   if ~isfolder(family_root)
      mkdir(family_root)
   end
   colocation = struct();
   if strlength(string(data_file)) > 0
      colocation.mar = struct('kind', 'point_met', 'staged', true, ...
         'source', 'mar', 'source_id', 'mar3.11', ...
         'data_files', {{char(data_file)}}, ...
         'sample_method', char(method));
   end
   manifest = struct( ...
      'dataset_family', char(family), ...
      'source_doi', '', ...
      'source_url', '', ...
      'source_version', '', ...
      'retrieval_date', '', ...
      'cases', struct( ...
      'case_id', char(case_id), ...
      'site_location', struct('lat_wgs84', lat, 'lon_wgs84', lon), ...
      'colocation', colocation), ...
      'skipped', []);
   fid = fopen(fullfile(family_root, 'manifest.json'), 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s\n', jsonencode(manifest));
   clear cleaner
end

function writeBareExistingData(data_file)
   %WRITEBAREEXISTINGDATA Create a tiny Data artifact without metadata.
   Data = timetable((datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:1))');
   Data.tair = [260; 261];

   ensureParent(data_file);
   save(data_file, 'Data')
end

function writeLegacyRacmoData(data_file)
   %WRITELEGACYRACMODATA Create a pre-canonical artifact with auxiliary payload.
   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:1))';
   Data = timetable(Time, [1; 2] * 1e-6, [3; 4], ...
      'VariableNames', {'precip', 'runoff'});
   auxiliary = struct('label', "preserve", 'values', [1, 2, 3]);
   ensureParent(data_file);
   save(data_file, 'Data', 'auxiliary')
end

function writeLegacyMerraData(data_file)
   %WRITELEGACYMERRADATA Create an explicitly positive-upward MERRA artifact.
   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:1))';
   Data = timetable(Time, [10; -5], [20; -4], [260; 261], ...
      'VariableNames', {'shf', 'lhf', 'tair'});
   artifact_metadata = struct( ...
      'merra_flux_sign_convention', 'positive_upward');
   Data.Properties.UserData = artifact_metadata;
   ensureParent(data_file);
   save(data_file, 'Data', 'artifact_metadata')
end

function writeBareModisPair(met_file, data_file)
   %WRITEBAREMODISPAIR Create paired artifacts sharing one location/time axis.
   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2012, 1, 2, 23, 0, 0, 'TimeZone', 'UTC'))';
   tair = 260 * ones(numel(Time), 1);
   Data = timetable(Time, tair);
   met = Data;
   ensureParent(met_file);
   ensureParent(data_file);
   save(met_file, 'met')
   save(data_file, 'Data')
end

function writeAllNanModisData(data_file, yyyy)
   %WRITEALLNANMODISDATA Create a temporary all-NaN uncovered MODIS artifact.
   Time = (datetime(yyyy, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:1))';
   Data = timetable(Time, [260; 261], nan(2, 1), ...
      'VariableNames', {'tair', 'modis'});
   ensureParent(data_file);
   save(data_file, 'Data')
end

function writeTinyRepairPairManifest( ...
      eval_root, case_id, lat, lon, met_file, data_file)
   %WRITETINYREPAIRPAIRMANIFEST Reference paired met/Data repair fixtures.
   family_root = fullfile(eval_root, 'promice');
   if ~isfolder(family_root)
      mkdir(family_root)
   end
   leg = struct('kind', 'point_met', 'staged', true, ...
      'source', 'mar', 'source_id', 'mar3.11', ...
      'met_files', {{char(met_file)}}, ...
      'data_files', {{char(data_file)}}, ...
      'sample_method', 'nearest');
   manifest = struct( ...
      'dataset_family', 'promice', ...
      'source_doi', '', 'source_url', '', 'source_version', '', ...
      'retrieval_date', '', ...
      'cases', struct( ...
      'case_id', char(case_id), ...
      'site_location', struct('lat_wgs84', lat, 'lon_wgs84', lon), ...
      'colocation', struct('mar', leg)), ...
      'skipped', []);
   fid = fopen(fullfile(family_root, 'manifest.json'), 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s\n', jsonencode(manifest));
   clear cleaner
end

function writeTinyTwoPointRepairManifest(eval_root, data_file1, data_file2)
   %WRITETINYTWOPOINTREPAIRMANIFEST Reference two distinct point artifacts.
   family_root = fullfile(eval_root, 'promice');
   if ~isfolder(family_root)
      mkdir(family_root)
   end
   cases = repmat(struct(), 2, 1);
   aliases = ["case1", "case2"];
   locations = [67.0, -48.0; 68.0, -47.0];
   files = [string(data_file1), string(data_file2)];
   for k = 1:2
      leg = struct('kind', 'point_met', 'staged', true, ...
         'source', 'mar', 'source_id', 'mar3.11', ...
         'data_files', {{char(files(k))}}, ...
         'sample_method', 'nearest');
      cases(k).case_id = char(aliases(k));
      cases(k).site_location = struct( ...
         'lat_wgs84', locations(k, 1), ...
         'lon_wgs84', locations(k, 2));
      cases(k).colocation = struct('mar', leg);
   end
   manifest = struct( ...
      'dataset_family', 'promice', ...
      'source_doi', '', 'source_url', '', 'source_version', '', ...
      'retrieval_date', '', 'cases', cases, 'skipped', []);
   fid = fopen(fullfile(family_root, 'manifest.json'), 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s\n', jsonencode(manifest));
   clear cleaner
end

function writeLegacyMarPair(met_file, data_file)
   %WRITELEGACYMARPAIR Create paired pre-native-daily MAR artifacts.
   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2012, 1, 2, 23, 0, 0, 'TimeZone', 'UTC'))';
   tair = 260 * ones(numel(Time), 1);
   runoff = 9 * ones(numel(Time), 1);
   smb = 8 * ones(numel(Time), 1);
   melt = (1:numel(Time))';
   Data = timetable(Time, tair, runoff, smb, melt);
   met = Data(:, ["tair", "runoff"]);
   artifact_metadata = struct('sample_method', 'nearest', ...
      'lat_wgs84', 67, 'lon_wgs84', -48);
   Data.Properties.UserData = artifact_metadata;
   met.Properties.UserData = artifact_metadata;
   ensureParent(met_file);
   ensureParent(data_file);
   save(met_file, 'met', 'artifact_metadata')
   save(data_file, 'Data', 'artifact_metadata')
end

function met = makeLegacyMarMet15m(met_file)
   %MAKELEGACYMARMET15M Expand the legacy hourly met payload without interpolation.
   loaded = load(met_file, 'met', 'artifact_metadata');
   source = loaded.met;
   Time = (source.Time(1):minutes(15):source.Time(end) + minutes(45))';
   values = repelem(source.Variables, 4, 1);
   met = array2timetable(values, RowTimes=Time, ...
      VariableNames=source.Properties.VariableNames);
   artifact_metadata = loaded.artifact_metadata;
   met.Properties.UserData = artifact_metadata;
   save(met_file, 'met', 'artifact_metadata')
end

function [Data, met] = writeLegacySignedMarPair(met_file, data_file)
   %WRITELEGACYSIGNEDMARPAIR Save paired RZ values under the old sign policy.
   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2012, 1, 2, 23, 0, 0, 'TimeZone', 'UTC'))';
   refreeze_deposition = [ ...
      -2e-4 * ones(24, 1); 3e-4 * ones(24, 1)];
   Data = timetable(Time, refreeze_deposition);
   met_time = (Time(1):minutes(15):Time(end) + minutes(45))';
   met = timetable(met_time, repelem(refreeze_deposition, 4), ...
      VariableNames="refreeze_deposition");
   artifact_metadata = struct('sample_method', 'nearest', ...
      'lat_wgs84', 67, 'lon_wgs84', -48, ...
      'mar_diagnostic_refreeze_deposition_sign', ...
      'positive_gain_tiny_negative_roundoff_allowed', ...
      'mar_diagnostic_refreeze_negative_tolerance_mwe_h', 1e-8);
   Data.Properties.UserData = artifact_metadata;
   met.Properties.UserData = artifact_metadata;
   ensureParent(met_file);
   ensureParent(data_file);
   save(met_file, 'met', 'artifact_metadata')
   save(data_file, 'Data', 'artifact_metadata')
end

function [Data, met, sentinel_payload] = ...
      writeOverlongMarLedgerPair(met_file, data_file)
   %WRITEOVERLONGMARLEDGERPAIR Save clipped payloads with full-2012 provenance.
   Time = (datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2012, 12, 31, 23, 0, 0, 'TimeZone', 'UTC'))';
   n = numel(Time);
   tair = 260 + zeros(n, 1);
   runoff = 1e-3 + zeros(n, 1);
   smb = -1e-3 + zeros(n, 1);
   melt = 2e-3 + zeros(n, 1);
   Data = timetable(Time, tair, runoff, smb, melt);

   % Generate production-shaped whole-calendar RU/SMB and ME/MEH ledgers, then
   % retain only a mid-day three-day window without altering that metadata.
   replacements = struct('runoff', runoff, 'smb', smb);
   [Data, artifact_metadata] = ...
      icemodel.forcing.helpers.applyMarDailyQualityControl( ...
      Data, replacements, sector=1);
   artifact_metadata = icemodel.forcing.helpers.marDiagnosticMetadata( ...
      Data, melt, artifact_metadata, sector=1);
   artifact_metadata.mar_qc_replaced_runoff_count = 17;
   artifact_metadata.mar_qc_replaced_smb_count = 23;
   artifact_metadata.source_files = "full-mar-2012.nc";
   artifact_metadata.sentinel_policy = "preserve";
   keep = Time >= datetime(2012, 1, 10, 12, 0, 0, 'TimeZone', 'UTC') ...
      & Time <= datetime(2012, 1, 12, 11, 0, 0, 'TimeZone', 'UTC');
   Data = Data(keep, :);
   Data.Properties.UserData = artifact_metadata;

   % Derived met repeats each hourly rate over four interval-start samples while
   % retaining the same overlong source ledger, matching the writer boundary.
   met_time = (Data.Time(1):minutes(15):Data.Time(end) + minutes(45))';
   values = repelem(Data{:, ["tair", "runoff"]}, 4, 1);
   met = array2timetable(values, RowTimes=met_time, ...
      VariableNames={'tair', 'runoff'});
   met.Properties.UserData = artifact_metadata;
   sentinel_payload = magic(3);
   ensureParent(met_file);
   ensureParent(data_file);
   save(met_file, 'met', 'artifact_metadata', 'sentinel_payload')
   save(data_file, 'Data', 'artifact_metadata', 'sentinel_payload')
end

function writeTinyMarDailySource(mar_dir, yyyy)
   %WRITETINYMARDailySOURCE Create a two-cell MAR RU/SMB sector fixture.
   if ~isfolder(mar_dir)
      mkdir(mar_dir)
   end
   filename = fullfile(mar_dir, sprintf('MARv3.11-test-%d.nc', yyyy));
   [LAT, LON] = ndgrid([67, 68], [-48, -47]);
   static = zeros(2, 2);
   surface = 4 * ones(2, 2);

   % Coordinate variables share the LON dimension names because marGridInfo
   % reconstructs the regular native grid from those named 1-D axes.
   nccreate(filename, 'x', 'Dimensions', {'x', 2});
   nccreate(filename, 'y', 'Dimensions', {'y', 2});
   nccreate(filename, 'LON', 'Dimensions', {'x', 2, 'y', 2});
   nccreate(filename, 'LAT', 'Dimensions', {'x', 2, 'y', 2});
   nccreate(filename, 'SH', 'Dimensions', {'x', 2, 'y', 2});
   nccreate(filename, 'SLO', 'Dimensions', {'x', 2, 'y', 2});
   nccreate(filename, 'SRF', 'Dimensions', {'x', 2, 'y', 2});
   ncwrite(filename, 'x', [0; 15]);
   ncwrite(filename, 'y', [0; 15]);
   ncwrite(filename, 'LON', LON);
   ncwrite(filename, 'LAT', LAT);
   ncwrite(filename, 'SH', static);
   ncwrite(filename, 'SLO', static);
   ncwrite(filename, 'SRF', surface);

   % Sector 1 is permanent ice and sector 2 is tundra. Spatially uniform
   % values keep the test focused on category selection and UTC-day expansion.
   dimensions = {'x', 2, 'y', 2, 'SECTOR', 2, 'TIME', 2};
   nccreate(filename, 'RU', 'Dimensions', dimensions);
   nccreate(filename, 'SMB', 'Dimensions', dimensions);
   runoff = zeros(2, 2, 2, 2);
   runoff(:, :, 1, 1) = 24;
   runoff(:, :, 1, 2) = 48;
   runoff(:, :, 2, 1) = 240;
   runoff(:, :, 2, 2) = 480;
   ncwrite(filename, 'RU', runoff);
   ncwrite(filename, 'SMB', -runoff);
   ncwriteatt(filename, 'RU', 'units', 'mmWE/day');
   ncwriteatt(filename, 'SMB', 'units', 'mmWE/day');
end

function writeTinyGeusModis(modis_dir, yyyy, daily_values)
   %WRITETINYGEUSMODIS Create a two-cell/two-day GEUS-compatible NetCDF.
   if ~isfolder(modis_dir)
      mkdir(modis_dir)
   end
   filename = fullfile(modis_dir, sprintf( ...
      'Greenland_Reflectivity_%d_5km_C6.nc', yyyy));
   [X, Y] = ndgrid([0; 5000], [0, 5000]);
   [lat, lon] = projinv( ...
      icemodel.forcing.helpers.geusModisProjection(), X, Y);
   nccreate(filename, 'lat', 'Dimensions', {'x', 2, 'y', 2});
   nccreate(filename, 'lon', 'Dimensions', {'x', 2, 'y', 2});
   nccreate(filename, 'albedo', ...
      'Dimensions', {'x', 2, 'y', 2, 'time', numel(daily_values)});
   ncwrite(filename, 'lat', lat);
   ncwrite(filename, 'lon', lon);
   values = repmat(reshape(daily_values, 1, 1, []), 2, 2, 1);
   ncwrite(filename, 'albedo', values);
end

function writeBareExistingFiles(met_file, data_file)
   %WRITEBAREEXISTINGFILES Create paired RCM artifacts without new metadata.
   Data = timetable((datetime(2012, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:1))');
   Data.tair = [260; 261];
   Data.ppt = [0; 0];
   met = Data(:, "tair");

   ensureParent(met_file);
   ensureParent(data_file);
   save(met_file, 'met')
   save(data_file, 'Data')
end

function label = missingPromiceLabel(~, ~) %#ok<STOUT>
   %MISSINGPROMICELABEL Simulate an absent station during label resolution.
   error('icemodel:forcing:readPromiceAws:stationNotFound', ...
      'missing station')
end

function state = missingPromiceSource() %#ok<STOUT>
   %MISSINGPROMICESOURCE Simulate a missing PROMICE source root.
   error('icemodel:forcing:readPromiceAws:sourceNotFound', ...
      'source missing')
end

function state = stageShouldNotRun() %#ok<STOUT>
   %STAGESHOULDNOTRUN Guard that label failures skip before staging.
   error('test:stageShouldNotRun', 'stage callback should not run')
end

function state = emptyProtocolWindow() %#ok<STOUT>
   %EMPTYPROTOCOLWINDOW Simulate a RetMIP protocol file with no window rows.
   error('icemodel:verification:importRetmip:emptyProtocolWindow', ...
      'no rows in requested window')
end

function ensureParent(filename)
   %ENSUREPARENT Create a file's parent directory if needed.
   folder = fileparts(filename);
   if ~isfolder(folder)
      mkdir(folder)
   end
end

function bytes = fileBytes(filename)
   %FILEBYTES Read one binary artifact for no-mutation assertions.
   fid = fopen(filename, 'r');
   cleanup = onCleanup(@() fclose(fid));
   bytes = fread(fid, Inf, '*uint8');
end

function entry = helperEntry(s)
   %HELPERENTRY Convert a helper state record into a manifest case entry.
   [forcing_sources, eval_sources] = ...
      icemodel.verification.setup.colocationSourceLists(s.colocation);
   case_values = { ...
      char(s.case_id)
      'firn_observational'
      char(s.site_id)
      char(s.site_name)
      'unknown'
      {'firn'}
      'none'
      s.site_location
      s.period
      s.evaluation_file_rel
      cellstr(forcing_sources(:))
      cellstr(eval_sources(:))
      cellstr(s.comparison_variables(:))
      s.observation_variables
      s.colocation
      'daily'
      'helper test case'};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(case_values);
end

function L = fixedClippedLeg(~, ~)
   %FIXEDCLIPPEDLEG Return one attainable 2012 source leg for ranking tests.
   L = struct('staged', true, 'years', 2012, ...
      'start', datetime(2012, 1, 1, 'TimeZone', 'UTC'), ...
      'end', datetime(2012, 12, 31, 23, 0, 0, 'TimeZone', 'UTC'), ...
      'reason', "");
end

function L = skippedLeg(~, src)
   %SKIPPEDLEG Return a no-coverage leg for a requested source.
   L = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
      'reason', "missing " + src);
end

function sites = tinySiteCatalog()
   %TINYSITECATALOG Return ordered source-site rows for selector tests.
   sites = struct('site_id', {"A", "B", "C"});
end

function manifest_file = writePriorFamilyManifest(testCase)
   %WRITEPRIORFAMILYMANIFEST Write one canonical case for prior-load tests.
   manifest_file = fullfile(testCase.TestData.tmp, 'prior-family.json');
   cases = struct('case_id', "case1");
   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "helper", "", "", "test", string(datetime('today')), cases, ...
      struct('site', {}, 'reason', {}));
   icemodel.verification.setup.writeManifest(manifest_file, manifest);
end

function touch(filename)
   %TOUCH Create a tiny placeholder file.
   fid = fopen(filename, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, 'placeholder\n');
   clear cleaner
end

function state = malformedStage()
   %MALFORMEDSTAGE Throw a non-missing staging error.
   state = helperState();
   if isfield(state, 'case_id')
      error('test:malformedStage', 'malformed staging artifact')
   end
end
