function tests = test_verification_source_cache
   %TEST_VERIFICATION_SOURCE_CACHE Verify shared source-cache fetch contracts.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   testCase.TestData.tmp = tempname;
   icemodel.helpers.ensureDirExists(testCase.TestData.tmp);
end

function teardownOnce(testCase)
   if exist(testCase.TestData.tmp, 'dir') == 7
      rmdir(testCase.TestData.tmp, 's');
   end
end

function test_fetchEsmSnowmip_strict_errors_on_empty_cache(testCase)
   % strict=true must error with a stable error id when no source
   % files are present.
   verifyError(testCase, ...
      @() icemodel.verification.setup.fetchEsmSnowmip( ...
         cache_dir=string(testCase.TestData.tmp), ...
         strict=true, silent=true), ...
      'icemodel:verification:fetchEsmSnowmip:missingSources');
end

function test_fetchEsmSnowmip_nonstrict_returns_path(testCase)
   % strict=false returns the cache path even when the cache is
   % empty, so callers can decide how to handle the missing-data
   % state without try/catch.
   [src, status] = icemodel.verification.setup.fetchEsmSnowmip( ...
      cache_dir=string(testCase.TestData.tmp), ...
      strict=false, silent=true);
   verifyClass(testCase, src, 'string');
   verifyTrue(testCase, exist(src, 'dir') == 7);
   verifyNotEmpty(testCase, status);
   verifyTrue(testCase, all(~[status.present]));
end

function test_fetchEsmSnowmip_creates_missing_directory(testCase)
   % If the cache directory does not exist yet, the helper creates
   % it (so the user can drop files in afterwards).
   subdir = fullfile(testCase.TestData.tmp, 'esm_snowmip_fresh');
   verifyTrue(testCase, exist(subdir, 'dir') ~= 7);
   icemodel.verification.setup.fetchEsmSnowmip( ...
      cache_dir=string(subdir), strict=false, silent=true);
   verifyTrue(testCase, exist(subdir, 'dir') == 7);
end

function test_fetchLaughTests_strict_errors_on_empty_cache(testCase)
   verifyError(testCase, ...
      @() icemodel.verification.setup.fetchLaughTests( ...
         cache_dir=string(testCase.TestData.tmp), ...
         strict=true, silent=true), ...
      'icemodel:verification:fetchLaughTests:missingSources');
end

function test_fetchLaughTests_nonstrict_returns_path(testCase)
   [src, status] = icemodel.verification.setup.fetchLaughTests( ...
      cache_dir=string(testCase.TestData.tmp), ...
      strict=false, silent=true);
   verifyClass(testCase, src, 'string');
   verifyNotEmpty(testCase, status);
end

function test_fetchers_can_validate_without_creating_cache(testCase)
   % All manual-cache fetchers expose the same non-mutating resolution option.
   esm = fullfile(testCase.TestData.tmp, 'esm-no-create');
   laugh = fullfile(testCase.TestData.tmp, 'laugh-no-create');
   icemodel.verification.setup.fetchEsmSnowmip( ...
      cache_dir=esm, strict=false, silent=true, create_cache_dir=false);
   icemodel.verification.setup.fetchLaughTests( ...
      cache_dir=laugh, strict=false, silent=true, create_cache_dir=false);
   verifyFalse(testCase, isfolder(esm));
   verifyFalse(testCase, isfolder(laugh));
end

function test_empty_fetch_selections_are_side_effect_free(testCase)
   % Empty selectors keep each public status schema but never create caches.
   retmip = fullfile(testCase.TestData.tmp, 'retmip-empty');
   imau = fullfile(testCase.TestData.tmp, 'imau-empty');
   gcnet = fullfile(testCase.TestData.tmp, 'gcnet-empty');
   esm = fullfile(testCase.TestData.tmp, 'esm-empty');
   cache_paths = string({retmip, imau, gcnet, esm});

   [retmip_dir, retmip_status] = ...
      icemodel.verification.setup.fetchRetmip( ...
      cache_dir=retmip, products=strings(1, 0), silent=true);
   [imau_dir, imau_status] = icemodel.verification.setup.fetchImau( ...
      cache_dir=imau, products=strings(1, 0), silent=true);
   [gcnet_dir, gcnet_status] = icemodel.verification.setup.fetchGcnet( ...
      cache_dir=gcnet, products=strings(1, 0), silent=true);
   [esm_dir, esm_status] = icemodel.verification.setup.fetchEsmSnowmip( ...
      cache_dir=esm, stations=strings(1, 0), silent=true);

   verifyEqual(testCase, string([retmip_dir, imau_dir, gcnet_dir, esm_dir]), ...
      cache_paths);
   verifySize(testCase, retmip_status, [1, 0]);
   verifySize(testCase, imau_status, [1, 0]);
   verifySize(testCase, gcnet_status, [1, 0]);
   verifySize(testCase, esm_status, [1, 0]);
   verifyEqual(testCase, string(fieldnames(retmip_status)), [ ...
      "product"; "doi"; "landing_url"; "present"; "cache_dir"; ...
      "found_files"; "dataverse_api_url"]);
   verifyEqual(testCase, string(fieldnames(imau_status)), [ ...
      "product"; "doi"; "landing_url"; "present"; "cache_dir"; ...
      "found_files"]);
   verifyFalse(testCase, any(isfolder(cache_paths)));
end

function test_fetch_product_validation_precedes_cache_creation(testCase)
   % Canonical product registries reject typos before filesystem side effects.
   retmip = fullfile(testCase.TestData.tmp, 'retmip-invalid');
   imau = fullfile(testCase.TestData.tmp, 'imau-invalid');

   verifyError(testCase, @() icemodel.verification.setup.fetchRetmip( ...
      cache_dir=retmip, products="unknown", silent=true), ...
      'icemodel:verification:fetchRetmip:unknownProduct');
   verifyError(testCase, @() icemodel.verification.setup.fetchImau( ...
      cache_dir=imau, products="unknown", silent=true), ...
      'icemodel:verification:fetchImau:unknownProduct');

   verifyFalse(testCase, isfolder(retmip));
   verifyFalse(testCase, isfolder(imau));
end

function test_registry_fetchers_preserve_family_specific_success_rows(testCase)
   % Shared registry control flow must retain family-specific file matching and
   % the richer RetMIP Dataverse provenance field on successful cache probes.
   retmip = fullfile(testCase.TestData.tmp, 'retmip-present');
   imau = fullfile(testCase.TestData.tmp, 'imau-present');
   mkdir(fullfile(retmip, 'forcing'));
   mkdir(fullfile(imau, 'hourly'));
   writelines("fixture", fullfile(retmip, 'forcing', 'protocol.tab'));
   writelines("fixture", fullfile(imau, 'hourly', 'S21.tab'));

   [retmip_dir, retmip_status] = ...
      icemodel.verification.setup.fetchRetmip( ...
      cache_dir=retmip, products="forcing", strict=true, silent=true);
   [imau_dir, imau_status] = icemodel.verification.setup.fetchImau( ...
      cache_dir=imau, products="hourly", strict=true, silent=true);

   verifyEqual(testCase, string([retmip_dir, imau_dir]), ...
      string({retmip, imau}));
   verifyTrue(testCase, retmip_status.present);
   verifyTrue(testCase, imau_status.present);
   verifyEqual(testCase, string(retmip_status.product), "forcing");
   verifyEqual(testCase, string(imau_status.product), "hourly");
   verifyTrue(testCase, contains(retmip_status.dataverse_api_url, ...
      "10.22008/FK2/GZ3CSN"));
   verifyEqual(testCase, string(imau_status.doi), ...
      "10.1594/PANGAEA.971647");
end

function test_optional_manual_cache_fetches_are_quiet(testCase)
   % Optional missing caches return status without printing retrieval banners.
   retmip = @() icemodel.verification.setup.fetchRetmip( ...
      cache_dir=fullfile(testCase.TestData.tmp, 'retmip-quiet'), ...
      strict=false, silent=true); %#ok<NASGU>
   imau = @() icemodel.verification.setup.fetchImau( ...
      cache_dir=fullfile(testCase.TestData.tmp, 'imau-quiet'), ...
      strict=false, silent=true); %#ok<NASGU>

   verifyEqual(testCase, string(evalc('retmip();')), "");
   verifyEqual(testCase, string(evalc('imau();')), "");
end

function test_strict_manual_cache_fetches_print_doi_guidance(testCase)
   % Required missing caches fail after printing actionable DOI provenance.
   retmip = @() icemodel.verification.setup.fetchRetmip( ...
      cache_dir=fullfile(testCase.TestData.tmp, 'retmip-strict'), ...
      strict=true, silent=false); %#ok<NASGU>
   imau = @() icemodel.verification.setup.fetchImau( ...
      cache_dir=fullfile(testCase.TestData.tmp, 'imau-strict'), ...
      strict=true, silent=false); %#ok<NASGU>

   retmip_output = evalc([ ...
      'testCase.verifyError(retmip, ' ...
      '''icemodel:verification:fetchRetmip:missingSources'');']);
   imau_output = evalc([ ...
      'testCase.verifyError(imau, ' ...
      '''icemodel:verification:fetchImau:missingSources'');']);

   verifyTrue(testCase, contains(retmip_output, "10.22008/FK2/GZ3CSN"));
   verifyTrue(testCase, contains(retmip_output, "10.22008/FK2/CVPUJL"));
   verifyTrue(testCase, contains(imau_output, "10.1594/PANGAEA.971647"));
   verifyTrue(testCase, contains(imau_output, "10.1594/PANGAEA.970127"));
end

function test_gcnet_product_registry_drives_spec_and_validation(testCase)
   % Every canonical selector resolves metadata; unknown input fails before mkdir.
   products = icemodel.verification.setup.gcnetProductNames();
   for product = products
      spec = icemodel.verification.setup.gcnetProductSpec(product);
      verifyEqual(testCase, string(spec.product), product);
      verifyNotEmpty(testCase, string(spec.doi));
      verifyNotEmpty(testCase, string(spec.station_suffixes));
      verifyEqual(testCase, string(spec.data_patterns), ...
         ["*.nc", fullfile("**", "*.nc")]);
      % Simulated state files add a numeric bin id; observed products name the
      % complete file directly, and the registry must preserve that distinction.
      if product == "simulated_firn"
         verifyEqual(testCase, string(spec.station_file_mode), "numbered");
      else
         verifyEqual(testCase, string(spec.station_file_mode), "exact");
      end
   end

   cache_dir = fullfile(testCase.TestData.tmp, 'gcnet-invalid');
   verifyError(testCase, @() icemodel.verification.setup.fetchGcnet( ...
      cache_dir=cache_dir, products="unknown", silent=true), ...
      'icemodel:validators:mustBeGcnetProductSelection');
   verifyFalse(testCase, isfolder(cache_dir));
end

function test_atomic_importers_accept_omitted_scalar_source_dir(testCase)
   % Dry-run atomic imports resolve their default scalar source path themselves.
   eval_root = fullfile(testCase.TestData.tmp, 'dry-eval');
   esm = icemodel.verification.setup.importEsmSnowmip( ...
      dry_run=true, evaluation_data_root=eval_root, case_ids="cdp");
   laugh = icemodel.verification.setup.importLaughTests( ...
      dry_run=true, evaluation_data_root=eval_root);
   verifyEqual(testCase, string(esm.dataset_family), "esm_snowmip");
   verifyEqual(testCase, string(laugh.dataset_family), "laugh_tests");
   verifyFalse(testCase, isfolder(eval_root));
end

function test_atomic_importer_rejects_invalid_met_cadence_before_writes(testCase)
   % Shared setup wrappers accept only the public 15m or explicit-native token.
   output_root = fullfile(testCase.TestData.tmp, 'bad-esm-cadence');
   verifyError(testCase, @() ...
      icemodel.verification.setup.importEsmSnowmip( ...
      "", case_ids="cdp", output_root=output_root, dry_run=true, ...
      dt_out="1h"), 'MATLAB:validators:mustBeMember');
   verifyFalse(testCase, isfolder(output_root));
end

function test_buildLaughTestsArtifacts_rejects_unknown_case(testCase)
   % The extracted builder owns case-specific dispatch and a stable error id.
   verifyError(testCase, @() ...
      icemodel.verification.setup.buildLaughTestsArtifacts( ...
      string(testCase.TestData.tmp), "unknown"), ...
      'icemodel:verification:buildLaughTestsArtifacts:unsupportedCase');
end
