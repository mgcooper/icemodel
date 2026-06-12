function tests = test_verification_source_cache
   %TEST_VERIFICATION_SOURCE_CACHE Verify the ESM-SnowMIP and Laugh-Tests
   %  source-cache helpers (icemodel.verification.setup.fetch*).
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
   src = icemodel.verification.setup.fetchEsmSnowmip( ...
      cache_dir=string(testCase.TestData.tmp), ...
      strict=false, silent=true);
   verifyClass(testCase, src, 'string');
   verifyTrue(testCase, exist(src, 'dir') == 7);
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
   src = icemodel.verification.setup.fetchLaughTests( ...
      cache_dir=string(testCase.TestData.tmp), ...
      strict=false, silent=true);
   verifyClass(testCase, src, 'string');
end
