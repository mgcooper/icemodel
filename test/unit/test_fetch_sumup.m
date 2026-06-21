function tests = test_fetch_sumup
   %TEST_FETCH_SUMUP Verify fetchSumup's missing-cache messaging and contract.
   %
   % SUMup is access-gated (NASA Earthdata, NSIDC G02288), so the real data
   % is never on disk in CI. This suite exercises fetchSumup against an empty
   % cache: strict=true must error with the missing-sources id, strict=false
   % must return the cache dir, and the helper must create the cache dir so a
   % user can drop files into a path that already exists. Mirrors the
   % fetchEsmSnowmip retrieval-banner contract.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.cache = tempname;
end

function teardown(testCase)
   if isfolder(testCase.TestData.cache)
      rmdir(testCase.TestData.cache, 's')
   end
   clear testCase.TestData.cleanup
end

function test_strict_errors_on_empty_cache(testCase)
   % strict=true (default) must raise the missing-sources error on an empty
   % cache rather than returning silently.

   testCase.verifyError(@() icemodel.verification.setup.fetchSumup( ...
      cache_dir=testCase.TestData.cache, strict=true, silent=true), ...
      'icemodel:verification:fetchSumup:missingSources');
end

function test_nonstrict_returns_cache_dir(testCase)
   % strict=false must return the (partial) cache dir instead of erroring.

   source_dir = icemodel.verification.setup.fetchSumup( ...
      cache_dir=testCase.TestData.cache, strict=false, silent=true);

   testCase.verifyEqual(string(source_dir), string(testCase.TestData.cache));
end

function test_creates_cache_directory(testCase)
   % The helper must create the cache directory so the retrieval banner points
   % at a path that already exists.

   testCase.assertFalse(isfolder(testCase.TestData.cache));
   icemodel.verification.setup.fetchSumup( ...
      cache_dir=testCase.TestData.cache, strict=false, silent=true);
   testCase.verifyTrue(isfolder(testCase.TestData.cache));
end

function test_missing_cache_prints_retrieval_instructions(testCase)
   % The non-silent retrieval banner must name the SUMup DOI, NSIDC URL, and
   % the expected 2025 release filenames so a developer can act on it.

   output = evalc(['icemodel.verification.setup.fetchSumup(' ...
      'cache_dir="' testCase.TestData.cache '", strict=false, silent=false);']);

   testCase.verifyTrue(contains(output, "SUMup firn source cache incomplete"));
   testCase.verifyTrue(contains(output, "10.18739/A2M61BR5M"));
   testCase.verifyTrue(contains(output, "nsidc.org/data/g02288"));
   testCase.verifyTrue(contains(output, "SUMup_2025_density_greenland.nc"));
end

function test_partial_cache_still_flags_missing_variable(testCase)
   % A cache with only one variable file must still flag the other required
   % SUMup variable groups as missing.

   mkdir(testCase.TestData.cache)
   % Drop a density file only; SMB + temperature remain missing. The real
   % 2025 release ships NetCDF per ice sheet (SUMup_2025_<var>_greenland.nc).
   fid = fopen(fullfile(testCase.TestData.cache, ...
      'SUMup_2025_density_greenland.nc'), 'w');
   fclose(fid);

   output = evalc(['icemodel.verification.setup.fetchSumup(' ...
      'cache_dir="' testCase.TestData.cache '", strict=false, silent=false);']);

   testCase.verifyTrue(contains(output, "SMB"));
   testCase.verifyTrue(contains(output, "temperature"));
   testCase.verifyFalse(contains(output, "*density*greenland*"));
end

function test_full_cache_passes(testCase)
   % A cache holding all three Greenland variable files satisfies the check;
   % strict=true returns the cache dir without erroring.

   mkdir(testCase.TestData.cache)
   for v = ["density", "SMB", "temperature"]
      fid = fopen(fullfile(testCase.TestData.cache, ...
         char("SUMup_2025_" + v + "_greenland.nc")), 'w');
      fclose(fid);
   end

   source_dir = icemodel.verification.setup.fetchSumup( ...
      cache_dir=testCase.TestData.cache, strict=true, silent=true);
   testCase.verifyEqual(string(source_dir), string(testCase.TestData.cache));
end
