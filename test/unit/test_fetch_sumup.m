function tests = test_fetch_sumup
   %TEST_FETCH_SUMUP Verify fetchSumup's missing-cache messaging and contract.
   %
   % SUMup is access-gated (NASA Earthdata, NSIDC G02288), so the real data
   % is never on disk in CI. This suite exercises fetchSumup against an empty
   % cache: strict=true must error with the missing-sources id, strict=false
   % must return the cache dir, and the helper must create the cache dir so a
   % user can drop files into a path that already exists. Empty and invalid
   % variable selectors also exercise the shared non-mutating fetch contract.
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
   % A strict failure must print actionable guidance even when silent=true.

   output = evalc(['testCase.verifyError(@() ' ...
      'icemodel.verification.setup.fetchSumup(' ...
      'cache_dir="' testCase.TestData.cache '", strict=true, silent=true), ' ...
      '''icemodel:verification:fetchSumup:missingSources'');']);

   testCase.verifyTrue(contains(output, "SUMup firn source cache incomplete"));
   testCase.verifyTrue(contains(output, "10.18739/A2M61BR5M"));
end

function test_nonstrict_returns_cache_dir(testCase)
   % A silent optional probe returns partial status without printing guidance.

   output = evalc(['[source_dir, status] = ' ...
      'icemodel.verification.setup.fetchSumup(' ...
      'cache_dir="' testCase.TestData.cache '", strict=false, silent=true);']);

   testCase.verifyEqual(string(source_dir), string(testCase.TestData.cache));
   testCase.verifyNotEmpty(status);
   testCase.verifyTrue(all(~[status.present]));
   testCase.verifyEmpty(strtrim(output));
end

function test_blank_cache_dir_uses_default_root(testCase)
   % Shared fetch option plumbing treats blank cache_dir the same as omitted.
   resolved = icemodel.verification.setup.resolveFetchCacheDir("", ...
      testCase.TestData.cache);

   testCase.verifyEqual(string(resolved), string(testCase.TestData.cache));
end

function test_finish_fetch_status_accepts_column_missing_list(testCase)
   % Missing collectors may return column vectors; the shared status adapter
   % keeps the fetch completion contract independent of row/column shape.
   status = icemodel.verification.setup.fetchMissingStatus(["density"; "SMB"]);
   source_dir = icemodel.verification.setup.finishFetchStatus( ...
      string(testCase.TestData.cache), status, strict=false, silent=true);

   testCase.verifyEqual(string(source_dir), string(testCase.TestData.cache));
end

function test_finish_fetch_status_strict_without_callback_errors(testCase)
   % A strict shared fetch check must fail even when a caller has no
   % dataset-specific error callback.
   status = icemodel.verification.setup.fetchMissingStatus("density");

   testCase.verifyError(@() icemodel.verification.setup.finishFetchStatus( ...
      string(testCase.TestData.cache), status, strict=true, silent=true), ...
      'icemodel:verification:finishFetchStatus:missingSources');
end

function test_creates_cache_directory(testCase)
   % The helper must create the cache directory so the retrieval banner points
   % at a path that already exists.

   testCase.assertFalse(isfolder(testCase.TestData.cache));
   icemodel.verification.setup.fetchSumup( ...
      cache_dir=testCase.TestData.cache, strict=false, silent=true);
   testCase.verifyTrue(isfolder(testCase.TestData.cache));
end

function test_can_validate_without_creating_cache_directory(testCase)
   % create_cache_dir=false keeps dry-run/source-light callers non-mutating.
   testCase.assertFalse(isfolder(testCase.TestData.cache));
   icemodel.verification.setup.fetchSumup( ...
      cache_dir=testCase.TestData.cache, strict=false, silent=true, ...
      create_cache_dir=false);
   testCase.verifyFalse(isfolder(testCase.TestData.cache));
end

function test_empty_variable_selection_is_side_effect_free(testCase)
   % An empty batch returns the resolved path and status without cache creation.
   [source_dir, status] = icemodel.verification.setup.fetchSumup( ...
      cache_dir=testCase.TestData.cache, variables=strings(1, 0), ...
      silent=true);

   testCase.verifyEqual(string(source_dir), string(testCase.TestData.cache));
   testCase.verifySize(status, [1, 0]);
   testCase.verifyEqual(string(fieldnames(status)), ["product"; "present"]);
   testCase.verifyFalse(isfolder(testCase.TestData.cache));
end

function test_unknown_variable_rejects_before_cache_creation(testCase)
   % The canonical variable registry rejects typos before mkdir or globbing.
   testCase.verifyError(@() icemodel.verification.setup.fetchSumup( ...
      cache_dir=testCase.TestData.cache, variables="denisty", silent=true), ...
      'icemodel:verification:fetchSumup:unknownVariable');

   testCase.verifyFalse(isfolder(testCase.TestData.cache));
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
