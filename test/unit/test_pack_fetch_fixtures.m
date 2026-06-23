function tests = test_pack_fetch_fixtures
   %TEST_PACK_FETCH_FIXTURES Roundtrip the fixture pack/fetch release tooling.
   %
   % Exercises icemodel.verification.setup.packFixtures + fetchFixtures against a
   % synthetic temporary demo-data root (no network, no committed bundle):
   %   - packFixtures bundles the fixture .mat set + a SHA-256 manifest;
   %   - fetchFixtures verifies an intact tree against that manifest;
   %   - extracting the archive into an empty root restores the files and
   %     re-verifies (sha256 matches);
   %   - a tampered / missing fixture is detected (mismatched / missing);
   %   - with no manifest the helper is a presence-only no-op, and empty roots
   %     surface the download banner.
   % Mirrors the fetchSumup / fetchEsmSnowmip verify-cache contract.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;

   % Build a self-contained fake demo-data root holding the three fixture-data
   % subtrees the real tree has, so the roundtrip never touches the committed
   % fixtures and never writes into the repo.
   root = string(tempname);
   testCase.TestData.root = root;
   makeFakeFixture(fullfile(root, "eval", "fam", "case1", "observations.mat"));
   makeFakeFixture(fullfile(root, "eval", "fam", "case2", "reference.mat"));
   makeFakeFixture(fullfile(root, "input", "met", "met_a_b_2016_1hr.mat"));
   makeFakeFixture(fullfile(root, "input", "userdata", "a_b_2016.mat"));

   % A lean companion that must NOT be bundled (manifests stay committed).
   mkfile(fullfile(root, "eval", "fam", "manifest.json"), '{"k":1}');

   % Gitignored staging dir for the bundle, isolated per test.
   testCase.TestData.staging = string(tempname);
   testCase.TestData.version = "v0.0.0-test";
end

function teardown(testCase)
   for f = ["root", "staging"]
      p = testCase.TestData.(f);
      if isfolder(p)
         rmdir(p, 's')
      end
   end
   clear testCase.TestData.cleanup
end

function test_filelist_excludes_manifest_and_lists_data(testCase)
   % The single-source file list enumerates the .mat data only, not the lean
   % committed manifest.json companion.
   returned = icemodel.verification.setup.fixtureFileList( ...
      root=testCase.TestData.root);

   expected = sort([ ...
      "eval/fam/case1/observations.mat"; ...
      "eval/fam/case2/reference.mat"; ...
      "input/met/met_a_b_2016_1hr.mat"; ...
      "input/userdata/a_b_2016.mat"]);
   testCase.verifyEqual(returned, expected);
end

function test_pack_writes_archive_and_manifest(testCase)
   % packFixtures produces both the .tar.gz and the SHA-256 MANIFEST.json in the
   % gitignored staging dir, and reports a saving.
   result = icemodel.verification.setup.packFixtures( ...
      testCase.TestData.version, root=testCase.TestData.root, ...
      staging_dir=testCase.TestData.staging, silent=true);

   testCase.verifyTrue(isfile(result.archive_file));
   testCase.verifyTrue(isfile(result.manifest_file));
   testCase.verifyEqual(result.num_files, 4);
   testCase.verifyEqual(result.manifest.version, testCase.TestData.version);
   testCase.verifyTrue(result.source_bytes > 0);
   % Each manifest entry carries a 64-char hex sha256.
   for k = 1:numel(result.manifest.files)
      testCase.verifyEqual(strlength(result.manifest.files(k).sha256), 64);
   end
end

function test_pack_refuses_clobber_without_overwrite(testCase)
   % A second pack of the same version errors unless overwrite=true.
   icemodel.verification.setup.packFixtures(testCase.TestData.version, ...
      root=testCase.TestData.root, staging_dir=testCase.TestData.staging, ...
      silent=true);

   testCase.verifyError(@() icemodel.verification.setup.packFixtures( ...
      testCase.TestData.version, root=testCase.TestData.root, ...
      staging_dir=testCase.TestData.staging, silent=true), ...
      'icemodel:verification:packFixtures:exists');

   % overwrite=true succeeds.
   testCase.verifyWarningFree(@() icemodel.verification.setup.packFixtures( ...
      testCase.TestData.version, root=testCase.TestData.root, ...
      staging_dir=testCase.TestData.staging, overwrite=true, silent=true));
end

function test_fetch_verifies_intact_tree_against_manifest(testCase)
   % An intact demo tree verifies "verified" against its own bundle manifest.
   result = icemodel.verification.setup.packFixtures(testCase.TestData.version, ...
      root=testCase.TestData.root, staging_dir=testCase.TestData.staging, ...
      silent=true);

   returned = icemodel.verification.setup.fetchFixtures( ...
      testCase.TestData.version, root=testCase.TestData.root, ...
      manifest=result.manifest_file, strict=true, silent=true);

   testCase.verifyTrue(returned.ok);
   testCase.verifyEqual(returned.mode, "verified");
   testCase.verifyEmpty(returned.missing);
   testCase.verifyEmpty(returned.mismatched);
end

function test_roundtrip_extract_into_empty_root_matches_sha256(testCase)
   % Pack from the fixture root, extract the archive into a fresh empty root,
   % and verify every file's SHA-256 matches the manifest. This is the
   % network-free analogue of "download the asset and provision".
   result = icemodel.verification.setup.packFixtures(testCase.TestData.version, ...
      root=testCase.TestData.root, staging_dir=testCase.TestData.staging, ...
      silent=true);

   empty_root = string(tempname);
   cleaner = onCleanup(@() rmIfPresent(empty_root));

   returned = icemodel.verification.setup.fetchFixtures(testCase.TestData.version, ...
      root=empty_root, archive=result.archive_file, extract=true, ...
      manifest=result.manifest_file, strict=true, silent=true);

   testCase.verifyTrue(returned.ok);
   testCase.verifyEqual(returned.mode, "verified");
   % The extracted files exist with the original relative layout.
   testCase.verifyTrue(isfile(fullfile(empty_root, ...
      "input", "met", "met_a_b_2016_1hr.mat")));
end

function test_fetch_detects_tampered_fixture(testCase)
   % Mutating a fixture after packing must be caught as a SHA-256 mismatch.
   result = icemodel.verification.setup.packFixtures(testCase.TestData.version, ...
      root=testCase.TestData.root, staging_dir=testCase.TestData.staging, ...
      silent=true);

   % Append bytes to one fixture so its hash changes.
   target = fullfile(testCase.TestData.root, "input", "met", ...
      "met_a_b_2016_1hr.mat");
   fid = fopen(target, 'a');
   fwrite(fid, uint8(1:64), 'uint8');
   fclose(fid);

   returned = icemodel.verification.setup.fetchFixtures(testCase.TestData.version, ...
      root=testCase.TestData.root, manifest=result.manifest_file, ...
      strict=false, silent=true);

   testCase.verifyFalse(returned.ok);
   testCase.verifyEqual(returned.mode, "missing");
   testCase.verifyTrue(any(returned.mismatched == "input/met/met_a_b_2016_1hr.mat"));
end

function test_fetch_strict_errors_on_missing_fixture(testCase)
   % Deleting a manifest-listed fixture must error under strict verification.
   result = icemodel.verification.setup.packFixtures(testCase.TestData.version, ...
      root=testCase.TestData.root, staging_dir=testCase.TestData.staging, ...
      silent=true);

   delete(fullfile(testCase.TestData.root, "input", "userdata", "a_b_2016.mat"));

   testCase.verifyError(@() icemodel.verification.setup.fetchFixtures( ...
      testCase.TestData.version, root=testCase.TestData.root, ...
      manifest=result.manifest_file, strict=true, silent=true), ...
      'icemodel:verification:fetchFixtures:incompleteFixtures');
end

function test_fetch_no_manifest_is_present_noop(testCase)
   % With no manifest and the fixtures present, fetchFixtures is a no-op: this
   % is the non-breaking bootstrap path used today.
   returned = icemodel.verification.setup.fetchFixtures( ...
      root=testCase.TestData.root, strict=true, silent=true);

   testCase.verifyTrue(returned.ok);
   testCase.verifyEqual(returned.mode, "present");
end

function test_fetch_empty_root_prints_download_banner(testCase)
   % An empty root with no manifest surfaces the actionable download banner and
   % errors under strict.
   empty_root = string(tempname);
   mkdir(empty_root)
   cleaner = onCleanup(@() rmIfPresent(empty_root));

   output = evalc(['try, icemodel.verification.setup.fetchFixtures(' ...
      '"v0.0.0-test", root="' char(empty_root) '", strict=false, ' ...
      'silent=false); catch, end']);

   testCase.verifyTrue(contains(output, "demo fixture data incomplete"));
   testCase.verifyTrue(contains(output, "gh release download"));
   testCase.verifyTrue(contains(output, "icemodel-fixtures-v0.0.0-test.tar.gz"));
end

%% Local fixture helpers
function makeFakeFixture(pathname)
   %MAKEFAKEFIXTURE Save a tiny but real .mat so dir()/sha256 see real bytes.
   folder = fileparts(pathname);
   if ~isfolder(folder)
      mkdir(folder)
   end
   payload = rand(8);
   save(char(pathname), 'payload');
end

function mkfile(pathname, text)
   %MKFILE Write a small text companion file.
   folder = fileparts(pathname);
   if ~isfolder(folder)
      mkdir(folder)
   end
   fid = fopen(pathname, 'w');
   fwrite(fid, text, 'char');
   fclose(fid);
end

function rmIfPresent(pathname)
   %RMIFPRESENT Remove a temp directory tree if it exists.
   if isfolder(pathname)
      rmdir(pathname, 's')
   end
end
