function tests = test_pack_fetch_fixtures
   %TEST_PACK_FETCH_FIXTURES Validate v1.1 release-data pack/fetch transactions.
   %
   % Uses tiny temporary capability trees and local archives. No test touches
   % canonical data, publishes an asset, or reaches an external network.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   %SETUP Install repo paths and build one deterministic synthetic data root.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.root = canonicalTempname();
   testCase.TestData.staging = canonicalTempname();
   mkdir(testCase.TestData.root)
   mkdir(testCase.TestData.staging)

   % Give formal-core two parents so rollback can be forced after one promote.
   mkfile(fullfile(testCase.TestData.root, "a", "formal1.mat"), "formal one");
   mkfile(fullfile(testCase.TestData.root, "b", "formal2.mat"), "formal two");
   mkfile(fullfile(testCase.TestData.root, "eval", "showcase", ...
      "manifest.json"), '{"family":"showcase"}');
   mkfile(fullfile(testCase.TestData.root, "eval", "showcase", ...
      "reference.mat"), "showcase reference");
   mkfile(fullfile(testCase.TestData.root, "forcing", "source.nc"), ...
      "optional forcing");

   % The source manifest has complete file metadata and deliberately unstaged
   % archive metadata; packFixtures fills the latter in its output manifest.
   testCase.TestData.manifest = fullfile(testCase.TestData.staging, ...
      "source-manifest.json");
   writeManifest(testCase.TestData.manifest, ...
      sourceManifest(testCase.TestData.root));
   testCase.TestData.version = "v1.1";
end

function teardown(testCase)
   %TEARDOWN Remove only this test's temporary trees and restore configuration.
   for name = ["root", "staging"]
      pathname = testCase.TestData.(name);
      if isfolder(pathname)
         rmdir(pathname, 's')
      end
   end
   clear testCase.TestData.cleanup
end

function test_tracked_manifest_declares_v11_boundary(testCase)
   % The tracked manifest names all three release archives and required rows.
   [files, selection] = icemodel.verification.setup.fixtureFileList();
   testCase.verifyEqual(selection.version, "v1.1");
   testCase.verifyEqual(selection.capabilities, ...
      ["formal-core"; "verification-showcase"]);
   testCase.verifyEqual(string({selection.archives.name})', [ ...
      "icemodel-v1.1-formal-core.tar.gz"; ...
      "icemodel-v1.1-verification-showcase.tar.gz"]);
   testCase.verifyTrue(all([selection.archives.bytes] > 0));
   testCase.verifyTrue(all(strlength( ...
      string({selection.archives.sha256})) == 64));
   testCase.verifyEqual(numel(files), 11);

   % Optional forcing is a complete, separately selectable source-fixture set.
   [optional_files, optional] = ...
      icemodel.verification.setup.fixtureFileList( ...
      capabilities="forcing-integration");
   testCase.verifyNumElements(optional_files, 141);
   testCase.verifyTrue(all(startsWith(optional_files, "forcing/")));
   testCase.verifyFalse(any(contains(optional_files, ".DS_Store")));
   testCase.verifyFalse(optional.archives.required);
   testCase.verifyEqual(string(optional.archives.name), ...
      "icemodel-v1.1-forcing-integration.tar.gz");
   testCase.verifyGreaterThan(optional.archives.bytes, 0);
   testCase.verifyEqual(strlength(string(optional.archives.sha256)), 64);
   testCase.verifyFalse(any([optional.files.required]));
end

function test_capability_order_is_canonical_for_pack_and_fetch(testCase)
   % Caller order and duplicates must not alter manifest or archive mapping.
   requested = [ ...
      "forcing-integration", "formal-core", "verification-showcase", ...
      "forcing-integration"];
   packed = icemodel.verification.setup.packFixtures("v1.1", ...
      capabilities=requested, root=testCase.TestData.root, ...
      manifest=testCase.TestData.manifest, ...
      staging_dir=testCase.TestData.staging, silent=true);

   % Packing emits one canonical capability/archive order.
   expected = [ ...
      "formal-core"; "verification-showcase"; "forcing-integration"];
   testCase.verifyEqual(packed.capabilities, expected);
   testCase.verifyEqual(string({packed.manifest.archives.capability})', ...
      expected);

   % Fetching the matching packed vector remains valid even when the caller
   % repeats the original noncanonical request.
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   fetched = icemodel.verification.setup.fetchFixtures("v1.1", ...
      capabilities=requested, root=destination, ...
      manifest=packed.manifest_file, archive=packed.archive_files, silent=true);
   testCase.verifyTrue(fetched.ok);
   testCase.verifyEqual(fetched.capabilities, expected);
   clear cleaner
end

function test_manifest_rejects_noncanonical_file_order(testCase)
   % Selection must not silently rewrite authoritative manifest file ordering.
   manifest = sourceManifest(testCase.TestData.root);
   manifest.files = flip(manifest.files);
   invalid = fullfile(testCase.TestData.staging, ...
      "noncanonical-file-order.json");
   writeManifest(invalid, manifest);

   testCase.verifyError(@() ...
      icemodel.verification.setup.fixtureFileList( ...
      root=testCase.TestData.root, manifest=invalid), ...
      'icemodel:verification:fixtureFileList:invalidManifest');
end

function test_filelist_default_and_explicit_selection(testCase)
   % Default selection returns only required rows; optional is explicit.
   files = icemodel.verification.setup.fixtureFileList( ...
      root=testCase.TestData.root, manifest=testCase.TestData.manifest);
   testCase.verifyEqual(files, sort([ ...
      "a/formal1.mat"; "b/formal2.mat"; ...
      "eval/showcase/manifest.json"; ...
      "eval/showcase/reference.mat"]));

   optional = icemodel.verification.setup.fixtureFileList( ...
      capabilities="forcing-integration", root=testCase.TestData.root, ...
      manifest=testCase.TestData.manifest);
   testCase.verifyEqual(optional, "forcing/source.nc");
end

function test_relative_posix_normalizes_windows_paths(testCase)
   % A single Windows separator normalizes identically on every test host.
   backslash = string(char(92));
   root = "C:" + backslash + "fixtures" + backslash;
   pathname = root + "input" + backslash + "met" + backslash + "file.mat";
   testCase.verifyEqual( ...
      icemodel.verification.setup.fixtureRelativePosix(root, pathname), ...
      "input/met/file.mat");
end

function test_manifest_rejects_unsafe_and_inconsistent_rows(testCase)
   % Manifest-controlled paths, names, tags, and capability policy are trusted input.
   manifest = sourceManifest(testCase.TestData.root);
   manifest.files(1).path = "../escape.mat";
   writeManifest(testCase.TestData.manifest, manifest);
   testCase.verifyError(@() syntheticFileList(testCase), ...
      'icemodel:verification:fixtureFileList:unsafePath');

   % A single Windows separator must not disguise parent traversal.
   manifest = sourceManifest(testCase.TestData.root);
   manifest.files(1).path = ".." + string(char(92)) + "escape.mat";
   writeManifest(testCase.TestData.manifest, manifest);
   testCase.verifyError(@() syntheticFileList(testCase), ...
      'icemodel:verification:fixtureFileList:unsafePath');

   manifest = sourceManifest(testCase.TestData.root);
   manifest.archives(1).name = "../formal-core.tar.gz";
   writeManifest(testCase.TestData.manifest, manifest);
   testCase.verifyError(@() syntheticFileList(testCase), ...
      'icemodel:verification:fixtureFileList:invalidManifest');

   manifest = sourceManifest(testCase.TestData.root);
   manifest.version = "../v1.1";
   writeManifest(testCase.TestData.manifest, manifest);
   testCase.verifyError(@() syntheticFileList(testCase), ...
      'icemodel:verification:fixtureFileList:invalidManifest');

   manifest = sourceManifest(testCase.TestData.root);
   manifest.files(1).required = false;
   writeManifest(testCase.TestData.manifest, manifest);
   testCase.verifyError(@() syntheticFileList(testCase), ...
      'icemodel:verification:fixtureFileList:invalidManifest');
end

function test_pack_writes_one_archive_per_capability(testCase)
   % Packing required capabilities produces two archives and one release manifest.
   result = packRequired(testCase);
   testCase.verifyEqual(numel(result.archive_files), 2);
   testCase.verifyTrue(all(isfile(result.archive_files)));
   testCase.verifyTrue(isfile(result.manifest_file));
   testCase.verifyEqual(string({result.manifest.archives.name})', [ ...
      "icemodel-v1.1-formal-core.tar.gz"; ...
      "icemodel-v1.1-verification-showcase.tar.gz"]);
   testCase.verifyTrue(all([result.manifest.archives.bytes] > 0));
   testCase.verifyTrue(all(strlength( ...
      string({result.manifest.archives.sha256})) == 64));
end

function test_required_archives_are_byte_deterministic(testCase)
   % Repeated packing of unchanged mandatory sources must be byte-identical.
   staging_one = canonicalTempname();
   staging_two = canonicalTempname();
   cleaner = onCleanup(@() cleanupTrees([staging_one; staging_two]));
   first = icemodel.verification.setup.packFixtures("v1.1", ...
      root=testCase.TestData.root, manifest=testCase.TestData.manifest, ...
      staging_dir=staging_one, silent=true);
   pause(2.1)
   second = icemodel.verification.setup.packFixtures("v1.1", ...
      root=testCase.TestData.root, manifest=testCase.TestData.manifest, ...
      staging_dir=staging_two, silent=true);

   first_hashes = arrayfun( ...
      @icemodel.verification.setup.fileSha256, first.archive_files);
   second_hashes = arrayfun( ...
      @icemodel.verification.setup.fileSha256, second.archive_files);
   for k = 1:numel(first.archive_files)
      testCase.verifyEqual(readBytes(first.archive_files(k)), ...
         readBytes(second.archive_files(k)));
   end
   testCase.verifyEqual(first_hashes, second_hashes);
   testCase.verifyEqual(readBytes(first.manifest_file), ...
      readBytes(second.manifest_file));
   testCase.verifyEqual( ...
      icemodel.verification.setup.fileSha256(first.manifest_file), ...
      icemodel.verification.setup.fileSha256(second.manifest_file));
end

function test_pack_rejects_source_drift_and_clobber(testCase)
   % A changed source cannot be silently blessed into a release artifact.
   mkfile(fullfile(testCase.TestData.root, "a", "formal1.mat"), "changed");
   testCase.verifyError(@() packRequired(testCase), ...
      'icemodel:verification:packFixtures:sourceMismatch');

   % Restoring manifest metadata allows one pack, but a second requires overwrite.
   writeManifest(testCase.TestData.manifest, ...
      sourceManifest(testCase.TestData.root));
   packRequired(testCase);
   testCase.verifyError(@() packRequired(testCase), ...
      'icemodel:verification:packFixtures:exists');
end

function test_multi_capability_pack_failure_leaves_no_partial_outputs(testCase)
   % A later archive failure must leave no final artifact and permit plain retry.
   writeManifest(testCase.TestData.manifest, ...
      lateArchiveFailureManifest(testCase.TestData.root));
   before = treeSnapshot(testCase.TestData.staging);
   caught = [];
   try
      packRequired(testCase);
   catch cause
      caught = cause;
   end
   testCase.verifyNotEmpty(caught);
   testCase.verifyEqual(treeSnapshot(testCase.TestData.staging), before);

   % With no partial artifact left behind, overwrite=false retry must succeed.
   writeManifest(testCase.TestData.manifest, ...
      sourceManifest(testCase.TestData.root));
   result = packRequired(testCase);
   testCase.verifyTrue(all(isfile(result.archive_files)));
   testCase.verifyTrue(isfile(result.manifest_file));
end

function test_multi_capability_pack_failure_preserves_existing_outputs(testCase)
   % overwrite=true must not disturb prior artifacts when a later build fails.
   packRequired(testCase);
   mkfile(fullfile(testCase.TestData.root, "a", "formal1.mat"), ...
      "replacement formal one");
   writeManifest(testCase.TestData.manifest, ...
      lateArchiveFailureManifest(testCase.TestData.root));
   before = treeSnapshot(testCase.TestData.staging);
   caught = [];
   try
      icemodel.verification.setup.packFixtures("v1.1", ...
         root=testCase.TestData.root, manifest=testCase.TestData.manifest, ...
         staging_dir=testCase.TestData.staging, overwrite=true, silent=true);
   catch cause
      caught = cause;
   end
   testCase.verifyNotEmpty(caught);
   testCase.verifyEqual(treeSnapshot(testCase.TestData.staging), before);
end

function test_pack_rejects_directory_artifact_collision(testCase)
   % overwrite=true must not replace a directory at a release artifact path.
   conflict = fullfile(testCase.TestData.staging, ...
      "icemodel-v1.1-formal-core.tar.gz");
   mkdir(conflict)
   mkfile(fullfile(conflict, "unrelated.txt"), "keep me");
   before = treeSnapshot(testCase.TestData.staging);

   testCase.verifyError(@() ...
      icemodel.verification.setup.packFixtures("v1.1", ...
      capabilities="formal-core", root=testCase.TestData.root, ...
      manifest=testCase.TestData.manifest, ...
      staging_dir=testCase.TestData.staging, overwrite=true, silent=true), ...
      'icemodel:verification:packFixtures:destinationTypeConflict');
   testCase.verifyEqual(treeSnapshot(testCase.TestData.staging), before);
end

function test_pack_rejects_manifest_source_symlink(testCase)
   % Packing must reject a declared link even when its target bytes match.
   source = fullfile(testCase.TestData.root, "a", "formal1.mat");
   outside = canonicalTempname();
   mkdir(outside)
   target = fullfile(outside, "formal1.mat");
   mkfile(target, "formal one");
   delete(source)
   makeSymbolicLink(target, source);
   cleaner = onCleanup(@() cleanupSymlinkFixture(source, outside));

   testCase.verifyError(@() packCapability(testCase, "formal-core"), ...
      'icemodel:verification:packFixtures:symlinkSource');
   archive = fullfile(testCase.TestData.staging, ...
      "icemodel-v1.1-formal-core.tar.gz");
   testCase.verifyFalse(isfile(archive));
end

function test_pack_rejects_symlinked_source_ancestor(testCase)
   % A linked ancestor below the system temp root must not redirect source reads.
   source_root = testCase.TestData.root;
   source_parent = string(fileparts(source_root));
   [~, source_name] = fileparts(source_root);
   alias_parent = canonicalTempname();
   linked_parent = fullfile(alias_parent, "source-parent");
   staging = canonicalTempname();
   mkdir(alias_parent)
   makeSymbolicLink(source_parent, linked_parent);
   cleaner = onCleanup(@() cleanupLinkedParent( ...
      linked_parent, [alias_parent; staging]));

   testCase.verifyError(@() ...
      icemodel.verification.setup.packFixtures("v1.1", ...
      capabilities="formal-core", ...
      root=fullfile(linked_parent, string(source_name)), ...
      manifest=testCase.TestData.manifest, staging_dir=staging, silent=true), ...
      'icemodel:verification:packFixtures:symlinkSource');
end

function test_pack_rejects_symlinked_staging_paths(testCase)
   % Staging at or below a caller-created link must not write outside its root.
   outside = canonicalTempname();
   alias_parent = canonicalTempname();
   link = fullfile(alias_parent, "staging-link");
   mkdir(outside)
   mkdir(alias_parent)
   makeSymbolicLink(outside, link);
   cleaner = onCleanup(@() cleanupLinkedParent( ...
      link, [alias_parent; outside]));

   % Exercise both a linked staging root and a linked ancestor of a new root.
   staging_paths = [link; fullfile(link, "nested")];
   for k = 1:numel(staging_paths)
      testCase.verifyError(@() ...
         icemodel.verification.setup.packFixtures("v1.1", ...
         capabilities="formal-core", root=testCase.TestData.root, ...
         manifest=testCase.TestData.manifest, ...
         staging_dir=staging_paths(k), silent=true), ...
         'icemodel:verification:packFixtures:symlinkDestination');
   end
   expected_outputs = [ ...
      "icemodel-v1.1-formal-core.tar.gz"; ...
      "icemodel-v1.1-data-manifest.json"];
   for k = 1:numel(expected_outputs)
      testCase.verifyFalse(isfile(fullfile(outside, expected_outputs(k))));
      testCase.verifyFalse(isfile(fullfile( ...
         outside, "nested", expected_outputs(k))));
   end
end

function test_default_call_uses_local_multi_archive_override(testCase)
   % Default downloading defers to supplied local archives for offline use.
   packed = packRequired(testCase);
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   result = icemodel.verification.setup.fetchFixtures("v1.1", ...
      root=destination, manifest=packed.manifest_file, ...
      archive=packed.archive_files, silent=true);
   testCase.verifyTrue(result.ok);
   testCase.verifyEqual(result.mode, "provisioned");
   testCase.verifyTrue(isfile(fullfile(destination, "a", "formal1.mat")));
   testCase.verifyTrue(isfile(fullfile(destination, "eval", "showcase", ...
      "reference.mat")));
end

function test_scalar_archive_requires_one_capability(testCase)
   % A scalar local archive must never be guessed across two capabilities.
   packed = packRequired(testCase);
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   testCase.verifyError(@() ...
      icemodel.verification.setup.fetchFixtures("v1.1", ...
      root=destination, manifest=packed.manifest_file, ...
      archive=packed.archive_files(1), silent=true), ...
      'icemodel:verification:fetchFixtures:archiveCapabilityMismatch');
end

function test_optional_capability_roundtrip_is_explicit(testCase)
   % Optional forcing packs and provisions only when selected by name.
   packed = icemodel.verification.setup.packFixtures("v1.1", ...
      capabilities="forcing-integration", root=testCase.TestData.root, ...
      manifest=testCase.TestData.manifest, ...
      staging_dir=testCase.TestData.staging, silent=true);
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   result = icemodel.verification.setup.fetchFixtures("v1.1", ...
      capabilities="forcing-integration", root=destination, ...
      manifest=packed.manifest_file, archive=packed.archive_file, silent=true);
   testCase.verifyTrue(result.ok);
   testCase.verifyTrue(isfile(fullfile(destination, "forcing", "source.nc")));
end

function test_explicit_no_download_prints_exact_command(testCase)
   % Missing data with download=false returns one copy-paste explicit command.
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   output = evalc(['result = icemodel.verification.setup.fetchFixtures(' ...
      '"v1.1", root=destination, manifest=testCase.TestData.manifest, ' ...
      'download=false, strict=false, silent=false);']);
   exact = "icemodel.verification.setup.fetchFixtures(""v1.1"", " ...
      + "capabilities=[""formal-core"",""verification-showcase""], root=" ...
      + quoteMatlabLiteral(destination) + ", manifest=" ...
      + quoteMatlabLiteral(testCase.TestData.manifest) + ", download=true)";
   testCase.verifyFalse(result.ok);
   testCase.verifyEqual(result.mode, "missing");
   testCase.verifyEqual(result.command, exact);
   testCase.verifyTrue(contains(string(output), exact));
   testCase.verifyFalse(isfolder(destination));
end

function test_explicit_no_download_errors_without_network(testCase)
   % Strict network-free verification errors before resolving a release URL.
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   testCase.verifyError(@() ...
      icemodel.verification.setup.fetchFixtures("v1.1", ...
      root=destination, manifest=testCase.TestData.manifest, ...
      download=false, release_url="not-a-url", silent=true), ...
      'icemodel:verification:fetchFixtures:incompleteFixtures');
   testCase.verifyFalse(isfolder(destination));
end

function test_fetch_command_override_matrix(testCase)
   % Exact commands cover default, individual, and combined override branches.
   quote = string(char(34));
   prefix = "icemodel.verification.setup.fetchFixtures(" + quote + "v1.1" ...
      + quote + ", capabilities=[" + quote + "formal-core" + quote + "]";
   default_exact = prefix + ", download=true)";
   custom_root = "root" + quote + "value";
   custom_manifest = "manifest" + quote + "value.json";
   custom_release = "https://fixtures.example/v1.1";
   custom_repo = "owner/project";
   root_literal = quote + "root" + quote + quote + "value" + quote;
   manifest_literal = quote + "manifest" + quote + quote ...
      + "value.json" + quote;

   testCase.verifyEqual( ...
      icemodel.verification.setup.fixtureFetchCommand( ...
      "v1.1", "formal-core"), default_exact);
   testCase.verifyEqual( ...
      icemodel.verification.setup.fixtureFetchCommand( ...
      "v1.1", "formal-core", root=custom_root), ...
      prefix + ", root=" + root_literal + ", download=true)");
   testCase.verifyEqual( ...
      icemodel.verification.setup.fixtureFetchCommand( ...
      "v1.1", "formal-core", manifest=custom_manifest), ...
      prefix + ", manifest=" + manifest_literal + ", download=true)");
   testCase.verifyEqual( ...
      icemodel.verification.setup.fixtureFetchCommand( ...
      "v1.1", "formal-core", root=custom_root, manifest=custom_manifest), ...
      prefix + ", root=" + root_literal + ", manifest=" ...
      + manifest_literal + ", download=true)");

   % Release URL and repository overrides remain independently copy-pasteable.
   testCase.verifyEqual( ...
      icemodel.verification.setup.fixtureFetchCommand( ...
      "v1.1", "formal-core", release_url=custom_release), ...
      prefix + ", release_url=" + quote + custom_release + quote ...
      + ", download=true)");
   testCase.verifyEqual( ...
      icemodel.verification.setup.fixtureFetchCommand( ...
      "v1.1", "formal-core", repo=custom_repo), ...
      prefix + ", repo=" + quote + custom_repo + quote + ", download=true)");
end

function test_missing_command_preserves_individual_overrides(testCase)
   % Root-only and manifest-only checks must repair the tree each call inspected.
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   root_result = icemodel.verification.setup.fetchFixtures("v1.1", ...
      capabilities="formal-core", root=destination, ...
      download=false, strict=false, silent=true);
   root_exact = "icemodel.verification.setup.fetchFixtures(""v1.1"", " ...
      + "capabilities=[""formal-core""], root=" ...
      + quoteMatlabLiteral(destination) + ", download=true)";
   testCase.verifyEqual(root_result.command, root_exact);

   % Point a custom manifest at unique missing paths under the default root.
   manifest = sourceManifest(testCase.TestData.root);
   prefix = "__fixture_command_missing_" + string(java.util.UUID.randomUUID());
   for k = 1:numel(manifest.files)
      manifest.files(k).path = prefix + "/" + string(manifest.files(k).path);
   end
   writeManifest(testCase.TestData.manifest, manifest);
   manifest_result = icemodel.verification.setup.fetchFixtures("v1.1", ...
      capabilities="formal-core", manifest=testCase.TestData.manifest, ...
      download=false, strict=false, silent=true);
   manifest_exact = "icemodel.verification.setup.fetchFixtures(""v1.1"", " ...
      + "capabilities=[""formal-core""], manifest=" ...
      + quoteMatlabLiteral(testCase.TestData.manifest) + ", download=true)";
   testCase.verifyEqual(manifest_result.command, manifest_exact);

   % Network-free inspection also preserves the selected release endpoint.
   release_url = "https://fixtures.example/releases/v1.1";
   repo = "owner/project";
   release_result = icemodel.verification.setup.fetchFixtures("v1.1", ...
      capabilities="formal-core", root=destination, ...
      manifest=testCase.TestData.manifest, download=false, ...
      release_url=release_url, repo=repo, strict=false, silent=true);
   release_exact = "icemodel.verification.setup.fetchFixtures(""v1.1"", " ...
      + "capabilities=[""formal-core""], root=" ...
      + quoteMatlabLiteral(destination) + ", manifest=" ...
      + quoteMatlabLiteral(testCase.TestData.manifest) + ", release_url=" ...
      + quoteMatlabLiteral(release_url) + ", repo=" ...
      + quoteMatlabLiteral(repo) + ", download=true)";
   testCase.verifyEqual(release_result.command, release_exact);
end

function test_default_download_failure_is_stable(testCase)
   % The default call reaches the download branch and fails invalid URLs clearly.
   packed = packRequired(testCase);
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   testCase.verifyError(@() ...
      icemodel.verification.setup.fetchFixtures("v1.1", ...
      root=destination, manifest=packed.manifest_file, ...
      release_url="not-a-url", silent=true), ...
      'icemodel:verification:fetchFixtures:downloadFailed');
end

function test_default_manifest_reaches_download_validation(testCase)
   % Mandatory default archive metadata must pass before default download fails.
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   testCase.verifyError(@() ...
      icemodel.verification.setup.fetchFixtures("v1.1", ...
      root=destination, release_url="not-a-url", silent=true), ...
      'icemodel:verification:fetchFixtures:downloadFailed');
end

function test_idempotence_and_unrelated_file_preservation(testCase)
   % A second verified call does not rewrite files; unrelated data survives.
   packed = packCapability(testCase, "formal-core");
   destination = canonicalTempname();
   mkdir(destination)
   cleaner = onCleanup(@() removeTree(destination));
   mkfile(fullfile(destination, "keep.txt"), "unrelated");
   first = icemodel.verification.setup.fetchFixtures("v1.1", ...
      capabilities="formal-core", root=destination, ...
      manifest=packed.manifest_file, archive=packed.archive_file, silent=true);
   target = fullfile(destination, "a", "formal1.mat");
   timestamp = java.io.File(char(target)).lastModified();
   second = icemodel.verification.setup.fetchFixtures("v1.1", ...
      capabilities="formal-core", root=destination, ...
      manifest=packed.manifest_file, download=false, silent=true);
   testCase.verifyEqual(first.mode, "provisioned");
   testCase.verifyEqual(second.mode, "verified");
   testCase.verifyEqual(java.io.File(char(target)).lastModified(), timestamp);
   testCase.verifyEqual(string(fileread(fullfile(destination, "keep.txt"))), ...
      "unrelated");
end

function test_archive_hash_mismatch_is_rejected_before_extract(testCase)
   % Altering archive bytes without its manifest hash must stop before extraction.
   packed = packCapability(testCase, "formal-core");
   flipFirstByte(packed.archive_file);
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   testCase.verifyError(@() ...
      icemodel.verification.setup.fetchFixtures("v1.1", ...
      capabilities="formal-core", root=destination, ...
      manifest=packed.manifest_file, archive=packed.archive_file, silent=true), ...
      'icemodel:verification:fetchFixtures:archiveHashMismatch');
   testCase.verifyFalse(isfolder(destination));
end

function test_extracted_file_hash_mismatch_is_rejected(testCase)
   % A rehashed archive with changed member bytes still fails per-file hashes.
   packed = packCapability(testCase, "formal-core");
   mkfile(fullfile(testCase.TestData.root, "a", "formal1.mat"), "tampered");
   replaceArchive(packed.archive_file, testCase.TestData.root, ...
      ["a/formal1.mat"; "b/formal2.mat"]);
   refreshArchiveMetadata(packed.manifest_file, "formal-core", ...
      packed.archive_file);
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   testCase.verifyError(@() ...
      icemodel.verification.setup.fetchFixtures("v1.1", ...
      capabilities="formal-core", root=destination, ...
      manifest=packed.manifest_file, archive=packed.archive_file, silent=true), ...
      'icemodel:verification:fetchFixtures:fileHashMismatch');
end

function test_missing_and_extra_members_are_rejected(testCase)
   % A complete manifest cannot be paired with an archive missing one file.
   packed = packCapability(testCase, "formal-core");
   replaceArchive(packed.archive_file, testCase.TestData.root, ...
      "a/formal1.mat");
   refreshArchiveMetadata(packed.manifest_file, "formal-core", ...
      packed.archive_file);
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   testCase.verifyError(@() fetchFormal(destination, packed), ...
      'icemodel:verification:fetchFixtures:missingMember');

   % Exact file headers without required ancestor headers are also incomplete.
   missing_directories = rawArchiveMembers(testCase.TestData.staging, ...
      ["a/formal1.mat"; "b/formal2.mat"], ['0'; '0'], "missing-directories");
   refreshArchiveMetadata(packed.manifest_file, "formal-core", ...
      missing_directories);
   packed.archive_file = missing_directories;
   testCase.verifyError(@() fetchFormal(destination, packed), ...
      'icemodel:verification:fetchFixtures:missingMember');

   % An undeclared regular file is rejected before extraction as well.
   mkfile(fullfile(testCase.TestData.root, "extra.bin"), "extra");
   replaceArchive(packed.archive_file, testCase.TestData.root, ...
      ["a/formal1.mat"; "b/formal2.mat"; "extra.bin"]);
   refreshArchiveMetadata(packed.manifest_file, "formal-core", ...
      packed.archive_file);
   testCase.verifyError(@() fetchFormal(destination, packed), ...
      'icemodel:verification:fetchFixtures:extraMember');
end

function test_unsafe_paths_and_links_are_rejected(testCase)
   % Parent traversal is rejected from raw headers before untar runs.
   packed = packCapability(testCase, "formal-core");
   unsafe = rawArchive(testCase.TestData.staging, "../escape.mat", '0');
   refreshArchiveMetadata(packed.manifest_file, "formal-core", unsafe);
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   packed.archive_file = unsafe;
   testCase.verifyError(@() fetchFormal(destination, packed), ...
      'icemodel:verification:fetchFixtures:unsafeMember');

   % Windows-style parent traversal is unsafe even on a POSIX test host.
   windows_unsafe = rawArchive(testCase.TestData.staging, ...
      ".." + string(char(92)) + "escape.mat", '0');
   refreshArchiveMetadata(packed.manifest_file, "formal-core", windows_unsafe);
   packed.archive_file = windows_unsafe;
   testCase.verifyError(@() fetchFormal(destination, packed), ...
      'icemodel:verification:fetchFixtures:unsafeMember');

   % Symbolic-link headers are forbidden even when their member path is declared.
   link_archive = rawArchive(testCase.TestData.staging, ...
      "a/formal1.mat", '2');
   refreshArchiveMetadata(packed.manifest_file, "formal-core", link_archive);
   packed.archive_file = link_archive;
   testCase.verifyError(@() fetchFormal(destination, packed), ...
      'icemodel:verification:fetchFixtures:unsafeMember');
end

function test_failed_promotion_restores_prior_data(testCase)
   % A parent-file conflict after the first move forces rollback of that move.
   packed = packCapability(testCase, "formal-core");
   destination = canonicalTempname();
   mkdir(destination)
   cleaner = onCleanup(@() removeTree(destination));
   mkfile(fullfile(destination, "a", "formal1.mat"), "prior formal one");
   mkfile(fullfile(destination, "b"), "blocking parent file");

   did_throw = false;
   try
      fetchFormal(destination, packed);
   catch
      did_throw = true;
   end
   testCase.verifyTrue(did_throw);
   testCase.verifyEqual(string(fileread(fullfile(destination, ...
      "a", "formal1.mat"))), "prior formal one");
   testCase.verifyEqual(string(fileread(fullfile(destination, "b"))), ...
      "blocking parent file");
end

function test_later_capability_failure_rolls_back_entire_call(testCase)
   % A showcase promotion failure must undo the earlier formal-core promotion.
   packed = packRequired(testCase);
   destination = canonicalTempname();
   mkdir(destination)
   cleaner = onCleanup(@() removeTree(destination));
   mkfile(fullfile(destination, "a", "formal1.mat"), "prior formal one");
   mkfile(fullfile(destination, "b", "formal2.mat"), "prior formal two");
   mkfile(fullfile(destination, "eval"), "blocking showcase parent");
   before = treeSnapshot(destination);

   % Both archives validate; failure occurs only after formal-core is promoted
   % and verification-showcase encounters the blocking eval parent file.
   caught = [];
   try
      icemodel.verification.setup.fetchFixtures("v1.1", ...
         root=destination, manifest=packed.manifest_file, ...
         archive=packed.archive_files, silent=true);
   catch cause
      caught = cause;
   end
   testCase.verifyNotEmpty(caught);
   testCase.verifyEqual(string(caught.identifier), ...
      "icemodel:verification:fetchFixtures:promotionFailed");
   testCase.verifyEqual(treeSnapshot(destination), before);
end

function test_destination_symlink_cannot_write_outside_root(testCase)
   % Promotion preflight must reject a linked parent before moving any file.
   packed = packCapability(testCase, "formal-core");
   destination = canonicalTempname();
   outside = canonicalTempname();
   mkdir(destination)
   mkdir(outside)
   link = fullfile(destination, "a");
   makeSymbolicLink(outside, link);
   cleaner = onCleanup(@() cleanupDestinationSymlink( ...
      link, destination, outside));

   testCase.verifyError(@() fetchFormal(destination, packed), ...
      'icemodel:verification:fetchFixtures:symlinkDestination');
   testCase.verifyFalse(isfile(fullfile(outside, "formal1.mat")));
   testCase.verifyFalse(isfile(fullfile(destination, "b", "formal2.mat")));
end

function test_file_destination_directory_conflict_preserves_children(testCase)
   % A declared file cannot replace a directory containing unrelated data.
   packed = packCapability(testCase, "formal-core");
   destination = canonicalTempname();
   conflict = fullfile(destination, "a", "formal1.mat");
   mkdir(conflict)
   mkfile(fullfile(conflict, "unrelated.txt"), "keep me");
   before = treeSnapshot(destination);
   cleaner = onCleanup(@() removeTree(destination));

   testCase.verifyError(@() fetchFormal(destination, packed), ...
      'icemodel:verification:fetchFixtures:destinationTypeConflict');
   testCase.verifyEqual(treeSnapshot(destination), before);
   testCase.verifyEqual(string(fileread( ...
      fullfile(conflict, "unrelated.txt"))), "keep me");
   clear cleaner
end

function test_symlinked_destination_root_is_rejected(testCase)
   % The selected root itself must not redirect promotion outside its boundary.
   packed = packCapability(testCase, "formal-core");
   destination = canonicalTempname();
   outside = canonicalTempname();
   mkdir(outside)
   makeSymbolicLink(outside, destination);
   cleaner = onCleanup(@() cleanupSymlinkFixture(destination, outside));

   testCase.verifyError(@() fetchFormal(destination, packed), ...
      'icemodel:verification:fetchFixtures:symlinkDestination');
   testCase.verifyFalse(isfile(fullfile(outside, "a", "formal1.mat")));
   testCase.verifyFalse(isfile(fullfile(outside, "b", "formal2.mat")));
   clear cleaner
end

function test_symlinked_destination_ancestor_is_rejected(testCase)
   % A caller-created ancestor must not redirect a prospective destination root.
   packed = packCapability(testCase, "formal-core");
   real_parent = canonicalTempname();
   alias_parent = canonicalTempname();
   linked_parent = fullfile(alias_parent, "destination-parent");
   mkdir(real_parent)
   mkdir(alias_parent)
   makeSymbolicLink(real_parent, linked_parent);
   cleaner = onCleanup(@() cleanupLinkedParent( ...
      linked_parent, [alias_parent; real_parent]));
   destination = fullfile(linked_parent, "new-root");

   testCase.verifyError(@() fetchFormal(destination, packed), ...
      'icemodel:verification:fetchFixtures:symlinkDestination');
   testCase.verifyFalse(isfolder(fullfile(real_parent, "new-root")));
end

function test_known_tmp_alias_allows_safe_pack_and_fetch(testCase)
   % The system /tmp alias is safe when no caller-created descendant is linked.
   source = aliasTempname();
   staging = aliasTempname();
   destination = aliasTempname();
   mkdir(source)
   mkfile(fullfile(source, "a", "formal1.mat"), "formal one");
   mkfile(fullfile(source, "b", "formal2.mat"), "formal two");
   cleaner = onCleanup(@() cleanupTrees([source; staging; destination]));

   packed = icemodel.verification.setup.packFixtures("v1.1", ...
      capabilities="formal-core", root=source, ...
      manifest=testCase.TestData.manifest, staging_dir=staging, silent=true);
   result = fetchFormal(destination, packed);
   testCase.verifyTrue(result.ok);
   testCase.verifyTrue(isfile(fullfile(destination, "a", "formal1.mat")));
end

function test_matching_installed_symlink_is_not_verified(testCase)
   % Matching outside bytes must not let the verified fast path accept a link.
   packed = packCapability(testCase, "formal-core");
   destination = canonicalTempname();
   outside = canonicalTempname();
   mkdir(outside)
   fetchFormal(destination, packed);
   link = fullfile(destination, "a", "formal1.mat");
   outside_file = fullfile(outside, "formal1.mat");
   [ok, message] = movefile(link, outside_file);
   assert(ok, message)
   makeSymbolicLink(outside_file, link);
   cleaner = onCleanup(@() cleanupDestinationSymlink( ...
      link, destination, outside));

   testCase.verifyEqual( ...
      icemodel.verification.setup.fileSha256(link), ...
      icemodel.verification.setup.fileSha256(outside_file));
   testCase.verifyError(@() ...
      icemodel.verification.setup.fetchFixtures("v1.1", ...
      capabilities="formal-core", root=destination, ...
      manifest=packed.manifest_file, silent=true), ...
      'icemodel:verification:fetchFixtures:symlinkDestination');
end

function test_installed_mismatch_is_reported_without_rewrite(testCase)
   % Verify-only mode reports a tampered installed path without changing it.
   packed = packCapability(testCase, "formal-core");
   destination = canonicalTempname();
   cleaner = onCleanup(@() removeTree(destination));
   fetchFormal(destination, packed);
   target = fullfile(destination, "a", "formal1.mat");
   mkfile(target, "installed tamper");
   result = icemodel.verification.setup.fetchFixtures("v1.1", ...
      capabilities="formal-core", root=destination, ...
      manifest=packed.manifest_file, download=false, ...
      strict=false, silent=true);
   testCase.verifyFalse(result.ok);
   testCase.verifyEqual(result.mismatched, "a/formal1.mat");
   testCase.verifyEqual(string(fileread(target)), "installed tamper");
   mismatch_exact = "icemodel.verification.setup.fetchFixtures(""v1.1"", " ...
      + "capabilities=[""formal-core""], root=" ...
      + quoteMatlabLiteral(destination) + ", manifest=" ...
      + quoteMatlabLiteral(packed.manifest_file) + ", download=true)";
   testCase.verifyEqual(result.command, mismatch_exact);
end

%% Packing and fetching helpers
function result = packRequired(testCase)
   %PACKREQUIRED Pack both required capabilities from the synthetic root.
   result = icemodel.verification.setup.packFixtures("v1.1", ...
      root=testCase.TestData.root, manifest=testCase.TestData.manifest, ...
      staging_dir=testCase.TestData.staging, silent=true);
end

function files = syntheticFileList(testCase)
   %SYNTHETICFILELIST Read the per-test source manifest through the public gate.
   files = icemodel.verification.setup.fixtureFileList( ...
      root=testCase.TestData.root, manifest=testCase.TestData.manifest);
end

function result = packCapability(testCase, capability)
   %PACKCAPABILITY Pack one named capability for isolated trust-boundary tests.
   result = icemodel.verification.setup.packFixtures("v1.1", ...
      capabilities=capability, root=testCase.TestData.root, ...
      manifest=testCase.TestData.manifest, ...
      staging_dir=testCase.TestData.staging, silent=true);
end

function result = fetchFormal(destination, packed)
   %FETCHFORMAL Provision one local formal-core archive.
   result = icemodel.verification.setup.fetchFixtures("v1.1", ...
      capabilities="formal-core", root=destination, ...
      manifest=packed.manifest_file, archive=packed.archive_file, silent=true);
end

%% Manifest helpers
function manifest = sourceManifest(root)
   %SOURCEMANIFEST Build complete synthetic file rows and unstaged archives.
   archives = [ ...
      archiveRow("formal-core", "icemodel-v1.1-formal-core.tar.gz", true); ...
      archiveRow("verification-showcase", ...
         "icemodel-v1.1-verification-showcase.tar.gz", true); ...
      archiveRow("forcing-integration", ...
         "icemodel-v1.1-forcing-integration.tar.gz", false)];
   definitions = { ...
      "formal-core", "a/formal1.mat", true; ...
      "formal-core", "b/formal2.mat", true; ...
      "verification-showcase", "eval/showcase/manifest.json", true; ...
      "verification-showcase", "eval/showcase/reference.mat", true; ...
      "forcing-integration", "forcing/source.nc", false};
   files = repmat(fileRow("", "", false, root), size(definitions, 1), 1);
   for k = 1:size(definitions, 1)
      files(k) = fileRow(definitions{k, 1}, definitions{k, 2}, ...
         definitions{k, 3}, root);
   end
   manifest = struct('schema_version', 1, 'version', "v1.1", ...
      'archives', archives, 'files', files);
end

function manifest = lateArchiveFailureManifest(root)
   %LATEARCHIVEFAILUREMANIFEST Force capability two to fail after one archive.
   manifest = sourceManifest(root);
   overlong_name = string(repmat('x', 1, 300)) + ".tar.gz";
   manifest.archives(2).name = overlong_name;
end

function row = archiveRow(capability, name, required)
   %ARCHIVEROW Build one pre-pack archive declaration.
   row = struct('capability', capability, 'name', name, ...
      'required', required, 'bytes', 0, 'sha256', "");
end

function row = fileRow(capability, relpath, required, root)
   %FILEROW Build one immutable file declaration from current bytes.
   if strlength(relpath) == 0
      row = struct('capability', "", 'path', "", 'required', false, ...
         'bytes', 0, 'sha256', "");
      return
   end
   pathname = underRoot(root, relpath);
   info = dir(char(pathname));
   row = struct('capability', capability, 'path', relpath, ...
      'required', required, 'bytes', info.bytes, ...
      'sha256', icemodel.verification.setup.fileSha256(pathname));
end

function writeManifest(pathname, manifest)
   %WRITEMANIFEST Write pretty JSON while preserving one-row array shapes.
   payload = manifest;
   payload.archives = num2cell(payload.archives);
   payload.files = num2cell(payload.files);
   fid = fopen(pathname, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fwrite(fid, jsonencode(payload, PrettyPrint=true), 'char');
end

function refreshArchiveMetadata(manifest_file, capability, archive_file)
   %REFRESHARCHIVEMETADATA Trust a test archive so deeper gates are exercised.
   manifest = jsondecode(fileread(manifest_file));
   index = find(string({manifest.archives.capability}) == capability, 1);
   info = dir(char(archive_file));
   manifest.archives(index).bytes = info.bytes;
   manifest.archives(index).sha256 = ...
      icemodel.verification.setup.fileSha256(archive_file);
   writeManifest(manifest_file, manifest);
end

%% Archive mutation helpers
function replaceArchive(archive_file, root, paths)
   %REPLACEARCHIVE Rebuild a test archive with an explicit member list.
   if isfile(archive_file)
      delete(archive_file)
   end
   paths = reshape(string(paths), [], 1);
   if ismac
      members = [ancestorDirectories(paths); paths];
      quoted_members = strings(numel(members), 1);
      for k = 1:numel(members)
         quoted_members(k) = quoteTestPath(members(k));
      end
      command = "COPYFILE_DISABLE=1 /usr/bin/tar --format ustar -czf " ...
         + quoteTestPath(archive_file) + " -C " + quoteTestPath(root) ...
         + " --no-recursion -- " + strjoin(quoted_members, " ");
      [status, message] = system(char(command));
      assert(status == 0, message)
      return
   end
   tar(char(archive_file), cellstr(paths), char(root));
end

function directories = ancestorDirectories(files)
   %ANCESTORDIRECTORIES Derive the exact directory headers for test archives.
   directories = strings(sum(count(files, "/")), 1);
   cursor = 0;
   for k = 1:numel(files)
      parts = split(files(k), "/");
      for n = 1:numel(parts) - 1
         cursor = cursor + 1;
         directories(cursor) = strjoin(parts(1:n), "/");
      end
   end
   directories = unique(directories);
end

function quoted = quoteTestPath(pathname)
   %QUOTETESTPATH Quote one path for a fixed native-tar test command.
   quote = char(39);
   escaped = strrep(char(pathname), quote, [quote '"' quote '"' quote]);
   quoted = string([quote escaped quote]);
end

function pathname = rawArchive(staging, member_name, type)
   %RAWARCHIVE Write one minimal gzipped tar header for pre-extraction rejection.
   pathname = rawArchiveMembers(staging, string(member_name), type, ...
      "raw-" + string(type));
end

function pathname = rawArchiveMembers(staging, member_names, types, label)
   %RAWARCHIVEMEMBERS Write zero-length tar members with no implicit directories.
   tar_file = fullfile(staging, string(label) + ".tar");
   if isfile(tar_file)
      delete(tar_file)
   end
   fid = fopen(tar_file, 'w');
   cleaner = onCleanup(@() fclose(fid));
   member_names = reshape(string(member_names), [], 1);
   types = reshape(types, [], 1);
   for k = 1:numel(member_names)
      header = zeros(1, 512, 'uint8');
      name_bytes = uint8(char(member_names(k)));
      header(1:numel(name_bytes)) = name_bytes;
      header(125:135) = uint8(sprintf('%011o', 0));
      header(136) = 0;
      header(157) = uint8(types(k));
      fwrite(fid, header, 'uint8');
   end
   fwrite(fid, zeros(1, 1024, 'uint8'), 'uint8');
   clear cleaner
   zipped = gzip(char(tar_file), char(staging));
   pathname = string(zipped{1});
end

function flipFirstByte(pathname)
   %FLIPFIRSTBYTE Corrupt an archive hash without changing its byte size.
   fid = fopen(pathname, 'r+');
   byte = fread(fid, 1, '*uint8');
   fseek(fid, 0, 'bof');
   fwrite(fid, bitxor(byte, uint8(1)), 'uint8');
   fclose(fid);
end

function bytes = readBytes(pathname)
   %READBYTES Read exact binary file bytes for byte-for-byte comparisons.
   fid = fopen(pathname, 'r');
   assert(fid >= 0, 'Cannot open test file: %s', pathname)
   cleaner = onCleanup(@() fclose(fid));
   bytes = fread(fid, Inf, '*uint8');
end

%% Filesystem helpers
function pathname = canonicalTempname()
   %CANONICALTEMPNAME Avoid platform temp aliases in ancestor-safety tests.
   pathname = string(java.io.File(char(tempname)).getCanonicalPath());
end

function pathname = aliasTempname()
   %ALIASTEMPNAME Exercise /tmp lexical paths where that system root exists.
   parent = string(tempdir);
   if isfolder("/tmp")
      parent = "/tmp";
   end
   pathname = string(tempname(parent));
end

function cleanupTrees(pathnames)
   %CLEANUPTREES Remove each independent temporary tree.
   for k = 1:numel(pathnames)
      removeTree(pathnames(k));
   end
end

function literal = quoteMatlabLiteral(value)
   %QUOTEMATLABLITERAL Build an independent expected quoted string literal.
   quote = string(char(34));
   literal = quote + replace(string(value), quote, quote + quote) + quote;
end

function mkfile(pathname, text)
   %MKFILE Write a deterministic small binary/text test file.
   folder = string(fileparts(pathname));
   if ~isfolder(folder)
      mkdir(folder)
   end
   fid = fopen(pathname, 'w');
   fwrite(fid, char(text), 'char');
   fclose(fid);
end

function makeSymbolicLink(target, link)
   %MAKESYMBOLICLINK Create a real filesystem link without shell quoting.
   target_object = java.nio.file.Paths.get(char(target), ...
      javaArray('java.lang.String', 0));
   link_object = java.nio.file.Paths.get(char(link), ...
      javaArray('java.lang.String', 0));
   attributes = javaArray('java.nio.file.attribute.FileAttribute', 0);
   java.nio.file.Files.createSymbolicLink( ...
      link_object, target_object, attributes);
end

function cleanupSymlinkFixture(link, outside)
   %CLEANUPSYMLINKFIXTURE Remove a link before deleting its external target.
   deleteSymbolicLink(link);
   removeTree(outside);
end

function cleanupDestinationSymlink(link, destination, outside)
   %CLEANUPDESTINATIONSYMLINK Remove the link before either test tree.
   deleteSymbolicLink(link);
   removeTree(destination);
   removeTree(outside);
end

function cleanupLinkedParent(link, containers)
   %CLEANUPLINKEDPARENT Remove a link before deleting its independent containers.
   deleteSymbolicLink(link);
   cleanupTrees(containers);
end

function deleteSymbolicLink(pathname)
   %DELETESYMBOLICLINK Delete only the link entry, never its target.
   path_object = java.nio.file.Paths.get(char(pathname), ...
      javaArray('java.lang.String', 0));
   if java.nio.file.Files.isSymbolicLink(path_object)
      java.nio.file.Files.deleteIfExists(path_object);
   end
end

function pathname = underRoot(root, relpath)
   %UNDERROOT Convert a relative POSIX manifest path to a local path.
   parts = cellstr(split(string(relpath), "/"));
   pathname = string(fullfile(root, parts{:}));
end

function removeTree(pathname)
   %REMOVETREE Remove a test-created tree if it exists.
   if isfolder(pathname)
      rmdir(pathname, 's')
   end
end

function snapshot = treeSnapshot(root)
   %TREESNAPSHOT Capture relative types and file bytes for rollback assertions.
   hits = dir(char(fullfile(root, "**", "*")));
   snapshot = strings(numel(hits), 1);
   for k = 1:numel(hits)
      pathname = fullfile(string(hits(k).folder), string(hits(k).name));
      relpath = extractAfter(pathname, strlength(root) + 1);
      if hits(k).isdir
         snapshot(k) = "D " + relpath;
      else
         digest = icemodel.verification.setup.fileSha256(pathname);
         snapshot(k) = "F " + relpath + " " + hits(k).bytes + " " + digest;
      end
   end
   snapshot = sort(snapshot);
end
