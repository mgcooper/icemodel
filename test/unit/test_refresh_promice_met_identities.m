function tests = test_refresh_promice_met_identities
   %TEST_REFRESH_PROMICE_MET_IDENTITIES Verify metadata-only native met pinning.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install project dependencies and create one isolated paired data tree.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.root = string(tempname);
   testCase.TestData.met_root = fullfile( ...
      testCase.TestData.root, 'input', 'met');
   testCase.TestData.manifest_file = fullfile( ...
      testCase.TestData.root, 'eval', 'promice', 'manifest.json');
   mkdir(fullfile(testCase.TestData.met_root, 'promice'))
   mkdir(fileparts(testCase.TestData.manifest_file))
end

function teardown(testCase)
   % Remove the generated data tree after every test.
   if isfolder(testCase.TestData.root)
      rmdir(testCase.TestData.root, 's')
   end
   clear testCase.TestData.cleanup
end

function test_refresh_pins_declared_files_and_is_idempotent(testCase)
   % A successful refresh records exact relative paths, sizes, and hashes while
   % retaining unrelated manifest fields and leaving an identical replay stable.
   files = ["promice/met_one_promice_2000_15m.mat"; ...
      "promice/met_one_promice_2001_15m.mat"];
   absolute = writeMetFiles(testCase.TestData.met_root, files);
   manifest = fixtureManifest("one", files);
   manifest.note = "preserve me";
   icemodel.verification.setup.writeManifest( ...
      testCase.TestData.manifest_file, manifest);

   returned = icemodel.verification.setup.refreshPromiceMetIdentities( ...
      manifest_file=testCase.TestData.manifest_file);
   identities = returned.cases.colocation.promice.met_file_identities;
   testCase.verifyEqual(string({identities.file})', files);
   for n = 1:numel(files)
      info = dir(absolute(n));
      testCase.verifyEqual(identities(n).size_bytes, info.bytes);
      testCase.verifyEqual(string(identities(n).sha256), ...
         icemodel.verification.setup.fileSha256(absolute(n)));
   end
   testCase.verifyEqual(string(returned.note), "preserve me");

   first_payload = fileread(testCase.TestData.manifest_file);
   replay = icemodel.verification.setup.refreshPromiceMetIdentities( ...
      manifest_file=testCase.TestData.manifest_file);
   testCase.verifyEqual(fileread(testCase.TestData.manifest_file), first_payload);
   testCase.verifyEqual(replay, returned);
end

function test_refresh_clears_stale_identities_without_met_files(testCase)
   % A PROMICE leg with no declared met files receives the typed empty identity
   % set rather than retaining a stale fingerprint.
   manifest = fixtureManifest("empty", strings(0, 1));
   manifest.cases.colocation.promice.met_file_identities = struct( ...
      'file', "stale.mat", 'size_bytes', 1, 'sha256', "stale");
   icemodel.verification.setup.writeManifest( ...
      testCase.TestData.manifest_file, manifest);

   returned = icemodel.verification.setup.refreshPromiceMetIdentities( ...
      manifest_file=testCase.TestData.manifest_file);
   identities = returned.cases.colocation.promice.met_file_identities;
   testCase.verifyEmpty(identities);
   testCase.verifyEqual(fieldnames(identities), ...
      {'file'; 'size_bytes'; 'sha256'});
end

function test_refresh_rejects_missing_inputs(testCase)
   % Missing manifest and met-root failures occur before any manifest write.
   missing_manifest = fullfile(testCase.TestData.root, 'missing.json');
   testCase.verifyError(@() ...
      icemodel.verification.setup.refreshPromiceMetIdentities( ...
      manifest_file=missing_manifest), ...
      'icemodel:verification:refreshPromiceMetIdentities:missingManifest');

   manifest = fixtureManifest("one", "promice/met_one.mat");
   icemodel.verification.setup.writeManifest( ...
      testCase.TestData.manifest_file, manifest);
   testCase.verifyError(@() ...
      icemodel.verification.setup.refreshPromiceMetIdentities( ...
      manifest_file=testCase.TestData.manifest_file, ...
      met_root=fullfile(testCase.TestData.root, 'absent')), ...
      'icemodel:verification:refreshPromiceMetIdentities:missingMetRoot');
end

function test_refresh_rejects_invalid_or_ambiguous_manifest(testCase)
   % A missing case array and duplicate declared paths are rejected explicitly.
   icemodel.verification.setup.writeManifest( ...
      testCase.TestData.manifest_file, struct('dataset_family', "promice"));
   testCase.verifyError(@() ...
      icemodel.verification.setup.refreshPromiceMetIdentities( ...
      manifest_file=testCase.TestData.manifest_file), ...
      'icemodel:verification:refreshPromiceMetIdentities:invalidManifest');

   duplicate = fixtureManifest("one", ...
      ["promice/met_one.mat"; "promice/met_one.mat"]);
   icemodel.verification.setup.writeManifest( ...
      testCase.TestData.manifest_file, duplicate);
   testCase.verifyError(@() ...
      icemodel.verification.setup.refreshPromiceMetIdentities( ...
      manifest_file=testCase.TestData.manifest_file), ...
      'icemodel:verification:refreshPromiceMetIdentities:ambiguousMetFiles');
end

function test_refresh_rejects_missing_file_without_partial_write(testCase)
   % Every identity is assembled before publication, so one missing file leaves
   % the manifest bytes unchanged even when an earlier file was hashable.
   files = ["promice/met_present.mat"; "promice/met_missing.mat"];
   writeMetFiles(testCase.TestData.met_root, files(1));
   manifest = fixtureManifest("one", files);
   icemodel.verification.setup.writeManifest( ...
      testCase.TestData.manifest_file, manifest);
   before = fileread(testCase.TestData.manifest_file);

   testCase.verifyError(@() ...
      icemodel.verification.setup.refreshPromiceMetIdentities( ...
      manifest_file=testCase.TestData.manifest_file), ...
      'icemodel:verification:refreshPromiceMetIdentities:missingMetFile');
   testCase.verifyEqual(fileread(testCase.TestData.manifest_file), before);
end

function test_refresh_rejects_path_outside_met_root(testCase)
   % A relative traversal cannot authorize hashing an arbitrary file outside
   % the selected staged met root.
   outside = fullfile(testCase.TestData.root, 'outside.mat');
   payload = 1;
   save(outside, 'payload')
   manifest = fixtureManifest("one", "../../outside.mat");
   icemodel.verification.setup.writeManifest( ...
      testCase.TestData.manifest_file, manifest);

   testCase.verifyError(@() ...
      icemodel.verification.setup.refreshPromiceMetIdentities( ...
      manifest_file=testCase.TestData.manifest_file), ...
      'icemodel:verification:refreshPromiceMetIdentities:artifactOutsideRoot');
end

function manifest = fixtureManifest(case_id, met_files)
   %FIXTUREMANIFEST Build the smallest valid PROMICE producer manifest.
   leg = struct('met_files', reshape(string(met_files), [], 1));
   entry = struct('case_id', case_id, ...
      'colocation', struct('promice', leg));
   manifest = struct('dataset_family', "promice", 'cases', entry);
end

function absolute = writeMetFiles(met_root, relative)
   %WRITEMETFILES Create stable MAT payloads at declared fixture paths.
   relative = reshape(string(relative), [], 1);
   absolute = strings(size(relative));
   for n = 1:numel(relative)
      absolute(n) = fullfile(met_root, relative(n));
      folder = fileparts(absolute(n));
      if ~isfolder(folder)
         mkdir(folder)
      end
      payload = n;
      save(absolute(n), 'payload')
   end
end
