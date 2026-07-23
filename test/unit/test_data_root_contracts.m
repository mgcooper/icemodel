function tests = test_data_root_contracts
   %TEST_DATA_ROOT_CONTRACTS Verify scoped data-root selection and cleanup.

   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   %SETUPONCE Add project paths and retain a temporary test workspace.

   rootdir = icemodel.internal.fullpath();
   addpath(genpath(fullfile(rootdir, 'icemodel')))
   addpath(genpath(fullfile(rootdir, 'test')))
   testCase.TestData.rootdir = string(rootdir);
   testCase.TestData.tmp = string(tempname);
   mkdir(testCase.TestData.tmp)
end

function teardownOnce(testCase)
   %TEARDOWNONCE Remove the temporary test workspace.

   if isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
end

function test_bootstrap_data_root_precedence_and_cleanup(testCase)
   % DATA_ROOT must override the case and restore every raw config env value.

   [names, original] = snapshotConfig();
   restore = onCleanup(@() restoreConfig(names, original));
   path_before = path;

   % Install mixed caller state: custom raw paths/metadata plus one truly unset
   % config field. Cleanup must preserve this exact distinction.
   caller_values = mixedCallerConfig(names, "caller");
   setConfig(names, caller_values)
   testCase.verifyEqual(currentRawConfig(names), caller_values)
   data_root = fullfile(testCase.TestData.tmp, 'bootstrap-data');
   mkdir(fullfile(data_root, 'input'))
   mkdir(fullfile(data_root, 'eval'))

   % An invalid case is deliberately ignored because DATA_ROOT has precedence.
   [~, input_path, output_path, eval_path, cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment( ...
      data_root=data_root, icemodel_config_casename="ignored");
   testCase.verifyClass(cleanup, 'onCleanup')
   testCase.verifyEqual(input_path, string(fullfile(data_root, 'input')))
   testCase.verifyEqual(output_path, string(fullfile(data_root, 'output')))
   testCase.verifyEqual(eval_path, string(fullfile(data_root, 'eval')))
   testCase.verifyTrue(any(strcmp(strsplit(path, pathsep), data_root)))

   % Clearing cleanup must restore the exact caller config and MATLAB path.
   clear cleanup
   testCase.verifyEqual(currentRawConfig(names), caller_values)
   testCase.verifyEqual(path, path_before)
   clear restore
end

function test_bootstrap_cleanup_runs_after_error(testCase)
   % A thrown suite error must still restore the caller's scoped configuration.

   [names, original] = snapshotConfig();
   restore = onCleanup(@() restoreConfig(names, original));
   path_before = path;
   caller_values = mixedCallerConfig(names, "error-caller");
   setConfig(names, caller_values)
   testCase.verifyEqual(currentRawConfig(names), caller_values)
   data_root = fullfile(testCase.TestData.tmp, 'error-data');
   mkdir(fullfile(data_root, 'input'))
   mkdir(fullfile(data_root, 'eval'))

   % The local helper keeps cleanup in its stack frame, then fails deliberately.
   testCase.verifyError(@() failAfterBootstrap(data_root), ...
      'test:dataRoots:forcedFailure')
   testCase.verifyEqual(currentRawConfig(names), caller_values)
   testCase.verifyEqual(path, path_before)
   clear restore
end

function test_formal_classes_inherit_runner_data_root(testCase)
   % Nested class setup must preserve the candidate root selected by its runner.

   config_names = configNames();
   selector_names = ["ICEMODEL_TEST_DATA_ROOT"; ...
      "ICEMODEL_TEST_SMBMODEL"; "ICEMODEL_TEST_SITENAME"; ...
      "ICEMODEL_TEST_FORCINGS"; "ICEMODEL_TEST_USERDATA"; ...
      "ICEMODEL_TEST_USERVARS"; "ICEMODEL_TEST_SIMYEAR"; ...
      "ICEMODEL_TEST_SIMYEARS"; "ICEMODEL_TEST_N_SPINUP_YEARS"; ...
      "ICEMODEL_TEST_SOLVER"];
   names = [config_names; selector_names];
   original = currentRawConfig(names);
   restore = onCleanup(@() setConfig(names, original));

   % Emulate the outer formal runner with the complete 15-minute demo tree.
   data_root = fullfile(testCase.TestData.rootdir, 'demo', 'data');
   [~, ~, ~, ~, outer_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(data_root=data_root);
   testCase.verifyClass(outer_cleanup, 'onCleanup')
   setenv('ICEMODEL_TEST_DATA_ROOT', data_root)

   % Regression setup previously replaced this outer root with test/data.
   regression_case = IcemodelRegressionTest();
   regression_case.configureCases();
   testCase.verifyEqual(string(getenv('ICEMODEL_DATA_PATH')), data_root)
   testCase.verifyEqual(string(getenv('ICEMODEL_INPUT_PATH')), ...
      string(fullfile(data_root, 'input')))
   regression_case.restoreConfig();

   % Performance setup resolves a concrete option bundle, exposing the selected
   % input root and exact 15-minute forcing paths without timing the model.
   setenv('ICEMODEL_TEST_SMBMODEL', 'skinmodel')
   setenv('ICEMODEL_TEST_SITENAME', 'kanm')
   setenv('ICEMODEL_TEST_FORCINGS', 'kanm')
   setenv('ICEMODEL_TEST_USERDATA', '')
   setenv('ICEMODEL_TEST_USERVARS', '')
   setenv('ICEMODEL_TEST_SIMYEAR', '2016')
   setenv('ICEMODEL_TEST_SIMYEARS', '2015,2016')
   setenv('ICEMODEL_TEST_N_SPINUP_YEARS', '1')
   setenv('ICEMODEL_TEST_SOLVER', '1')
   perf_case = IcemodelPerfTest();
   perf_case.configureCases();
   testCase.verifyEqual(string(perf_case.opts.pathinput), ...
      string(fullfile(data_root, 'input')))
   testCase.verifyTrue(all(contains(string(perf_case.opts.metfname), ...
      '_15m.mat')))
   perf_case.restoreConfig();

   % Each nested cleanup restores the still-active outer runner selection.
   testCase.verifyEqual(string(getenv('ICEMODEL_DATA_PATH')), data_root)
   clear outer_cleanup regression_case perf_case
   clear restore
end

function test_resolver_precedence_and_legacy_alias(testCase)
   % Whole-root selection precedes leaves and the legacy alias remains supported.

   data_root = fullfile(testCase.TestData.tmp, 'selected');
   [eval_root, input_root] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      data_root=data_root, evaluation_data_root="ignored-eval", ...
      input_data_root="ignored-input", icemodel_config_casename="ignored");
   testCase.verifyEqual(eval_root, string(fullfile(data_root, 'eval')))
   testCase.verifyEqual(input_root, string(fullfile(data_root, 'input')))

   % OUTPUT_ROOT is an exact legacy alias, but combining aliases is ambiguous.
   [legacy_eval, legacy_input] = ...
      icemodel.verification.setup.resolveStagingRoots(output_root=data_root);
   testCase.verifyEqual(legacy_eval, eval_root)
   testCase.verifyEqual(legacy_input, input_root)
   testCase.verifyError(@() ...
      icemodel.verification.setup.resolveStagingRoots( ...
      data_root=data_root, output_root=data_root), ...
      'icemodel:verification:resolveStagingRoots:conflictingRoots')
end

function test_resolver_infers_siblings_in_both_directions(testCase)
   % Either explicit leaf must derive its peer from the same parent.

   parent = fullfile(testCase.TestData.tmp, 'paired');
   eval_root = fullfile(parent, 'eval');
   input_root = fullfile(parent, 'input');
   [from_eval, inferred_input] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      evaluation_data_root=eval_root);
   [inferred_eval, from_input] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      input_data_root=input_root);
   testCase.verifyEqual(from_eval, string(eval_root))
   testCase.verifyEqual(inferred_input, string(input_root))
   testCase.verifyEqual(inferred_eval, string(eval_root))
   testCase.verifyEqual(from_input, string(input_root))

   % Supplying both leaves keeps both exact paths even when parents differ.
   other_input = fullfile(testCase.TestData.tmp, 'other', 'input');
   [paired_eval, paired_input] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      evaluation_data_root=eval_root, input_data_root=other_input);
   testCase.verifyEqual(paired_eval, string(eval_root))
   testCase.verifyEqual(paired_input, string(other_input))
end

function test_resolver_cases_are_nonmutating(testCase)
   % Case resolution must not replace an interactive caller's active config.

   [names, original] = snapshotConfig();
   restore = onCleanup(@() restoreConfig(names, original));
   sentinels = "interactive-" + string(1:numel(names))';
   setConfig(names, cellstr(sentinels))

   % Resolve every public case and compare the exact owned roots.
   cases = ["test", "verification", "demo"];
   roots = [fullfile(testCase.TestData.rootdir, 'test', 'data'); ...
      fullfile(testCase.TestData.rootdir, 'data'); ...
      fullfile(testCase.TestData.rootdir, 'demo', 'data')];
   for n = 1:numel(cases)
      [eval_root, input_root] = ...
         icemodel.verification.setup.resolveStagingRoots( ...
         icemodel_config_casename=cases(n));
      testCase.verifyEqual(eval_root, string(fullfile(roots(n), 'eval')))
      testCase.verifyEqual(input_root, string(fullfile(roots(n), 'input')))
   end
   testCase.verifyEqual(currentConfig(names), sentinels)
   clear restore
end

function test_manifest_discovery_does_not_cross_roots(testCase)
   % Default discovery must never append demo/test fallback manifests.

   eval_root = string(fullfile(testCase.TestData.rootdir, 'data', 'eval'));
   files = icemodel.verification.helpers.familyManifestFiles();
   testCase.verifyTrue(all(startsWith(files, eval_root + filesep)))

   % An explicitly empty scope remains empty even when other repo roots exist.
   empty_root = fullfile(testCase.TestData.tmp, 'empty-eval');
   mkdir(empty_root)
   explicit = icemodel.verification.helpers.familyManifestFiles( ...
      evaluation_data_root=empty_root);
   testCase.verifyEmpty(explicit)
end

function test_family_manifest_helper_data_root_precedence_and_isolation(testCase)
   % The direct helper must honor whole-root precedence without fallback reads.

   selected_root = fullfile(testCase.TestData.tmp, 'manifest-selected');
   ignored_root = fullfile(testCase.TestData.tmp, 'manifest-ignored');
   selected_file = writeFamilyManifest(selected_root, 'selected_family');
   ignored_file = writeFamilyManifest(ignored_root, 'ignored_family');

   % DATA_ROOT owns discovery even when lower-precedence selectors conflict.
   files = icemodel.verification.helpers.familyManifestFiles( ...
      data_root=selected_root, ...
      evaluation_data_root=fullfile(ignored_root, 'eval'), ...
      icemodel_config_casename="not-a-real-case");
   testCase.verifyEqual(files, string(selected_file))
   testCase.verifyFalse(any(files == string(ignored_file)))

   % Existing explicit-leaf and nonmutating config selectors remain supported.
   leaf_files = icemodel.verification.helpers.familyManifestFiles( ...
      evaluation_data_root=fullfile(ignored_root, 'eval'));
   testCase.verifyEqual(leaf_files, string(ignored_file))
   % The minimal demo root intentionally has no evaluation manifests, and the
   % helper must not fall back to populated test or verification roots.
   testCase.verifyEmpty( ...
      icemodel.verification.helpers.familyManifestFiles( ...
      icemodel_config_casename="demo"))

   % An empty whole-root scope stays empty even when all repo roots are populated.
   empty_root = fullfile(testCase.TestData.tmp, 'manifest-empty');
   mkdir(fullfile(empty_root, 'eval'))
   testCase.verifyEmpty( ...
      icemodel.verification.helpers.familyManifestFiles(data_root=empty_root))
end

function test_public_readers_forward_data_root(testCase)
   % Empty scoped trees expose forwarding through every public reader by error.

   data_root = fullfile(testCase.TestData.tmp, 'empty-public');
   mkdir(fullfile(data_root, 'eval'))
   mkdir(fullfile(data_root, 'input'))
   expected_eval = string(fullfile(data_root, 'eval'));

   % Catalog discovery and the artifact audit stay inside the selected tree.
   testCase.verifyEmpty(icemodel.verification.listcases(data_root=data_root))
   report = icemodel.verification.auditArtifacts( ...
      data_root=data_root, families="laugh_tests", report_dir="");
   paths = string({report.findings.path});
   testCase.verifyTrue(any(startsWith(paths, expected_eval)))

   % Manifest-backed comparison, plotting, Colbeck, and runner paths must all
   % fail at the same selected eval root rather than falling back elsewhere.
   calls = { ...
      @() icemodel.verification.loadmanifest("colbeck1976", data_root=data_root)
      @() icemodel.verification.comparecase("colbeck1976", ...
         data_root=data_root, make_plot=false)
      @() icemodel.verification.plotcase("colbeck1976", ...
         data_root=data_root, visible="off")
      @() icemodel.verification.plotscatter("cdp", ...
         data_root=data_root, visible="off")};
   for n = 1:numel(calls)
      message = errorMessage(calls{n});
      testCase.verifySubstring(message, expected_eval)
   end

   % Public showcase entry points report one exact, non-downloading repair
   % command and restore the scoped configuration on failure.
   expected_command = expectedShowcaseFetchCommand(data_root);
   message = errorMessage(@() ...
      icemodel.verification.colbeck.compareSolutions( ...
      data_root=data_root, make_plot=false));
   testCase.verifySubstring(message, expected_command)

   [names, before] = snapshotConfig();
   message = errorMessage(@() run_snow_verification_suite(data_root=data_root));
   testCase.verifySubstring(message, expected_command)
   testCase.verifyEqual(currentConfig(names), string(before))
end

function test_suite_rejects_disjoint_explicit_leaf_roots(testCase)
   % Model execution cannot assign dependent config paths to two owner trees.

   eval_root = fullfile(testCase.TestData.tmp, 'eval-owner', 'eval');
   input_root = fullfile(testCase.TestData.tmp, 'input-owner', 'input');
   mkdir(eval_root)
   mkdir(input_root)
   [names, before] = snapshotConfig();

   % The read-only Colbeck path owns its public evaluation files independently
   % and therefore reports a repair command for the evaluation parent.
   message = errorMessage(@() ...
      icemodel.verification.colbeck.compareSolutions( ...
      evaluation_data_root=eval_root, input_data_root=input_root, ...
      make_plot=false));
   testCase.verifySubstring(message, ...
      expectedShowcaseFetchCommand(string(fileparts(eval_root))))

   testCase.verifyError(@() run_snow_verification_suite( ...
      evaluation_data_root=eval_root, input_data_root=input_root), ...
      ['icemodel:verification:runSnowVerificationSuite:', ...
      'disjointDataRoots'])
   testCase.verifyEqual(currentConfig(names), string(before))
end

function test_colbeck_rejects_mismatched_showcase_bytes(testCase)
   % Existing but modified public files must fail hash verification before use.
   source_root = string(fullfile(testCase.TestData.rootdir, 'test', 'data'));
   data_root = fullfile(testCase.TestData.tmp, 'tampered-showcase');
   files = icemodel.verification.setup.fixtureFileList( ...
      capabilities="verification-showcase", root=source_root);
   source_paths = strings(numel(files), 1);
   for k = 1:numel(files)
      source_paths(k) = fixturePath(source_root, files(k));
   end

   % The public archive is intentionally ignored; a clean checkout skips this
   % integration fixture while the synthetic provisioning tests remain active.
   testCase.assumeTrue(all(isfile(source_paths)), ...
      "Public showcase is not provisioned. Run: " ...
      + expectedShowcaseFetchCommand(source_root))
   for k = 1:numel(files)
      source = source_paths(k);
      destination = fixturePath(data_root, files(k));
      folder = string(fileparts(destination));
      if ~isfolder(folder)
         mkdir(folder)
      end
      [ok, message] = copyfile(source, destination);
      assert(ok, message)
   end

   % Appending one byte preserves path completeness while invalidating the hash.
   target = fixturePath(data_root, files(1));
   fid = fopen(target, 'a');
   assert(fid >= 0, 'Cannot open showcase fixture for tampering.')
   cleaner = onCleanup(@() fclose(fid));
   fwrite(fid, uint8(0), 'uint8');
   clear cleaner

   message = errorMessage(@() ...
      icemodel.verification.colbeck.compareSolutions( ...
      data_root=data_root, make_plot=false));
   testCase.verifySubstring(message, expectedShowcaseFetchCommand(data_root))
end

function failAfterBootstrap(data_root)
   %FAILAFTERBOOTSTRAP Exercise cleanup while unwinding a thrown suite error.

   [~, ~, ~, ~, cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(data_root=data_root);
   assert(isa(cleanup, 'onCleanup'))
   error('test:dataRoots:forcedFailure', 'forced cleanup failure')
end

function filename = writeFamilyManifest(data_root, family)
   %WRITEFAMILYMANIFEST Create one minimal direct-discovery fixture.

   family_root = fullfile(data_root, 'eval', family);
   mkdir(family_root)
   filename = fullfile(family_root, 'manifest.json');

   % The discovery helper only requires the declared filename; valid empty JSON
   % keeps the fixture harmless if inspected during a failed assertion.
   fid = fopen(filename, 'w');
   cleanup = onCleanup(@() fclose(fid));
   fprintf(fid, '{}\n');
   clear cleanup
end

function message = errorMessage(call)
   %ERRORMESSAGE Return the message from one expected failing call.

   try
      call();
      message = "";
   catch err
      message = string(err.message);
   end
end

function command = expectedShowcaseFetchCommand(data_root)
   %EXPECTEDSHOWCASEFETCHCOMMAND Build the canonical public-data repair command.
   manifest = string(fullfile(icemodel.getpath('test'), 'assets', ...
      'icemodel-v1.1-data-manifest.json'));
   command = icemodel.verification.setup.fixtureFetchCommand( ...
      "v1.1", "verification-showcase", ...
      root=data_root, manifest=manifest);
end

function pathname = fixturePath(root, relpath)
   %FIXTUREPATH Convert one validated POSIX manifest path to a local path.
   parts = cellstr(split(string(relpath), "/"));
   pathname = string(fullfile(root, parts{:}));
end

function [names, values] = snapshotConfig()
   %SNAPSHOTCONFIG Capture every environment value written by config.

   names = configNames();
   values = currentRawConfig(names);
end

function names = configNames()
   %CONFIGNAMES Derive the complete environment schema from config itself.

   names = string(fieldnames(icemodel.config('setenv', false)));
end

function values = currentConfig(names)
   %CURRENTCONFIG Read raw config values without applying defaults.

   values = strings(numel(names), 1);
   for n = 1:numel(names)
      values(n) = string(getenv(names(n)));
   end
end

function values = currentRawConfig(names)
   %CURRENTRAWCONFIG Read exact raw values, including genuinely unset fields.

   values = cell(numel(names), 1);
   for n = 1:numel(names)
      values{n} = getenv(names(n));
   end
end

function values = mixedCallerConfig(names, prefix)
   %MIXEDCALLERCONFIG Build custom caller metadata with one unset field.

   values = cell(numel(names), 1);
   for n = 1:numel(names)
      values{n} = char(prefix + "-" + names(n));
   end

   % Exercise metadata independently from path restoration and leave contact
   % genuinely unset so cleanup must replay an empty raw value.
   version_index = find(names == "ICEMODEL_VERSION", 1);
   reference_index = find(names == "ICEMODEL_REFERENCE", 1);
   contact_index = find(names == "ICEMODEL_CONTACT", 1);
   assert(~isempty(version_index))
   assert(~isempty(reference_index))
   assert(~isempty(contact_index))
   values{version_index} = char(prefix + "-custom-version");
   values{reference_index} = char(prefix + "-custom-reference");
   values{contact_index} = '';
end

function setConfig(names, values)
   %SETCONFIG Install the supplied raw config values.

   for n = 1:numel(names)
      setenv(names(n), values{n})
   end
end

function restoreConfig(names, values)
   %RESTORECONFIG Restore raw config values after a test.

   setConfig(names, values)
end
