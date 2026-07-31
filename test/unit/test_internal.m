function tests = test_internal
   %TEST_INTERNAL Test the toolbox internal functions.
   tests = functiontests(localfunctions);
end

function setup(testCase) %#ok<INUSD>
   %SETUP Reset the fixture scratch folder before each internal helper test.

end

function teardown(testCase) %#ok<INUSD>
   %TEARDOWN Remove the fixture scratch folder after each internal helper test.

end

function test_ispathinside_containment_contract(testCase)
   %TEST_ISPATHINSIDE_CONTAINMENT_CONTRACT Canonical containment predicate.
   % The shared helper must accept root-equality and true descendants,
   % reject siblings whose names merely share a prefix and dot-dot
   % escapes, resolve symlinked roots, and stay well defined for
   % nonexistent leaves under an existing root.
   root = fullfile(tempdir, 'ipi_root');
   sibling = fullfile(tempdir, 'ipi_root2');
   mkdir(root);
   mkdir(sibling);
   cleaner = onCleanup(@() cellfun(@(d) rmdir(d, 's'), {root, sibling}));

   % Root equality and a direct child are inside.
   testCase.verifyTrue(icemodel.internal.isPathInside(root, root));
   child = fullfile(root, 'a.mat');
   testCase.verifyTrue(icemodel.internal.isPathInside(child, root));
   % A sibling sharing the root's name prefix is outside.
   testCase.verifyFalse(icemodel.internal.isPathInside( ...
      fullfile(sibling, 'a.mat'), root));
   % A dot-dot escape resolves outside the root.
   testCase.verifyFalse(icemodel.internal.isPathInside( ...
      fullfile(root, '..', 'ipi_root2', 'a.mat'), root));
   % A nonexistent leaf under the root still resolves inside.
   testCase.verifyTrue(icemodel.internal.isPathInside( ...
      fullfile(root, 'missing', 'leaf.mat'), root));
   % A symlinked alias of the root resolves to the same canonical root.
   link = fullfile(tempdir, 'ipi_link');
   if isfile(link) || isfolder(link)
      delete(link);
   end
   status = system(sprintf('ln -s "%s" "%s"', root, link));
   if status == 0
      link_cleaner = onCleanup(@() delete(link));
      testCase.verifyTrue(icemodel.internal.isPathInside( ...
         fullfile(link, 'a.mat'), root));
   end
end

function test_basepath(testCase)
   %TEST_BASEPATH Verify internal.basepath path composition.
   modelpath = fullfile(icemodel.internal.fullpath(), 'icemodel');

   testCase.verifyTrue(isfolder(modelpath), ...
      'Expected icemodel/ folder to exist in the top level project folder.');
end

function test_functionSignatures(testCase)
   % Validate the toolbox-level signatures.
   T = validateFunctionSignaturesJSON(fullfile( ...
      icemodel.internal.fullpath(), 'icemodel', 'functionSignatures.json'));
   T_test = validateFunctionSignaturesJSON(fullfile( ...
      icemodel.internal.fullpath(), 'test', 'functionSignatures.json'));

   testCase.verifyEmpty(T, ...
      'The functionSignatures.json file contains invalid entries.');
   testCase.verifyEmpty(T_test, ...
      'The test/functionSignatures.json file contains invalid entries.');
end

function test_buildpath(testCase)
   %TEST_BUILDPATH Verify internal.fullpath path composition.
   demofilepath = icemodel.internal.fullpath('demo');
   [parentpath, foldername] = fileparts(demofilepath);

   testCase.verifyEqual(parentpath, icemodel.internal.fullpath(), ...
      'Expected fullpath() to return the full path to the demo/ folder.');

   testCase.verifyEqual(foldername, 'demo', ...
      'Expected fullpath() to return foldername "demo".');
end

function test_config_cases_resolve_owned_data_roots(testCase)
   %TEST_CONFIG_CASES_RESOLVE_OWNED_DATA_ROOTS Verify the three scoped cases.

   % Resolve without mutating the test process so the case mapping itself is the
   % only behavior under test.
   demo_cfg = icemodel.config('casename', 'demo', 'setenv', false);
   test_cfg = icemodel.config('casename', 'test', 'setenv', false);
   verification_cfg = icemodel.config( ...
      'casename', 'verification', 'setenv', false);

   % Each case must own a distinct, centralized data root.
   testCase.verifyEqual(string(demo_cfg.ICEMODEL_DATA_PATH), ...
      string(icemodel.internal.fullpath('demo', 'data')));
   testCase.verifyEqual(string(test_cfg.ICEMODEL_DATA_PATH), ...
      string(icemodel.internal.fullpath('test', 'data')));
   testCase.verifyEqual(string(verification_cfg.ICEMODEL_DATA_PATH), ...
      string(icemodel.internal.fullpath('data')));
end

function test_getpath_demo_resolves_demo_root(testCase)
   %TEST_GETPATH_DEMO_RESOLVES_DEMO_ROOT Verify the demo path getter.

   expected = icemodel.internal.fullpath('demo');
   returned = icemodel.getpath('demo');
   spectral = fullfile(returned, 'data', 'input', 'spectral');

   testCase.verifyEqual(returned, expected);
   testCase.verifyEqual(spectral, ...
      fullfile(expected, 'data', 'input', 'spectral'));
end

% function test_docpath(testCase)
%    docfilename = icemodel.internal.fullpath('icemodel_gettingStarted');
%    [foldername, filename, fileext] = fileparts(docfilename);
%
%    testCase.verifyEqual(foldername, icemodel.internal.fullpath('doc', 'html'), ...
%       'Expected docpath() to return doc/html folder.');
%
%    testCase.verifyEqual(filename, 'icemodel_gettingStarted', ...
%       'Expected docpath() to return icemodel_gettingStarted.html file.');
%
%    testCase.verifyEqual(fileext, '.html', ...
%       'Expected docpath() to return a .html file.');
% end

function test_version(testCase)
   %TEST_VERSION Verify the CFF-backed default and process-local override.
   cleanup = onCleanup(@() icemodel.internal.version('reset'));

   % A reset must resolve the persisted version instead of a hard-coded value.
   version = icemodel.internal.version('reset');
   testCase.verifyEqual(version, '1.1');
   testCase.verifyTrue(ischar(version));

   % A character-row override must remain active until the next reset.
   testCase.verifyEqual(icemodel.internal.version('9.8.7-test'), ...
      '9.8.7-test');
   testCase.verifyEqual(icemodel.internal.version(), '9.8.7-test');

   % Reset must discard the override and reread the citation source.
   testCase.verifyEqual(icemodel.internal.version('reset'), '1.1');
end

function test_version_failedResetClearsOverride(testCase)
   %TEST_VERSION_FAILEDRESETCLEARSOVERRIDE Reject stale state after failure.

   % Shadow the two version functions from a disposable release-like layout so
   % the failure test never removes the repository's live CITATION.cff.
   fixture_root = string(tempname);
   fixture_source = fullfile(fixture_root, "icemodel");
   fixture_internal = fullfile( ...
      fixture_source, "+icemodel", "+internal");
   mkdir(fixture_internal);
   fixture_cleanup = onCleanup(@() rmdir(fixture_root, 's'));
   copyfile(fullfile(icemodel.internal.fullpath("icemodel"), ...
      "+icemodel", "+internal", "version.m"), fixture_internal);
   copyfile(fullfile(icemodel.internal.fullpath("icemodel"), ...
      "+icemodel", "+internal", "readCffVersion.m"), fixture_internal);
   cffpath = fullfile(fixture_root, "CITATION.cff");
   writeText(cffpath, "version: 1.1");

   original_path = path;
   path_cleanup = onCleanup(@() restoreVersionPath(original_path));
   addpath(fixture_source, '-begin');
   clear('icemodel.internal.version', 'icemodel.internal.readCffVersion')
   testCase.assertEqual( ...
      string(which('icemodel.internal.version')), ...
      fullfile(fixture_internal, "version.m"));

   % Once reset fails, a subsequent lookup must retry the CFF and fail again
   % instead of returning the previous process-local override.
   testCase.verifyEqual(icemodel.internal.version('reset'), '1.1');
   testCase.verifyEqual(icemodel.internal.version('9.8.7-test'), ...
      '9.8.7-test');
   delete(cffpath);
   testCase.verifyError(@() icemodel.internal.version('reset'), ...
      'icemodel:internal:readCffVersion:FileNotFound');
   testCase.verifyError(@() icemodel.internal.version(), ...
      'icemodel:internal:readCffVersion:FileNotFound');
   clear path_cleanup
end

function test_releaseMetadata_resetsVersionForRelativeCff(testCase)
   %TEST_RELEASEMETADATA_RESETSVERSIONFORRELATIVECFF Reset alias-backed cache.

   % Shadow the release functions in a disposable repository layout so this
   % test can replace its canonical CFF without touching the live checkout.
   fixture_root = string(tempname);
   fixture_source = fullfile(fixture_root, "icemodel");
   fixture_internal = fullfile( ...
      fixture_source, "+icemodel", "+internal");
   mkdir(fixture_internal);
   fixture_cleanup = onCleanup(@() rmdir(fixture_root, 's'));
   source_internal = fullfile(icemodel.internal.fullpath("icemodel"), ...
      "+icemodel", "+internal");
   files = ["version.m", "readCffVersion.m", "releaseMetadata.m", ...
      "fullpath.m"];
   for k = 1:numel(files)
      copied = copyfile(fullfile(source_internal, files(k)), fixture_internal);
      testCase.assertTrue(copied);
   end
   fixture_private = fullfile(fixture_internal, "private");
   copied = copyfile(fullfile(source_internal, "private"), fixture_private);
   testCase.assertTrue(copied);
   copied = copyfile(icemodel.internal.fullpath("CITATION.cff"), fixture_root);
   testCase.assertTrue(copied);

   % A relative spelling of the canonical path must still clear an override
   % after prepare installs the validated CFF.
   original_path = path;
   path_cleanup = onCleanup(@() restoreVersionPath(original_path));
   original_folder = pwd;
   folder_cleanup = onCleanup(@() cd(original_folder));
   addpath(fixture_source, '-begin');
   cd(fixture_root)
   clear('icemodel.internal.fullpath', 'icemodel.internal.version', ...
      'icemodel.internal.readCffVersion', ...
      'icemodel.internal.releaseMetadata')
   testCase.assertEqual( ...
      string(which('icemodel.internal.releaseMetadata')), ...
      fullfile(fixture_internal, "releaseMetadata.m"));
   icemodel.internal.version('9.8.7-test');
   icemodel.internal.releaseMetadata("prepare", version="1.1", ...
      date_released="2026-07-30", cff_file=fullfile(".", "CITATION.cff"), ...
      preflight=@successfulPreflight, validator=@validateTestCff);
   testCase.verifyEqual(icemodel.internal.version(), '1.1');
   clear folder_cleanup path_cleanup
end

function test_readCffVersion(testCase)
   %TEST_READCFFVERSION Verify CFF scalar parsing and actionable failures.
   cffpath = string(tempname) + ".cff";
   cleanup = onCleanup(@() delete(cffpath));

   % The canonical unquoted scalar may include a trailing YAML comment.
   writeText(cffpath, sprintf( ...
      'cff-version: 1.2.0\nversion: 2.3.4 # candidate\n'));
   testCase.verifyEqual(icemodel.internal.readCffVersion(char(cffpath)), ...
      '2.3.4');

   % Optional key spacing remains valid, and matching quotes preserve a hash
   % plus YAML's doubled single-quote escape.
   writeText(cffpath, sprintf( ...
      'cff-version: 1.2.0\nversion : ''2 # 4''''beta'' # candidate\n'));
   testCase.verifyEqual(icemodel.internal.readCffVersion(cffpath), ...
      '2 # 4''beta');
   writeText(cffpath, sprintf( ...
      'cff-version: 1.2.0\nversion: "2 # 5" # candidate\n'));
   testCase.verifyEqual(icemodel.internal.readCffVersion(cffpath), '2 # 5');
   writeText(cffpath, 'version: "null"');
   testCase.verifyEqual(icemodel.internal.readCffVersion(cffpath), 'null');
   writeText(cffpath, sprintf( ...
      'metadata:\n  version: 9.9\nversion: 2.5\n'));
   testCase.verifyEqual(icemodel.internal.readCffVersion(cffpath), '2.5');

   % A missing, duplicate, empty, or partially quoted version must fail with a
   % specific contract error instead of silently choosing a value.
   writeText(cffpath, 'cff-version: 1.2.0');
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(cffpath), ...
      'icemodel:internal:readCffVersion:MissingVersion');
   writeText(cffpath, sprintf('metadata:\n  version: 9.9\n'));
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(cffpath), ...
      'icemodel:internal:readCffVersion:MissingVersion');
   writeText(cffpath, 'version:2.5');
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(cffpath), ...
      'icemodel:internal:readCffVersion:MissingVersion');
   writeText(cffpath, sprintf('version: 2.5\nversion: 2.6\n'));
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(cffpath), ...
      'icemodel:internal:readCffVersion:DuplicateVersion');
   writeText(cffpath, 'version:');
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(cffpath), ...
      'icemodel:internal:readCffVersion:InvalidVersion');
   writeText(cffpath, 'version: # missing');
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(cffpath), ...
      'icemodel:internal:readCffVersion:InvalidVersion');
   writeText(cffpath, 'version: null');
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(cffpath), ...
      'icemodel:internal:readCffVersion:InvalidVersion');
   writeText(cffpath, 'version: ~');
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(cffpath), ...
      'icemodel:internal:readCffVersion:InvalidVersion');
   writeText(cffpath, 'version: ''2.5');
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(cffpath), ...
      'icemodel:internal:readCffVersion:InvalidVersion');
   writeText(cffpath, 'version: "2.5');
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(cffpath), ...
      'icemodel:internal:readCffVersion:InvalidVersion');

   % A nonexistent path must report the file contract directly.
   missing = string(tempname) + ".cff";
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(missing), ...
      'icemodel:internal:readCffVersion:FileNotFound');

   % Nonscalar strings and nontext values must fail before file access.
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(["a", "b"]), ...
      'icemodel:internal:readCffVersion:InvalidPath');
   testCase.verifyError( ...
      @() icemodel.internal.readCffVersion(42), ...
      'icemodel:internal:readCffVersion:InvalidPath');
end

function test_releaseMetadata_prepare(testCase)
   %TEST_RELEASEMETADATA_PREPARE Verify staged CFF updates and rollback.

   cffpath = copyCffFixture(testCase);
   spaced_text = regexprep(fileread(cffpath), ...
      '(?m)^(version|date-released):', '$1 :');
   writeText(cffpath, spaced_text);
   result = icemodel.internal.releaseMetadata("prepare", ...
      version="1.2", date_released="2026-08-01", cff_file=cffpath, ...
      preflight=@successfulPreflight, validator=@validateTestCff);

   % Prepare must update both persisted scalars and report the ready state.
   testCase.verifyEqual(result.status, "ready");
   testCase.verifyEqual(result.version, "1.2");
   testCase.verifyEqual( ...
      icemodel.internal.readCffVersion(cffpath), '1.2');
   testCase.verifySubstring(fileread(cffpath), ...
      "date-released: '2026-08-01'");

   % A rejected preflight or validator must leave the live CFF untouched.
   before = fileread(cffpath);
   testCase.verifyError(@() icemodel.internal.releaseMetadata("prepare", ...
      version="1.3", date_released="2026-08-02", cff_file=cffpath, ...
      preflight=@failedPreflight, validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:preconditionFailed');
   testCase.verifyEqual(fileread(cffpath), before);
   testCase.verifyError(@() icemodel.internal.releaseMetadata("prepare", ...
      version="1.3", date_released="2026-08-02", cff_file=cffpath, ...
      preflight=@successfulPreflight, validator=@rejectValidation), ...
      'test:releaseMetadata:validationRejected');
   testCase.verifyEqual(fileread(cffpath), before);
end

function test_releaseMetadata_defaultPreflight(testCase)
   %TEST_RELEASEMETADATA_DEFAULTPREFLIGHT Verify Git release-state checks.

   project_dir = string(tempname) + " $ICEMODEL_RELEASE_QUOTE 'safe'";
   mkdir(project_dir);
   cleanup = onCleanup(@() rmdir(project_dir, 's'));
   cffpath = fullfile(project_dir, "CITATION.cff");
   copied = copyfile(icemodel.internal.fullpath("CITATION.cff"), cffpath);
   testCase.assertTrue(copied);
   runTestGit(project_dir, "init -q");
   runTestGit(project_dir, "checkout -qb dev");
   runTestGit(project_dir, "add CITATION.cff");
   runTestGit(project_dir, ...
      "-c user.name=Test -c user.email=test@example.invalid " ...
      + "commit -qm initial");
   remote_dir = string(tempname) + " $ICEMODEL_RELEASE_REMOTE 'safe'";
   runTestGit(project_dir, ...
      "init -q --bare " + quoteTestShellPath(remote_dir));
   remote_cleanup = onCleanup(@() rmdir(remote_dir, 's'));
   runTestGit(project_dir, ...
      "remote add origin " + quoteTestShellPath(remote_dir));

   % A clean dev branch and unused tag pass even when the path contains shell
   % expansion characters.
   original_folder = pwd;
   folder_cleanup = onCleanup(@() cd(original_folder));
   cd(project_dir)
   ready = icemodel.internal.releaseMetadata("prepare", ...
      version="1.2", date_released="2026-08-01", ...
      cff_file="CITATION.cff", validator=@validateTestCff);
   clear folder_cleanup
   testCase.verifyEqual(ready.status, "ready");
   runTestGit(project_dir, "add CITATION.cff");
   runTestGit(project_dir, ...
      "-c user.name=Test -c user.email=test@example.invalid " ...
      + "commit -qm prepared");

   % Existing tags, untracked files, and non-dev branches each block prepare.
   runTestGit(project_dir, "tag v1.3");
   testCase.verifyError(@() icemodel.internal.releaseMetadata("prepare", ...
      version="1.3", date_released="2026-08-02", cff_file=cffpath, ...
      validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:tagExists');
   runTestGit(project_dir, "tag -d v1.3");
   runTestGit(project_dir, "tag v1.4");
   runTestGit(project_dir, "push -q origin refs/tags/v1.4");
   runTestGit(project_dir, "tag -d v1.4");
   testCase.verifyError(@() icemodel.internal.releaseMetadata("prepare", ...
      version="1.4", date_released="2026-08-02", cff_file=cffpath, ...
      validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:tagExists');
   writeText(fullfile(project_dir, "untracked.txt"), "not releasable");
   testCase.verifyError(@() icemodel.internal.releaseMetadata("prepare", ...
      version="1.5", date_released="2026-08-02", cff_file=cffpath, ...
      validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:trackedChanges');
   delete(fullfile(project_dir, "untracked.txt"));
   runTestGit(project_dir, "checkout -qb main");
   testCase.verifyError(@() icemodel.internal.releaseMetadata("prepare", ...
      version="1.5", date_released="2026-08-02", cff_file=cffpath, ...
      validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:notDevBranch');
end

function test_releaseMetadata_observe(testCase)
   %TEST_RELEASEMETADATA_OBSERVE Verify ready and bounded pending results.

   cffpath = copyCffFixture(testCase);
   ready = runSingleObservation(cffpath, @successfulReleaseFetcher);
   testCase.verifyEqual(ready.status, "ready");
   testCase.verifyTrue(ready.github_release_seen);
   testCase.verifyTrue(ready.zenodo_record_seen);
   testCase.verifyEqual(ready.version_doi, ...
      "10.5281/zenodo.99999999");
   testCase.verifyEqual(ready.attempts, 1);

   % Even a successful final response remains pending if it reaches the hard
   % deadline before the command can report ready.
   clock_calls = 0;
   at_deadline = icemodel.internal.releaseMetadata("observe", ...
      version="1.1", cff_file=cffpath, ...
      fetcher=@successfulReleaseFetcher, timeout_seconds=1, ...
      waiter=@failWait, clock=@deadlineClock);
   testCase.verifyEqual(at_deadline.status, "pending");
   testCase.verifyTrue(at_deadline.github_release_seen);
   testCase.verifyTrue(at_deadline.zenodo_record_seen);

   % A deterministic clock proves one bounded wait before the command returns
   % resumable pending state.
   clock_value = 0;
   wait_count = 0;
   pending = icemodel.internal.releaseMetadata("observe", ...
      version="1.1", cff_file=cffpath, ...
      fetcher=@pendingReleaseFetcher, timeout_seconds=2, ...
      poll_interval_seconds=1, waiter=@advanceWait, clock=@advanceClock);
   testCase.verifyEqual(pending.status, "pending");
   testCase.verifyEqual(pending.attempts, 1);
   testCase.verifyEqual(wait_count, 1);
   testCase.verifyFalse(pending.github_release_seen);
   testCase.verifyFalse(pending.zenodo_record_seen);

   function value = advanceClock()
      %ADVANCECLOCK Return the deterministic current polling time.
      value = clock_value;
   end

   function advanceWait(delay)
      %ADVANCEWAIT Record the single bounded wait requested by observe mode.
      testCase.verifyEqual(delay, 1);
      wait_count = wait_count + 1;
      clock_value = 2;
   end

   function value = deadlineClock()
      %DEADLINECLOCK Reach the deadline only after both public responses.
      clock_calls = clock_calls + 1;
      value = double(clock_calls >= 4);
   end

   function failWait(~)
      %FAILWAIT Fail if a deadline result unexpectedly requests a wait.
      error('test:releaseMetadata:unexpectedWait', ...
         'observe unexpectedly requested a wait')
   end
end

function test_releaseMetadata_observeFailures(testCase)
   %TEST_RELEASEMETADATA_OBSERVEFAILURES Verify public-response error states.

   cffpath = copyCffFixture(testCase);

   % A public draft and an HTTP 404 are valid pending states.
   draft = runSingleObservation(cffpath, @draftReleaseFetcher);
   testCase.verifyEqual(draft.status, "pending");
   prerelease = runSingleObservation( ...
      cffpath, @prereleaseReleaseFetcher);
   testCase.verifyEqual(prerelease.status, "pending");
   not_found = runSingleObservation(cffpath, @notFoundReleaseFetcher);
   testCase.verifyEqual(not_found.status, "pending");

   % Malformed, mismatched, and failed responses must never be called pending.
   malformed = runSingleObservation(cffpath, @malformedReleaseFetcher);
   testCase.verifyEqual(malformed.status, "error");
   testCase.verifyEqual(malformed.error_identifier, ...
      "icemodel:internal:releaseMetadata:githubMalformed");
   wrong_tag = runSingleObservation(cffpath, @wrongTagReleaseFetcher);
   testCase.verifyEqual(wrong_tag.status, "error");
   testCase.verifyEqual(wrong_tag.error_identifier, ...
      "icemodel:internal:releaseMetadata:githubTagMismatch");
   failed = runSingleObservation(cffpath, @failedReleaseFetcher);
   testCase.verifyEqual(failed.status, "error");
   testCase.verifyEqual(failed.error_identifier, ...
      "test:releaseMetadata:httpFailure");
   wrong_concept = runSingleObservation( ...
      cffpath, @wrongConceptReleaseFetcher);
   testCase.verifyEqual(wrong_concept.status, "error");
   testCase.verifyEqual(wrong_concept.error_identifier, ...
      "icemodel:internal:releaseMetadata:zenodoConceptMismatch");
end

function test_releaseMetadata_finalize(testCase)
   %TEST_RELEASEMETADATA_FINALIZE Verify DOI lineage and idempotent insertion.

   cffpath = copyCffFixture(testCase);
   version_doi = "10.5281/zenodo.99999999";
   first = icemodel.internal.releaseMetadata("finalize", ...
      version="1.1", version_doi=version_doi, cff_file=cffpath, ...
      fetcher=@successfulReleaseFetcher);
   testCase.verifyFalse(first.added_concept);
   testCase.verifyTrue(first.added_version);
   once = fileread(cffpath);
   testCase.verifyEqual(count(once, version_doi), 1);

   % Repeating finalize must validate but leave the CFF bytes unchanged.
   second = icemodel.internal.releaseMetadata("finalize", ...
      version="1.1", version_doi=version_doi, cff_file=cffpath, ...
      fetcher=@successfulReleaseFetcher);
   testCase.verifyFalse(second.added_concept);
   testCase.verifyFalse(second.added_version);
   testCase.verifyEqual(fileread(cffpath), once);

   % A CFF without an identifiers block receives both CFF-form DOI entries.
   bare_cff = copyCffFixture(testCase);
   bare_text = regexprep(fileread(bare_cff), ...
      '(?ms)^identifiers:.*?(?=^repository-code:)', '');
   writeText(bare_cff, bare_text);
   added = icemodel.internal.releaseMetadata("finalize", ...
      version="1.1", version_doi=version_doi, cff_file=bare_cff, ...
      fetcher=@successfulReleaseFetcher);
   testCase.verifyTrue(added.added_concept);
   testCase.verifyTrue(added.added_version);
   final_text = fileread(bare_cff);
   testCase.verifyEqual(count(final_text, ...
      "10.5281/zenodo.11539329"), 1);
   testCase.verifyEqual(count(final_text, version_doi), 1);

   % An identifiers block at EOF without a newline still receives a separated
   % sequence item and remains parseable.
   eof_cff = string(tempname) + ".cff";
   eof_cleanup = onCleanup(@() delete(eof_cff));
   canonical = fileread(icemodel.internal.fullpath("CITATION.cff"));
   identifier_block = regexp(canonical, ...
      '(?ms)^identifiers:.*?(?=^repository-code:)', 'match', 'once');
   without_identifiers = regexprep(canonical, ...
      '(?ms)^identifiers:.*?(?=^repository-code:)', '');
   eof_text = [regexprep(without_identifiers, '[\r\n]+$', ''), ...
      newline, regexprep(identifier_block, '[\r\n]+$', '')];
   writeText(eof_cff, eof_text);
   eof_result = icemodel.internal.releaseMetadata("finalize", ...
      version_doi=version_doi, cff_file=eof_cff, ...
      fetcher=@successfulReleaseFetcher);
   testCase.verifyTrue(eof_result.added_version);
   testCase.verifySubstring(fileread(eof_cff), ...
      sprintf('11539329\n  - description:'));
end

function test_releaseMetadata_rejectsInvalidInput(testCase)
   %TEST_RELEASEMETADATA_REJECTSINVALIDINPUT Verify mutation boundaries.

   cffpath = copyCffFixture(testCase);
   testCase.verifyError(@() icemodel.internal.releaseMetadata("prepare", ...
      cff_file=cffpath, preflight=@successfulPreflight, ...
      validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:versionRequired');
   testCase.verifyError(@() icemodel.internal.releaseMetadata("prepare", ...
      version="release-one", cff_file=cffpath, ...
      preflight=@successfulPreflight, validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:versionInvalid');
   testCase.verifyError(@() icemodel.internal.releaseMetadata("prepare", ...
      version="1.2", date_released="2026-02-30", cff_file=cffpath, ...
      preflight=@successfulPreflight, validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:dateInvalid');
   testCase.verifyError(@() icemodel.internal.releaseMetadata("observe", ...
      version="1.1", repository="not-a-repository", cff_file=cffpath, ...
      fetcher=@successfulReleaseFetcher, timeout_seconds=1), ...
      'icemodel:internal:releaseMetadata:repositoryInvalid');
   testCase.verifyError(@() icemodel.internal.releaseMetadata("observe", ...
      version="2.0", cff_file=cffpath, ...
      fetcher=@successfulReleaseFetcher, timeout_seconds=1), ...
      'icemodel:internal:releaseMetadata:cffVersionMismatch');
   testCase.verifyError(@() icemodel.internal.releaseMetadata("finalize", ...
      version="1.1", version_doi="not-a-doi", cff_file=cffpath, ...
      fetcher=@successfulReleaseFetcher, validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:doiInvalid');
   testCase.verifyError(@() icemodel.internal.releaseMetadata("finalize", ...
      version="1.1", version_doi="10.5281/zenodo.11539329", ...
      cff_file=cffpath, fetcher=@successfulReleaseFetcher, ...
      validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:versionDoiIsConceptDoi');
   testCase.verifyError(@() icemodel.internal.releaseMetadata("finalize", ...
      version="1.1", version_doi="10.5281/zenodo.88888888", ...
      cff_file=cffpath, fetcher=@successfulReleaseFetcher, ...
      validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:versionDoiNotFound');
   testCase.verifyError(@() icemodel.internal.releaseMetadata("finalize", ...
      version="1.1", version_doi="10.5281/zenodo.99999999", ...
      cff_file=cffpath, fetcher=@wrongVersionReleaseFetcher, ...
      validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:zenodoVersionMismatch');
   testCase.verifyError(@() icemodel.internal.releaseMetadata("finalize", ...
      version="1.1", version_doi="10.5281/zenodo.99999999", ...
      cff_file=cffpath, fetcher=@wrongConceptReleaseFetcher, ...
      validator=@validateTestCff), ...
      'icemodel:internal:releaseMetadata:zenodoConceptMismatch');

   % The real validator command must reject malformed staged YAML without
   % changing the live CFF; the fixture path also exercises shell quoting.
   invalid_cff = copyCffFixture(testCase);
   invalid_text = string(fileread(invalid_cff)) + newline + "malformed: [";
   writeText(invalid_cff, invalid_text);
   before = fileread(invalid_cff);
   testCase.verifyError(@() icemodel.internal.releaseMetadata("prepare", ...
      version="1.2", date_released="2026-08-01", cff_file=invalid_cff, ...
      preflight=@successfulPreflight), ...
      'icemodel:internal:releaseMetadata:cffInvalid');
   testCase.verifyEqual(fileread(invalid_cff), before);

   % On Windows, reject percent-expanded command paths before validation can
   % inspect or replace the wrong file.
   if ispc
      percent_dir = string(tempname) + " %TEMP%";
      mkdir(percent_dir);
      percent_cleanup = onCleanup(@() rmdir(percent_dir, 's'));
      percent_cff = fullfile(percent_dir, "CITATION.cff");
      copied = copyfile(icemodel.internal.fullpath("CITATION.cff"), ...
         percent_cff);
      testCase.assertTrue(copied);
      before = fileread(percent_cff);
      testCase.verifyError(@() icemodel.internal.releaseMetadata("prepare", ...
         version="1.2", date_released="2026-08-01", ...
         cff_file=percent_cff, preflight=@successfulPreflight), ...
         'icemodel:internal:releaseMetadata:unsafeWindowsPath');
      testCase.verifyEqual(fileread(percent_cff), before);
   end
end

function cffpath = copyCffFixture(testCase)
   %COPYCFFFIXTURE Copy the repository CFF into one disposable test file.

   % Let the test framework remove the shell-sensitive path even after a
   % verification failure.
   cffpath = string(tempname) + " $ICEMODEL_CFF 'safe'.cff";
   copied = copyfile(icemodel.internal.fullpath("CITATION.cff"), cffpath);
   assert(copied == 1)
   testCase.addTeardown(@delete, cffpath);
end

function ok = successfulPreflight(~, ~)
   %SUCCESSFULPREFLIGHT Accept one mocked local release preflight.
   ok = true;
end

function ok = failedPreflight(~, ~)
   %FAILEDPREFLIGHT Reject one mocked local release preflight.
   ok = false;
end

function validateTestCff(filename)
   %VALIDATETESTCFF Exercise the required staged CFF version contract.
   icemodel.internal.readCffVersion(filename);
end

function rejectValidation(~)
   %REJECTVALIDATION Simulate a schema validator failure.
   error('test:releaseMetadata:validationRejected', ...
      'mock validator rejected the staged CFF')
end

function result = runSingleObservation(cffpath, fetcher)
   %RUNSINGLEOBSERVATION Execute one fetch attempt before the deadline.

   clock_value = 0;
   result = icemodel.internal.releaseMetadata("observe", ...
      version="1.1", cff_file=cffpath, fetcher=fetcher, ...
      timeout_seconds=1, poll_interval_seconds=1, ...
      waiter=@finishWait, clock=@currentTime);

   function value = currentTime()
      %CURRENTTIME Return the test's controlled wall-clock value.
      value = clock_value;
   end

   function finishWait(~)
      %FINISHWAIT Advance directly to the observation deadline.
      clock_value = 1;
   end
end

function response = successfulReleaseFetcher(url, ~)
   %SUCCESSFULRELEASEFETCHER Return matching public GitHub or Zenodo state.

   if contains(url, "api.github.com")
      response = struct("tag_name", "v1.1", "draft", false);
   else
      response = zenodoTestResponse("v1.1", ...
         "10.5281/zenodo.99999999", ...
         "10.5281/zenodo.11539329", "11539329");
   end
end

function response = pendingReleaseFetcher(url, ~)
   %PENDINGRELEASEFETCHER Return valid but incomplete public release state.

   if contains(url, "api.github.com")
      response = struct("not_found", true);
   else
      response = struct("hits", struct("hits", []));
   end
end

function response = draftReleaseFetcher(url, ~)
   %DRAFTRELEASEFETCHER Return a matching draft plus no Zenodo version.

   if contains(url, "api.github.com")
      response = struct("tag_name", "v1.1", "draft", true);
   else
      response = struct("hits", struct("hits", []));
   end
end

function response = prereleaseReleaseFetcher(url, ~)
   %PRERELEASERELEASEFETCHER Return a prerelease plus no Zenodo version.

   if contains(url, "api.github.com")
      response = struct("tag_name", "v1.1", "draft", false, ...
         "prerelease", true);
   else
      response = struct("hits", struct("hits", []));
   end
end

function response = notFoundReleaseFetcher(url, ~)
   %NOTFOUNDRELEASEFETCHER Simulate a GitHub 404 and empty Zenodo query.

   if contains(url, "api.github.com")
      error('test:releaseMetadata:HTTP404', 'HTTP 404 release not found')
   end
   response = struct("hits", struct("hits", []));
end

function response = malformedReleaseFetcher(url, ~)
   %MALFORMEDRELEASEFETCHER Return a malformed GitHub response.

   if contains(url, "api.github.com")
      response = struct("name", "missing tag");
   else
      response = struct("hits", struct("hits", []));
   end
end

function response = wrongTagReleaseFetcher(url, ~)
   %WRONGTAGRELEASEFETCHER Return an unexpected GitHub tag.

   if contains(url, "api.github.com")
      response = struct("tag_name", "v2.0", "draft", false);
   else
      response = struct("hits", struct("hits", []));
   end
end

function response = failedReleaseFetcher(~, ~)
   %FAILEDRELEASEFETCHER Simulate a non-404 public HTTP failure.

   % Build and throw one stable exception so the fetcher retains its required
   % output signature without a dead assignment.
   response = MException( ...
      'test:releaseMetadata:httpFailure', 'HTTP 503 unavailable');
   throw(response)
end

function response = wrongConceptReleaseFetcher(url, ~)
   %WRONGCONCEPTRELEASEFETCHER Return a version under another concept.

   if contains(url, "api.github.com")
      response = struct("tag_name", "v1.1", "draft", false);
   else
      response = zenodoTestResponse("v1.1", ...
         "10.5281/zenodo.99999999", ...
         "10.5281/zenodo.77777777", "77777777");
   end
end

function response = wrongVersionReleaseFetcher(~, ~)
   %WRONGVERSIONRELEASEFETCHER Return the requested DOI with another version.
   response = zenodoTestResponse("v2.0", ...
      "10.5281/zenodo.99999999", ...
      "10.5281/zenodo.11539329", "11539329");
end

function response = zenodoTestResponse(version, doi, concept_doi, concept_id)
   %ZENODOTESTRESPONSE Build one public-record-shaped Zenodo response.

   record = struct("doi", doi, "conceptdoi", concept_doi, ...
      "conceptrecid", concept_id, ...
      "metadata", struct("version", version));
   response = struct("hits", struct("hits", record));
end

function runTestGit(project_dir, arguments)
   %RUNTESTGIT Run one checked Git command inside a disposable fixture.

   command = "git -C " + quoteTestShellPath(project_dir) + " " + arguments;
   [status, output] = system(command);
   assert(status == 0, '%s', output)
end

function quoted = quoteTestShellPath(value)
   %QUOTETESTSHELLPATH Protect the disposable Git fixture path.

   value = char(string(value));
   if ispc
      quoted = string([char(34), value, char(34)]);
      return
   end
   embedded_quote = char([39, 34, 39, 34, 39]);
   value = strrep(value, char(39), embedded_quote);
   quoted = string([char(39), value, char(39)]);
end

function restoreVersionPath(original_path)
   %RESTOREVERSIONPATH Restore package resolution after the shadowed test.

   path(original_path);
   clear('icemodel.internal.fullpath', 'icemodel.internal.version', ...
      'icemodel.internal.readCffVersion', ...
      'icemodel.internal.releaseMetadata')
   icemodel.internal.version('reset');
end

function writeText(filename, text)
   %WRITETEXT Write an R2017a-compatible temporary text fixture.

   % Close the fixture deterministically even when a parser assertion fails.
   file_id = fopen(char(filename), 'w');
   cleanup = onCleanup(@() fclose(file_id));
   fprintf(file_id, '%s', text);
end
