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
   %TEST_VERSION Verify version resolution and parsing helpers.
   version = icemodel.internal.version();
   testCase.verifyTrue(ischar(version))
end
