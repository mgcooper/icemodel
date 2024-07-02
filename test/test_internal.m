function tests = test_internal
   %TEST_INTERNAL Test the toolbox internal functions.
   tests = functiontests(localfunctions);
end

function setup(testCase) %#ok<INUSD>

end

function teardown(testCase) %#ok<INUSD>

end

function test_basepath(testCase)
   modelpath = fullfile(icemodel.internal.fullpath(), 'icemodel');
   
   testCase.verifyTrue(isfolder(modelpath), ...
      'Expected icemodel/ folder to exist in the top level project folder.');
end

function test_functionSignatures(testCase)
   % Validate the JSON and get the table
   T = validateFunctionSignaturesJSON(fullfile( ...
      icemodel.internal.fullpath(), 'icemodel', 'functionSignatures.json'));

   % Check if the table is empty
   testCase.verifyEmpty(T, ...
      'The functionSignatures.json file contains invalid entries.');
end

function test_buildpath(testCase)
   demofilepath = icemodel.internal.fullpath('demo');
   [parentpath, foldername] = fileparts(demofilepath);
   
   testCase.verifyEqual(parentpath, icemodel.internal.fullpath(), ...
      'Expected fullpath() to return the full path to the demo/ folder.');
   
   testCase.verifyEqual(foldername, 'demo', ...
      'Expected fullpath() to return foldername "demo".');
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
   version = icemodel.internal.version();
   testCase.verifyTrue(ischar(version))
end
