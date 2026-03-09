classdef IcemodelPerfTest < matlab.perftest.TestCase
   %ICEMODELPERFTEST Formal runtime benchmark for one model case.

   properties
      opts
   end

   methods (TestMethodSetup)
      function configureCase(testCase)
         % Resolve the configured formal test case before timing begins.
         rootdir = fileparts(fileparts(mfilename('fullpath')));
         addpath(fullfile(rootdir, 'test'));
         test.helpers.configureModelPaths(rootdir);
         testCase.opts = buildCaseOpts();
      end
   end

   methods (Test)
      function testCoreRuntime(testCase)
         % Time only the core model call. Setup and reporting stay outside the
         % measured region.
         testCase.startMeasuring
         runModel(testCase.opts);
         testCase.stopMeasuring
      end
   end
end

function opts = buildCaseOpts()
   % Read one perf case definition from environment variables and resolve it
   % through the formal case builder.
   sitename = getenvOrDefault('ICEMODEL_TEST_SITENAME', 'kanm');
   forcings = getenvOrDefault('ICEMODEL_TEST_FORCINGS', sitename);
   userdata = getenvOrDefault('ICEMODEL_TEST_USERDATA', '');
   uservars = getenvOrDefault('ICEMODEL_TEST_USERVARS', '');
   simyear = str2double(getenvOrDefault('ICEMODEL_TEST_SIMYEAR', '2016'));
   solver_mode = str2double(getenvOrDefault('ICEMODEL_TEST_SOLVER_MODE', '2'));
   smbmodel = getenvOrDefault('ICEMODEL_TEST_SMBMODEL', 'icemodel');
   c = struct('smbmodel', string(smbmodel), 'sitename', string(sitename), ...
      'forcings', string(forcings), 'userdata', string(userdata), ...
      'uservars', string(uservars), 'simyear', simyear, ...
      'solver_mode', solver_mode);
   opts = test.helpers.buildFormalCaseOpts(c);
end

function runModel(opts)
   % Dispatch to the requested model kernel.
   switch opts.smbmodel
      case 'icemodel'
         icemodel(opts);
      case 'skinmodel'
         skinmodel(opts);
      otherwise
         error('unsupported smbmodel: %s', opts.smbmodel)
   end
end

function s = getenvOrDefault(name, default)
   s = getenv(name);
   if isempty(s)
      s = default;
   end
end
