classdef IcemodelPerfTest < matlab.perftest.TestCase
   %ICEMODELPERFTEST Formal runtime benchmark for one model case.

   properties
      opts
   end

   methods (TestMethodSetup)
      function configureCase(testCase)
         % Resolve the configured formal test case before timing begins.
         rootdir = icemodel.internal.fullpath();
         addpath(fullfile(rootdir, 'test'));
         icemodel.test.helpers.configureModelPaths(rootdir);
         testCase.opts = buildCaseOpts();
      end
   end

   methods (Test)
      function testCoreRuntime(testCase)
         % Time only the core model call. Setup and reporting stay outside the
         % measured region.
         testCase.startMeasuring
         icemodel.test.helpers.runSmbModel(testCase.opts);
         testCase.stopMeasuring
      end
   end
end

function opts = buildCaseOpts()
   % Read one perf case definition from environment variables and resolve it
   % through the formal case builder. Perf opts always exclude spinup years so
   % the timed region reflects only the analyzed output years.
   smbmodel = getenvRequired('ICEMODEL_TEST_SMBMODEL');
   sitename = getenvRequired('ICEMODEL_TEST_SITENAME');
   forcings = getenvRequired('ICEMODEL_TEST_FORCINGS');
   userdata = getenv('ICEMODEL_TEST_USERDATA');
   uservars = getenv('ICEMODEL_TEST_USERVARS');
   simyear = str2double(getenvRequired('ICEMODEL_TEST_SIMYEAR'));
   simyears = getenv('ICEMODEL_TEST_SIMYEARS');
   n_spinup_years = str2double(getenvRequired('ICEMODEL_TEST_N_SPINUP_YEARS'));
   solver = str2double(getenvRequired('ICEMODEL_TEST_SOLVER'));

   c = struct('smbmodel', string(smbmodel), 'sitename', string(sitename), ...
      'forcings', string(forcings), 'userdata', string(userdata), ...
      'uservars', string(uservars), 'simyear', simyear, 'solver', solver);

   if ~isempty(simyears)
      c.simyears = parseSimyears(simyears);
      c.n_spinup_years = n_spinup_years;
   end
   opts = icemodel.test.helpers.setModelOptsForCase(c, include_spinup=false);
end

function s = getenvRequired(name)
   s = getenv(name);
   if isempty(s)
      error('missing required perf env var: %s', name)
   end
end

function years = parseSimyears(s)
   parts = split(string(s), ',');
   years = str2double(parts);
   if any(~isfinite(years))
      error('invalid ICEMODEL_TEST_SIMYEARS: %s', s)
   end
   years = years(:).';
end
