classdef IcemodelPerfTest < matlab.perftest.TestCase
   %ICEMODELPERFTEST Formal runtime benchmark for one model case.

   properties
      opts
      env_cleanup
   end

   methods (TestMethodSetup)
      function configureCases(testCase)
         % Install the canonical test config before timing begins so direct
         % class runs and runner-based runs see the same environment.
         [~, ~, ~, ~, testCase.env_cleanup] = ...
            icemodel.test.helpers.bootstrapTestEnvironment();
         testCase.opts = buildCaseOpts();
      end
   end

   methods (TestMethodTeardown)
      function restoreConfig(testCase)
         % Release the setup cleanup handle after each timed case.
         testCase.env_cleanup = [];
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
   %BUILDCASEOPTS Resolve one formal perf case from environment variables.

   % Read one perf case definition from environment variables and resolve it
   % through the shared formal case builder.
   smbmodel = getenvRequired('ICEMODEL_TEST_SMBMODEL');
   sitename = getenvRequired('ICEMODEL_TEST_SITENAME');
   forcings = getenvRequired('ICEMODEL_TEST_FORCINGS');
   userdata = getenv('ICEMODEL_TEST_USERDATA');
   uservars = getenv('ICEMODEL_TEST_USERVARS');
   simyear = str2double(getenvRequired('ICEMODEL_TEST_SIMYEAR'));
   simyears = getenv('ICEMODEL_TEST_SIMYEARS');
   n_spinup_years = str2double(getenvRequired('ICEMODEL_TEST_N_SPINUP_YEARS'));
   solver = str2double(getenvRequired('ICEMODEL_TEST_SOLVER'));
   % Rebuild the formal case struct and hand it to the shared opts helper.
   c = struct('smbmodel', string(smbmodel), 'sitename', string(sitename), ...
      'forcings', string(forcings), 'userdata', string(userdata), ...
      'uservars', string(uservars), 'simyear', simyear, 'solver', solver);

   if ~isempty(simyears)
      c.simyears = parseSimyears(simyears);
      c.n_spinup_years = n_spinup_years;
   end
   opts = icemodel.test.helpers.setModelOptsForCase(c);
end

function s = getenvRequired(name)
   %GETENVREQUIRED Read one required case-selector env var or error cleanly.
   %
   % The test bootstrap installs config paths, but the runner still owns the
   % concrete ICEMODEL_TEST_* case selector variables. Error immediately when
   % a direct class run or a runner regression forgets to install one.

   s = getenv(name);
   if isempty(s)
      error('missing required perf env var: %s', name)
   end
end

function years = parseSimyears(s)
   %PARSESIMYEARS Parse a comma-delimited SIMYEARS env var.

   parts = split(string(s), ',');
   years = str2double(parts);
   if any(~isfinite(years))
      error('invalid ICEMODEL_TEST_SIMYEARS: %s', s)
   end
   years = years(:).';
end
