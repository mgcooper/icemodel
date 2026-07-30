function tests = test_demo_contracts
   %TEST_DEMO_CONTRACTS Exercise the three self-contained 15-minute demos.

   tests = functiontests(localfunctions);
end

function test_three_demo_option_bundles_use_15_minute_tree(testCase)
   % All public demo option bundles must run from only the scoped demo tree.

   % Install a caller config first so the demo scope and its cleanup are both
   % observable rather than relying on the process defaults.
   names = configNames();
   original = currentRawConfig(names);
   restore = onCleanup(@() setConfig(names, original));
   caller_root = string(tempname);
   mkdir(fullfile(caller_root, 'input'))
   mkdir(fullfile(caller_root, 'eval'))
   icemodel.config('ICEMODEL_DATA_PATH', caller_root);
   caller_values = currentRawConfig(names);
   caller_values{names == "ICEMODEL_VERSION"} = 'demo-caller-version';
   caller_values{names == "ICEMODEL_REFERENCE"} = 'demo-caller-reference';
   caller_values{names == "ICEMODEL_CONTACT"} = '';
   setConfig(names, caller_values)
   path_before = path;

   % Bootstrap the demo case and verify its exact tracked root.
   [rootdir, input_path, ~, ~, cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment( ...
      icemodel_config_casename="demo");
   testCase.verifyClass(cleanup, 'onCleanup')
   demo_root = string(fullfile(rootdir, 'demo', 'data'));
   testCase.verifyEqual(input_path, fullfile(demo_root, 'input'))

   % The primary forcing carries 15-minute KANM data and inline MODIS albedo.
   primary_file = fullfile(input_path, 'met', ...
      'met_kanm_kanm_2016_15m.mat');
   primary = load(primary_file, 'met');
   testCase.verifyEqual(seconds(median(diff(primary.met.Time))), 900)
   modis_index = find(strcmpi( ...
      primary.met.Properties.VariableNames, 'modis'), 1);
   testCase.verifyNotEmpty(modis_index)
   modis_name = primary.met.Properties.VariableNames{modis_index};
   testCase.verifyTrue(any(isfinite(primary.met.(modis_name))))

   % The external MERRA2 swap is the exact 2016 15-minute slice.
   merra2_file = fullfile(input_path, 'met', 'merra2', ...
      'met_kanm_merra2_20160101_20161231_15m.mat');
   merra2 = load(merra2_file, 'met');
   testCase.verifyEqual(height(merra2.met), 35136)
   testCase.verifyEqual(seconds(median(diff(merra2.met.Time))), 900)
   testCase.verifyTrue(all(year(merra2.met.Time) == 2016))

   % The tracked demo boundary contains only the two public forcing files and
   % four shared spectral tables; verification and formal data live elsewhere.
   entries = dir(fullfile(demo_root, '**', '*'));
   entries = entries(~[entries.isdir]);
   entries = entries(string({entries.name}) ~= ".DS_Store");
   actual = strings(numel(entries), 1);
   for n = 1:numel(entries)
      pathname = fullfile(entries(n).folder, entries(n).name);
      actual(n) = icemodel.verification.setup.fixtureRelativePosix( ...
         demo_root, pathname);
   end
   expected = [ ...
      "input/met/merra2/met_kanm_merra2_20160101_20161231_15m.mat"
      "input/met/met_kanm_kanm_2016_15m.mat"
      "input/spectral/kabs.mat"
      "input/spectral/kice.mat"
      "input/spectral/mie.mat"
      "input/spectral/solar.mat"];
   testCase.verifyEqual(sort(actual), sort(expected))

   % Run the three option bundles. Routine output is captured so the unit result
   % remains concise; a missing or wrongly selected asset still fails the call.
   sources = ["", "modis", "merra2"];
   variables = ["", "albedo", "tair"]; %#ok<NASGU> Used by the captured command below.
   for n = 1:numel(sources)
      args = {"saveflag", false, "sitename", "kanm", ...
         "forcings", "kanm", "smbmodel", "skinmodel", ...
         "simyears", 2016, "backupflag", false}; %#ok<NASGU> Used by evalc below.
      if strlength(sources(n)) > 0
         routine_output = evalc( ...
            ['[ice1, ice2, met, opts] = icemodel.run.point(args{:}, ' ...
            '"userdata", sources(n), "uservars", variables(n));']);
      else
         routine_output = evalc( ...
            '[ice1, ice2, met, opts] = icemodel.run.point(args{:});');
      end
      testCase.verifyClass(routine_output, 'char')
      testCase.verifyFalse(isempty(ice1))
      testCase.verifyFalse(isempty(ice2))
      testCase.verifyFalse(isempty(met))
      testCase.verifyEqual(opts.dt, 900)
      testCase.verifyTrue(all(startsWith(string(opts.metfname), input_path)))
      testCase.verifyTrue(all(contains(string(opts.metfname), '_15m.mat')))
   end

   % The former one-hour primary forcing must not remain as a hidden fallback.
   one_hour = fullfile(input_path, 'met', 'met_kanm_kanm_2016_1hr.mat');
   testCase.verifyFalse(isfile(one_hour))

   % Leaving the demo scope restores the interactive caller configuration.
   clear cleanup
   testCase.verifyEqual(currentRawConfig(names), caller_values)
   testCase.verifyEqual(path, path_before)
   clear restore
   rmdir(caller_root, 's')
end

function names = configNames()
   %CONFIGNAMES Return the complete environment-backed config field list.

   names = string(fieldnames(icemodel.config('setenv', false)));
end

function values = currentRawConfig(names)
   %CURRENTRAWCONFIG Snapshot exact raw values without applying defaults.

   values = cell(numel(names), 1);
   for n = 1:numel(names)
      values{n} = getenv(names(n));
   end
end

function setConfig(names, values)
   %SETCONFIG Restore raw environment-backed config values.

   for n = 1:numel(names)
      if iscell(values)
         value = values{n};
      else
         value = values(n);
      end
      setenv(names(n), value)
   end
end
