function results = run_test_bootstrap(kwargs)
   %RUN_TEST_BOOTSTRAP Full formal-suite setup and refresh workflow.
   %
   % Use this function when you want one command to:
   %  1. optionally clean old formal artifacts and managed baselines
   %  2. rebuild the rolling regression baseline(s)
   %  3. snapshot release regression baseline(s)
   %  4. rebuild the rolling perf baseline(s)
   %  5. snapshot release perf baseline(s)
   %  6. run smoke/full regression compares against rolling and release
   %  7. run smoke/full perf compares against rolling and release
   %
   % This is the highest-level orchestration entrypoint in test/. It is mainly
   % for first-time setup, full refreshes, and larger development resets.
   %
   %  results = run_test_bootstrap()
   %  results = run_test_bootstrap(baseline_tag="v1.1")
   %  results = run_test_bootstrap(baseline_tag="v1.1", smbmodel="skinmodel")
   %  results = run_test_bootstrap(baseline_tag="v1.1", smbmodel="all", ...
   %     clean_artifacts=true, clean_baselines=true)
   %
   % Use the lower-level tools instead when you only need one part of the
   % workflow:
   %  - build_regression_baseline(...) : accept a new rolling or versioned
   %    regression baseline without running compare suites
   %  - build_perf_baseline(...)       : accept a new rolling or versioned
   %    perf baseline without running compare suites
   %  - snapshot_regression_baseline(...) / snapshot_perf_baseline(...)
   %    : freeze release baselines from the current rolling baselines
   %  - run_regression_suite(...) / run_perf_suite(...)
   %    : compare current code against existing baselines and write artifacts
   %
   % Cleanup behavior:
   %  - clean_artifacts=true empties test/artifacts/ completely
   %  - clean_baselines=true deletes all baseline .mat files in
   %    test/baselines/
   %  - backup_before_clean=true first creates zip backups in test/backups/
   %
   % Notes:
   %  - This bootstrap does not rebuild test/references/runoff_reference.mat.
   %  - smbmodel="all" is virtual: it rebuilds per-model files for each formal
   %    model and runs the union of those cases.
   %  - Compare runs are read-only with respect to baselines. Baseline updates
   %    happen through the explicit build/snapshot actions above.

   arguments (Input)
      kwargs.baseline_tag (1, :) string = "v1.1"
      kwargs.smbmodel (1, :) string {mustBeMember(kwargs.smbmodel, ...
         ["all", "icemodel", "skinmodel"])} = "all"
      kwargs.clean_artifacts (1, 1) logical = false
      kwargs.clean_baselines (1, 1) logical = false
      kwargs.backup_before_clean (1, 1) logical = true
   end
   [baseline_tag, smbmodel, clean_artifacts, clean_baselines, ...
      backup_before_clean] = deal( ...
      kwargs.baseline_tag, kwargs.smbmodel, kwargs.clean_artifacts, ...
      kwargs.clean_baselines, kwargs.backup_before_clean);

   % Ensure test/ and dependencies/ are on path
   rootdir = fileparts(mfilename('fullpath'));
   addpath(rootdir)
   addpath(fullfile(fileparts(rootdir), 'icemodel', 'dependencies'))

   % Set the directory paths
   baselinesdir = fullfile(rootdir, 'baselines');
   artifactsdir = fullfile(rootdir, 'artifacts');
   backupzipdir = fullfile(rootdir, 'backups');

   % Create the YYYYMMDD-HHMMSS run_name
   [~, ~, run_name] = test.helpers.resolveRunStamp(string.empty());

   % Backup and delete existing artifacts and baselines
   if clean_artifacts && exist(artifactsdir, 'dir') == 7
      if backup_before_clean && hasContents(artifactsdir)
         backupToFolder(artifactsdir, backupzipdir);
      end
      removeFolder(artifactsdir);
   end
   if clean_baselines
      if backup_before_clean && hasContents(baselinesdir)
         backupToFolder(baselinesdir, backupzipdir);
      end
      removeMatFiles(baselinesdir);
   end

   % Generate new baselines and artifacts covering the canonical workflow.
   results = struct();
   results.run_name = run_name;
   results.baseline_tag = baseline_tag;
   results.smbmodel = smbmodel;

   % Generate the canonical list of all formal suite cases and run them.
   cases = test.helpers.getFormalTestSuiteCases();
   for i = 1:height(cases)
      c = cases(i, :);
      results.(c.result_field) = runStep(c, baseline_tag, smbmodel, run_name);
   end

   results.regression_passed = collectPassFlags(results, cases, "regression");
   results.perf_passed = collectPassFlags(results, cases, "perf");
end

function out = runStep(c, baseline_tag, smbmodel, run_name)
   switch c.suite
      case "regression"
         switch c.action
            case "build"
               % Accept the current modeled outputs as the rolling baseline.
               out = build_regression_baseline( ...
                  baseline="rolling", tier=c.tier, smbmodel=smbmodel);

            case "snapshot"
               % Freeze the current rolling baseline into a versioned release.
               out = snapshot_regression_baseline( ...
                  baseline_tag=baseline_tag, smbmodel=smbmodel, ...
                  overwrite=true);

            case "run"
               % Compare current outputs against the requested baseline.
               out = run_regression_suite( ...
                  tier=c.tier, smbmodel=smbmodel, ...
                  baseline=resolveBaseline(c.baseline_mode, baseline_tag), ...
                  run_name=run_name);

            otherwise
               error('unsupported regression bootstrap action: %s', c.action)
         end

      case "perf"
         switch c.action
            case "build"
               % Accept the current runtimes as the rolling baseline.
               out = build_perf_baseline( ...
                  baseline="rolling", tier=c.tier, smbmodel=smbmodel);

            case "snapshot"
               % Freeze the current rolling baseline into a versioned release.
               out = snapshot_perf_baseline( ...
                  baseline_tag=baseline_tag, smbmodel=smbmodel, ...
                  overwrite=true);

            case "run"
               % Compare current runtimes against the requested baseline.
               out = run_perf_suite( ...
                  tier=c.tier, smbmodel=smbmodel, ...
                  baseline=resolveBaseline(c.baseline_mode, baseline_tag), ...
                  run_name=run_name);

            otherwise
               error('unsupported perf bootstrap action: %s', c.action)
         end

      otherwise
         error('unsupported bootstrap suite: %s', c.suite)
   end
end

function passed = collectPassFlags(results, cases, suite)
   run_cases = cases(cases.suite == suite & cases.action == "run", :);
   passed = struct();

   for i = 1:height(run_cases)
      key = run_cases.result_field(i);
      field = erase(key, suite + "_");
      if suite == "regression"
         passed.(field) = all([results.(key).Passed]);
      else
         passed.(field) = results.(key).passed;
      end
   end
end

function baseline = resolveBaseline(baseline_mode, baseline_tag)
   if baseline_mode == "rolling"
      baseline = "rolling";
   else
      baseline = baseline_tag;
   end
end

function tf = hasContents(folder)
   listing = dir(folder);
   tf = any(~ismember(string({listing.name}), [".", ".."]));
end

function removeFolder(folder)
   listing = dir(folder);
   listing = listing(~ismember(string({listing.name}), [".", ".."]));
   for i = 1:numel(listing)
      pathname = fullfile(listing(i).folder, listing(i).name);
      if listing(i).isdir
         rmdir(pathname, 's');
      else
         delete(pathname);
      end
   end
end

function removeMatFiles(folder)
   files = dir(fullfile(folder, '*.mat'));
   for i = 1:numel(files)
      delete(fullfile(files(i).folder, files(i).name));
   end
end

function backupToFolder(source, backupdir)
   if exist(backupdir, 'dir') ~= 7
      mkdir(backupdir);
   end
   pathname = backupfile(source, true, true);
   [~, name, ext] = fileparts(pathname);
   movefile(pathname, fullfile(backupdir, [name ext]));
end
