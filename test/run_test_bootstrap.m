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
   %  results = run_test_bootstrap(baseline_tag="v1.1", smbmodel="icemodel", ...
   %     solver=2)
   %  results = run_test_bootstrap(baseline_tag="v1.1", smbmodel="icemodel", ...
   %     solver=[1 3])
   %  results = run_test_bootstrap(simyear=2017, smoke_sites="kanm", ...
   %     full_sites=["kanm"; "kanl"])
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
   %  - SMOKE_SITES and FULL_SITES are advanced overrides for the site lists
   %    used by each formal tier.
   %  - The optional solver filter accepts any subset of [1 2 3].
   %  - Compare runs are read-only with respect to baselines. Baseline updates
   %    happen through the explicit build/snapshot actions above.

   arguments (Input)

      kwargs.baseline_tag (1, :) string ...
         = "v1.1"

      kwargs.smbmodel (1, :) string ...
         {icemodel.validators.mustBeTestSmbmodelSelector(kwargs.smbmodel)} ...
         = "all"

      kwargs.solver {icemodel.validators.mustBeSolverFilter(kwargs.solver)} ...
         = []

      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} ...
         = 2016

      kwargs.smoke_sites string ...
         = "kanm"

      kwargs.full_sites string ...
         = ["kanm"; "kanl"]

      kwargs.clean_artifacts (1, 1) logical ...
         = false

      kwargs.clean_baselines (1, 1) logical ...
         = false

      kwargs.backup_before_clean (1, 1) logical ...
         = true
   end

   % Deal out arguments.
   [baseline_tag, smbmodel, solver, simyear, smoke_sites, full_sites, ...
      clean_artifacts, clean_baselines, backup_before_clean] = deal( ...
      kwargs.baseline_tag, kwargs.smbmodel, kwargs.solver, kwargs.simyear, ...
      reshape(kwargs.smoke_sites, [], 1), reshape(kwargs.full_sites, [], 1), ...
      kwargs.clean_artifacts, kwargs.clean_baselines, kwargs.backup_before_clean);

   % Bootstrap the broad source/test trees once for the full workflow.
   % Keep the cleanup handle in scope so the caller's config is restored
   % when this entrypoint returns.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Set the directory paths
   rootdir = icemodel.getpath('test');
   baselinesdir = fullfile(rootdir, 'baselines');
   artifactsdir = fullfile(rootdir, 'artifacts');
   backupzipdir = fullfile(rootdir, 'backups');

   % Create the YYYYMMDD-HHMMSS run_name
   [~, ~, run_name] = icemodel.test.helpers.resolveRunStamp(string.empty());

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
   results.solver = solver;
   results.simyear = simyear;
   results.smoke_sites = smoke_sites;
   results.full_sites = full_sites;

   % Generate the canonical list of all formal suite cases and run them.
   cases = icemodel.test.helpers.getFormalTestSuiteCases();
   for i = 1:height(cases)
      c = cases(i, :);
      results.(c.result_field) = runStep(c, baseline_tag, smbmodel, solver, ...
         simyear, smoke_sites, full_sites, run_name);
   end

   results.regression_passed = collectPassFlags(results, cases, "regression");
   results.perf_passed = collectPassFlags(results, cases, "perf");
end

function out = runStep(c, baseline_tag, smbmodel, solver, simyear, ...
      smoke_sites, full_sites, run_name)
   %RUNSTEP Execute one suite lifecycle step from the bootstrap matrix.
   switch c.suite
      case "regression"
         switch c.action
            case "build"
               % Accept the current modeled outputs as the rolling baseline.
               out = build_regression_baseline( ...
                  baseline="rolling", tier=c.tier, smbmodel=smbmodel, ...
                  solver=solver, simyear=simyear, smoke_sites=smoke_sites, ...
                  full_sites=full_sites);

            case "snapshot"
               % Freeze the current rolling baseline into a versioned release.
               out = snapshot_regression_baseline( ...
                  baseline_tag=baseline_tag, smbmodel=smbmodel, ...
                  overwrite=true);

            case "run"
               % Compare current outputs against the requested baseline.
               out = run_regression_suite( ...
                  tier=c.tier, smbmodel=smbmodel, solver=solver, ...
                  simyear=simyear, smoke_sites=smoke_sites, ...
                  full_sites=full_sites, ...
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
                  baseline="rolling", tier=c.tier, smbmodel=smbmodel, ...
                  solver=solver, simyear=simyear, ...
                  smoke_sites=smoke_sites, full_sites=full_sites);

            case "snapshot"
               % Freeze the current rolling baseline into a versioned release.
               out = snapshot_perf_baseline( ...
                  baseline_tag=baseline_tag, smbmodel=smbmodel, ...
                  overwrite=true);

            case "run"
               % Compare current runtimes against the requested baseline.
               out = run_perf_suite( ...
                  tier=c.tier, smbmodel=smbmodel, solver=solver, ...
                  simyear=simyear, smoke_sites=smoke_sites, ...
                  full_sites=full_sites, ...
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
   %COLLECTPASSFLAGS Extract pass/fail flags from one suite result struct.
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
   %RESOLVEBASELINE Resolve the baseline selector for one bootstrap step.
   if baseline_mode == "rolling"
      baseline = "rolling";
   else
      baseline = baseline_tag;
   end
end

function tf = hasContents(folder)
   %HASCONTENTS Return true when a folder exists and is non-empty.
   listing = dir(folder);
   tf = any(~ismember(string({listing.name}), [".", ".."]));
end

function removeFolder(folder)
   %REMOVEFOLDER Remove one folder if it exists.
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
   %REMOVEMATFILES Remove top-level MAT files from a folder.
   files = dir(fullfile(folder, '*.mat'));
   for i = 1:numel(files)
      delete(fullfile(files(i).folder, files(i).name));
   end
end

function backupToFolder(source, backupdir)
   %BACKUPTOFOLDER Copy one file or folder into a backup location.
   if exist(backupdir, 'dir') ~= 7
      mkdir(backupdir);
   end
   pathname = backupfile(source, true, true);
   [~, name, ext] = fileparts(pathname);
   movefile(pathname, fullfile(backupdir, [name ext]));
end
