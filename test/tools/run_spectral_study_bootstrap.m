function results = run_spectral_study_bootstrap(kwargs)
   %RUN_SPECTRAL_STUDY_BOOTSTRAP Run the spectral study workflow.
   %
   %  results = run_spectral_study_bootstrap()
   %  results = run_spectral_study_bootstrap( ...
   %     include_regression_acceptance=true, include_perf_acceptance=true)
   %  results = run_spectral_study_bootstrap( ...
   %     include_regression_acceptance=true, include_perf_acceptance=true, ...
   %     rebuild_regression_baseline=true, rebuild_perf_baseline=true)
   %
   % Use this tool when you want one entrypoint for the current spectral study
   % workflow. By default it runs the study diagnostics and saves their outputs
   % under:
   %
   %  /Users/mattcooper/MATLAB/projects/icemodel/test/benchmarks/reports/
   %  /Users/mattcooper/MATLAB/projects/icemodel/test/benchmarks/figures/
   %
   % The expensive formal acceptance and baseline rebuild steps are opt-in.
   % Those steps use the requested SPECTRAL_VARIANT for the formal suite.

   arguments
      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} ...
         = 2016

      kwargs.solver (1, 1) double {mustBeMember(kwargs.solver, [1 2 3])} ...
         = 2

      kwargs.smbmodel (1, :) string ...
         {icemodel.validators.mustBeTestSmbmodelSelector(kwargs.smbmodel)} ...
         = "all"

      kwargs.smoke_sites string ...
         = "kanm"

      kwargs.full_sites string ...
         = ["kanm"; "kanl"]

      kwargs.spectral_variant (1, :) string ...
         {mustBeMember(kwargs.spectral_variant, ...
         ["inlined", "functions", "lookup"])} ...
         = "lookup"

      kwargs.n_direct_runs (1, 1) double {mustBeInteger, mustBePositive} ...
         = 1

      kwargs.include_density_floor (1, 1) logical ...
         = true

      kwargs.include_spectral_perf (1, 1) logical ...
         = true

      kwargs.include_figures (1, 1) logical ...
         = true

      kwargs.include_regression_acceptance (1, 1) logical ...
         = false

      kwargs.include_perf_acceptance (1, 1) logical ...
         = false

      kwargs.include_benchmarks (1, 1) logical ...
         = false

      kwargs.rebuild_regression_baseline (1, 1) logical ...
         = false

      kwargs.rebuild_perf_baseline (1, 1) logical ...
         = false

      kwargs.output_root (1, :) string ...
         = ""

      kwargs.run_name (1, :) string ...
         = ""
   end

   % Deal out arguments.
   [simyear, solver, smbmodel, smoke_sites, full_sites, spectral_variant, ...
      n_direct_runs, include_density_floor, include_spectral_perf, ...
      include_figures, include_regression_acceptance, ...
      include_perf_acceptance, include_benchmarks, ...
      rebuild_regression_baseline, rebuild_perf_baseline, output_root, ...
      run_name] = deal(kwargs.simyear, kwargs.solver, kwargs.smbmodel, ...
      reshape(kwargs.smoke_sites, [], 1), reshape(kwargs.full_sites, [], 1), ...
      kwargs.spectral_variant, kwargs.n_direct_runs, ...
      kwargs.include_density_floor, kwargs.include_spectral_perf, ...
      kwargs.include_figures, kwargs.include_regression_acceptance, ...
      kwargs.include_perf_acceptance, kwargs.include_benchmarks, ...
      kwargs.rebuild_regression_baseline, kwargs.rebuild_perf_baseline, ...
      kwargs.output_root, kwargs.run_name);

   % Install the canonical suite config once for the whole workflow.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Resolve the run stamp and the benchmark-owned output folders.
   [run_date, run_id, run_name] = ...
      icemodel.test.helpers.resolveRunStamp(run_name);
   if isblanktext(output_root)
      output_root = fullfile(icemodel.getpath('test'), 'benchmarks');
   end
   report_root = fullfile(output_root, 'reports');
   figure_root = fullfile(output_root, 'figures');
   spectral_report_dir = fullfile(report_root, 'spectral');
   spectral_figure_dir = fullfile(figure_root, 'spectral', char(run_name));

   % Create the output folders once before running any diagnostics.
   ensureDir(spectral_report_dir);
   ensureDir(spectral_figure_dir);

   % Record the resolved bootstrap configuration in the returned struct.
   results = struct();
   results.run_date = run_date;
   results.run_id = run_id;
   results.run_name = run_name;
   results.simyear = simyear;
   results.solver = solver;
   results.smbmodel = smbmodel;
   results.smoke_sites = smoke_sites;
   results.full_sites = full_sites;
   results.spectral_variant = spectral_variant;
   results.output_root = string(output_root);
   results.report_root = string(report_root);
   results.figure_root = string(figure_root);
   results.spectral_report_dir = string(spectral_report_dir);
   results.spectral_figure_dir = string(spectral_figure_dir);

   % Run the spectral density-floor audit when requested.
   if include_density_floor
      pathname = fullfile(spectral_report_dir, ...
         sprintf('spectral_density_floor_%s.mat', char(run_name)));
      results.spectral_density_floor = summarize_spectral_density_floor( ...
         simyear=simyear, solver=solver, smoke_site=smoke_sites(1), ...
         output_file=pathname);
      results.spectral_density_floor_file = string(pathname);
   end

   % Run the spectral kernel/direct-model study and save its report.
   if include_spectral_perf
      pathname = fullfile(spectral_report_dir, ...
         sprintf('spectral_perf_%s.mat', char(run_name)));
      results.spectral_perf = summarize_spectral_perf( ...
         simyear=simyear, smoke_site=smoke_sites(1), solver=solver, ...
         include_formal_perf=false, n_direct_runs=n_direct_runs, ...
         output_file=pathname);
      results.spectral_perf_file = string(pathname);
   end

   % Export the spectral profile figures when requested.
   if include_figures
      results.spectral_figures = plot_spectral_variant_profiles( ...
         simyear=simyear, output_dir=spectral_figure_dir);
      results.spectral_figures_dir = string(spectral_figure_dir);
   end

   % Run the 2-year smoke regression acceptance with the requested variant.
   if include_regression_acceptance
      results.regression_acceptance = run_regression_suite( ...
         tier="smoke", smbmodel=smbmodel, solver=solver, simyear=simyear, ...
         smoke_sites=smoke_sites, full_sites=full_sites, ...
         baseline="rolling", spectral_variant=spectral_variant, ...
         run_name=run_name);
   end

   % Run the 2-year smoke perf acceptance with the requested variant.
   if include_perf_acceptance
      results.perf_acceptance = run_perf_suite( ...
         tier="smoke", smbmodel=smbmodel, solver=solver, simyear=simyear, ...
         smoke_sites=smoke_sites, full_sites=full_sites, ...
         baseline="rolling", n_runs=1, include_benchmarks=include_benchmarks, ...
         spectral_variant=spectral_variant, run_name=run_name);
   end

   % Optionally rebuild the rolling regression baseline with the study variant.
   if rebuild_regression_baseline
      results.regression_baseline = build_regression_baseline( ...
         baseline="rolling", smbmodel=smbmodel, solver=solver, ...
         simyear=simyear, smoke_sites=smoke_sites, full_sites=full_sites, ...
         spectral_variant=spectral_variant);
   end

   % Optionally rebuild the rolling perf baseline with the study variant.
   if rebuild_perf_baseline
      results.perf_baseline = build_perf_baseline( ...
         baseline="rolling", smbmodel=smbmodel, solver=solver, ...
         simyear=simyear, smoke_sites=smoke_sites, full_sites=full_sites, ...
         include_benchmarks=include_benchmarks, ...
         spectral_variant=spectral_variant);
   end

   % Save one bootstrap summary so the run can be reproduced from one file.
   bootstrap_file = fullfile(spectral_report_dir, ...
      sprintf('spectral_study_bootstrap_%s.mat', char(run_name)));
   save(bootstrap_file, 'results');
   results.bootstrap_file = string(bootstrap_file);

   % Print the key saved locations for interactive use.
   fprintf('bootstrap_file: %s\n', bootstrap_file);
   if isfield(results, 'spectral_perf_file')
      fprintf('spectral_perf_file: %s\n', results.spectral_perf_file);
   end
   if isfield(results, 'spectral_density_floor_file')
      fprintf('spectral_density_floor_file: %s\n', ...
         results.spectral_density_floor_file);
   end
   if isfield(results, 'spectral_figures_dir')
      fprintf('spectral_figures_dir: %s\n', results.spectral_figures_dir);
   end
end

function ensureDir(pathname)
   %ENSUREDIR Create one folder when it does not already exist.

   if exist(pathname, 'dir') ~= 7
      mkdir(pathname);
   end
end
