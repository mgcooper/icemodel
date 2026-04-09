function results = run_thf_cases(options)
   %RUN_THF_CASES Run KANM/KANL THF validation cases and save results.
   %
   %  results = run_thf_cases()
   %  results = run_thf_cases(case_idx=1:4)
   %  results = run_thf_cases(save_results=false)
   %
   % Runs all eight THF validation cases (2 sites × 2 year windows × 2
   % schemes) and saves the full run data to
   %   test/interactive/data/thf_validation_results.mat
   % for later plotting by plot_thf_cases.
   %
   % solver=1 and seb_solver=2 are fixed for every case so comparisons
   % isolate the turbulent-flux scheme rather than the outer coupling mode.

   arguments
      options.case_idx (1,:) double = []
      options.save_results (1,1) logical = true
   end

   icemodel.config('casename', 'demo');

   % Build all eight cases: 2 sites × 2 year windows × 2 schemes.
   cases = icemodel.test.helpers.buildThfValidationCases( ...
      spinup_sites=["kanm", "kanl"]);
   if ~isempty(options.case_idx)
      cases = cases(options.case_idx);
   end

   rows = repmat(init_result_row(), numel(cases), 1);
   run_data = cell(numel(cases), 1);

   for n = 1:numel(cases)
      c = cases(n);
      fprintf('\n[%d/%d] %s\n', n, numel(cases), c.label);
      [ice1, met, opts, metrics, runtime_seconds] = ...
         icemodel.test.helpers.runThfValidationCase(c);
      rows(n) = build_result_row(c, opts, metrics, runtime_seconds);
      run_data{n} = struct('ice1', ice1, 'met', met, 'c', c);
   end

   % Save to data/ (git-ignored).
   datadir = fullfile(fileparts(mfilename('fullpath')), 'data');
   if ~exist(datadir, 'dir')
      mkdir(datadir);
   end
   datafile = fullfile(datadir, 'thf_validation_results.mat');
   timestamp = datetime('now');

   if options.save_results
      save(datafile, 'cases', 'run_data', 'rows', 'timestamp');
      fprintf('\nSaved: %s\n', datafile);
   end

   % Summarize results.
   results = struct2table(rows, 'AsArray', true);
   disp(results(:, {'label', 'site', 'scheme', 'output_years', ...
      'gof_shf_bias', 'gof_shf_rmse', 'gof_shf_nse', ...
      'gof_lhf_bias', 'gof_lhf_rmse', 'gof_lhf_nse', ...
      'gof_netr_bias', 'gof_netr_rmse', 'gof_netr_nse', ...
      'stability_n_Tice_not_converged', 'stability_n_Tsfc_not_converged', ...
      'closure_seb_rmse', 'runtime_seconds'}));
end

% =========================================================================

function row = init_result_row()
   %INIT_RESULT_ROW Initialize one result row for the summary table.

   row = struct( ...
      'label', "", ...
      'site', "", ...
      'scheme', "", ...
      'simyears', "", ...
      'output_years', "", ...
      'n_spinup_years', 0, ...
      'gof_shf_bias', nan, ...
      'gof_shf_rmse', nan, ...
      'gof_shf_nse', nan, ...
      'gof_lhf_bias', nan, ...
      'gof_lhf_rmse', nan, ...
      'gof_lhf_nse', nan, ...
      'gof_netr_bias', nan, ...
      'gof_netr_rmse', nan, ...
      'gof_netr_nse', nan, ...
      'closure_seb_rmse', nan, ...
      'stability_n_Tice_not_converged', nan, ...
      'stability_n_Tsfc_not_converged', nan, ...
      'runtime_seconds', nan);
end

function row = build_result_row(c, opts, metrics, runtime_seconds)
   %BUILD_RESULT_ROW Flatten case metadata and metrics into one summary row.

   row = init_result_row();
   row.label = c.label;
   row.site = c.site;
   row.scheme = c.scheme;
   row.simyears = join(string(c.simyears), "-");
   row.output_years = join(string(opts.output_years), "-");
   row.n_spinup_years = c.n_spinup_years;
   row.gof_shf_bias = metrics.gof_shf_bias;
   row.gof_shf_rmse = metrics.gof_shf_rmse;
   row.gof_shf_nse = metrics.gof_shf_nse;
   row.gof_lhf_bias = metrics.gof_lhf_bias;
   row.gof_lhf_rmse = metrics.gof_lhf_rmse;
   row.gof_lhf_nse = metrics.gof_lhf_nse;
   row.gof_netr_bias = metrics.gof_netr_bias;
   row.gof_netr_rmse = metrics.gof_netr_rmse;
   row.gof_netr_nse = metrics.gof_netr_nse;
   row.closure_seb_rmse = metrics.closure_seb_rmse;
   row.stability_n_Tice_not_converged = metrics.stability_n_Tice_not_converged;
   row.stability_n_Tsfc_not_converged = metrics.stability_n_Tsfc_not_converged;
   row.runtime_seconds = runtime_seconds;
end
