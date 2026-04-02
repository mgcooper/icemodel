function results = validate_thf_cases(make_plots, case_idx)
%VALIDATE_THF_CASES Run focused KANM/KANL THF validation cases.
%
%  results = validate_thf_cases()
%  results = validate_thf_cases(true)
%  results = validate_thf_cases(false, 1:2)
%
% This interactive helper runs the current turbulent-flux validation cases
% discussed in Bead `icemodel-3sm`:
%   - KANM 2016, bulk_richardson
%   - KANM 2016, bulk_mo
%   - KANM 2015-2016 with 2015 spinup, bulk_richardson
%   - KANM 2015-2016 with 2015 spinup, bulk_mo
%   - KANL 2016, bulk_richardson
%   - KANL 2016, bulk_mo
%
% The helper forces `solver=1`, `seb_solver=2` for every case so the
% comparison isolates the turbulent-flux scheme rather than the outer
% surface/subsurface coupling mode. It returns a compact table of SHF/LHF
% and net-radiation GOF metrics based on the processed met contract.

if nargin < 1
   make_plots = false;
end
if nargin < 2 || isempty(case_idx)
   case_idx = [];
end

cases = validation_cases();
if ~isempty(case_idx)
   cases = cases(case_idx);
end
rows = repmat(init_result_row(), numel(cases), 1);

% Use the repo-local demo data tree, which contains the KANM/KANL point
% forcing files needed by this validation pass.
icemodel.config('casename', 'demo');

for i = 1:numel(cases)
   c = cases(i);
   fprintf('\n[%d/%d] %s\n', i, numel(cases), c.label);

   tic
   [ice1, met, opts] = run_validation_case(c);
   runtime_seconds = toc;

   metrics = icemodel.test.helpers.summarizeIce1Metrics(ice1, met);
   rows(i) = build_result_row(c, opts, metrics, runtime_seconds);

   if make_plots
      plot_validation_case(ice1, met, c.label);
   end
end

results = struct2table(rows, 'AsArray', true);
disp(results(:, {'label', 'site', 'scheme', 'output_years', ...
   'gof_shf_bias', 'gof_shf_rmse', 'gof_shf_nse', ...
   'gof_lhf_bias', 'gof_lhf_rmse', 'gof_lhf_nse', ...
   'gof_netr_bias', 'gof_netr_rmse', 'gof_netr_nse', ...
   'stability_n_Tice_not_converged', 'stability_n_Tsfc_not_converged', ...
   'closure_seb_rmse', 'runtime_seconds'}));
end

function cases = validation_cases()
   %VALIDATION_CASES Return the current focused THF validation cases.

   cases = [ ...
      struct('label', "KANM 2016 bulk_richardson", ...
         'site', "kanm", 'forcings', "kanm", 'simyears', 2016, ...
         'n_spinup_years', 0, 'scheme', "bulk_richardson"), ...
      struct('label', "KANM 2016 bulk_mo", ...
         'site', "kanm", 'forcings', "kanm", 'simyears', 2016, ...
         'n_spinup_years', 0, 'scheme', "bulk_mo"), ...
      struct('label', "KANM 2015-2016 spinup bulk_richardson", ...
         'site', "kanm", 'forcings', "kanm", 'simyears', [2015 2016], ...
         'n_spinup_years', 1, 'scheme', "bulk_richardson"), ...
      struct('label', "KANM 2015-2016 spinup bulk_mo", ...
         'site', "kanm", 'forcings', "kanm", 'simyears', [2015 2016], ...
         'n_spinup_years', 1, 'scheme', "bulk_mo"), ...
      struct('label', "KANL 2016 bulk_richardson", ...
         'site', "kanl", 'forcings', "kanl", 'simyears', 2016, ...
         'n_spinup_years', 0, 'scheme', "bulk_richardson"), ...
      struct('label', "KANL 2016 bulk_mo", ...
         'site', "kanl", 'forcings', "kanl", 'simyears', 2016, ...
         'n_spinup_years', 0, 'scheme', "bulk_mo")];
end

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

function [ice1, met, opts] = run_validation_case(c)
   %RUN_VALIDATION_CASE Configure and run one THF validation case.

   opts = icemodel.setopts('icemodel', c.site, c.simyears, c.forcings, ...
      [], [], [], false, false, ...
      'n_spinup_years', c.n_spinup_years, ...
      'solver', 1, ...
      'seb_solver', 2, ...
      'turbulent_flux_scheme', c.scheme);

   [ice1, ice2, opts] = icemodel(opts);
   [ice1, ~, met] = icemodel.postprocess(ice1, ice2, opts, opts.output_years);
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

function plot_validation_case(ice1, met, label)
   %PLOT_VALIDATION_CASE Show the standard overlay and 1:1 scatter panels.

   icemodel.plot.enbal(ice1, met);
   title(label);

   figure('Name', label + " scatter");
   tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
   plot_one_to_one(ice1, met, "netr");
   plot_one_to_one(ice1, met, "shf");
   plot_one_to_one(ice1, met, "lhf");
end

function plot_one_to_one(ice1, met, varname)
   %PLOT_ONE_TO_ONE Compare one processed met/model series against 1:1.

   model = timetable(ice1.Time, ice1.(char(varname)), 'VariableNames', {'model'});
   ref = timetable(met.Time, met.(char(varname)), 'VariableNames', {'ref'});
   TT = synchronize(model, ref, 'intersection');
   keep = isfinite(TT.model) & isfinite(TT.ref);

   nexttile
   scatter(TT.model(keep), TT.ref(keep), 8, 'filled');
   hold on
   lo = min([TT.model(keep); TT.ref(keep)]);
   hi = max([TT.model(keep); TT.ref(keep)]);
   plot([lo hi], [lo hi], '-', 'LineWidth', 1.2);
   axis equal
   xlim([lo hi]);
   ylim([lo hi]);
   xlabel('icemodel');
   ylabel('forcing');
   title(upper(char(varname)));
end
