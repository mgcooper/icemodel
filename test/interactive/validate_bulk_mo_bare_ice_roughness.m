function results = validate_bulk_mo_bare_ice_roughness(z0_ice_values, make_plots, case_idx)
%VALIDATE_BULK_MO_BARE_ICE_ROUGHNESS Sweep z0_ice on real KANM/KANL cases.
%
%  results = validate_bulk_mo_bare_ice_roughness()
%  results = validate_bulk_mo_bare_ice_roughness([0.001 0.003 0.01])
%  results = validate_bulk_mo_bare_ice_roughness([], true)
%  results = validate_bulk_mo_bare_ice_roughness([], false, 1:2)
%
% This helper keeps the full-year THF validation contract and focuses only
% on the `monin_obukhov` bare-ice momentum roughness question. It uses the shared
% real-case case-builder plus explicit `z0_ice` overrides so site-specific
% calibration can be compared without baking any hidden site policy into the
% production defaults. The current cases are:
%   - KANM 2016
%   - KANM 2015-2016 with one spinup year retained only in the runtime
%   - KANL 2016
%   - KANL 2015-2016 with one spinup year retained only in the runtime

if nargin < 1 || isempty(z0_ice_values)
   z0_ice_values = [0.001 0.003 0.01];
end
if nargin < 2
   make_plots = false;
end
if nargin < 3 || isempty(case_idx)
   case_idx = [];
end

cases = icemodel.test.helpers.buildThfValidationCases( ...
   schemes="monin_obukhov", ...
   spinup_sites=["kanm" "kanl"]);
if ~isempty(case_idx)
   cases = cases(case_idx);
end
rows = repmat(init_result_row(), numel(cases) * numel(z0_ice_values), 1);

% Use the repo-local demo data tree, which contains the KANM/KANL point
% forcing files needed by this validation pass.
icemodel.config('casename', 'demo');

k = 0;
for i = 1:numel(cases)
   c = cases(i);
   for j = 1:numel(z0_ice_values)
      z0_ice = z0_ice_values(j);
      label = sprintf('%s z0_ice=%.4f', c.label, z0_ice);
      fprintf('\n[%d/%d] %s\n', k + 1, numel(rows), label);

      overrides = struct( ...
         'z0_ice', z0_ice, ...
         'testname', sprintf('monin_obukhov_rough_%s_%s', c.site, ...
         strrep(num2str(z0_ice, '%.4f'), '.', 'p')));
      [ice1, met, opts, metrics, runtime_seconds] = ...
         icemodel.test.helpers.runThfValidationCase(c, overrides=overrides);

      k = k + 1;
      rows(k) = build_result_row(c, opts, metrics, runtime_seconds, z0_ice);

      if make_plots
         icemodel.plot.enbal(ice1, met);
         title(label);
      end
   end
end

results = struct2table(rows, 'AsArray', true);
disp(results(:, {'label', 'site', 'output_years', 'z0_ice', ...
   'gof_shf_bias', 'gof_shf_rmse', 'gof_shf_nse', ...
   'gof_lhf_bias', 'gof_lhf_rmse', 'gof_lhf_nse', ...
   'gof_netr_bias', 'gof_netr_rmse', 'gof_netr_nse', ...
   'stability_n_Tice_not_converged', 'stability_n_Tsfc_not_converged', ...
   'closure_seb_rmse', 'runtime_seconds'}));
end

function row = init_result_row()
   %INIT_RESULT_ROW Initialize one result row for the roughness summary.

   row = struct( ...
      'label', "", ...
      'site', "", ...
      'simyears', "", ...
      'output_years', "", ...
      'n_spinup_years', 0, ...
      'z0_ice', nan, ...
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

function row = build_result_row(c, opts, metrics, runtime_seconds, z0_ice)
   %BUILD_RESULT_ROW Flatten one roughness-sensitivity result.

   row = init_result_row();
   row.label = c.label;
   row.site = c.site;
   row.simyears = join(string(c.simyears), "-");
   row.output_years = join(string(opts.output_years), "-");
   row.n_spinup_years = c.n_spinup_years;
   row.z0_ice = z0_ice;
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
