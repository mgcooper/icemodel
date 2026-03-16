function S = summarizeIce1Metrics(ice1, met, refrow)
%SUMMARIZEICE1METRICS Extract formal regression metrics from output and refs.
%
%  S = icemodel.test.helpers.summarizeIce1Metrics(ice1)
%  S = icemodel.test.helpers.summarizeIce1Metrics(ice1, met)
%  S = icemodel.test.helpers.summarizeIce1Metrics(ice1, met, refrow)
%
% Inputs:
%  ice1  - postprocessed model output timetable or struct
%  met   - processed met timetable used as the external energy-balance
%          reference (optional)
%  refrow - one-row runoff reference table/struct for the evaluation window
%           and catchment-scale comparisons (optional)
%
% Output fields include:
%  runoff_final, melt_final, runoff_eval, melt_eval,
%  mean_Tice_numiter, max_Tice_numiter, n_not_converged,
%  stability_n_Tsfc_not_converged,
%  closure_seb_mae, closure_seb_rmse, closure_seb_max_abs,
%  gof_* energy-balance metrics, and evaluation-period runoff GOF metrics.

   if nargin < 2
      met = [];
   end
   if nargin < 3
      refrow = [];
   end

   S = initMetrics();

   if istimetable(ice1)
      vn = string(ice1.Properties.VariableNames);
   elseif isstruct(ice1)
      vn = string(fieldnames(ice1));
   else
      error('ice1 must be timetable or struct')
   end

   S.runoff_final = extractFinal(ice1, vn, "runoff");
   S.melt_final = extractFinal(ice1, vn, "melt");

   [S.mean_Tice_numiter, S.max_Tice_numiter] = summarizeIterCounts(ice1, vn);
   S.n_not_converged = countFailures(ice1, vn, "Tice_converged");
   S.stability_mean_Tice_numiter = S.mean_Tice_numiter;
   S.stability_max_Tice_numiter = S.max_Tice_numiter;
   S.stability_n_Tice_not_converged = S.n_not_converged;
   S.stability_n_Tsfc_not_converged = countFailures(ice1, vn, "Tsfc_converged");

   seb_resid = extractPreferredSeries(ice1, ["Qbal", "balance"]);
   if ~isempty(seb_resid)
      [S.closure_seb_mae, S.closure_seb_rmse, S.closure_seb_max_abs] = ...
         summarizeResidual(seb_resid);
   end

   if ~isempty(met)
      S = addEnergyBalanceGOF(S, ice1, met);
   end

   if ~isempty(refrow)
      S = addEvaluationWindowMetrics(S, ice1, refrow);
   end
end

function S = initMetrics()

   S = struct( ...
      'runoff_final', nan, ...
      'melt_final', nan, ...
      'runoff_eval', nan, ...
      'melt_eval', nan, ...
      'mean_Tice_numiter', nan, ...
      'max_Tice_numiter', nan, ...
      'n_not_converged', nan, ...
      'stability_mean_Tice_numiter', nan, ...
      'stability_max_Tice_numiter', nan, ...
      'stability_n_Tice_not_converged', nan, ...
      'stability_n_Tsfc_not_converged', nan, ...
      'closure_seb_mae', nan, ...
      'closure_seb_rmse', nan, ...
      'closure_seb_max_abs', nan, ...
      'gof_tsfc_bias', nan, ...
      'gof_tsfc_rmse', nan, ...
      'gof_tsfc_nse', nan, ...
      'gof_netr_bias', nan, ...
      'gof_netr_rmse', nan, ...
      'gof_netr_nse', nan, ...
      'gof_lhf_bias', nan, ...
      'gof_lhf_rmse', nan, ...
      'gof_lhf_nse', nan, ...
      'gof_shf_bias', nan, ...
      'gof_shf_rmse', nan, ...
      'gof_shf_nse', nan, ...
      'obs_eval_m3', nan, ...
      'mar_eval_m3', nan, ...
      'merra_eval_m3', nan, ...
      'racmo_eval_m3', nan, ...
      'icemodel_eval_m3', nan, ...
      'gof_obs_eval_diff_m3', nan, ...
      'gof_mar_eval_diff_m3', nan, ...
      'gof_merra_eval_diff_m3', nan, ...
      'gof_racmo_eval_diff_m3', nan, ...
      'gof_obs_eval_pct_diff', nan, ...
      'gof_mar_eval_pct_diff', nan, ...
      'gof_merra_eval_pct_diff', nan, ...
      'gof_racmo_eval_pct_diff', nan);
end

function x = extractFinal(ice1, vn, name)

   x = nan;
   if any(vn == name)
      series = getField(ice1, name);
      if ~isempty(series)
         series = series(:, 1);
         x = series(end);
      end
   end
end

function [mean_iter, max_iter] = summarizeIterCounts(ice1, vn)

   mean_iter = nan;
   max_iter = nan;
   if any(vn == "Tice_numiter")
      niter = getField(ice1, "Tice_numiter");
      niter = niter(isfinite(niter));
      if ~isempty(niter)
         mean_iter = mean(niter);
         max_iter = max(niter);
      end
   end
end

function nfail = countFailures(ice1, vn, name)

   nfail = nan;
   if any(vn == name)
      conv = logical(getField(ice1, name));
      nfail = sum(~conv);
   end
end

function [mae, rmse, maxabs] = summarizeResidual(x)

   x = x(isfinite(x));
   if isempty(x)
      mae = nan;
      rmse = nan;
      maxabs = nan;
      return
   end

   mae = mean(abs(x));
   rmse = sqrt(mean(x .^ 2));
   maxabs = max(abs(x));
end

function S = addEnergyBalanceGOF(S, ice1, met)

   metrics = {
      "tsfc", "gof_tsfc";
      "netr", "gof_netr";
      "lhf",  "gof_lhf";
      "shf",  "gof_shf"
      };

   for i = 1:size(metrics, 1)
      varname = metrics{i, 1};
      prefix = metrics{i, 2};
      [bias, rmse, nse] = computeSeriesGOF(ice1, met, varname);
      S.(char(prefix + "_bias")) = bias;
      S.(char(prefix + "_rmse")) = rmse;
      S.(char(prefix + "_nse")) = nse;
   end
end

function S = addEvaluationWindowMetrics(S, ice1, refrow)

   if istable(refrow)
      if height(refrow) == 0
         return
      end
      refrow = table2struct(refrow(1, :));
   end

   if isempty(refrow) || ~isfield(refrow, 't1') || ~isfield(refrow, 't2')
      return
   end

   t1 = refrow.t1;
   t2 = refrow.t2;

   S.runoff_eval = icemodel.test.helpers.computeCumulativeWindowDelta( ...
      ice1, "runoff", t1, t2);
   S.melt_eval = icemodel.test.helpers.computeCumulativeWindowDelta( ...
      ice1, "melt", t1, t2);

   if isfield(refrow, 'area_med_m2') && isfinite(refrow.area_med_m2)
      S.icemodel_eval_m3 = icemodel.test.helpers.computeCatchmentRunoffFinal( ...
         ice1, refrow.area_med_m2, t1, t2);
   end

   S.obs_eval_m3 = getRefScalar(refrow, 'obs_final_m3');
   S.mar_eval_m3 = getRefScalar(refrow, 'mar_final_m3');
   S.merra_eval_m3 = getRefScalar(refrow, 'merra_final_m3');
   S.racmo_eval_m3 = getRefScalar(refrow, 'racmo_final_m3');

   [S.gof_obs_eval_diff_m3, S.gof_obs_eval_pct_diff] = ...
      computeDiffAndPct(S.icemodel_eval_m3, S.obs_eval_m3);
   [S.gof_mar_eval_diff_m3, S.gof_mar_eval_pct_diff] = ...
      computeDiffAndPct(S.icemodel_eval_m3, S.mar_eval_m3);
   [S.gof_merra_eval_diff_m3, S.gof_merra_eval_pct_diff] = ...
      computeDiffAndPct(S.icemodel_eval_m3, S.merra_eval_m3);
   [S.gof_racmo_eval_diff_m3, S.gof_racmo_eval_pct_diff] = ...
      computeDiffAndPct(S.icemodel_eval_m3, S.racmo_eval_m3);
end

function [bias, rmse, nse] = computeSeriesGOF(ice1, met, varname)

   bias = nan;
   rmse = nan;
   nse = nan;

   model = extractSeries(ice1, varname);
   if isempty(model) || ~ismember(char(varname), met.Properties.VariableNames)
      return
   end

   model_tt = timetable(getTimeVector(ice1), model, 'VariableNames', {'model'});
   ref = met.(char(varname));
   ref_tt = timetable(met.Time, ref, 'VariableNames', {'ref'});
   TT = synchronize(model_tt, ref_tt, 'intersection');

   mask = isfinite(TT.model) & isfinite(TT.ref);
   if nnz(mask) < 2
      return
   end

   model = TT.model(mask);
   ref = TT.ref(mask);
   resid = model - ref;

   bias = mean(resid);
   rmse = sqrt(mean(resid .^ 2));

   denom = sum((ref - mean(ref)) .^ 2);
   if denom > 0
      nse = 1 - sum(resid .^ 2) / denom;
   end
end

function [diff_value, pct_diff] = computeDiffAndPct(model_value, ref_value)

   diff_value = nan;
   pct_diff = nan;

   if isfinite(model_value) && isfinite(ref_value)
      diff_value = model_value - ref_value;
      if ref_value ~= 0
         pct_diff = 100 * diff_value / ref_value;
      end
   end
end

function x = getRefScalar(refrow, name)

   if isfield(refrow, name)
      x = refrow.(name);
   else
      x = nan;
   end
end

function x = extractPreferredSeries(ice1, names)

   x = [];
   for i = 1:numel(names)
      x = extractSeries(ice1, names(i));
      if ~isempty(x)
         return
      end
   end
end

function x = extractSeries(ice1, name)

   x = [];
   name = char(name);
   if istimetable(ice1)
      if ~ismember(name, ice1.Properties.VariableNames)
         return
      end
      x = ice1.(name);
   elseif isstruct(ice1) && isfield(ice1, name)
      x = ice1.(name);
   end

   if isempty(x) || ~isnumeric(x)
      x = [];
      return
   end
   if size(x, 2) > 1
      x = x(:, 1);
   end
   x = x(:);
end

function timevec = getTimeVector(ice1)

   if istimetable(ice1)
      timevec = ice1.Time;
   elseif isstruct(ice1) && isfield(ice1, 'Time')
      timevec = ice1.Time;
   else
      error('ice1 must provide Time to compute time-series GOF')
   end
end

function x = getField(s, name)

   x = s.(char(name));
end
