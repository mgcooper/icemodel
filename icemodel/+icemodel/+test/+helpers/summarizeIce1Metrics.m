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
   %INITMETRICS Initialize the scalar regression metric struct.

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
   %EXTRACTFINAL Extract the retained final value of one named series.

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
   %SUMMARIZEITERCOUNTS Summarize solver iteration counts from ICE1.

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
   %COUNTFAILURES Count failed convergence flags in ICE1.

   nfail = nan;
   if any(vn == name)
      conv = logical(getField(ice1, name));
      nfail = sum(~conv);
   end
end

function [mae, rmse, maxabs] = summarizeResidual(x)
   %SUMMARIZERESIDUAL Summarize one residual series with scalar norms.

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
   %ADDENERGYBALANCEGOF Add energy-balance goodness-of-fit diagnostics.

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
   %ADDEVALUATIONWINDOWMETRICS Add evaluation-window runoff diagnostics.

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

   S.runoff_eval = computeCumulativeWindowDelta(ice1, "runoff", t1, t2);
   S.melt_eval = computeCumulativeWindowDelta(ice1, "melt", t1, t2);

   if isfield(refrow, 'area_med_m2') && isfinite(refrow.area_med_m2)
      S.icemodel_eval_m3 = computeCatchmentRunoffFinal( ...
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

function x = computeCumulativeWindowDelta(ice1, varname, t1, t2, timevec)
   %COMPUTECUMULATIVEWINDOWDELTA Return cumulative-variable change over a window.

   if nargin < 5 || isempty(timevec)
      if istimetable(ice1)
         timevec = ice1.Time;
      elseif isstruct(ice1) && isfield(ice1, 'Time')
         timevec = ice1.Time;
      else
         error('time vector is required when ice1 does not contain Time')
      end
   end

   x = nan;
   series = getSeries(ice1, varname);
   series = series(:);
   timevec = timevec(:);

   if numel(series) ~= numel(timevec)
      error('%s and time vector length mismatch', varname)
   end
   if ~isdatetime(timevec)
      error('timevec must be datetime')
   end

   inc = [0; diff(series, 1, 1)];
   keep = isbetween(timevec, t1, t2);
   if ~any(keep)
      return
   end

   inc_window = inc(keep);
   inc_window(~isfinite(inc_window)) = 0;
   if ~isempty(inc_window)
      x = sum(inc_window, 1);
   end
end

function final_m3 = computeCatchmentRunoffFinal(ice1, area_med_m2, ...
      t1, t2, timevec)
   %COMPUTECATCHMENTRUNOFFFINAL Convert point runoff to catchment cumulative.

   if nargin < 5 || isempty(timevec)
      if istimetable(ice1)
         timevec = ice1.Time;
      elseif isstruct(ice1) && isfield(ice1, 'Time')
         timevec = ice1.Time;
      else
         error('time vector is required when ice1 does not contain Time')
      end
   end

   runoff = getSeries(ice1, "runoff");
   runoff = runoff(:);
   timevec = timevec(:);

   if numel(runoff) ~= numel(timevec)
      error('runoff and time vector length mismatch')
   end
   if ~isdatetime(timevec)
      error('timevec must be datetime')
   end

   r_inc = [0; diff(runoff, 1, 1)];
   keep = isbetween(timevec, t1, t2);
   if ~any(keep)
      final_m3 = nan;
      return
   end

   r_inc_window = r_inc(keep);
   r_inc_window(~isfinite(r_inc_window)) = 0;
   catchment_cum_m3 = area_med_m2 * cumsum(r_inc_window, 1);
   if isempty(catchment_cum_m3)
      final_m3 = nan;
   else
      final_m3 = catchment_cum_m3(end);
   end
end

function x = getSeries(ice1, varname)
   %GETSERIES Extract one named cumulative series from ICE1.

   name = char(varname);
   if istimetable(ice1)
      if ~ismember(name, ice1.Properties.VariableNames)
         error('ice1 timetable missing %s variable', varname)
      end
      x = ice1.(name);
   elseif isstruct(ice1) && isfield(ice1, name)
      x = ice1.(name);
   else
      error('ice1 must provide %s', varname)
   end

   if ~isnumeric(x)
      error('%s must be numeric', varname)
   end
   if size(x, 2) > 1
      x = x(:, 1);
   end
end

function [bias, rmse, nse] = computeSeriesGOF(ice1, met, varname)
   %COMPUTESERIESGOF Compute bias, RMSE, and NSE for one output series.

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
   %COMPUTEDIFFANDPCT Compute absolute and percent deltas against a reference.

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
   %GETREFSCALAR Read one scalar field from a runoff-reference row.

   if isfield(refrow, name)
      x = refrow.(name);
   else
      x = nan;
   end
end

function x = extractPreferredSeries(ice1, names)
   %EXTRACTPREFERREDSERIES Return the first available named series.

   x = [];
   for i = 1:numel(names)
      x = extractSeries(ice1, names(i));
      if ~isempty(x)
         return
      end
   end
end

function x = extractSeries(ice1, name)
   %EXTRACTSERIES Extract one series from an ICE1 timetable or struct.

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
   %GETTIMEVECTOR Return the time vector associated with ICE1.

   if istimetable(ice1)
      timevec = ice1.Time;
   elseif isstruct(ice1) && isfield(ice1, 'Time')
      timevec = ice1.Time;
   else
      error('ice1 must provide Time to compute time-series GOF')
   end
end

function x = getField(s, name)
   %GETFIELD Read one named field/variable from a struct or timetable.

   x = s.(char(name));
end
