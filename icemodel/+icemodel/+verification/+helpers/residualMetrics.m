function metrics = residualMetrics(model, observation)
   %RESIDUALMETRICS Bias, MAE, RMSE, max error, and NSE for one paired series.
   %
   %  metrics = icemodel.verification.helpers.residualMetrics(model, obs)
   %
   % Defines how this repository scores a modeled series against an observed
   % one. The diagnostics, the regression summaries, and the ablation
   % evaluation all score through it.
   %
   % Pairs where either side is not finite are dropped. A single remaining
   % pair still gives an exact bias, MAE, RMSE, and max error; only NSE needs
   % spread, so it alone is NaN there. NSE is also NaN when the observations
   % have zero variance, since the denominator is zero and the skill score is
   % undefined rather than infinite.
   %
   % NSE here is the standard form, with the OBSERVED mean in the denominator.
   % matfunclib's nashsutcliffe defaults to a modified form using the model
   % mean, so the two are not interchangeable.
   %
   % Inputs
   %  model       - modeled values
   %  observation - observed values, same size as model
   %
   % Outputs
   %  metrics - struct with fields bias, mae, rmse, max_abs_error, nse, and
   %            n_pairs

   model = model(:);
   observation = observation(:);

   metrics = struct('bias', NaN, 'mae', NaN, 'rmse', NaN, ...
      'max_abs_error', NaN, 'nse', NaN, 'n_pairs', 0);

   paired = isfinite(model) & isfinite(observation);
   metrics.n_pairs = nnz(paired);
   if metrics.n_pairs < 1
      return
   end

   model = model(paired);
   observation = observation(paired);
   residual = model - observation;

   metrics.bias = mean(residual);
   metrics.mae = mean(abs(residual));
   metrics.rmse = sqrt(mean(residual .^ 2));
   metrics.max_abs_error = max(abs(residual));

   % One pair, or zero observed variance, leaves NSE undefined not infinite.
   denominator = sum((observation - mean(observation)) .^ 2);
   if metrics.n_pairs > 1 && denominator > 0
      metrics.nse = 1 - sum(residual .^ 2) / denominator;
   end
end
