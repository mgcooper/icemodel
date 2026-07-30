function gate = admissionGate(channel, metrics, baseline_rmse, kwargs)
   %ADMISSIONGATE Apply the approved per-variable admission thresholds.
   %
   %  gate = icemodel.forcing.reconstruct.admissionGate("lwd", ...
   %     metrics.overall, baseline_rmse)
   %
   % Role
   %  Method-admission decision for one channel and one metric row (gap-fill
   %  policy): a method is admitted where, on held-out
   %  draws, |bias| stays within the instrument-class cap, RMSE improves on
   %  the best available baseline by at least the required margin, and no
   %  physical-bound violation occurred. The engine consumes the emitted
   %  admit/deny row per stratum; a denied stratum stays missing rather
   %  than receiving a weaker fill.
   %
   % Inputs
   %  channel : canonical channel name (bias caps are per channel).
   %  metrics : one-row metric table from validationMetrics (overall or one
   %     stratum row).
   %  baseline_rmse : held-out RMSE of the policy baseline on the exact
   %     finite candidate support (persistence for the ≤6 h bucket,
   %     station day-of-year climatology otherwise). NaN disables the test ONLY
   %     when the baseline itself produced no finite prediction, and the
   %     gate records that condition.
   %
   % Name-value
   %  rmse_improvement : required fractional improvement (default 0.10).
   %  min_variability_ratio, max_variability_ratio : accepted range for
   %     predicted/observed standard deviation when held-out truth has
   %     measurable spread. This prevents a low-RMSE mean/climatology
   %     estimate from winning by suppressing weather variability. Missing
   %     within-gap spread evidence is denied; proven zero spread is exempt.
    %  min_coverage : minimum reconstructed fraction of the drawn samples
   %     (default 0.10). This is a usefulness floor, not a completeness
   %     requirement: the orchestrator composes ordered methods and
   %     cascades uncovered leftovers to the next tier, so a method that
   %     skillfully covers part of a stratum (e.g. a donor whose record
   %     ends before the target's) is admitted for the samples it covers
    %     — support-held coarse-cadence donors reach only fractional
    %     coverage on a finer target axis while beating the climatology
    %     baseline severalfold where they do cover.
    %  metrics must contain finite provenance_accounting equal to one; an
    %     absent accounting result is a failed gate, never an opt-out.
   %
   % Returns
   %  gate : struct — admit (logical), reasons (string column, empty when
   %     admitted), bias_cap, rmse_improvement_required, baseline_rmse.
   %
   % See also: icemodel.forcing.reconstruct.validationMetrics

   arguments
      channel (1, 1) string
      metrics (1, :) table
      baseline_rmse (1, 1) double
      kwargs.rmse_improvement (1, 1) double {mustBeNonnegative} = ...
         icemodel.forcing.reconstruct.setopts().rmse_improvement
      kwargs.min_variability_ratio (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().min_variability_ratio
      kwargs.max_variability_ratio (1, 1) double {mustBePositive} = ...
         icemodel.forcing.reconstruct.setopts().max_variability_ratio
      kwargs.min_coverage (1, 1) double = ...
         icemodel.forcing.reconstruct.setopts().min_coverage
      kwargs.typical_magnitude (1, 1) double = NaN
   end

   bias_cap = biasCap(channel);
   % Degenerate metric rows (NaN bias/rmse/coverage) must be denied, not
   % silently admitted through NaN comparisons.
   if ~isfinite(metrics.bias) || ~isfinite(metrics.coverage) || ...
         (~isfinite(metrics.rmse) && isfinite(baseline_rmse))
      gate = struct('admit', false, ...
         'reasons', "nonfinite metrics row", 'bias_cap', bias_cap, ...
         'rmse_improvement_required', kwargs.rmse_improvement, ...
         'baseline_rmse', baseline_rmse, ...
         'baseline_available', isfinite(baseline_rmse));
      return
   end
   if channel == "wspd" && isfinite(kwargs.typical_magnitude)
      % Policy: wspd bias cap is 1 m/s OR 10% — the relative alternative
      % applies when the caller supplies the stratum's typical wind speed.
      bias_cap = max(bias_cap, 0.10 * kwargs.typical_magnitude);
   end
   reasons = strings(0, 1);

   % Each failed criterion is its own recorded reason so the report can say
   % exactly why a stratum stays missing.
   if abs(metrics.bias) > bias_cap
      reasons = [reasons; sprintf("bias %.3g exceeds cap %.3g", ...
         metrics.bias, bias_cap)];
   end
   if isfinite(baseline_rmse)
      required = (1 - kwargs.rmse_improvement) * baseline_rmse;
      if ~(metrics.rmse <= required)
         reasons = [reasons; sprintf( ...
            "rmse %.3g not <= %.3g (>=%.0f%% better than baseline %.3g)", ...
            metrics.rmse, required, 100 * kwargs.rmse_improvement, ...
            baseline_rmse)];
      end
   end
   % A NaN baseline means the baseline itself produced no finite prediction;
   % the improvement test is disabled and baseline_available records it.
   if metrics.bound_violations > 0
      reasons = [reasons; sprintf("%d physical-bound violation(s)", ...
         metrics.bound_violations)];
   end
   if metrics.coverage < kwargs.min_coverage
      reasons = [reasons; sprintf("coverage %.2f below %.2f", ...
         metrics.coverage, kwargs.min_coverage)];
   end
    if ~ismember('provenance_accounting', metrics.Properties.VariableNames) ...
          || ~isfinite(metrics.provenance_accounting)
       reasons = [reasons; "provenance accounting unavailable"];
    elseif metrics.provenance_accounting ~= 1
       reasons = [reasons; sprintf( ...
          "provenance accounting %.2f is not 1.00", ...
          metrics.provenance_accounting)];
   end
   if ~ismember('within_gap_observed_spread', ...
         metrics.Properties.VariableNames) ...
         || ~isfinite(metrics.within_gap_observed_spread)
      reasons = [reasons; "within-gap variability evidence unavailable"];
   elseif metrics.within_gap_observed_spread > 0
      if ~isfinite(metrics.variability_ratio)
         reasons = [reasons; "variability ratio unavailable"];
      elseif metrics.variability_ratio < kwargs.min_variability_ratio ...
            || metrics.variability_ratio > kwargs.max_variability_ratio
         reasons = [reasons; sprintf( ...
            "variability ratio %.2f outside [%.2f, %.2f]", ...
            metrics.variability_ratio, kwargs.min_variability_ratio, ...
            kwargs.max_variability_ratio)];
      end
   end

   gate = struct('admit', isempty(reasons), 'reasons', reasons, ...
      'bias_cap', bias_cap, ...
      'rmse_improvement_required', kwargs.rmse_improvement, ...
      'baseline_rmse', baseline_rmse, ...
      'baseline_available', isfinite(baseline_rmse));
end

function cap = biasCap(channel)
   %BIASCAP Instrument-class absolute bias caps (policy admission gates).
   switch channel
      case "tair"
         cap = 0.5;      % K
      case "rh"
         cap = 5;        % percent
      case "wspd"
         cap = 1;        % m/s
      case "psfc"
         cap = 100;      % Pa
      case {"swd", "swu"}
         cap = 20;       % W/m2 (daily-mean sense)
      case "lwd"
         cap = 15;       % W/m2
      case "albedo"
         cap = 0.05;     % fraction
      case "ppt"
         cap = Inf;      % precip admits on RMSE/bounds only (no native truth)
      otherwise
         error('icemodel:reconstruct:admissionGate:unknownChannel', ...
            'no approved bias cap for channel: %s', channel)
   end
end
