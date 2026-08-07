function [compare_metric, passed, reason] = ...
      formalRegressionMetricEvidence(actual, baseline, row, varname, selector)
   %FORMALREGRESSIONMETRICEVIDENCE Check one saved metric is comparable.
   %
   %  [compare_metric, passed, reason] = ...
   %     icemodel.test.helpers.formalRegressionMetricEvidence( ...
   %     actual, baseline, row, varname, selector)
   %
   % Rolling baselines use the current complete metric schema and therefore
   % fail closed. Frozen release baselines may omit metrics added after their
   % release, but every metric they do persist must remain finite.

   policy = icemodel.test.helpers.formalBaselinePolicy(selector);
   names = string(baseline.Properties.VariableNames);
   persisted = ismember(string(varname), names);
   if ~persisted
      compare_metric = false;
      passed = policy.baseline_type ~= "rolling";
      reason = "baseline metric is absent: " + string(varname);
      return
   end

   % A persisted metric is evidence only when both sides are finite scalars.
   expected = baseline.(varname)(row);
   compare_metric = true;
   passed = isscalar(actual) && isscalar(expected) ...
      && isfinite(actual) && isfinite(expected);
   reason = "nonfinite regression evidence for metric " + string(varname);
end
