function tests = test_performance_gate
   %TEST_PERFORMANCE_GATE Test two-sided runtime review gates.

   tests = functiontests(localfunctions);
end

function test_runtime_inside_or_on_band_passes(testCase)
   % In-band measurements, including boundaries, are accepted.

   [passed, lower, upper, reason] = ...
      icemodel.test.helpers.performanceGate(10, 10, 0.2);
   verifyTrue(testCase, passed)
   verifyEqual(testCase, [lower upper], [8 12], AbsTol=eps)
   verifyEqual(testCase, reason, "")
   verifyTrue(testCase, ...
      icemodel.test.helpers.performanceGate(lower, 10, 0.2))
   verifyTrue(testCase, ...
      icemodel.test.helpers.performanceGate(upper, 10, 0.2))
end

function test_suspicious_speedup_fails_with_reason(testCase)
   % A runtime below the accepted band can expose an inflated baseline.

   [passed, ~, ~, reason] = ...
      icemodel.test.helpers.performanceGate(7, 10, 0.2);
   verifyFalse(testCase, passed)
   verifySubstring(testCase, reason, "below the accepted")
end

function test_slowdown_fails_with_reason(testCase)
   % A runtime above the accepted band remains a normal regression failure.

   [passed, ~, ~, reason] = ...
      icemodel.test.helpers.performanceGate(13, 10, 0.2);
   verifyFalse(testCase, passed)
   verifySubstring(testCase, reason, "above the accepted")
end
