function tests = test_coupler_acceleration
   %TEST_COUPLER_ACCELERATION Cover the shared coupler iterate accelerator.
   %
   % All three couplers run the same Picard loop through these helpers, so the
   % history bookkeeping and the Aitken-then-secant selection are tested here
   % rather than through a full model run.
   tests = functiontests(localfunctions);
end

function test_initialized_history_starts_empty(testCase)
   % NaN means no history, so the first iteration falls back
   % to plain relaxation.
   returned = icemodel.couplers.initialize_coupler_history();

   testCase.verifyEqual(sort(fieldnames(returned)), ...
      sort({'T_sfc_1'; 'T_sfc_2'; 'T_sfc_prev'; 'res_prev'}))
   testCase.verifyTrue(all(structfun(@isnan, returned)))
end

function test_first_call_relaxes_toward_the_new_iterate(testCase)
   % With no history, neither Aitken nor the secant can act, so the result is
   % the relaxed Picard step T_sfc_old + alpha * residual.
   hist = icemodel.couplers.initialize_coupler_history();
   T_sfc_old = 270.0;
   T_sfc_new = 272.0;
   alpha = 0.5;

   [returned, hist] = icemodel.couplers.accelerate_coupler_iterate( ...
      hist, T_sfc_old, T_sfc_new, alpha, 10.0, true);

   expected = T_sfc_old + alpha * (T_sfc_new - T_sfc_old);
   testCase.verifyEqual(returned, expected, AbsTol=1e-12)

   % The history must now carry this iteration for the next one.
   testCase.verifyEqual(hist.T_sfc_1, T_sfc_new)
   testCase.verifyEqual(hist.T_sfc_prev, T_sfc_old)
   testCase.verifyEqual(hist.res_prev, T_sfc_new - T_sfc_old, AbsTol=1e-12)
end

function test_history_advances_by_one_iteration_each_call(testCase)
   % Aitken needs the last two iterates and the secant needs the last
   % (iterate, residual) pair, so both histories must shift every call.
   hist = icemodel.couplers.initialize_coupler_history();
   [~, hist] = icemodel.couplers.accelerate_coupler_iterate( ...
      hist, 270.0, 271.0, 0.5, 10.0, true);
   first_iterate = hist.T_sfc_1;
   [~, hist] = icemodel.couplers.accelerate_coupler_iterate( ...
      hist, 271.0, 271.5, 0.5, 10.0, true);

   testCase.verifyEqual(hist.T_sfc_2, first_iterate)
   testCase.verifyEqual(hist.T_sfc_1, 271.5)
   testCase.verifyEqual(hist.T_sfc_prev, 271.0)
end

function test_bracketing_residuals_take_a_secant_step(testCase)
   % Once two residuals bracket a root, the secant replaces the relaxed step.
   % Seed the history directly so the bracket exists on this call.
   hist = icemodel.couplers.initialize_coupler_history();
   hist.T_sfc_prev = 270.0;
   hist.res_prev = 2.0;

   T_sfc_old = 274.0;
   T_sfc_new = 272.0;
   returned = icemodel.couplers.accelerate_coupler_iterate( ...
      hist, T_sfc_old, T_sfc_new, 0.5, 100.0, true);

   % The secant root of (270, +2) and (274, -2) is 272, which is not the
   % relaxed step, so this shows that the secant stage ran.
   relaxed = T_sfc_old + 0.5 * (T_sfc_new - T_sfc_old);
   testCase.verifyEqual(returned, 272.0, AbsTol=1e-9)
   testCase.verifyNotEqual(round(returned, 6), round(relaxed, 6))
end

function test_aitken_extrapolates_when_no_bracket_exists(testCase)
   % The only assertion that the Aitken stage changes the iterate at all.
   % icemodel.numerics.aitkenscalar takes (T_sfc_2, T_sfc_1, T_sfc_new,
   % fallback); a swapped argument order there degrades all three couplers to
   % plain relaxation, and no other test in this file would fail.
   % Seed three iterates with same-sign residuals so no bracket exists and
   % the secant stage returns the Aitken extrapolation as its fallback.
   hist = icemodel.couplers.initialize_coupler_history();
   hist.T_sfc_2 = 271.0;
   hist.T_sfc_1 = 272.0;

   % Same sign as the residual this call produces, so secantscalar finds no
   % bracket and passes the Aitken value straight through.
   hist.T_sfc_prev = 270.0;
   hist.res_prev = -1.0;

   T_sfc_old = 272.5;
   T_sfc_new = 272.25;
   cpl_alpha = 0.5;
   returned = icemodel.couplers.accelerate_coupler_iterate( ...
      hist, T_sfc_old, T_sfc_new, cpl_alpha, 100.0, true);

   % Aitken on (T_sfc_2, T_sfc_1, T_sfc_new) = (271, 272, 272.25): the
   % denominator is (272.25 - 272) - (272 - 271) = -0.75 and the extrapolation
   % is 271 - (272 - 271)^2 / -0.75.
   expected = icemodel.numerics.aitkenscalar(hist.T_sfc_2, hist.T_sfc_1, ...
      T_sfc_new, T_sfc_old + cpl_alpha * (T_sfc_new - T_sfc_old), 100.0, true);
   relaxed = T_sfc_old + cpl_alpha * (T_sfc_new - T_sfc_old);
   testCase.verifyEqual(returned, expected, AbsTol=1e-12)
   testCase.verifyNotEqual(round(returned, 9), round(relaxed, 9))
end

function test_disabling_acceleration_leaves_the_relaxed_step(testCase)
   % cpl_aitken false must disable both stages, so a caller can fall back to
   % plain relaxation.
   hist = icemodel.couplers.initialize_coupler_history();
   hist.T_sfc_prev = 270.0;
   hist.res_prev = 2.0;
   hist.T_sfc_1 = 273.0;
   hist.T_sfc_2 = 271.0;

   T_sfc_old = 274.0;
   T_sfc_new = 272.0;
   returned = icemodel.couplers.accelerate_coupler_iterate( ...
      hist, T_sfc_old, T_sfc_new, 0.5, 100.0, false);

   testCase.verifyEqual(returned, T_sfc_old + 0.5 * (T_sfc_new - T_sfc_old), ...
      AbsTol=1e-12)
end
