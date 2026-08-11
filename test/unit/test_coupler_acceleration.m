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
      sort({'Ts_1'; 'Ts_2'; 'Ts_prev'; 'res_prev'}))
   testCase.verifyTrue(all(structfun(@isnan, returned)))
end

function test_first_call_relaxes_toward_the_new_iterate(testCase)
   % With no history, neither Aitken nor the secant can act, so the result is
   % the relaxed Picard step Ts_old + alpha * residual.
   hist = icemodel.couplers.initialize_coupler_history();
   Ts_old = 270.0;
   Ts_new = 272.0;
   alpha = 0.5;

   [returned, hist] = icemodel.couplers.accelerate_coupler_iterate( ...
      hist, Ts_old, Ts_new, alpha, 10.0, true);

   expected = Ts_old + alpha * (Ts_new - Ts_old);
   testCase.verifyEqual(returned, expected, AbsTol=1e-12)

   % The history must now carry this iteration for the next one.
   testCase.verifyEqual(hist.Ts_1, Ts_new)
   testCase.verifyEqual(hist.Ts_prev, Ts_old)
   testCase.verifyEqual(hist.res_prev, Ts_new - Ts_old, AbsTol=1e-12)
end

function test_history_advances_by_one_iteration_each_call(testCase)
   % Aitken needs the last two iterates and the secant needs the last
   % (iterate, residual) pair, so both histories must shift every call.
   hist = icemodel.couplers.initialize_coupler_history();
   [~, hist] = icemodel.couplers.accelerate_coupler_iterate( ...
      hist, 270.0, 271.0, 0.5, 10.0, true);
   first_iterate = hist.Ts_1;
   [~, hist] = icemodel.couplers.accelerate_coupler_iterate( ...
      hist, 271.0, 271.5, 0.5, 10.0, true);

   testCase.verifyEqual(hist.Ts_2, first_iterate)
   testCase.verifyEqual(hist.Ts_1, 271.5)
   testCase.verifyEqual(hist.Ts_prev, 271.0)
end

function test_bracketing_residuals_take_a_secant_step(testCase)
   % Once two residuals bracket a root, the secant replaces the relaxed step.
   % Seed the history directly so the bracket exists on this call.
   hist = icemodel.couplers.initialize_coupler_history();
   hist.Ts_prev = 270.0;
   hist.res_prev = 2.0;

   Ts_old = 274.0;
   Ts_new = 272.0;
   returned = icemodel.couplers.accelerate_coupler_iterate( ...
      hist, Ts_old, Ts_new, 0.5, 100.0, true);

   % The secant root of (270, +2) and (274, -2) is 272, which is not the
   % relaxed step, so this shows that the secant stage ran.
   relaxed = Ts_old + 0.5 * (Ts_new - Ts_old);
   testCase.verifyEqual(returned, 272.0, AbsTol=1e-9)
   testCase.verifyNotEqual(round(returned, 6), round(relaxed, 6))
end

function test_aitken_extrapolates_when_no_bracket_exists(testCase)
   % The only assertion that the Aitken stage changes the iterate at all.
   % icemodel.numerics.aitkenscalar takes (Ts_2, Ts_1, Ts_new, fallback); a
   % swapped argument order there degrades all three couplers to plain
   % relaxation, and no other test in this file would fail.
   % Seed three iterates with same-sign residuals so no bracket exists and
   % the secant stage returns the Aitken extrapolation as its fallback.
   hist = icemodel.couplers.initialize_coupler_history();
   hist.Ts_2 = 271.0;
   hist.Ts_1 = 272.0;

   % Same sign as the residual this call produces, so secantscalar finds no
   % bracket and passes the Aitken value straight through.
   hist.Ts_prev = 270.0;
   hist.res_prev = -1.0;

   Ts_old = 272.5;
   Ts_new = 272.25;
   cpl_alpha = 0.5;
   returned = icemodel.couplers.accelerate_coupler_iterate( ...
      hist, Ts_old, Ts_new, cpl_alpha, 100.0, true);

   % Aitken on (Ts_2, Ts_1, Ts_new) = (271, 272, 272.25): the denominator is
   % (272.25 - 272) - (272 - 271) = -0.75 and the extrapolation is
   % 271 - (272 - 271)^2 / -0.75.
   expected = icemodel.numerics.aitkenscalar(hist.Ts_2, hist.Ts_1, Ts_new, ...
      Ts_old + cpl_alpha * (Ts_new - Ts_old), 100.0, true);
   relaxed = Ts_old + cpl_alpha * (Ts_new - Ts_old);
   testCase.verifyEqual(returned, expected, AbsTol=1e-12)
   testCase.verifyNotEqual(round(returned, 9), round(relaxed, 9))
end

function test_disabling_acceleration_leaves_the_relaxed_step(testCase)
   % cpl_aitken false must disable both stages, so a caller can fall back to
   % plain relaxation.
   hist = icemodel.couplers.initialize_coupler_history();
   hist.Ts_prev = 270.0;
   hist.res_prev = 2.0;
   hist.Ts_1 = 273.0;
   hist.Ts_2 = 271.0;

   Ts_old = 274.0;
   Ts_new = 272.0;
   returned = icemodel.couplers.accelerate_coupler_iterate( ...
      hist, Ts_old, Ts_new, 0.5, 100.0, false);

   testCase.verifyEqual(returned, Ts_old + 0.5 * (Ts_new - Ts_old), ...
      AbsTol=1e-12)
end
