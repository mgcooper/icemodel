function tests = test_core_thermo_and_timestep
   %TEST_CORE_THERMO_AND_TIMESTEP Verify local thermo and substep contracts.
   tests = functiontests(localfunctions);
end

function test_vappress_derivative_matches_finite_difference(testCase)
   % Check the analytic vapor-pressure derivative against a centered finite
   % difference so downstream linearizations have a stable local reference.

   Tf = icemodel.physicalConstant('Tf');
   T = Tf - 8.0;
   [es, des_dT] = VAPPRESS(T, false);
   h = 1e-3;
   fd = (VAPPRESS(T + h, false) - VAPPRESS(T - h, false)) / (2 * h);

   testCase.verifyGreaterThan(es, 0);
   testCase.verifyEqual(des_dT, fd, 'RelTol', 1e-6);
end

function test_tdewpoint_inverts_vappress(testCase)
   % TDEWPOINT should invert VAPPRESS: es(Tdew) should equal ea = es(T)*rh/100.

   Tf = icemodel.physicalConstant('Tf');
   T = Tf - 5.0;
   rh = 75;

   Tdew = TDEWPOINT(T, rh, false);
   ea = VAPPRESS(T, false) * rh / 100;
   es_at_Tdew = VAPPRESS(Tdew, false);

   testCase.verifyEqual(es_at_Tdew, ea, 'RelTol', 1e-8);
   testCase.verifyLessThan(Tdew, T, ...
      'Dew point should be below air temperature for rh < 100');
end

function test_loadmetdata_respects_liqflag(testCase)
   % LOADMETDATA should switch between ice and liquid saturation relations
   % when the caller flips the liquid-phase flag.

   [met, ~] = icemodel.test.fixtures.makeSyntheticMetFile(2016, 'nsteps', 1);
   met.De = 0.01;

   [~, ~, ~, ~, ~, ~, ~, ea_ice] = LOADMETDATA(met, 1, false);
   [~, ~, ~, ~, ~, ~, ~, ea_liq] = LOADMETDATA(met, 1, true);

   testCase.verifyGreaterThan(ea_liq, ea_ice);
end

function test_metsub_interpolates_midstep_for_enabled_option(testCase)
   % With substep interpolation enabled, METSUB should evaluate the forcing
   % state halfway through the parent step.

   opts.met_substep_interp = true;
   Tf = icemodel.physicalConstant('Tf');

   [tair, swd, lwd, albedo, wspd, rh, psfc, ppt, tppt, De, ea] = ...
      METSUB(opts, 1, 1800, 3600, true, Tf, [270; 272], [100; 200], ...
      [220; 260], [0.60; 0.70], [3; 5], [70; 90], [78000; 79000], ...
      [0; 2], [268; 269], [0.01; 0.02]);

   testCase.verifyEqual(tair, 271);
   testCase.verifyEqual(swd, 150);
   testCase.verifyEqual(lwd, 240);
   testCase.verifyEqual(albedo, 0.65, 'AbsTol', 1e-12);
   testCase.verifyEqual(wspd, 4);
   testCase.verifyEqual(rh, 80);
   testCase.verifyEqual(psfc, 78500);
   testCase.verifyEqual(ppt, 1);
   testCase.verifyEqual(tppt, 268.5);
   testCase.verifyEqual(De, 0.015);
   testCase.verifyEqual(ea, VAPPRESS(271, true) * 0.80, 'RelTol', 1e-12);
end

function test_metsub_returns_step_value_when_interp_disabled(testCase)
   % With interpolation disabled, METSUB should hold the parent-step value
   % and still compute a consistent vapor pressure from that state.

   opts.met_substep_interp = false;
   Tf = icemodel.physicalConstant('Tf');

   [tair, swd, lwd, albedo, wspd, rh, psfc, ppt, tppt, De, ea] = ...
      METSUB(opts, 1, 1800, 3600, false, Tf, [270; 272], [100; 200], ...
      [220; 260], [0.60; 0.70], [3; 5], [70; 90], [78000; 79000], ...
      [0; 2], [268; 269], [0.01; 0.02]);

   testCase.verifyEqual([tair swd lwd albedo wspd rh psfc ppt tppt De], ...
      [270 100 220 0.60 3 70 78000 0 268 0.01], 'AbsTol', 1e-12);
   testCase.verifyEqual(ea, VAPPRESS(270, false) * 0.70, 'RelTol', 1e-12);
end

function test_sfcflin_matches_surface_terms_at_linearization_point(testCase)
   % The linearized surface flux coefficients should reproduce the full
   % surface residual at the temperature used for the linearization.

   [cv_air, SB, roLs, Tf, epsilon] = icemodel.physicalConstant( ...
      'cv_air', 'SB', 'roLs', 'Tf', 'epsilon');
   emiss = icemodel.parameterLookup('emiss');
   [De, scoef] = WINDCOEF(4.0, 0.001, 2.0, 3.0);

   Ta = Tf - 10.0;
   Ts = Tf - 5.0;
   Qsi = 200.0;
   Qli = 240.0;
   albedo = 0.6;
   wspd = 4.0;
   Pa = 78000.0;
   chi = 1.0;
   liqflag = false;
   ea = 0.8 * VAPPRESS(Ta, liqflag);

   [Sc, Sp] = SFCFLIN(Ta, Qsi, Qli, albedo, wspd, Pa, De, ea, cv_air, ...
      emiss, SB, roLs, scoef, chi, Tf, Ts, liqflag);

   es = VAPPRESS(Ts, liqflag);
   S = STABLEFN(Ta, Ts, wspd, scoef);
   F = emiss * (Qli - SB * Ts ^ 4) ...
      + chi * Qsi * (1 - albedo) ...
      + cv_air * De * (Ta - Ts) * S ...
      + roLs * De * epsilon / Pa * (ea - es) * S;

   testCase.verifyEqual(Sc + Sp * Ts, F, 'RelTol', 1e-10);
end
