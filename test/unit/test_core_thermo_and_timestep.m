function tests = test_core_thermo_and_timestep
   %TEST_CORE_THERMO_AND_TIMESTEP Verify local thermo and substep contracts.
   tests = functiontests(localfunctions);
end

function test_vappress_derivative_matches_finite_difference(testCase)
   % Check the analytic vapor-pressure derivative against a centered finite
   % difference so downstream linearizations have a stable local reference.

   Tf = icemodel.physicalConstant('Tf');
   T = Tf - 8.0;
   [es, des_dT] = icemodel.vapor.saturation_vapor_pressure(T, false);
   h = 1e-3;
   fd = (icemodel.vapor.saturation_vapor_pressure(T + h, false) ...
      - icemodel.vapor.saturation_vapor_pressure(T - h, false)) / (2 * h);

   testCase.verifyGreaterThan(es, 0);
   testCase.verifyEqual(des_dT, fd, 'RelTol', 1e-6);
end

function test_tdewpoint_inverts_vappress(testCase)
   % icemodel.vapor.dew_point_temperature should invert icemodel.vapor.saturation_vapor_pressure: es(Tdew) should equal ea = es(T)*rh/100.

   Tf = icemodel.physicalConstant('Tf');
   T = Tf - 5.0;
   rh = 75;

   Tdew = icemodel.vapor.dew_point_temperature(T, rh, false);
   ea = icemodel.vapor.saturation_vapor_pressure(T, false) * rh / 100;
   es_at_Tdew = icemodel.vapor.saturation_vapor_pressure(Tdew, false);

   testCase.verifyEqual(es_at_Tdew, ea, 'RelTol', 1e-8);
   testCase.verifyLessThan(Tdew, T, ...
      'Dew point should be below air temperature for rh < 100');
end

function test_metsub_interpolates_midstep_for_enabled_option(testCase)
   % With substep interpolation enabled, getsubstepforcings should evaluate
   % the forcing state halfway through the parent step.

   opts.met_substep_interp = true;
   Tf = icemodel.physicalConstant('Tf');

   [tair, swd, lwd, albedo, wspd, rh, psfc, ppt, tppt, De, ea] = ...
      icemodel.timestepping.getsubstepforcings( ...
      opts, 1, 1800, 3600, true, Tf, [270; 272], [100; 200], ...
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
   testCase.verifyEqual(ea, icemodel.vapor.saturation_vapor_pressure(271, true) * 0.80, 'RelTol', 1e-12);
end

function test_metsub_returns_step_value_when_interp_disabled(testCase)
   % With interpolation disabled, getsubstepforcings should hold the
   % parent-step value and still compute a consistent vapor pressure.

   opts.met_substep_interp = false;
   Tf = icemodel.physicalConstant('Tf');

   [tair, swd, lwd, albedo, wspd, rh, psfc, ppt, tppt, De, ea] = ...
      icemodel.timestepping.getsubstepforcings( ...
      opts, 1, 1800, 3600, false, Tf, [270; 272], [100; 200], ...
      [220; 260], [0.60; 0.70], [3; 5], [70; 90], [78000; 79000], ...
      [0; 2], [268; 269], [0.01; 0.02]);

   testCase.verifyEqual([tair swd lwd albedo wspd rh psfc ppt tppt De], ...
      [270 100 220 0.60 3 70 78000 0 268 0.01], 'AbsTol', 1e-12);
   testCase.verifyEqual(ea, icemodel.vapor.saturation_vapor_pressure(270, false) * 0.70, 'RelTol', 1e-12);
end

function test_sfcflin_matches_surface_terms_at_linearization_point(testCase)
   % The linearized surface flux coefficients should reproduce the full
   % surface residual at the temperature used for the linearization.

   [cv_air, cv_liq, SB, roLs, Tf, epsilon] = icemodel.physicalConstant( ...
      'cv_air', 'cv_liq', 'SB', 'roLs', 'Tf', 'epsilon');
   emiss = icemodel.parameterLookup('emiss');
   wspd = 4.0;
   z0_bulk = 0.001;
   z_tair = 2.0;
   z_wind = 3.0;
   [De, br_coefs] = ...
      icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
      wspd, z0_bulk, z_tair, z_wind);

   tair = Tf - 10.0;
   T_sfc = Tf - 5.0;
   Qsi = 200.0;
   Qli = 240.0;
   albedo = 0.6;
   ppt = 2e-4;
   tppt = tair + 4.0;
   psfc = 78000.0;
   chi = 1.0;
   liqflag = false;
   ea_atm = 0.8 * icemodel.vapor.saturation_vapor_pressure(tair, liqflag);

   [Sc, Sp] = ...
      icemodel.surface.turbulence.bulk_richardson.surface_flux_linearization( ...
      T_sfc, tair, Qsi, Qli, albedo, wspd, ppt, tppt, psfc, De, ea_atm, ...
      br_coefs, roLs, liqflag, chi);

   es_sfc = icemodel.vapor.saturation_vapor_pressure(T_sfc, liqflag);
   stability = icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
      T_sfc, tair, wspd, br_coefs);
   F = emiss * (Qli - SB * T_sfc ^ 4) ...
      + chi * Qsi * (1 - albedo) ...
      + cv_air * De * (tair - T_sfc) * stability ...
      + roLs * De * epsilon / psfc * (ea_atm - es_sfc) * stability ...
      + icemodel.surface.advective_heat_flux(ppt, tppt, cv_liq);

   testCase.verifyEqual(Sc + Sp * T_sfc, F, 'RelTol', 1e-10);
end
