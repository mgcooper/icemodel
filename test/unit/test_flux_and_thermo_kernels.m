function tests = test_flux_and_thermo_kernels
   %TEST_FLUX_AND_THERMO_KERNELS Verify local thermo and surface-flux kernels.
   tests = functiontests(localfunctions);
end

function test_pressure_decreases_with_elevation(testCase)
   % Pressure should monotonically decrease with elevation in the helper.

   sea_level = PRESSURE(0);
   summit = PRESSURE(1500);

   testCase.verifyGreaterThan(sea_level, summit);
   testCase.verifyEqual(sea_level, 101300.0, 'AbsTol', 1e-12);
end

function test_longwave_fluxes_have_expected_signs(testCase)
   % Longwave-in should remain positive while longwave-out stays negative
   % under the project sign convention.

   SB = icemodel.physicalConstant('SB');
   Qli_dry = LONGIN(263.15, 150, SB);
   Qli_moist = LONGIN(263.15, 300, SB);
   Qle = LONGOUT(263.15, 0.98, SB);

   testCase.verifyGreaterThan(Qli_dry, 0);
   testCase.verifyGreaterThan(Qli_moist, Qli_dry);
   testCase.verifyLessThan(Qle, 0);
end

function test_turbulent_flux_kernels_follow_sign_convention(testCase)
   % Positive sensible and latent fluxes should correspond to transfer
   % toward the surface under the standard kernel sign convention.

   Qh = SENSIBLE(1e-3, 1.2, 270, 268, 1200);
   Qe = LATENT(1e-3, 1.2, 350, 300, 2.5e9, 0.622, 78000);

   testCase.verifyGreaterThan(Qh, 0);
   testCase.verifyGreaterThan(Qe, 0);
   testCase.verifyEqual(SENSIBLE(1e-3, 1.0, 270, 270, 1200), 0, ...
      'AbsTol', 1e-12);
end

function test_windcoef_and_stablefn_cover_neutral_stable_and_unstable(testCase)
   % The Monin-Obukhov stability helper should span neutral, stable, and
   % unstable branches with the expected ordering.

   [De, scoef] = WINDCOEF(4.0, 1e-3, 3.0, 3.0);
   S_neutral = STABLEFN(268, 268, 4.0, scoef);
   S_stable = STABLEFN(268, 265, 4.0, scoef);
   S_unstable = STABLEFN(268, 271, 4.0, scoef);

   testCase.verifyGreaterThan(De, 0);
   testCase.verifyEqual(S_neutral, 1.0, 'AbsTol', 5e-3);
   testCase.verifyLessThan(S_stable, 1.0);
   testCase.verifyGreaterThan(S_unstable, 1.0);
end

function test_vappress2rh_recovers_saturation_for_ice_and_water(testCase)
   % The RH conversion should map each saturation vapor pressure back to
   % roughly 100 percent for both phase relations.

   T = 268.15;
   rh_ice = VAPPRESS2RH(VAPPRESS2(T, false), T, false);
   rh_liq = VAPPRESS2RH(VAPPRESS2(T, true), T, true);

   testCase.verifyEqual(rh_ice, 100, 'AbsTol', 0.1);
   testCase.verifyEqual(rh_liq, 100, 'AbsTol', 0.1);
end

function test_getkvapor_matches_vaporheat_dry_branch(testCase)
   % GETKVAPOR is the dry-branch shortcut used by VAPORHEAT and should
   % match that branch exactly.

   [Ls, Rv, Tf] = icemodel.physicalConstant('Ls', 'Rv', 'Tf');
   T = [260; 265; 270];
   f_ice = 0.95 * ones(size(T));
   f_liq = zeros(size(T));

   [~, ~, k_vap_ref] = VAPORHEAT(T, f_ice, f_liq, Tf, Rv, Ls);
   k_vap = GETKVAPOR(T, Ls, Rv, Tf);

   testCase.verifyEqual(k_vap, k_vap_ref, 'RelTol', 1e-10);
end

function test_getgamma_stays_positive_and_responds_to_liquid(testCase)
   % Effective conductivity should stay positive and increase as the same
   % state gains more liquid water.

   [ro_ice, k_liq, Ls, Rv, Tf] = icemodel.physicalConstant( ...
      'ro_ice', 'k_liq', 'Ls', 'Rv', 'Tf');
   T = 268.15 * ones(3, 1);
   f_ice = [0.7; 0.7; 0.7];
   f_liq_dry = [0.00; 0.01; 0.02];
   f_liq_wet = [0.05; 0.08; 0.10];

   k_dry = GETGAMMA(T, f_ice, f_liq_dry, ro_ice, k_liq, Ls, Rv, Tf);
   k_wet = GETGAMMA(T, f_ice, f_liq_wet, ro_ice, k_liq, Ls, Rv, Tf);

   testCase.verifyGreaterThan(min(k_dry), 0);
   testCase.verifyGreaterThan(min(k_wet), 0);
   testCase.verifyGreaterThan(mean(k_wet), mean(k_dry));
end

function test_updateState_matches_component_kernels(testCase)
   % UPDATESTATE should stay consistent with the lower-level thermo helpers
   % it wraps into one state-update call.

   [ro_ice, ro_liq, ro_air, cv_ice, cv_liq, k_liq, roLf, Ls, Rv, Tf, ...
      fcp] = icemodel.physicalConstant('ro_ice', 'ro_liq', 'ro_air', ...
      'cv_ice', 'cv_liq', 'k_liq', 'roLf', 'Ls', 'Rv', 'Tf', 'fcp');

   T = [266; 267; 268];
   f_ice = [0.90; 0.88; 0.85];
   f_liq = [0.01; 0.02; 0.03];
   f_wat = f_liq + f_ice * ro_ice / ro_liq;

   [H, k_eff, dHdT, dLdT, drovdT, ro_vap] = UPDATESTATE(T, f_ice, ...
      f_liq, f_wat, ro_ice, ro_liq, ro_air, cv_ice, cv_liq, k_liq, roLf, ...
      Ls, Rv, Tf, fcp);

   [ro_vap_ref, drovdT_ref, k_vap_ref] = VAPORHEAT(T, f_ice, f_liq, Tf, ...
      Rv, Ls);
   k_eff_ref = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, k_vap_ref);
   H_ref = TOTALHEAT(T, f_ice, f_liq, cv_ice, cv_liq, roLf, ...
      Ls * ro_vap_ref, Tf);
   dLdT_ref = FREEZECURVE(T, f_ice, f_liq, ro_ice, ro_liq, fcp, Tf);

   testCase.verifyEqual(ro_vap, ro_vap_ref, 'RelTol', 1e-12);
   testCase.verifyEqual(drovdT, drovdT_ref, 'RelTol', 1e-12);
   testCase.verifyEqual(k_eff, k_eff_ref, 'RelTol', 1e-12);
   testCase.verifyEqual(H, H_ref, 'RelTol', 1e-12);
   testCase.verifyEqual(dHdT, cv_ice * f_ice + cv_liq * f_liq, ...
      'RelTol', 1e-12);
   testCase.verifyEqual(dLdT, dLdT_ref, 'RelTol', 1e-12);
end

function test_solvewb_returns_air_temperature_at_saturation(testCase)
   % Wet-bulb temperature should collapse to air temperature at saturation
   % and fall below it for drier air.

   [Ls, cp_air, Pa] = icemodel.physicalConstant('Ls', 'cp_air', 'P0');
   Ta = 268.15;

   [Tw_sat, ok_sat] = SOLVEWB(Ta, 100, Ls, cp_air, Pa, false);
   [Tw_dry, ok_dry] = SOLVEWB(Ta, 60, Ls, cp_air, Pa, false);

   testCase.verifyTrue(ok_sat);
   testCase.verifyTrue(ok_dry);
   testCase.verifyEqual(Tw_sat, Ta, 'AbsTol', 0.3);
   testCase.verifyLessThan(Tw_dry, Ta);
end
