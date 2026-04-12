function tests = test_flux_and_thermo_kernels
   %TEST_FLUX_AND_THERMO_KERNELS Verify local thermo and surface-flux kernels.
   tests = functiontests(localfunctions);
end

function test_pressure_decreases_with_elevation(testCase)
   % Pressure should monotonically decrease with elevation in the helper.

   sea_level = icemodel.surface.atmospheric_pressure_from_elevation(0);
   summit = icemodel.surface.atmospheric_pressure_from_elevation(1500);

   testCase.verifyGreaterThan(sea_level, summit);
   testCase.verifyEqual(sea_level, 101300.0, 'AbsTol', 1e-12);
end

function test_longwave_fluxes_have_expected_signs(testCase)
   % Longwave-in should remain positive while net longwave becomes less
   % negative as incoming longwave increases under the project sign convention.

   Qli_dry = 150;
   Qli_moist = 300;
   Qln_out = icemodel.surface.net_longwave_radiation(263.15, 0.0);
   Qln_dry = icemodel.surface.net_longwave_radiation(263.15, Qli_dry);
   Qln_moist = icemodel.surface.net_longwave_radiation(263.15, Qli_moist);

   testCase.verifyLessThan(Qln_out, 0);
   testCase.verifyGreaterThan(Qln_dry, Qln_out);
   testCase.verifyGreaterThan(Qln_moist, Qln_dry);
end

function test_turbulent_flux_kernels_follow_sign_convention(testCase)
   % Positive sensible and latent fluxes should correspond to transfer
   % toward the surface under the standard kernel sign convention.

   T_sfc = 268;
   tair = 270;
   De = 1e-3;
   stability = 1.2;
   es_sfc = 300;
   ea_atm = 350;
   psfc = 78000;
   roL = 2.5e9;

   Qh = icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux( ...
      T_sfc, tair, De, stability);
   Qe = icemodel.surface.turbulence.bulk_richardson.latent_heat_flux( ...
      es_sfc, ea_atm, De, stability, psfc, roL);

   testCase.verifyGreaterThan(Qh, 0);
   testCase.verifyGreaterThan(Qe, 0);
   testCase.verifyEqual( ...
      icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux( ...
      tair, tair, De, 1.0), 0, 'AbsTol', 1e-12);
end

function test_windcoef_and_stablefn_cover_neutral_stable_and_unstable(testCase)
   % The Monin-Obukhov stability helper should span neutral, stable, and
   % unstable branches with the expected ordering.

   T_sfc_neutral = 268;
   tair = 268;
   T_sfc_stable = 265;
   T_sfc_unstable = 271;
   wspd = 4.0;
   z0_bulk = 1e-3;
   z_tair = 3.0;
   z_wind = 3.0;

   [De, br_coefs] = ...
      icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
      wspd, z0_bulk, z_tair, z_wind);

   S_neutral = icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
      T_sfc_neutral, tair, wspd, br_coefs);

   S_stable = icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
      T_sfc_stable, tair, wspd, br_coefs);

   S_unstable = icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
      T_sfc_unstable, tair, wspd, br_coefs);

   testCase.verifyGreaterThan(De, 0);
   testCase.verifyEqual(S_neutral, 1.0, 'AbsTol', 5e-3);
   testCase.verifyLessThan(S_stable, 1.0);
   testCase.verifyGreaterThan(S_unstable, 1.0);
end

function test_stablefn_derivative_matches_finite_difference(testCase)
   % stability_factor should return a Ts derivative consistent with finite
   % differences in stable, unstable, and near-neutral regimes.

   tair = 268.15;
   wspd = 4.0;

   % Aerodynamics
   z0_bulk = 1e-3;
   z_tair = 3.0;
   z_wind = 3.0;
   [~, br_coefs] = ...
      icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
      4.0, z0_bulk, z_tair, z_wind);

   % Temperature perturbation
   h = 1e-6;

   % Three T_sfc test cases
   T_sfc_cases = [tair - 2.0, tair, tair + 2.0];

   for n = 1:numel(T_sfc_cases)

      T_sfc = T_sfc_cases(n);

      [~, dSdTs] = ...
         icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
         T_sfc, tair, wspd, br_coefs);

      S_plus = icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
         T_sfc + h, tair, wspd, br_coefs);

      S_minus = icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
         T_sfc - h, tair, wspd, br_coefs);

      dSdTs_fd = (S_plus - S_minus) / (2 * h);

      testCase.verifyEqual(dSdTs, dSdTs_fd, 'RelTol', 1e-5, ...
         sprintf('Derivative mismatch at T_sfc = %.3f K', T_sfc));
   end
end

function test_stablefn_neutral_blend_matches_endpoint_formulas(testCase)
   % The near-neutral branch should linearly blend the stable and unstable
   % endpoint formulas across the configured transition width.

   dT0 = icemodel.parameterLookup( ...
      'thf_bulk_richardson_neutral_transition_width');

   tair = 268.15;
   wspd = 4.0;

   % Bulk roughness and observation heights
   z0_bulk = 1e-3;
   z_tair = 3.0;
   z_wind = 3.0;

   % Bulk Richardson coefficients
   [~, br_coefs] = ...
      icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
      wspd, z0_bulk, z_tair, z_wind);

   % Louis bulk Richardson stability coefficients
   B1 = br_coefs(2) / (tair * wspd ^ 2);
   B2 = br_coefs(3) / (sqrt(tair) * wspd);

   % Stability at the transition
   S_stable_edge = (1 + 0.5 * B1 * dT0) ^ -2;
   S_unstable_edge = 1 + B1 * dT0 / (1 + B2 * sqrt(dT0));
   S_expected = 0.5 * (S_stable_edge + S_unstable_edge);
   dSdT_expected = (S_unstable_edge - S_stable_edge) / (2 * dT0);

   [S_neutral, dSdT_neutral] = ...
      icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
      tair, tair, wspd, br_coefs);

   testCase.verifyEqual(S_neutral, S_expected, 'RelTol', 1e-12);
   testCase.verifyEqual(dSdT_neutral, dSdT_expected, 'RelTol', 1e-12);
end

function test_vappress2rh_recovers_saturation_for_ice_and_water(testCase)
   % The RH conversion should map each saturation vapor pressure back to
   % roughly 100 percent for both phase relations.

   T = 268.15;
   rh_ice = icemodel.vapor.relative_humidity_from_vapor_pressure( ...
      icemodel.vapor.saturation_vapor_pressure(T, false), T, false);
   rh_liq = icemodel.vapor.relative_humidity_from_vapor_pressure( ...
      icemodel.vapor.saturation_vapor_pressure(T, true), T, true);

   testCase.verifyEqual(rh_ice, 100, 'AbsTol', 0.1);
   testCase.verifyEqual(rh_liq, 100, 'AbsTol', 0.1);
end

function test_vapork_matches_vapordensity_times_diffusivity(testCase)
   % icemodel.vapor.vapor_thermal_diffusion_coefficient should equal
   % Ls * De * dro_vapdT for dry cells and Lv * De * dro_vapdT for wet cells.

   [Ls, Lv] = icemodel.physicalConstant('Ls', 'Lv');
   T = [260; 265; 270];
   f_liq = zeros(size(T));

   [~, dro_vapdT] = icemodel.vapor.saturation_vapor_density(T, f_liq);
   De = icemodel.vapor.vapor_diffusivity(T);
   k_vap_manual = Ls * De .* dro_vapdT;
   k_vap = icemodel.vapor.vapor_thermal_diffusion_coefficient(T, f_liq);

   testCase.verifyEqual(k_vap, k_vap_manual, 'RelTol', 1e-10);

   % Test wet cells use Lv
   De_wet = icemodel.vapor.vapor_diffusivity(T);
   f_liq_wet = 0.05 * ones(size(T));

   k_vap_wet = icemodel.vapor.vapor_thermal_diffusion_coefficient(T, f_liq_wet);

   [~, dro_vapdT_wet] = icemodel.vapor.saturation_vapor_density(T, f_liq_wet);

   k_vap_wet_manual = Lv * De_wet .* dro_vapdT_wet;

   testCase.verifyEqual(k_vap_wet, k_vap_wet_manual, 'RelTol', 1e-10);
end

function test_bulk_thermalk_stays_positive_and_responds_to_liquid(testCase)
   % Effective conductivity should stay positive and increase as the same
   % state gains more liquid water.

   T = 268.15 * ones(3, 1);
   f_ice = [0.7; 0.7; 0.7];
   f_liq_dry = [0.00; 0.01; 0.02];
   f_liq_wet = [0.05; 0.08; 0.10];

   k_dry = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq_dry);
   k_wet = icemodel.column.bulk_thermal_conductivity(T, f_ice, f_liq_wet);

   testCase.verifyGreaterThan(min(k_dry), 0);
   testCase.verifyGreaterThan(min(k_wet), 0);
   testCase.verifyGreaterThan(mean(k_wet), mean(k_dry));
end

function test_updateState_matches_component_kernels(testCase)
   % icemodel.column.updatestate should stay consistent with the lower-level
   % thermo helpers it wraps into one state-update call.

   [cv_ice, cv_liq] = icemodel.physicalConstant('cv_ice', 'cv_liq');

   T = [266; 267; 268];
   f_ice = [0.90; 0.88; 0.85];
   f_liq = [0.01; 0.02; 0.03];
   f_wat = icemodel.water_fraction(f_ice, f_liq);

   [H, k_eff, dHdT, dLdT, drovdT, ro_vap] = icemodel.column.updatestate( ...
      T, f_ice, f_liq, f_wat);

   [ro_vap_ref, drovdT_ref] = icemodel.vapor.saturation_vapor_density( ...
      T, f_liq);

   k_vap_ref = icemodel.vapor.vapor_thermal_diffusion_coefficient( ...
      T, f_liq, drovdT_ref);

   k_eff_ref = icemodel.column.bulk_thermal_conductivity( ...
      T, f_ice, f_liq, k_vap_ref);

   H_ref = icemodel.column.total_enthalpy( ...
      T, f_ice, f_liq, f_wat, ro_vap_ref);

   dLdT_ref = icemodel.column.liquid_fraction_derivative( ...
      T, f_ice, f_liq);

   testCase.verifyEqual(ro_vap, ro_vap_ref, 'RelTol', 1e-12);
   testCase.verifyEqual(drovdT, drovdT_ref, 'RelTol', 1e-12);
   testCase.verifyEqual(k_eff, k_eff_ref, 'RelTol', 1e-12);
   testCase.verifyEqual(H, H_ref, 'RelTol', 1e-12);
   testCase.verifyEqual(dHdT, cv_ice * f_ice + cv_liq * f_liq, 'RelTol', 1e-12);
   testCase.verifyEqual(dLdT, dLdT_ref, 'RelTol', 1e-12);
end

function test_twetbulb_returns_air_temperature_at_saturation(testCase)
   % Wet-bulb temperature should collapse to air temperature at saturation
   % and fall below it for drier air.

   Pa = icemodel.physicalConstant('P0');
   Ta = 268.15;
   rh_sat = 100;
   rh_dry = 60;

   [Tw_sat, ok_sat] = icemodel.vapor.wet_bulb_temperature(Ta, rh_sat, Pa, false);
   [Tw_dry, ok_dry] = icemodel.vapor.wet_bulb_temperature(Ta, rh_dry, Pa, false);

   testCase.verifyTrue(ok_sat);
   testCase.verifyTrue(ok_dry);
   testCase.verifyEqual(Tw_sat, Ta, 'AbsTol', 0.3);
   testCase.verifyLessThan(Tw_dry, Ta);
end

function test_vaporinit_coefficients_consistent(testCase)
   % icemodel.vapor.initialize_vapor_model output should reproduce the reference
   % vapor pressure at the triple point when plugged into the Rankine-Kirchhoff
   % formula.

   [Tf, es0] = icemodel.physicalConstant('Tf', 'es0');
   [al, bl, cl, ai, bi, ci] = icemodel.vapor.initialize_vapor_model();

   es_liq = al * exp(bl / Tf) * Tf ^ cl;
   es_ice = ai * exp(bi / Tf) * Tf ^ ci;

   % Ambaum formula is a curve fit; es at the triple point is close but not
   % identical to the reference es0. Both phases should agree to ~0.1%.
   testCase.verifyEqual(es_liq, es0, 'RelTol', 1e-3);
   testCase.verifyEqual(es_ice, es0, 'RelTol', 1e-3);

   % Liquid and ice curves must agree at the triple point.
   testCase.verifyEqual(es_liq, es_ice, 'RelTol', 1e-10);
end

function test_ambaum_buck_agreement(testCase)
   % Ambaum and Buck saturation vapor pressure should agree within 1% over
   % the temperature range 230-273 K (ice regime).

   T = (230:273)';
   Tf = 273.16;

   es_ambaum = icemodel.vapor.saturation_vapor_pressure(T, false);
   es_buck = icemodel.kernels.buckVaporModel(T, Tf, false);

   testCase.verifyEqual(es_ambaum, es_buck, 'RelTol', 0.01);
end

function test_parameterLookup_returns_expected_fields(testCase)
   % parameterLookup should return a struct with all canonical fields when
   % called with 'all'.

   params = icemodel.parameterLookup('all');

   testCase.verifyTrue(isstruct(params));
   expected = {'al','bl','cl','ai','bi','ci','nd','De0','emiss','fcp', ...
      'thf_bulk_richardson_eta', 'thf_bulk_richardson_psi', ...
      'thf_bulk_richardson_neutral_transition_width', 'thf_bulk_iter_max', ...
      'thf_bulk_holtslag_aa', 'thf_bulk_dyer_gamma', ...
      'thf_bulk_andreas_ch1', 'thf_bulk_smeets_a0'};
   for k = 1:numel(expected)
      testCase.verifyTrue(isfield(params, expected{k}), ...
         sprintf('Missing field: %s', expected{k}));
   end
end

function test_atmosphericVaporPressure_matches_direct_contract(testCase)
   % The centralized helper should preserve the current ea contract.

   Ta = 268.15;
   rh = 72.0;
   ea_ref = icemodel.vapor.saturation_vapor_pressure(Ta, false) * rh / 100;
   ea = icemodel.surface.atmospheric_vapor_pressure(Ta, rh, false);

   testCase.verifyEqual(ea, ea_ref, 'RelTol', 1e-12);
end

function test_firn_thermalk_positive_and_density_dependent(testCase)
   % icemodel.column.firn_thermal_conductivity should produce positive
   % conductivity that increases with ice fraction (density).

   T = [255; 265; 270];
   f_ice = [0.95; 0.80; 0.60];

   k_firn = icemodel.column.firn_thermal_conductivity(T, f_ice);

   testCase.verifyGreaterThan(min(k_firn), 0);
   testCase.verifyGreaterThan(k_firn(1), k_firn(3), ...
      'Higher ice fraction should yield higher thermal conductivity');
end

function test_ambaum_derivative_chain_consistency(testCase)
   % Verify that all analytical derivatives in the Ambaum/Romps chain
   % match centered finite differences to within RelTol 1e-6.
   %
   % This tests the complete chain: es -> des_dT -> d2es_dT2,
   % ro_vap -> dro_vapdT -> d2ro_vapdT2, and cross-function agreement
   % between icemodel.vapor.saturation_vapor_pressure,
   % icemodel.vapor.saturation_vapor_density, and
   % icemodel.vapor.vapor_thermal_diffusion_coefficient.

   Ls = icemodel.physicalConstant('Ls');

   T = (235:0.5:273)';
   h = 1e-4;  % finite difference step [K]

   %%% Vapor pressure derivatives via icemodel.vapor.saturation_vapor_pressure

   [~, des_dT, d2es_dT2] = icemodel.vapor.saturation_vapor_pressure(T, false);
   es_p = icemodel.vapor.saturation_vapor_pressure(T + h, false);
   es_m = icemodel.vapor.saturation_vapor_pressure(T - h, false);
   [~, des_dT_p] = icemodel.vapor.saturation_vapor_pressure(T + h, false);
   [~, des_dT_m] = icemodel.vapor.saturation_vapor_pressure(T - h, false);

   % 1. des_dT: numerical first derivative of es
   des_dT_num = (es_p - es_m) / (2 * h);
   testCase.verifyEqual(des_dT, des_dT_num, 'RelTol', 1e-6, ...
      'des_dT analytical vs finite difference');

   % 2. d2es_dT2: numerical second derivative of es
   d2es_dT2_num = (des_dT_p - des_dT_m) / (2 * h);
   testCase.verifyEqual(d2es_dT2, d2es_dT2_num, 'RelTol', 1e-6, ...
      'd2es_dT2 analytical vs finite difference');

   %%% Vapor density derivatives via icemodel.vapor.saturation_vapor_density

   f_liq = zeros(size(T));
   [~, dro_vapdT, d2ro_vapdT2] = icemodel.vapor.saturation_vapor_density(T, f_liq);
   ro_vap_p = icemodel.vapor.saturation_vapor_density(T + h, f_liq);
   ro_vap_m = icemodel.vapor.saturation_vapor_density(T - h, f_liq);
   [~, dro_vapdT_p] = icemodel.vapor.saturation_vapor_density(T + h, f_liq);
   [~, dro_vapdT_m] = icemodel.vapor.saturation_vapor_density(T - h, f_liq);

   % 3. dro_vapdT: numerical first derivative of ro_vap
   dro_vapdT_num = (ro_vap_p - ro_vap_m) / (2 * h);
   testCase.verifyEqual(dro_vapdT, dro_vapdT_num, 'RelTol', 1e-6, ...
      'dro_vapdT analytical vs finite difference');

   % 4. d2ro_vapdT2: numerical second derivative of ro_vap
   d2ro_vapdT2_num = (dro_vapdT_p - dro_vapdT_m) / (2 * h);
   testCase.verifyEqual(d2ro_vapdT2, d2ro_vapdT2_num, 'RelTol', 1e-6, ...
      'd2ro_vapdT2 analytical vs finite difference');

   %%% Cross-function agreement

   % 5. icemodel.vapor.vapor_thermal_diffusion_coefficient reuse path should
   % match its internal icemodel.vapor.saturation_vapor_density path
   k_vap_vk = icemodel.vapor.vapor_thermal_diffusion_coefficient( ...
      T, f_liq);
   k_vap_reuse = icemodel.vapor.vapor_thermal_diffusion_coefficient( ...
      T, f_liq, dro_vapdT);
   testCase.verifyEqual(k_vap_reuse, k_vap_vk, 'RelTol', 1e-12, ...
      ['icemodel.vapor.vapor_thermal_diffusion_coefficient ' ...
      'reuse path vs internal dro_vapdT path']);

   % 6. icemodel.vapor.vapor_thermal_diffusion_coefficient k_vap matches Ls * De
   % * dro_vapdT for dry cells
   De = icemodel.vapor.vapor_diffusivity(T);
   k_vap_manual = Ls * De .* dro_vapdT;
   testCase.verifyEqual(k_vap_vk, k_vap_manual, 'RelTol', 1e-10, ...
      ['icemodel.vapor.vapor_thermal_diffusion_coefficient ' ...
      'k_vap vs manual Ls*De*dro_vapdT']);
end
