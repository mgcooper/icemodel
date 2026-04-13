function tests = test_turbulent_flux_schemes
   %TEST_TURBULENT_FLUX_SCHEMES Verify the shared THF helper layer.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Build one synthetic workspace that can drive helper-level and
   % operational turbulent-flux tests.

   testCase.TestData.workspace = ...
      icemodel.test.fixtures.makeSyntheticWorkspace( ...
      2016, configure=true, nsteps=12, dt_seconds=900);
end

function teardown(testCase)
   % Remove the shared synthetic workspace after the file-level tests end.

   icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      testCase.TestData.workspace);
end

function test_surface_roughness_length_switches_on_snow_depth(testCase)
   % Snow depth, not bulk density alone, should decide snow versus ice.

   snow_depth = [0.0, 0.2, 0.2];
   ro_sfc = [450.0, 350.0, 700.0];
   z0_ice = [0.02, 0.02, 0.02];
   z0_snow_low_density = [1e-3, 1e-3, 1e-3];
   z0_snow_high_density = [5e-4, 5e-4, 5e-4];

   z0m_ice = icemodel.surface.surface_roughness_length(snow_depth(1), ...
      ro_sfc(1), z0_ice(1), z0_snow_low_density(1), z0_snow_high_density(1));

   z0m_snow_low = icemodel.surface.surface_roughness_length(snow_depth(2), ...
      ro_sfc(2), z0_ice(2), z0_snow_low_density(2), z0_snow_high_density(2));

   z0m_snow_high = icemodel.surface.surface_roughness_length(snow_depth(3), ...
      ro_sfc(3), z0_ice(3), z0_snow_low_density(3), z0_snow_high_density(3));

   testCase.verifyEqual(z0m_ice, 0.02, 'AbsTol', 1e-12);
   testCase.verifyEqual(z0m_snow_low, 1e-3, 'AbsTol', 1e-12);
   testCase.verifyEqual(z0m_snow_high, 5e-4, 'AbsTol', 1e-12);
end

function test_surface_roughness_length_vector_diags_use_external_loop(testCase)
   % surface_roughness_length is intentionally scalar. Diagnostic callers
   % that want a vector of states should loop externally.

   snow_depth = [0.0; 0.2; 0.2];
   ro_sfc = [450.0; 350.0; 700.0];

   z0m = NaN(size(snow_depth));
   for n = 1:numel(snow_depth)
      z0m(n) = icemodel.surface.surface_roughness_length( ...
         snow_depth(n), ro_sfc(n), 0.02, 1e-3, 5e-4);
   end

   testCase.verifyEqual(z0m, [0.02; 1e-3; 5e-4], 'AbsTol', 1e-12);
end

function test_forcing_snow_depth_hook_switches_to_snow_roughness(testCase)
   % The optional forcing snow-depth hook should switch the bulk-MO
   % roughness path without requiring a full snow model.

   workspace = testCase.TestData.workspace;

   % Create synthetic snow depth data
   snow_depth = 0.08 + zeros(workspace.nsteps, 1);

   % Create a synthetic met file
   icemodel.test.fixtures.makeSyntheticMetFile(2016, ...
      'sitename', workspace.sitename, ...
      'forcings', workspace.forcings, ...
      'nsteps', workspace.nsteps, ...
      'dt_seconds', workspace.dt_seconds, ...
      'snow_depth', snow_depth, ...
      'metdir', workspace.metdir);

   % Create synthetic model opts
   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      workspace, 'icemodel', 2016, solver=1, seb_solver=2, ...
      turbulent_flux_scheme='monin_obukhov', ...
      use_forcing_snow_depth_for_thf=true, ...
      z0_ice=0.02, z0_snow_low_density=0.0015, ...
      testname='forcing_snow_depth_hook');

   % Load the met data
   met = icemodel.loadmet(opts);

   % Load snow depth via icemodel.surface.initialize_surface_forcings
   [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, forcing_snow_depth] = ...
      icemodel.surface.initialize_surface_forcings(opts);

   snow_depth_scalar = icemodel.surface.resolve_forcing_snow_depth( ...
      forcing_snow_depth, 1, opts.use_forcing_snow_depth_for_thf);

   % Set the momentum roughness length
   z0m = icemodel.surface.surface_roughness_length( ...
      snow_depth_scalar, 350.0, opts.z0_ice, ...
      opts.z0_snow_low_density, opts.z0_snow_high_density);

   testCase.verifyEqual(met.snow_depth(1), snow_depth(1), 'AbsTol', 1e-12);
   testCase.verifyEqual(snow_depth_scalar, snow_depth(1), 'AbsTol', 1e-12);
   testCase.verifyEqual(z0m, opts.z0_snow_low_density, 'AbsTol', 1e-12);
end

function test_bulk_richardson_dispatch_matches_legacy_flux_kernels(testCase)
   % The dispatcher should reproduce the legacy sensible/latent formulas
   % when the configured scheme is bulk_richardson.

   s = icemodel.test.fixtures.makeSyntheticColumnState( ...
      testCase.TestData.workspace, 'icemodel', solver=1, ...
      testname='bulk_richardson_dispatch');

   T_sfc = icemodel.surface.physical_surface_temperature(s.Ts);

   stability = icemodel.surface.turbulence.bulk_richardson.stability_factor( ...
      T_sfc, s.tair, s.wspd, s.br_coefs);

   es_sfc = icemodel.vapor.saturation_vapor_pressure(T_sfc, s.liqflag);

   Qe_ref = icemodel.surface.turbulence.bulk_richardson.latent_heat_flux( ...
      es_sfc, s.ea_atm, s.De, stability, s.psfc, s.ro_air_Lv);

   Qh_ref = icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux( ...
      T_sfc, s.tair, s.De, stability);

   [Qe, Qh, diag] = icemodel.surface.diagnose_turbulent_heat_fluxes(T_sfc, ...
      s.tair, s.wspd, s.psfc, s.ea_atm, s.De, s.br_coefs, s.ro_sfc, ...
      s.snow_depth, s.ro_air_Lv, s.liqflag, s.opts);

   testCase.verifyEqual(Qe, Qe_ref, 'RelTol', 1e-12);
   testCase.verifyEqual(Qh, Qh_ref, 'RelTol', 1e-12);
   testCase.verifyEqual(diag.scheme, 'bulk_richardson');
end

function test_bulk_richardson_scalar_exchange_weakens_rough_ice_fluxes(testCase)
   % The scalar-exchange diagnostic experiment should keep production fluxes
   % unchanged while showing weaker scalar transfer over rough ice.

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=1, ...
      z0_bulk=0.05, testname='bulk_richardson_scalar_experiment');

   [roLs, epsilon, cv_air] = icemodel.physicalConstant( ...
      'roLs', 'epsilon', 'cv_air');

   tair = 268.15;
   T_sfc = 263.15;
   wspd = 6.0;
   psfc = 90000.0;
   rh = 75.0;
   ea_atm = icemodel.surface.atmospheric_vapor_pressure(tair, rh, false);

   [De, br_coefs] = ...
      icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
      wspd, opts.z0_bulk, opts.z_tair, opts.z_wind);

   ro_sfc = 650.0;
   snow_depth = 0.0;
   liqflag = false;

   [Qe, Qh, diag] = icemodel.surface.diagnose_turbulent_heat_fluxes(T_sfc, ...
      tair, wspd, psfc, ea_atm, De, br_coefs, ro_sfc, snow_depth, roLs, ...
      liqflag, opts);

   Qh_expected = cv_air * De * diag.stability_factor * (tair - T_sfc);
   Qe_expected = roLs * De * diag.stability_factor ...
      * (epsilon / psfc * (ea_atm - diag.es_sfc));

   testCase.verifyTrue(diag.scalar_exchange_active);
   testCase.verifyGreaterThan(diag.scalar_exchange_z_obs, diag.z0m);
   testCase.verifyLessThan(diag.z0h, diag.z0m);
   testCase.verifyLessThan(diag.z0q, diag.z0m);
   testCase.verifyLessThan(diag.scalar_exchange_De_h, De);
   testCase.verifyLessThan(diag.scalar_exchange_De_e, De);

   testCase.verifyEqual(Qh, Qh_expected, 'RelTol', 1e-12);
   testCase.verifyEqual(Qe, Qe_expected, 'RelTol', 1e-12);
   testCase.verifyLessThan(abs(diag.scalar_exchange_Qh), abs(Qh));
   testCase.verifyLessThan(abs(diag.scalar_exchange_Qe), abs(Qe));
end

function test_bulk_richardson_scalar_exchange_noop_for_calm_air(testCase)
   % The scalar-exchange experiment should degrade gracefully when the
   % current De contract provides no usable aerodynamic signal.

   roLs = icemodel.physicalConstant('roLs');

   tair = 268.15;
   T_sfc = 263.15;
   wspd = 0.0;
   psfc = 90000.0;
   De = 0.0;
   ea_atm = 300.0;
   br_coefs = [NaN, NaN, NaN];
   liqflag = false;
   z0_bulk = 0.05;

   [Qe, Qh, diag] = ...
      icemodel.surface.turbulence.bulk_richardson.turbulent_heat_flux( ...
      T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, roLs, liqflag, z0_bulk);

   testCase.verifyFalse(diag.scalar_exchange_active);
   testCase.verifyEqual(diag.scalar_exchange_De_h, 0.0, 'AbsTol', 1e-12);
   testCase.verifyEqual(diag.scalar_exchange_De_e, 0.0, 'AbsTol', 1e-12);
   testCase.verifyEqual(diag.scalar_exchange_Qh, Qh, 'AbsTol', 1e-12);
   testCase.verifyEqual(diag.scalar_exchange_Qe, Qe, 'AbsTol', 1e-12);
end

function test_bulk_mo_finite_and_weaker_scalar_exchange_for_rough_ice(testCase)
   % Rough bare ice should not force scalar exchange to increase as much as
   % the legacy single-z0 formulation at the same momentum roughness.

   opts_liston = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=1, z0_bulk=0.05, ...
      testname='rough_bulk_richardson');

   opts_vanas = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=1, ...
      seb_solver=2, turbulent_flux_scheme='monin_obukhov', z0_bulk=0.05, ...
      z0_ice=0.05, testname='rough_bulk_mo');

   roLs = icemodel.physicalConstant('roLs');

   % Setup required forcings
   tair = 268.15;
   T_sfc = 263.15;
   wspd = 6.0;
   psfc = 90000.0;
   rh = 75.0;
   ea_atm = icemodel.surface.atmospheric_vapor_pressure(tair, rh, false);
   [De, br_coefs] = ...
      icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
      wspd, 0.05, opts_liston.z_tair, opts_liston.z_wind);

   % Compute turbulent heat fluxes
   ro_sfc = 650.0;
   snow_depth = 0.0;
   liqflag = false;

   [Qe_liston, Qh_liston] = icemodel.surface.diagnose_turbulent_heat_fluxes( ...
      T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, ro_sfc, snow_depth, ...
      roLs, liqflag, opts_liston);
   [Qe_vanas, Qh_vanas, diag_vanas] = ...
      icemodel.surface.diagnose_turbulent_heat_fluxes(T_sfc, tair, wspd, ...
      psfc, ea_atm, De, br_coefs, ro_sfc, snow_depth, roLs, liqflag, ...
      opts_vanas);

   testCase.verifyTrue(all(isfinite([Qe_vanas, Qh_vanas, diag_vanas.L, ...
      diag_vanas.u_star, diag_vanas.z0h, diag_vanas.z0q])));
   testCase.verifyLessThan(abs(Qh_vanas), abs(Qh_liston));
   testCase.verifyLessThan(abs(Qe_vanas), abs(Qe_liston));
end

function test_bulk_mo_vector_diagnostics_use_external_loop(testCase)
   % The bulk-MO solver is intentionally scalar in production. Diagnostic
   % callers that want a vector of independent states should loop outside
   % the helper rather than relying on internal broadcasting.

   opts_vanas = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=1, ...
      seb_solver=2, turbulent_flux_scheme='monin_obukhov', ...
      testname='vector_bulk_mo');

   roLs = icemodel.physicalConstant('roLs');

   % Setup forcings
   tair = 268.15 * ones(3, 1);
   T_sfc = [263.15; 264.15; 265.15];
   wspd = [4.0; 5.0; 6.0];
   psfc = 90000.0 * ones(3, 1);
   ea_atm = icemodel.surface.atmospheric_vapor_pressure( ...
      tair, 75 * ones(3, 1), false);

   ro_sfc = 650.0 * ones(3, 1);
   snow_depth = zeros(3, 1);
   liqflag = false(3, 1);

   % Compute turbulent heat fluxes
   Qe = NaN(size(T_sfc));
   Qh = NaN(size(T_sfc));
   diag = cell(size(T_sfc));
   for n = 1:numel(T_sfc)
      [Qe(n), Qh(n), diag{n}] = ...
         icemodel.surface.turbulence.monin_obukhov.turbulent_heat_flux( ...
         T_sfc(n), tair(n), wspd(n), psfc(n), ea_atm(n), ro_sfc(n), ...
         snow_depth(n), roLs, liqflag(n), opts_vanas);
   end

   testCase.verifyEqual(size(Qe), size(T_sfc));
   testCase.verifyEqual(size(Qh), size(T_sfc));
   testCase.verifyEqual(size(diag), size(T_sfc));
   testCase.verifyTrue(all(isfinite(Qe)));
   testCase.verifyTrue(all(isfinite(Qh)));
   testCase.verifyTrue(all(cellfun(@(x) isstruct(x) && isfield(x, 'L'), diag)));
end

function test_monin_obukhov_respects_surface_phase_switch(testCase)
   % Surface saturation should change when liqflag switches at the same
   % atmospheric state.

   opts_vanas = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=1, seb_solver=2, ...
      turbulent_flux_scheme='monin_obukhov', testname='phase_bulk_mo');

   roLs = icemodel.physicalConstant('roLs');

   % Setup forcings
   tair = 273.15;
   T_sfc = 272.15;
   wspd = 5.0;
   psfc = 90000.0;
   ea_atm = icemodel.surface.atmospheric_vapor_pressure(tair, 90.0, false);
   [De, br_coefs] = ...
      icemodel.surface.turbulence.bulk_richardson.exchange_coefficients( ...
      wspd, opts_vanas.z0_bulk, opts_vanas.z_tair, opts_vanas.z_wind);
   ro_sfc = 600.0;
   snow_depth = 0.0;

   % Compute turbulent heat fluxes
   [Qe_ice, ~, diag_ice] = icemodel.surface.diagnose_turbulent_heat_fluxes( ...
      T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, ro_sfc, snow_depth, ...
      roLs, false, opts_vanas);
   [Qe_liq, ~, diag_liq] = icemodel.surface.diagnose_turbulent_heat_fluxes( ...
      T_sfc, tair, wspd, psfc, ea_atm, De, br_coefs, ro_sfc, snow_depth, ...
      roLs, true, opts_vanas);

   testCase.verifyNotEqual(diag_ice.es_sfc, diag_liq.es_sfc);
   testCase.verifyNotEqual(Qe_ice, Qe_liq);
end

function test_monin_obukhov_cold_state_remains_continuous(testCase)
   % The cold stable bulk-MO replay state that used to trip the 2015
   % spinup run should vary smoothly across small Ts perturbations.

   opts_vanas = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=1, ...
      seb_solver=2, turbulent_flux_scheme='monin_obukhov', ...
      testname='cold_bulk_mo_continuity');

   tair = 260.755;
   wspd = 1.41;
   psfc = 87416.0;
   ea_atm = 150.425;
   T_sfc_1 = 245.2210;
   T_sfc_2 = 245.2220;
   ro_sfc = 135.536437794510;
   snow_depth = 0.0;
   liqflag = false;
   ro_air_Lv = 3.664362e6;
   De = 0.0;
   br_coefs = [];

   [Qe1, Qh1, diag1] = icemodel.surface.diagnose_turbulent_heat_fluxes( ...
      T_sfc_1, tair, wspd, psfc, ea_atm, De, br_coefs, ro_sfc, snow_depth, ...
      ro_air_Lv, liqflag, opts_vanas);

   [Qe2, Qh2, diag2] = icemodel.surface.diagnose_turbulent_heat_fluxes( ...
      T_sfc_2, tair, wspd, psfc, ea_atm, De, br_coefs, ro_sfc, snow_depth, ...
      ro_air_Lv, liqflag, opts_vanas);

   testCase.verifyTrue(all(isfinite([Qe1, Qh1, Qe2, Qh2, diag1.L, diag2.L])));
   testCase.verifyLessThan(abs(Qh2 - Qh1), 1e-2);
   testCase.verifyLessThan(abs(Qe2 - Qe1), 1e-3);
end

function test_monin_obukhov_linearization_matches_finite_difference(testCase)
   % The bulk-MO Robin linearization should match a finite-difference slope
   % of the non-conductive surface flux.

   s = icemodel.test.fixtures.makeSyntheticColumnState( ...
      testCase.TestData.workspace, 'icemodel', solver=3, seb_solver=2, ...
      turbulent_flux_scheme='monin_obukhov', z0_ice=0.02, ...
      testname='bulk_mo_robin_linearization');

   [Fc, Fp, diag] = ...
      icemodel.surface.turbulence.monin_obukhov.surface_flux_linearization( ...
      s.Ts, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.De, s.ea_atm, s.chi, s.ro_air_Lv, s.liqflag, s.ro_sfc, ...
      s.snow_depth, s.opts);

   h = 1e-5;
   q_plus = surface_flux_bulk_mo(s.Ts + h, s);
   q_minus = surface_flux_bulk_mo(s.Ts - h, s);
   dq_fd = (q_plus - q_minus) / (2 * h);

   testCase.verifyTrue(all(isfinite([Fc, Fp, diag.q_surface, ...
      diag.dq_surface_dTs])));
   testCase.verifyEqual(Fc + Fp * s.Ts, diag.q_surface, 'RelTol', 1e-10);
   testCase.verifyEqual(Fp, dq_fd, 'RelTol', 5e-3);
end

function test_monin_obukhov_synthetic_icemodel_run_completes(testCase)
   % A synthetic icemodel run with the supported bulk-MO configuration
   % should complete and return finite Tsfc, Qh, and Qe outputs.

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=1, ...
      seb_solver=2, turbulent_flux_scheme='monin_obukhov', z0_ice=0.02, ...
      testname='icemodel_bulk_mo_run');

   [ice1, ~, opts_out] = icemodel(opts);

   testCase.verifyEqual(opts_out.turbulent_flux_scheme, 'monin_obukhov');
   testCase.verifyTrue(all(isfinite(ice1.Tsfc)));
   testCase.verifyTrue(all(isfinite(ice1.Qe)));
   testCase.verifyTrue(all(isfinite(ice1.Qh)));
end

function test_monin_obukhov_synthetic_icemodel_robin_run_completes(testCase)
   % A synthetic icemodel run with the Robin bulk-MO configuration should
   % complete and return finite surface state and turbulent flux outputs.

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=3, ...
      seb_solver=2, turbulent_flux_scheme='monin_obukhov', z0_ice=0.02, ...
      testname='icemodel_bulk_mo_robin_run');

   [ice1, ~, opts_out] = icemodel(opts);

   testCase.verifyEqual(opts_out.turbulent_flux_scheme, 'monin_obukhov');
   testCase.verifyEqual(opts_out.solver, 3);
   testCase.verifyTrue(all(isfinite(ice1.Tsfc)));
   testCase.verifyTrue(all(isfinite(ice1.Qe)));
   testCase.verifyTrue(all(isfinite(ice1.Qh)));
end

function q_surface = surface_flux_bulk_mo(T_sfc, s)
   %SURFACE_FLUX_BULK_MO Evaluate the non-conductive surface flux at T_sfc.

   [Qe, Qh] = icemodel.surface.diagnose_turbulent_heat_fluxes(T_sfc, ...
      s.tair, s.wspd, s.psfc, s.ea_atm, s.De, s.br_coefs, s.ro_sfc, ...
      s.snow_depth, s.ro_air_Lv, s.liqflag, s.opts);
   Qa = icemodel.surface.advective_heat_flux(s.ppt, s.tppt, s.cv_liq);
   Qsn = icemodel.surface.net_shortwave_radiation(s.swd, s.albedo, s.chi);
   Qln = icemodel.surface.net_longwave_radiation(T_sfc, s.lwd);
   q_surface = icemodel.surface.evaluate_surface_energy_balance( ...
      Qsn, Qln, Qh, Qe, 0.0, Qa, 0.0);
end

function test_vapor_tendency_uses_physical_surface_temperature(testCase)
   % Potential vapor tendency should evaluate Qe using the physical
   % diagnosed surface temperature even if the solver-internal Ts exceeds Tf.

   s = icemodel.test.fixtures.makeSyntheticColumnState( ...
      testCase.TestData.workspace, 'icemodel', solver=1, seb_solver=2, ...
      turbulent_flux_scheme='monin_obukhov', z0_ice=0.02, ...
      testname='potential_vapor_tendency');
   Lv = icemodel.physicalConstant('Lv');

   Ts_solver = s.Tf + 3.0;

   [d_pevp, pevp, Qe, Ts_phys] = ...
      icemodel.surface.potential_surface_vapor_tendency(Ts_solver, s.tair, ...
      s.wspd, s.psfc, s.ea_atm, s.De, s.br_coefs, s.snow_depth, s.f_ice(1), ...
      s.f_liq(1), s.opts.dt, s.dz(1), s.ro_air_Lv, s.liqflag, s.opts);

   [Qe_phys, ~] = icemodel.surface.diagnose_turbulent_heat_fluxes( ...
      icemodel.surface.physical_surface_temperature(Ts_solver), ...
      s.tair, s.wspd, s.psfc, s.ea_atm, s.De, s.br_coefs, s.ro_sfc, ...
      s.snow_depth, s.ro_air_Lv, s.liqflag, s.opts);

   [Qe_uncapped, ~] = icemodel.surface.diagnose_turbulent_heat_fluxes( ...
      Ts_solver, s.tair, s.wspd, s.psfc, s.ea_atm, s.De, s.br_coefs, ...
      s.ro_sfc, s.snow_depth, s.ro_air_Lv, s.liqflag, s.opts);

   testCase.verifyEqual(Ts_phys, s.Tf, 'AbsTol', 1e-12);
   testCase.verifyEqual(Qe, Qe_phys, 'RelTol', 1e-12);
   testCase.verifyNotEqual(Qe, Qe_uncapped);
   testCase.verifyEqual(pevp, Qe / (Lv * s.ro_liq), 'RelTol', 1e-12);
   testCase.verifyEqual(d_pevp, pevp * s.opts.dt / s.dz(1), 'RelTol', 1e-12);
end

function test_thf_debug_dump_reuses_diag_contract(testCase)
   % The THF failure dump should preserve the same diag-style diagnostics
   % returned interactively so failed states can be replayed later.

   workspace = testCase.TestData.workspace;
   s = icemodel.test.fixtures.makeSyntheticColumnState(workspace, ...
      'icemodel', solver=1, seb_solver=2, ...
      turbulent_flux_scheme='monin_obukhov', z0_ice=0.02, ...
      testname='thf_debug_dump');

   debug_file = [tempname, '.mat'];
   old_debug_file = getenv('ICEMODEL_DEBUG_THF_FILE');
   cleanup = onCleanup(@() restoreThfDebugEnv(old_debug_file, debug_file));
   setenv('ICEMODEL_DEBUG_THF_FILE', debug_file);

   icemodel.surface.dump_turbulent_heat_flux_debug_state( ...
      'unit_test', s.Ts, s.Ts - 0.5, s.tair, s.swd, s.lwd, s.albedo, ...
      s.wspd, s.ppt, s.tppt, s.psfc, s.De, s.ea_atm, s.chi, s.ro_air_Lv, s.br_coefs, ...
      s.liqflag, s.T, s.k_eff, s.dz, s.ro_sfc, s.snow_depth, s.opts);

   loaded = load(debug_file, 'debug_state');
   debug_state = loaded.debug_state;

   testCase.verifyEqual(debug_state.reason, 'unit_test');
   testCase.verifyEqual(debug_state.scheme, 'monin_obukhov');
   testCase.verifyEqual(debug_state.thf_at_Ts.diag.scheme, 'monin_obukhov');
   testCase.verifyTrue(isfield(debug_state.seb, 'dfdT_complex_step_at_Ta'));
   testCase.verifyTrue(isfield(debug_state.seb, 'fsearchzero_ok'));
end

function restoreThfDebugEnv(old_debug_file, debug_file)
   %RESTORETHFDEBUGENV Restore the caller's THF debug env var.

   if isempty(old_debug_file)
      setenv('ICEMODEL_DEBUG_THF_FILE', '');
   else
      setenv('ICEMODEL_DEBUG_THF_FILE', old_debug_file);
   end

   if exist(debug_file, 'file') == 2
      delete(debug_file);
   end
end
