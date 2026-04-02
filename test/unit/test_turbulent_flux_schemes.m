function tests = test_turbulent_flux_schemes
   %TEST_TURBULENT_FLUX_SCHEMES Verify the shared THF helper layer.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Build one synthetic workspace that can drive helper-level and
   % operational turbulent-flux tests.

   testCase.TestData.workspace = icemodel.test.fixtures.makeSyntheticWorkspace( ...
      2016, configure=true, nsteps=12, dt_seconds=900);
end

function teardown(testCase)
   % Remove the shared synthetic workspace after the file-level tests end.

   icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      testCase.TestData.workspace);
end

function test_surface_roughness_length_switches_on_snow_depth(testCase)
   % Snow depth, not bulk density alone, should decide snow versus ice.

   z0m_ice = icemodel.surface.surface_roughness_length(0.0, 450.0, 0.02, 1e-3, 5e-4);
   z0m_snow_low = icemodel.surface.surface_roughness_length(0.2, 350.0, 0.02, ...
      1e-3, 5e-4);
   z0m_snow_high = icemodel.surface.surface_roughness_length(0.2, 700.0, 0.02, ...
      1e-3, 5e-4);

   testCase.verifyEqual(z0m_ice, 0.02, 'AbsTol', 1e-12);
   testCase.verifyEqual(z0m_snow_low, 1e-3, 'AbsTol', 1e-12);
   testCase.verifyEqual(z0m_snow_high, 5e-4, 'AbsTol', 1e-12);
end

function test_bulk_richardson_dispatch_matches_legacy_flux_kernels(testCase)
   % The dispatcher should reproduce the legacy sensible/latent formulas
   % when the configured scheme is bulk_richardson.

   s = icemodel.test.fixtures.makeSyntheticColumnState( ...
      testCase.TestData.workspace, 'icemodel', solver=1, ...
      testname='bulk_richardson_dispatch');
   Ts = MELTTEMP(s.Ts, s.Tf);
   S = STABLEFN(s.tair, Ts, s.wspd, s.scoef);
   es_sfc = VAPPRESS(Ts, s.liqflag);
   epsilon = icemodel.physicalConstant('epsilon');
   Qe_ref = LATENT(s.De, S, s.ea, es_sfc, s.roL, epsilon, s.psfc);
   Qh_ref = SENSIBLE(s.De, S, s.tair, Ts, s.cv_air);

   [Qe, Qh, diag] = icemodel.surface.turbulent_heat_flux(s.tair, Ts, s.wspd, s.psfc, ...
      s.De, s.ea, s.cv_air, s.roL, s.scoef, s.liqflag, s.ro_sfc, ...
      s.snow_depth, s.opts);

   testCase.verifyEqual(Qe, Qe_ref, 'RelTol', 1e-12);
   testCase.verifyEqual(Qh, Qh_ref, 'RelTol', 1e-12);
   testCase.verifyEqual(diag.scheme, 'bulk_richardson');
end

function test_bulk_mo_finite_and_weaker_scalar_exchange_for_rough_ice(testCase)
   % Rough bare ice should not force scalar exchange to increase as much as
   % the legacy single-z0 formulation at the same momentum roughness.

   opts_liston = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=1, z0_bulk=0.05, ...
      testname='rough_bulk_richardson');
   opts_vanas = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=1, ...
      seb_solver=2, turbulent_flux_scheme='bulk_mo', z0_bulk=0.05, ...
      z0_ice=0.05, testname='rough_bulk_mo');

   Ta = 268.15;
   Ts = 263.15;
   wspd = 6.0;
   Pa = 90000.0;
   rh = 75.0;
   ea = icemodel.surface.atmospheric_vapor_pressure(Ta, rh, false);
   [cv_air, roLs] = icemodel.physicalConstant('cv_air', 'roLs');
   [De, scoef] = WINDCOEF(wspd, 0.05, opts_liston.z_tair, opts_liston.z_wind);

   [Qe_liston, Qh_liston] = icemodel.surface.turbulent_heat_flux(Ta, Ts, wspd, Pa, ...
      De, ea, cv_air, roLs, scoef, false, 650.0, 0.0, opts_liston);
   [Qe_vanas, Qh_vanas, diag_vanas] = icemodel.surface.turbulent_heat_flux(Ta, Ts, ...
      wspd, Pa, De, ea, cv_air, roLs, scoef, false, 650.0, 0.0, opts_vanas);

   testCase.verifyTrue(all(isfinite([Qe_vanas, Qh_vanas, diag_vanas.L, ...
      diag_vanas.u_star, diag_vanas.z0h, diag_vanas.z0q])));
   testCase.verifyLessThan(abs(Qh_vanas), abs(Qh_liston));
   testCase.verifyLessThan(abs(Qe_vanas), abs(Qe_liston));
end

function test_bulk_mo_respects_surface_phase_switch(testCase)
   % Surface saturation should change when liqflag switches at the same
   % atmospheric state.

   opts_vanas = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=1, ...
      seb_solver=2, turbulent_flux_scheme='bulk_mo', testname='phase_bulk_mo');

   Ta = 273.15;
   Ts = 272.15;
   wspd = 5.0;
   Pa = 90000.0;
   ea = icemodel.surface.atmospheric_vapor_pressure(Ta, 90.0, false);
   [cv_air, roLs] = icemodel.physicalConstant('cv_air', 'roLs');
   [De, scoef] = WINDCOEF(wspd, opts_vanas.z0_bulk, opts_vanas.z_tair, ...
      opts_vanas.z_wind);

   [Qe_ice, ~, diag_ice] = icemodel.surface.turbulent_heat_flux(Ta, Ts, wspd, Pa, ...
      De, ea, cv_air, roLs, scoef, false, 600.0, 0.0, opts_vanas);
   [Qe_liq, ~, diag_liq] = icemodel.surface.turbulent_heat_flux(Ta, Ts, wspd, Pa, ...
      De, ea, cv_air, roLs, scoef, true, 600.0, 0.0, opts_vanas);

   testCase.verifyNotEqual(diag_ice.es_sfc, diag_liq.es_sfc);
   testCase.verifyNotEqual(Qe_ice, Qe_liq);
end

function test_bulk_mo_cold_state_remains_continuous(testCase)
   % The cold stable bulk-MO replay state that used to trip the 2015
   % spinup run should vary smoothly across small Ts perturbations.

   opts_vanas = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=1, ...
      seb_solver=2, turbulent_flux_scheme='bulk_mo', ...
      testname='cold_bulk_mo_continuity');

   Ta = 260.755;
   wspd = 1.41;
   Pa = 87416.0;
   ea = 150.425;
   ro_sfc = 135.536437794510;
   snow_depth = 0.0;
   liqflag = false;
   roL = 3.664362e6;
   scoef = 6.170329907002;
   cv_air = icemodel.physicalConstant('cv_air');

   [Qe1, Qh1, diag1] = icemodel.surface.turbulent_heat_flux(Ta, 245.2210, ...
      wspd, Pa, 0.0, ea, cv_air, roL, scoef, liqflag, ro_sfc, ...
      snow_depth, opts_vanas);
   [Qe2, Qh2, diag2] = icemodel.surface.turbulent_heat_flux(Ta, 245.2220, ...
      wspd, Pa, 0.0, ea, cv_air, roL, scoef, liqflag, ro_sfc, ...
      snow_depth, opts_vanas);

   testCase.verifyTrue(all(isfinite([Qe1, Qh1, Qe2, Qh2, diag1.L, diag2.L])));
   testCase.verifyLessThan(abs(Qh2 - Qh1), 1e-2);
   testCase.verifyLessThan(abs(Qe2 - Qe1), 1e-3);
end

function test_bulk_mo_synthetic_icemodel_run_completes(testCase)
   % A synthetic icemodel run with the supported bulk-MO configuration
   % should complete and return finite Tsfc, Qh, and Qe outputs.

   opts = icemodel.test.helpers.buildSyntheticOpts( ...
      testCase.TestData.workspace, 'icemodel', 2016, solver=1, ...
      seb_solver=2, turbulent_flux_scheme='bulk_mo', z0_ice=0.02, ...
      testname='icemodel_bulk_mo_run');

   [ice1, ~, opts_out] = icemodel(opts);

   testCase.verifyEqual(opts_out.turbulent_flux_scheme, 'bulk_mo');
   testCase.verifyTrue(all(isfinite(ice1.Tsfc)));
   testCase.verifyTrue(all(isfinite(ice1.Qe)));
   testCase.verifyTrue(all(isfinite(ice1.Qh)));
end

function test_potential_surface_vapor_tendency_uses_physical_surface_temperature(testCase)
   % Potential vapor tendency should evaluate Qe using the physical
   % diagnosed surface temperature even if the solver-internal Ts exceeds Tf.

   s = icemodel.test.fixtures.makeSyntheticColumnState( ...
      testCase.TestData.workspace, 'icemodel', solver=1, seb_solver=2, ...
      turbulent_flux_scheme='bulk_mo', z0_ice=0.02, ...
      testname='potential_vapor_tendency');
   Lv = icemodel.physicalConstant('Lv');
   Ts_solver = s.Tf + 3.0;

   [d_pevp, pevp, Qe, Ts_phys] = ...
      icemodel.surface.potential_surface_vapor_tendency(Ts_solver, s.Tf, ...
      s.tair, s.wspd, s.psfc, s.De, s.ea, s.cv_air, s.roL, s.scoef, ...
      s.liqflag, s.ro_sfc, s.snow_depth, s.opts, Lv, s.ro_liq, ...
      s.opts.dt, s.dz(1));

   [Qe_phys, ~] = icemodel.surface.turbulent_heat_flux(s.tair, ...
      MELTTEMP(Ts_solver, s.Tf), s.wspd, s.psfc, s.De, s.ea, s.cv_air, ...
      s.roL, s.scoef, s.liqflag, s.ro_sfc, s.snow_depth, s.opts);
   [Qe_uncapped, ~] = icemodel.surface.turbulent_heat_flux(s.tair, ...
      Ts_solver, s.wspd, s.psfc, s.De, s.ea, s.cv_air, s.roL, s.scoef, ...
      s.liqflag, s.ro_sfc, s.snow_depth, s.opts);

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
      turbulent_flux_scheme='bulk_mo', z0_ice=0.02, ...
      testname='thf_debug_dump');

   debug_file = [tempname, '.mat'];
   old_debug_file = getenv('ICEMODEL_DEBUG_THF_FILE');
   cleanup = onCleanup(@() restoreThfDebugEnv(old_debug_file, debug_file)); %#ok<NASGU>
   setenv('ICEMODEL_DEBUG_THF_FILE', debug_file);

   icemodel.surface.dump_turbulent_heat_flux_debug_state( ...
      'unit_test', s.Ts, s.Ts - 0.5, s.tair, s.swd, s.lwd, s.albedo, ...
      s.wspd, s.ppt, s.tppt, s.psfc, s.De, s.ea, s.chi, s.roL, s.scoef, ...
      s.liqflag, s.T, s.k_eff, s.dz, s.ro_sfc, s.snow_depth, s.opts);

   loaded = load(debug_file, 'debug_state');
   debug_state = loaded.debug_state;

   testCase.verifyEqual(debug_state.reason, 'unit_test');
   testCase.verifyEqual(debug_state.scheme, 'bulk_mo');
   testCase.verifyEqual(debug_state.thf_at_Ts.diag.scheme, 'bulk_mo');
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
