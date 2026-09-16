function tests = test_surface_solver_kernels
   %TEST_SURFACE_SOLVER_KERNELS Verify local surface and column solver kernels.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Build paired skinmodel and icemodel synthetic columns so the surface
   % and coupled column solvers can be exercised on matched states.

   % Setup workspace test data
   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=24, dt_seconds=900);
   testCase.TestData.workspace = workspace;

   % Setup skinmodel test data
   testCase.TestData.skin = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'skinmodel', solver=1, testname='skin_kernel');

   % Setup icemodel test data
   testCase.TestData.ice = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'icemodel', solver=3, testname='ice_kernel');
end

function teardown(testCase)
   % Remove the shared synthetic columns after the file-level tests end.

   icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      testCase.TestData.workspace);
end

function test_sfcflux_derivative_matches_finite_difference(testCase)
   % numerical_surface_flux_linearization should return a derivative consistent
   % with a centered finite difference about the same surface state.

   s = testCase.TestData.skin;
   T_sfc = s.T_sfc;

   % Compute SEB residual and numerical derivative for the bulk_richardson
   % scheme.
   [Q_sfc, dQ_sfc_dT_sfc] = ...
      icemodel.surface.numerical_surface_flux( ...
      T_sfc, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, s.hv_atm, ...
      s.br_coefs, s.liqflag, s.chi, s.T_ice, s.k_eff, ...
      s.dz, s.ro_sfc, s.snow_depth, s.opts);

   % Recompute at a small positive perturbation for a centered difference check.
   h = 1e-5;
   Fplus = icemodel.surface.numerical_surface_flux( ...
      T_sfc + h, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, ...
      s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.T_ice, s.k_eff, s.dz, ...
      s.ro_sfc, s.snow_depth, s.opts);

   % Recompute at a small negative perturbation.
   Fminus = icemodel.surface.numerical_surface_flux( ...
      T_sfc - h, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, ...
      s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.T_ice, s.k_eff, s.dz, ...
      s.ro_sfc, s.snow_depth, s.opts);

   % Centered finite difference.
   dQ_sfc_dT_sfc_fd = (Fplus - Fminus) / (2 * h);

   % Verify the derivative is finite and matches the centered difference value.
   testCase.verifyTrue(isfinite(Q_sfc));
   testCase.verifyEqual(dQ_sfc_dT_sfc, dQ_sfc_dT_sfc_fd, 'RelTol', 2e-4);
end

function test_sfctemp_finds_small_surface_residual(testCase)
   % solve_surface_temperature should converge to a T_sfc where
   % numerical_surface_flux_linearization returns a small residual.

   s = testCase.TestData.skin;

   % Solve for T_sfc using the analytical Newton-Raphson solver. Pass s.T_sfc as
   % the initial guess along with the column state (k_eff, T_ice, dz).
   [T_sfc, ok] = icemodel.surface.solve_surface_temperature( ...
      s.T_sfc, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.ea_atm, s.H_h, s.H_e, s.br_coefs, s.liqflag, s.chi, s.T_ice, ...
      s.k_eff, s.dz);

   % Verify convergence using numerical_surface_flux_linearization, which
   % defines the residual that solve_surface_temperature minimizes.
   Q_sfc = icemodel.surface.numerical_surface_flux( ...
      T_sfc, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, s.hv_atm, ...
      s.br_coefs, s.liqflag, s.chi, s.T_ice, s.k_eff, ...
      s.dz, s.ro_sfc, s.snow_depth, s.opts);

   testCase.verifyTrue(ok);
   testCase.verifyLessThan(abs(Q_sfc), 1e-2);
end

function test_sfcflux_includes_precipitation_advection(testCase)
   % The explicit surface residual should include Qa so the analytic and
   % numerical SEB paths use the same forcing terms.

   s = testCase.TestData.skin;
   T_sfc = s.T_sfc;
   ppt = 2e-4;
   tppt = s.tair + 3.0;

   % Compute the SEB residual without Qa.
   Q_sfc_dry = icemodel.surface.numerical_surface_flux( ...
      T_sfc, s.tair, s.swd, s.lwd, s.albedo, s.wspd, 0.0, tppt, s.psfc, ...
      s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, s.hv_atm, ...
      s.br_coefs, s.liqflag, s.chi, s.T_ice, s.k_eff, ...
      s.dz, s.ro_sfc, s.snow_depth, s.opts);

   % Compute the SEB residual with Qa.
   Q_sfc_wet = icemodel.surface.numerical_surface_flux( ...
      T_sfc, s.tair, s.swd, s.lwd, s.albedo, s.wspd, ppt, tppt, s.psfc, ...
      s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, s.hv_atm, ...
      s.br_coefs, s.liqflag, s.chi, s.T_ice, s.k_eff, ...
      s.dz, s.ro_sfc, s.snow_depth, s.opts);

   % Verify the difference equals Qa.
   testCase.verifyEqual(Q_sfc_wet - Q_sfc_dry, ...
      icemodel.surface.advective_heat_flux(ppt, tppt, s.cv_liq), 'RelTol', 1e-12);
end

function test_sebsolve_converges_across_root_finders(testCase)
   % All standalone surface-energy-balance root-finder modes should converge
   % on the same synthetic forcing state.

   s = testCase.TestData.skin;

   % Solve for T_sfc using each seb_solver option and verify they are identical
   for seb_solver = 0:2
      opts_sv = s.opts;
      opts_sv.seb_solver = seb_solver;
      [T_sfc, ok] = icemodel.surface.solve_surface_energy_balance( ...
         s.T_sfc, s.tair, s.swd, s.lwd, s.albedo, s.wspd, ...
         s.ppt, s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, ...
         s.H_h, s.H_e, s.hv_atm, s.br_coefs, ...
         s.liqflag, s.chi, s.T_ice, s.k_eff, s.dz, s.ro_sfc, s.snow_depth, ...
         opts_sv);
      residual = icemodel.surface.surface_energy_balance_residual(T_sfc, ...
         s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
         s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, s.hv_atm, ...
         s.br_coefs, s.liqflag, s.chi, s.T_ice, s.k_eff, ...
         s.dz, s.ro_sfc, ...
         s.snow_depth, s.opts);

      testCase.verifyTrue(ok);
      testCase.verifyTrue(isfinite(T_sfc));
      testCase.verifyTrue(isreal(T_sfc));
      testCase.verifyLessThan(abs(residual), 1.0);
   end
end

function test_sebsolve_debug_dump_names_the_solver_variable(testCase)
   % A negative seb_solver limits the solve to one outer iteration, so the
   % solve fails on the synthetic state. With debug on, the failure dump
   % stores the solver id under seb_solver, the name of the variable.

   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   debug_file = fullfile(fixture.Folder, 'debug_sebsolve.mat');
   testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
      'ICEMODEL_DEBUG_SEBSOLVE_FILE', debug_file));

   s = testCase.TestData.skin;
   opts = s.opts;
   opts.debug = true;
   opts.seb_solver = -1;
   [~, ok] = icemodel.surface.solve_surface_energy_balance( ...
      s.T_sfc, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, ...
      s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.T_ice, s.k_eff, s.dz, ...
      s.ro_sfc, s.snow_depth, opts);

   testCase.verifyFalse(ok);
   loaded = load(debug_file, 'debug_state');
   testCase.verifyEqual(loaded.debug_state.seb_solver, 1);
   testCase.verifyFalse(isfield(loaded.debug_state, 'solver'));
end

function test_sebsolve_monin_obukhov_converges_with_numeric_derivative(testCase)
   % The new bulk-MO scheme should converge through the shared SEB residual
   % path under the supported solver=1, seb_solver=2 contract.

   workspace = testCase.TestData.workspace;
   s = icemodel.test.fixtures.makeSyntheticColumnState(workspace, ...
      'icemodel', solver=1, seb_solver=2, ...
      turbulent_flux_scheme='monin_obukhov', z0_ice=0.02, ...
      testname='ice_kernel_bulk_mo');

   % Solve for T_sfc
   [T_sfc, ok] = icemodel.surface.solve_surface_energy_balance( ...
      s.T_sfc, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, ...
      s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.T_ice, s.k_eff, s.dz, ...
      s.ro_sfc, s.snow_depth, s.opts);

   % Compute the residual
   residual = icemodel.surface.surface_energy_balance_residual( ...
      T_sfc, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, s.hv_atm, ...
      s.br_coefs, s.liqflag, s.chi, s.T_ice, s.k_eff, ...
      s.dz, s.ro_sfc, s.snow_depth, s.opts);

   % Verify valid solution and residual
   testCase.verifyTrue(ok);
   testCase.verifyTrue(isfinite(T_sfc));
   testCase.verifyTrue(isreal(T_sfc));
   testCase.verifyTrue(isfinite(residual));
   testCase.verifyLessThan(abs(residual), 1.0);
end

function test_one_sweep_dirichlet_uses_vapor_free_predictor(testCase)
   % With one coupler sweep, the pre-coupler predictor supplies the only
   % boundary used by the accepted state. The predictor must therefore use
   % the same vapor-free node conductivity as the column solver.
   %
   % The oracle independently replays the one-sweep sequence with a
   % vapor-free predictor. Including the node vapor term changes the trial
   % surface temperature and fails the exact comparison.

   s = testCase.TestData.ice;
   opts = s.opts;

   % Force one coupler sweep and reuse seb_solver as the inner solver mode,
   % matching the oracle replay below so both paths run the same solve.
   settings = s.settings;
   settings.solver = s.seb_solver;
   settings.cpl_maxiter = 1;

   % Count the predictor and inner-solver calls without modifying the kernel.
   % Restore the caller's profiler state because this test owns only the
   % temporary profiling interval below.
   prior = profile('status');
   restore_profiler = onCleanup(@() restoreProfiler(prior));
   profile off
   profile clear
   profile on

   [T_sfc, T_ice, f_ice, f_liq, k_eff, U_vap, ~, cpl_diag] = ...
      icemodel.couplers.solve_surface_column_dirichlet( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, ...
      s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, s.opts.f_res_pore_ice, settings, opts);
   profile off
   profile_data = profile('info');
   profile clear
   clear restore_profiler

   % One checkpoint evaluation plus one refresh after each solved Picard
   % iterate is the expected saturation-density call count.
   profile_names = string({profile_data.FunctionTable.FunctionName});
   saturation_rows = endsWith(profile_names, 'saturation_vapor_density');
   n_saturation_calls = sum( ...
      [profile_data.FunctionTable(saturation_rows).NumCalls]);
   testCase.verifyEqual(n_saturation_calls, cpl_diag.n_iters + 1);

   testCase.verifyTrue(cpl_diag.ok_seb && cpl_diag.ok_ieb && cpl_diag.ok_cpl);
   testCase.verifyTrue(isfinite(T_sfc) && all(isfinite([T_ice; f_ice; f_liq])));

   % The oracle: the same one-sweep sequence, with the predictor
   % conductivity built from the vapor-free form directly.
   k_pred = icemodel.column.bulk_thermal_conductivity( ...
      s.T_ice, s.f_ice, s.f_liq, 0);
   T_sfc_trial = icemodel.surface.solve_surface_energy_balance( ...
      s.T_sfc, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, ...
      s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.T_ice, k_pred, s.dz, ...
      s.ro_sfc, s.snow_depth, opts);
   [T_ice_x, f_ice_x, f_liq_x, k_x, U_vap_x, ~, ~, ~, ~, err_x] = ...
      icemodel.column.solve_column_enthalpy( ...
      T_sfc_trial, s.T_ice, s.f_ice, s.f_liq, 0.0, 0.0, s.Sc, s.Sp, s.dz, ...
      s.delz, s.fn, s.opts.dt, s.opts.f_res_pore_ice, settings);
   ro_sfc_x = icemodel.surface.surface_bulk_density( ...
      f_ice_x(1), f_liq_x(1));
   T_sfc_x = icemodel.surface.solve_surface_energy_balance( ...
      T_sfc_trial, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, ...
      s.hv_atm, s.br_coefs, s.liqflag, s.chi, T_ice_x, k_x, s.dz, ...
      ro_sfc_x, s.snow_depth, opts);

   testCase.verifyEqual(T_sfc, T_sfc_x, 'AbsTol', 0);
   testCase.verifyEqual(T_ice, T_ice_x, 'AbsTol', 0);
   testCase.verifyEqual(f_ice, f_ice_x, 'AbsTol', 0);
   testCase.verifyEqual(f_liq, f_liq_x, 'AbsTol', 0);
   testCase.verifyEqual(k_eff, k_x, 'AbsTol', 0);
   testCase.verifyEqual(U_vap, U_vap_x, 'AbsTol', 0);
   testCase.verifyTrue(isfinite(err_x));

   % Recompute the face flux independently from the accepted solver state.
   [ro_vap_x, dro_vapdT_x] = ...
      icemodel.vapor.saturation_vapor_density(T_ice_x, f_liq_x);
   [~, De_x] = icemodel.vapor.vapor_thermal_conductivity( ...
      T_ice_x, f_liq_x, dro_vapdT_x);
   [~, ~, ~, U_expected] = ...
      icemodel.column.vapor_transport_terms( ...
      T_ice_x, f_ice_x, f_liq_x, k_x, ro_vap_x, dro_vapdT_x, De_x, ...
      s.delz, s.fn, s.opts.f_res_pore_ice);
   testCase.verifyEqual(U_vap, U_expected, 'AbsTol', 0);
end

function test_iceebsolvedirichlet_converges_on_synthetic_column(testCase)
   % The coupled Dirichlet solver should converge on the synthetic ice
   % column, keep phase fractions bounded, and reduce the accepted
   % updated-state SEB residual relative to the current one-pass pattern.

   s = testCase.TestData.ice;

   % Solve the coupled surface-column energy balance. Reuse seb_solver as
   % the inner solver mode, matching the one-pass comparison built below.
   settings = s.settings;
   settings.solver = s.seb_solver;
   [T_sfc, T_ice, f_ice, f_liq, k_eff, ~, ~, cpl_diag] = ...
      icemodel.couplers.solve_surface_column_dirichlet( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, ...
      s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, s.opts.f_res_pore_ice, settings, s.opts);

   % Compute the coupled surface-column energy balance residual
   residual_coupled = icemodel.surface.surface_energy_balance_residual( ...
      T_sfc, ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, s.hv_atm, ...
      s.br_coefs, s.liqflag, s.chi, T_ice, k_eff, s.dz, ...
      s.ro_sfc, s.snow_depth, s.opts);

   % Solve for T_sfc using the surface solver to mimic a one-pass solve
   [T_sfc_old, ok_seb_old] = icemodel.surface.solve_surface_energy_balance( ...
      s.T_sfc, s.tair, s.swd, s.lwd, s.albedo, ...
      s.wspd, s.ppt, s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, ...
      s.nu_air, s.H_h, s.H_e, s.hv_atm, s.br_coefs, ...
      s.liqflag, s.chi, s.T_ice, s.k_eff, s.dz, s.ro_sfc, s.snow_depth, s.opts);

   % Solve for column T_ice using the subsurface column solver
   one_pass_settings = s.settings;
   one_pass_settings.solver = 1;
   one_pass_settings.debug = false;
   [T_ice_old, ~, ~, k_eff_old, ~, ~, ok_ieb_old] = ...
      icemodel.column.solve_column_enthalpy(T_sfc_old, s.T_ice, s.f_ice, ...
      s.f_liq, s.Fc, s.Fp, s.Sc, s.Sp, s.dz, s.delz, s.fn, s.opts.dt, ...
      s.opts.f_res_pore_ice, one_pass_settings);

   % Compute the one-pass residual
   residual_old = icemodel.surface.surface_energy_balance_residual( ...
      T_sfc_old, ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, s.hv_atm, ...
      s.br_coefs, s.liqflag, s.chi, T_ice_old, k_eff_old, s.dz, ...
      s.ro_sfc, s.snow_depth, s.opts);

   % Verify the coupled solutions are valid
   testCase.verifyTrue(cpl_diag.ok_seb);
   testCase.verifyTrue(cpl_diag.ok_ieb);
   testCase.verifyTrue(cpl_diag.ok_cpl);
   testCase.verifyTrue(isfinite(cpl_diag.n_iters));
   testCase.verifyTrue(isreal(T_sfc));
   testCase.verifyTrue(all(isfinite([T_sfc; T_ice; f_ice; f_liq; k_eff])));
   testCase.verifyLessThanOrEqual(max(f_ice + f_liq * s.ro_liq / s.ro_ice), ...
      1 + 1e-9);
   testCase.verifyLessThan(abs(residual_coupled), 1.0);

   % Verify the one-pass solutions are valid
   testCase.verifyTrue(ok_seb_old);
   testCase.verifyTrue(ok_ieb_old);

   % Verify the coupled solution residual is less than the one-pass residual
   testCase.verifyLessThanOrEqual(abs(residual_coupled), ...
      abs(residual_old) + 1e-6);
end

function test_surface_flux_diagnostics_match_energy_balance_residual(testCase)
   % The namespaced surface contracts should report the same residual as the
   % explicit ENBAL helper built from their returned component fluxes.

   s = testCase.TestData.skin;

   % Compute atmospheric vapor pressure
   rh = s.met.rh(s.metstep);
   ea_atm = icemodel.surface.atmospheric_vapor_pressure(s.tair, rh, s.liqflag);

   % Solve for T_sfc
   [T_sfc_raw, ok] = icemodel.surface.solve_surface_energy_balance(s.T_sfc, s.tair, ...
      s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ea_atm, ...
      s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, s.hv_atm, ...
      s.br_coefs, s.liqflag, s.chi, s.T_ice, s.k_eff, s.dz, s.ro_sfc, ...
      s.snow_depth, s.opts);

   % Cap T_sfc at the melt temp
   T_sfc = icemodel.surface.physical_surface_temperature(T_sfc_raw);

   % Diagnose the full surface energy balance term set.
   [Qe, Qh, Qc, Qsn, Qln, Qa, Qm, ~, balance] = ...
      icemodel.surface.diagnose_surface_energy_balance(T_sfc, s.tair, s.swd, ...
      s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ea_atm, ...
      s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, s.hv_atm, ...
      s.br_coefs, s.liqflag, s.chi, s.T_ice, ...
      s.k_eff, s.dz, s.ro_sfc, s.snow_depth, ...
      s.opts);

   % Verify the T_sfc solver worked
   testCase.verifyTrue(ok);

   % Verify the ea_atm calculation matches the test case setup
   testCase.verifyEqual(ea_atm, s.ea_atm, 'RelTol', 1e-12);

   % Verify the balance returned by diagnose_surface_energy_balance matches
   % the explicit SEB evaluation built from its diagnosed terms.
   testCase.verifyEqual(balance, icemodel.surface.evaluate_surface_energy_balance( ...
      Qsn, Qln, Qh, Qe, Qc, Qa, Qm), 'RelTol', 1e-12);

   % Verify the T_sfc solution matches the test case setup value.
   testCase.verifyLessThanOrEqual(T_sfc, s.Tf + 1e-12);
end

function test_skinsolve_returns_finite_bounded_state(testCase)
   % icemodel.column.solve_column_temperature should keep the skin column
   % finite and phase-bounded on the shared synthetic state.

   s = testCase.TestData.skin;

   % Solve for ice column temperature using the reduced-complexity solver
   settings = s.settings;
   settings.debug = false;
   [T_ice, f_ice, f_liq, k_eff, ok, iter] = ...
      icemodel.column.solve_column_temperature(s.T_sfc, s.T_ice, s.f_ice, ...
      s.f_liq, s.dz, s.delz, s.fn, s.opts.dt, settings);

   % Verify the solution is valid
   testCase.verifyTrue(ok);
   testCase.verifyGreaterThan(iter, 0);
   testCase.verifyTrue(all(isfinite(T_ice)));
   testCase.verifyTrue(all(isfinite(k_eff)));
   testCase.verifyLessThanOrEqual(max(f_ice + f_liq * s.ro_liq / s.ro_ice), ...
      1 + 1e-9);
end

function test_skinebsolve_converges_on_synthetic_column(testCase)
   % The coupled skin energy-balance solve should converge and keep the
   % phase fractions inside their physical bounds.

   s = testCase.TestData.skin;

   % Solve for T_sfc and T_ice using the coupled reduced-complexity solver
   [T_sfc, T_ice, f_ice, f_liq, k_eff, diag] = ...
      icemodel.couplers.solve_skin_surface_column( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, ...
      s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, s.settings, s.opts);

   % Verify the solution is valid
   testCase.verifyTrue(diag.ok_seb);
   testCase.verifyTrue(diag.ok_ieb);
   testCase.verifyTrue(diag.ok_cpl);
   testCase.verifyGreaterThan(diag.n_iters, 0);
   testCase.verifyEqual(sum(isfinite(diag.cpl_res_hist)), ...
      min(diag.cpl_iters, numel(diag.cpl_res_hist)));
   testCase.verifyEqual(diag.cpl_res_hist(end), diag.cpl_res, 'AbsTol', 0);
   testCase.verifyTrue(isreal(T_sfc));
   testCase.verifyTrue(all(isfinite([T_sfc; T_ice; k_eff])));
   testCase.verifyLessThanOrEqual(max(f_ice + f_liq * s.ro_liq / s.ro_ice), ...
      1 + 1e-9);
end

function test_iceenbal_and_iceebsolve_converge_on_synthetic_column(testCase)
   % The direct and coupled icemodel column solves should both converge on
   % the shared synthetic state and preserve phase bounds.

   s = testCase.TestData.ice;

   % Solve the ice column model using a direct call to
   % icemodel.column.solve_column_enthalpy (single sweep)
   solver = 1;
   settings = s.settings;
   settings.solver = solver;
   settings.debug = false;
   [T_ice_dir, f_ice_dir, f_liq_dir, k_eff_dir, ~, ~, ok_dir, iter_dir] = ...
      icemodel.column.solve_column_enthalpy(s.T_sfc, s.T_ice, s.f_ice, ...
      s.f_liq, s.Fc, s.Fp, s.Sc, s.Sp, s.dz, s.delz, s.fn, s.opts.dt, ...
      s.opts.f_res_pore_ice, settings);

   % Solve the ice column model using the fully coupled solver with robin bc.
   % s.settings.solver is already 3 (the ice fixture's build-time opts).
   [T_sfc, T_ice_rob, f_ice_rob, f_liq_rob, k_eff_rob, ~, ~, diag_rob] = ...
      icemodel.couplers.solve_surface_column_robin( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
      s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, ...
      s.H_e, s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, s.opts.f_res_pore_ice, s.settings, s.opts);

   % Verify the single-sweep direct solve is valid
   testCase.verifyTrue(ok_dir);
   testCase.verifyGreaterThanOrEqual(iter_dir, 0);
   testCase.verifyTrue(all(isfinite([T_ice_dir; k_eff_dir])));
   testCase.verifyLessThanOrEqual(max(f_ice_dir + f_liq_dir * ...
      s.ro_liq / s.ro_ice), 1 + 1e-9);

   % Verify the fully coupled robin bc solve is valid
   testCase.verifyTrue(diag_rob.ok_seb && diag_rob.ok_ieb && diag_rob.ok_cpl);
   testCase.verifyGreaterThan(diag_rob.n_iters, 0);
   testCase.verifyTrue(isfinite(T_sfc));
   testCase.verifyTrue(isreal(T_sfc));
   testCase.verifyTrue(all(isfinite([T_ice_rob; k_eff_rob])));
   testCase.verifyLessThanOrEqual(max(f_ice_rob + f_liq_rob * ...
      s.ro_liq / s.ro_ice), 1 + 1e-9);
end

function test_robin_adapter_returns_inner_accepted_vapor_flux(testCase)
   % A one-sweep Robin adapter must return the accepted inner solve's vapor
   % flux through its sixth output when solver is greater than one.

   s = testCase.TestData.ice;
   solver = 3;

   % Reproduce the coefficients that the adapter sends to its only inner
   % solve. The direct solve and adapter then start from the same state.
   [Fc, Fp] = icemodel.surface.surface_flux_linearization( ...
      s.T_sfc, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, ...
      s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, s.opts);
   inner_settings = s.settings;
   inner_settings.solver = solver;
   inner_settings.debug = s.opts.debug;
   [T_ice_inner, f_ice_inner, f_liq_inner, k_eff_inner, U_vap_inner, ~, ...
      ok_inner, ~, ~, ~] = icemodel.column.solve_column_enthalpy( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, Fc, Fp, s.Sc, s.Sp, s.dz, s.delz, ...
      s.fn, s.opts.dt, s.opts.f_res_pore_ice, inner_settings);

   settings = s.settings;
   settings.cpl_maxiter = 1;
   [~, T_ice_adapter, f_ice_adapter, f_liq_adapter, k_eff_adapter, ...
      U_vap_adapter, ~, diag_adapter] = ...
      icemodel.couplers.solve_surface_column_robin( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
      s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, ...
      s.H_e, s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, s.opts.f_res_pore_ice, settings, s.opts);

   testCase.verifyTrue(ok_inner && diag_adapter.ok_seb && ...
      diag_adapter.ok_ieb && diag_adapter.ok_cpl);
   testCase.verifyEqual(T_ice_adapter, T_ice_inner, 'AbsTol', 0);
   testCase.verifyEqual(f_ice_adapter, f_ice_inner, 'AbsTol', 0);
   testCase.verifyEqual(f_liq_adapter, f_liq_inner, 'AbsTol', 0);
   testCase.verifyEqual(k_eff_adapter, k_eff_inner, 'AbsTol', 0);
   testCase.verifyEqual(U_vap_adapter, U_vap_inner, 'AbsTol', 0);
end

function test_one_cell_column_preserves_solver_and_coupler_outputs(testCase)
   % The tenth solver output must stay a finite residual and the fifth must
   % stay the two-face vapor flux when a closed column has only one cell.

   s = testCase.TestData.ice;
   s.T_ice = s.T_ice(1);
   s.f_ice = s.f_ice(1);
   s.f_liq = s.f_liq(1);
   s.Sc = s.Sc(1);
   s.Sp = s.Sp(1);
   s.dz = s.dz(1);
   s.delz = 0.5 * s.dz * ones(2, 1);
   s.fn = [0; 1];

   % Exercise the direct output ordering that exposed the one-cell indexing.
   settings = s.settings;
   settings.solver = 1;
   settings.debug = s.opts.debug;
   [T_ice, f_ice, f_liq, k_eff, U_vap, ~, ok, ~, a1, err] = ...
      icemodel.column.solve_column_enthalpy( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, 0.0, 0.0, s.Sc, s.Sp, s.dz, ...
      s.delz, s.fn, s.opts.dt, s.opts.f_res_pore_ice, settings);
   testCase.verifyTrue(ok);
   testCase.verifyTrue(all(isfinite([T_ice; f_ice; f_liq; k_eff; a1; err])));
   testCase.verifyEqual(U_vap, zeros(2, 1), 'AbsTol', 0);

   % Both column couplers request output five and must therefore retain
   % the same one-cell behavior through the complete surface-column call.
   settings_dir = s.settings;
   settings_dir.solver = 1;
   settings_dir.cpl_maxiter = 1;
   [T_sfc_dir, T_ice_dir, ~, ~, ~, U_dir, ~, diag_dir] = ...
      icemodel.couplers.solve_surface_column_dirichlet( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
      s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, ...
      s.H_e, s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, s.opts.f_res_pore_ice, settings_dir, s.opts);
   testCase.verifyTrue(diag_dir.ok_seb && diag_dir.ok_ieb && diag_dir.ok_cpl);
   testCase.verifyTrue(all(isfinite([T_sfc_dir; T_ice_dir])));
   testCase.verifyEqual(U_dir, zeros(2, 1), 'AbsTol', 0);

   settings_rob = s.settings;
   settings_rob.solver = 3;
   settings_rob.cpl_maxiter = 1;
   [T_sfc_rob, T_ice_rob, ~, ~, ~, U_rob, ~, diag_rob] = ...
      icemodel.couplers.solve_surface_column_robin( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
      s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, ...
      s.H_e, s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, s.opts.f_res_pore_ice, settings_rob, s.opts);
   testCase.verifyTrue(diag_rob.ok_seb && diag_rob.ok_ieb && diag_rob.ok_cpl);
   testCase.verifyTrue(all(isfinite([T_sfc_rob; T_ice_rob])));
   testCase.verifyEqual(U_rob, zeros(2, 1), 'AbsTol', 0);

   % Both couplers return the same fields as the raw substep attempt.
   [~, step_diag] = icemodel.couplers.initialize_solver_diag();
   expected = fieldnames(step_diag);
   testCase.verifyEqual(fieldnames(diag_dir), expected);
   testCase.verifyEqual(fieldnames(diag_rob), expected);
   testCase.verifyEqual(size(diag_dir.cpl_res_hist), [16 1]);
   testCase.verifyEqual(size(diag_rob.cpl_res_hist), [16 1]);
end

function test_iceenbal_routes_exhausted_budget_to_substep_recovery(testCase)
   % A solve that uses its entire inner iteration budget must let the caller
   % restore the checkpoint and reduce dt instead of accepting that iterate.
   % With debug on, the failure dump stores the entering temperature under
   % T_ice_old, the name of the variable.

   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   debug_file = fullfile(fixture.Folder, 'debug_iceenbal.mat');
   testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
      'ICEMODEL_DEBUG_ICEENBAL_FILE', debug_file));

   s = testCase.TestData.ice;
   settings = s.settings;
   settings.solver = 1;
   settings.maxiter = 1;
   settings.debug = true;
   [~, ~, ~, ~, ~, ~, ok, iter] = icemodel.column.solve_column_enthalpy( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.Fc, s.Fp, s.Sc, s.Sp, s.dz, ...
      s.delz, s.fn, s.opts.dt, s.opts.f_res_pore_ice, settings);

   testCase.verifyFalse(ok);
   testCase.verifyEqual(iter, 1);
   loaded = load(debug_file, 'debug_state');
   testCase.verifyEqual(loaded.debug_state.reason, "maxiter_nonconvergence");
   testCase.verifyEqual(loaded.debug_state.T_ice_old, s.T_ice);
   testCase.verifyFalse(isfield(loaded.debug_state, 'T_old'));
end

function test_secantscalar_finds_root_after_acceleration(testCase)
   % Residual pairs remain valid even when evaluated states are not a
   % consecutive unaccelerated Picard sequence.

   [x_next, ok] = icemodel.numerics.secantscalar( ...
      1.25, 0.375, 2.25, -0.125, 2.125, 10, true);

   testCase.verifyTrue(ok);
   testCase.verifyEqual(x_next, 2.0, 'AbsTol', 10 * eps);
end

function test_secantscalar_uses_safeguarded_fallbacks(testCase)
   % Each invalid or unsafe secant state must retain the caller's relaxed
   % fallback instead of emitting an unusable accelerated temperature.

   cases = [ ...
      struct('xprev', 0.25, 'rprev', 0.875, 'x', 1.25, 'r', 0.375, ...
      'fallback', 1.625, 'jumpmax', 10, 'enabled', false); ...
      struct('xprev', nan, 'rprev', 0.875, 'x', 1.25, 'r', 0.375, ...
      'fallback', 1.625, 'jumpmax', 10, 'enabled', true); ...
      struct('xprev', 0.25, 'rprev', nan, 'x', 1.25, 'r', 0.375, ...
      'fallback', 1.625, 'jumpmax', 10, 'enabled', true); ...
      struct('xprev', 0.25, 'rprev', 0.875, 'x', nan, 'r', 0.375, ...
      'fallback', 1.625, 'jumpmax', 10, 'enabled', true); ...
      struct('xprev', 0.25, 'rprev', 0.875, 'x', 1.25, 'r', nan, ...
      'fallback', 1.625, 'jumpmax', 10, 'enabled', true); ...
      struct('xprev', 0.25, 'rprev', 0, 'x', 1.25, 'r', -0.375, ...
      'fallback', 1.625, 'jumpmax', 10, 'enabled', true); ...
      struct('xprev', 0.25, 'rprev', 0.375, 'x', 1.25, 'r', 0, ...
      'fallback', 1.625, 'jumpmax', 10, 'enabled', true); ...
      struct('xprev', 0.25, 'rprev', 0.375, 'x', 1.25, 'r', 0.375, ...
      'fallback', 1.625, 'jumpmax', 10, 'enabled', true); ...
      struct('xprev', 1, 'rprev', 1, 'x', 1 + eps, 'r', -1, ...
      'fallback', 3, 'jumpmax', 10, 'enabled', true)];

   for n = 1:numel(cases)
      c = cases(n);
      [x_next, ok] = icemodel.numerics.secantscalar(c.xprev, c.rprev, ...
         c.x, c.r, c.fallback, c.jumpmax, c.enabled);
      testCase.verifyFalse(ok, sprintf('case %d unexpectedly accelerated', n));
      testCase.verifyEqual(x_next, c.fallback, ...
         sprintf('case %d did not retain fallback', n));
   end
end

function test_secantscalar_keeps_bracket_when_fallback_is_outside(testCase)
   % Every endpoint orientation, exterior fallback side, and jump regime must
   % keep a sign-bracketed step inside the bracket and independent of fallback.

   endpoint_cases = [ ...
      struct('xprev', 0, 'rprev', 1, 'x', 2, 'r', -1); ...
      struct('xprev', 2, 'rprev', -1, 'x', 0, 'r', 1)];
   fallbacks = [-100, 100];
   jumpmax_cases = [realmax, 0.25];

   for n_endpoint = 1:numel(endpoint_cases)
      c = endpoint_cases(n_endpoint);
      for fallback = fallbacks
         for jumpmax = jumpmax_cases
            [x_next, ok] = icemodel.numerics.secantscalar( ...
               c.xprev, c.rprev, c.x, c.r, fallback, jumpmax, true);

            expected = c.x + sign(1 - c.x) * min(abs(1 - c.x), jumpmax);
            testCase.verifyTrue(ok);
            testCase.verifyEqual(x_next, expected, 'AbsTol', 10 * eps);
            testCase.verifyGreaterThan(x_next, 0);
            testCase.verifyLessThan(x_next, 2);
         end
      end
   end
end

function test_secantscalar_keeps_unsafe_interpolation_inside_bracket(testCase)
   % Ill-conditioned, overflowing, and endpoint-hugging interpolation should
   % take the finite interior midpoint rather than leave the root bracket.

   cases = [ ...
      struct('xprev', 0, 'rprev', realmin, 'x', 2, 'r', -realmin, ...
      'fallback', 1, 'expected', 1); ...
      struct('xprev', -realmax, 'rprev', -1, 'x', realmax, 'r', 1, ...
      'fallback', 0, 'expected', 0); ...
      struct('xprev', 0, 'rprev', 1, 'x', 2, 'r', -eps, ...
      'fallback', 1, 'expected', 1)];

   for n = 1:numel(cases)
      c = cases(n);
      [x_next, ok] = icemodel.numerics.secantscalar(c.xprev, c.rprev, ...
         c.x, c.r, c.fallback, realmax, true);
      testCase.verifyTrue(ok, sprintf('case %d rejected its bracket', n));
      testCase.verifyEqual(x_next, c.expected, ...
         sprintf('case %d did not use its bracket midpoint', n));
   end
end

function test_robin_debug_dump_handles_inner_failure(testCase)
   % The Robin debug path must preserve the failed solver flags instead of
   % masking the original inner-solve failure with an argument-count error.

   s = testCase.TestData.ice;
   debug_file = [tempname, '.mat'];
   old_debug_file = getenv('ICEMODEL_DEBUG_ICEEBSOLVE_FILE');
   cleanup = onCleanup(@() restoreIceEbDebugEnv( ...
      old_debug_file, debug_file));
   setenv('ICEMODEL_DEBUG_ICEEBSOLVE_FILE', debug_file)
   opts = s.opts;
   opts.debug = true;

   % A one-iteration inner budget deterministically exercises the failure
   % dump. The coupler reads debug from settings.debug, so the debug flag
   % needs the same override as opts.debug.
   settings = s.settings;
   settings.maxiter = 1;
   settings.debug = true;
   [~, ~, ~, ~, ~, ~, ~, diag] = ...
      icemodel.couplers.solve_surface_column_robin( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
      s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, ...
      s.H_e, s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, s.opts.f_res_pore_ice, settings, opts);

   % The returned flags and saved diagnostic must identify the same
   % failure. Dumps are sequence-numbered by a persistent counter that
   % survives across calls in one MATLAB session, so this single-event
   % probe discovers the one numbered file under its own fresh base path
   % instead of assuming the sequence starts at one.
   testCase.verifyTrue(diag.ok_seb);
   testCase.verifyFalse(diag.ok_ieb);
   testCase.verifyFalse(diag.ok_cpl);
   [dump_dir, dump_name, dump_ext] = fileparts(debug_file);
   dump_listing = dir(fullfile(dump_dir, ...
      sprintf('%s_*%s', dump_name, dump_ext)));
   testCase.assertNumElements(dump_listing, 1);
   dump_file = fullfile(dump_dir, dump_listing(1).name);
   loaded = load(dump_file, 'debug_state');
   testCase.verifyEqual(loaded.debug_state.reason, "iceenbal_failed");
   testCase.verifyFalse(loaded.debug_state.ok_ieb);
   testCase.verifyFalse(loaded.debug_state.ok_cpl);
   testCase.verifyEqual(size(loaded.debug_state.res_hist), [16 1]);
end

function test_dirichlet_debug_dump_keeps_coupling_residual_history(testCase)
   % A coupling failure dump must contain each evaluated residual.

   s = testCase.TestData.ice;
   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   debug_file = fullfile(fixture.Folder, 'ice-debug.mat');
   old_debug_file = getenv('ICEMODEL_DEBUG_ICEEBSOLVE_FILE');
   cleanup = onCleanup(@() restoreIceEbDebugEnv( ...
      old_debug_file, debug_file));
   setenv('ICEMODEL_DEBUG_ICEEBSOLVE_FILE', debug_file)

   % Negative tolerances force failure after two evaluated residuals.
   settings = s.settings;
   settings.solver = 1;
   settings.debug = true;
   settings.cpl_maxiter = 2;
   settings.cpl_Ts_tol = -1;
   settings.cpl_seb_tol = -1;
   opts = s.opts;
   opts.debug = true;
   [~, ~, ~, ~, ~, ~, ~, diag] = ...
      icemodel.couplers.solve_surface_column_dirichlet( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
      s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, ...
      s.H_e, s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, s.opts.f_res_pore_ice, settings, opts);

   dump_listing = dir(fullfile(fixture.Folder, 'ice-debug_dirichlet_*.mat'));
   testCase.assertNumElements(dump_listing, 1);
   loaded = load(fullfile(fixture.Folder, dump_listing(1).name), ...
      'debug_state');
   testCase.verifyEqual(loaded.debug_state.reason, "coupler_nonconvergence");
   testCase.verifyEqual(loaded.debug_state.res_hist, diag.cpl_res_hist);
   testCase.verifyEqual(sum(isfinite(diag.cpl_res_hist)), 2);
end

function test_robin_coupler_supports_monin_obukhov_on_synthetic_column(testCase)
   % The Robin coupler should accept the bulk-MO scheme and converge on the
   % synthetic ice state used by the shared solver tests.

   workspace = testCase.TestData.workspace;
   s = icemodel.test.fixtures.makeSyntheticColumnState(workspace, ...
      'icemodel', solver=3, seb_solver=2, turbulent_flux_scheme='monin_obukhov', ...
      z0_ice=0.02, testname='ice_kernel_robin_bulk_mo');

   % Solve the ice column model using the fully coupled solver with robin bc
   [T_sfc, T_ice_rob, f_ice_rob, f_liq_rob, k_eff_rob, ~, ~, diag_rob] = ...
      icemodel.couplers.solve_surface_column_robin( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.Sc, s.Sp, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
      s.tppt, s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, ...
      s.H_e, s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, s.opts.f_res_pore_ice, s.settings, s.opts);

   % Compute the residual
   residual = icemodel.surface.surface_energy_balance_residual(T_sfc, ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, s.hv_atm, ...
      s.br_coefs, s.liqflag, s.chi, T_ice_rob, k_eff_rob, ...
      s.dz, s.ro_sfc, s.snow_depth, s.opts);

   % Verify the solution is valid
   testCase.verifyTrue(diag_rob.ok_seb && diag_rob.ok_ieb && diag_rob.ok_cpl);
   testCase.verifyGreaterThan(diag_rob.n_iters, 0);
   testCase.verifyTrue(isreal(T_sfc));
   testCase.verifyTrue(all(isfinite([T_sfc; T_ice_rob; k_eff_rob; residual])));
   testCase.verifyLessThanOrEqual(max(f_ice_rob + f_liq_rob * ...
      s.ro_liq / s.ro_ice), 1 + 1e-9);
   testCase.verifyLessThan(abs(residual), 1.0);
end

function test_skin_debug_dump_keeps_coupling_residual_history(testCase)
   % A coupling failure dump must contain the evaluated residual history.

   s = testCase.TestData.skin;
   % Negative tolerances force an outer-coupling failure after two evaluated
   % residuals while leaving both inner solves unchanged.
   settings = s.settings;
   settings.debug = true;
   settings.cpl_maxiter = 2;
   settings.cpl_Ts_tol = -1;
   settings.cpl_seb_tol = -1;
   opts = s.opts;
   opts.debug = true;
   [diag, debug_state] = runSkinDebugFailure(testCase, s, settings, opts);

   testCase.verifyEqual(debug_state.reason, "coupler_nonconvergence");
   testCase.verifyEqual(debug_state.res_hist, diag.cpl_res_hist);
end

function test_skin_debug_dump_handles_column_failure(testCase)
   % A column-solve failure occurs before any coupling residual is evaluated.

   s = testCase.TestData.skin;
   settings = s.settings;
   settings.debug = true;
   settings.maxiter = 1;
   opts = s.opts;
   opts.debug = true;
   [diag, debug_state] = runSkinDebugFailure(testCase, s, settings, opts);

   testCase.verifyEqual(debug_state.reason, "skinsolve_failed");
   testCase.verifyEqual(debug_state.res_hist, diag.cpl_res_hist);
   testCase.verifyTrue(all(isnan(debug_state.res_hist)));
end

function test_skin_debug_dump_handles_surface_failure(testCase)
   % A surface-solve failure occurs before its residual can be recorded.

   s = testCase.TestData.skin;
   settings = s.settings;
   settings.debug = true;
   opts = s.opts;
   opts.debug = true;
   opts.seb_solver = -1;
   [diag, debug_state] = runSkinDebugFailure(testCase, s, settings, opts);

   testCase.verifyEqual(debug_state.reason, "sebsolve_failed");
   testCase.verifyEqual(debug_state.res_hist, diag.cpl_res_hist);
   testCase.verifyTrue(all(isnan(debug_state.res_hist)));
end

function [diag, debug_state] = runSkinDebugFailure(testCase, s, settings, opts)
   %RUNSKINDEBUGFAILURE Run one skin-coupler failure with an isolated dump.

   fixture = testCase.applyFixture( ...
      matlab.unittest.fixtures.TemporaryFolderFixture);
   debug_file = fullfile(fixture.Folder, 'skin-debug.mat');
   old_debug_file = getenv('ICEMODEL_DEBUG_SKINEBSOLVE_FILE');
   cleanup = onCleanup(@() restoreSkinEbDebugEnv(old_debug_file));
   setenv('ICEMODEL_DEBUG_SKINEBSOLVE_FILE', debug_file)

   [~, ~, ~, ~, ~, diag] = ...
      icemodel.couplers.solve_skin_surface_column( ...
      s.T_sfc, s.T_ice, s.f_ice, s.f_liq, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.ea_atm, s.ro_atm, s.cv_atm, s.nu_air, s.H_h, s.H_e, ...
      s.hv_atm, s.br_coefs, s.liqflag, s.chi, s.ro_sfc, ...
      s.snow_depth, settings, opts);

   dump_listing = dir(fullfile(fixture.Folder, 'skin-debug_*.mat'));
   testCase.assertNumElements(dump_listing, 1);
   loaded = load(fullfile(fixture.Folder, dump_listing(1).name), 'debug_state');
   debug_state = loaded.debug_state;
end

function restoreIceEbDebugEnv(old_debug_file, debug_file)
   %RESTOREICEEBDEBUGENV Restore the debug target and remove the test artifact.
   setenv('ICEMODEL_DEBUG_ICEEBSOLVE_FILE', old_debug_file)
   % The dump writer appends sequence numbers before the extension, so
   % remove every numbered file the probe produced.
   [dump_dir, dump_name, dump_ext] = fileparts(debug_file);
   numbered = dir(fullfile(dump_dir, [dump_name, '_*', dump_ext]));
   for k = 1:numel(numbered)
      delete(fullfile(numbered(k).folder, numbered(k).name))
   end
   if exist(debug_file, 'file') == 2
      delete(debug_file)
   end
end

function restoreSkinEbDebugEnv(old_debug_file)
   %RESTORESKINEBDEBUGENV Restore the debug target.
   setenv('ICEMODEL_DEBUG_SKINEBSOLVE_FILE', old_debug_file)
end

function restoreProfiler(prior)
   %RESTOREPROFILER Put the profiler back in the caller's on/off state.
   %
   % PROFILE('info') requires the temporary data collected by this test, so
   % only the caller's active state can be restored after that data is cleared.

   profile off
   profile clear
   if strcmp(prior.ProfilerStatus, 'on')
      profile on
   end
end
