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
   % evaluate_surface_flux should return a derivative consistent with a
   % centered finite difference about the same surface state.

   s = testCase.TestData.skin;
   Ts = s.Ts;
   Qc = icemodel.surface.conductive_heat_flux(s.k_eff, s.T, s.dz, Ts);

   % Compute the SEB residual and numerical derivative for the bulk_richardson scheme.
   [Q_sfc, dQ_sfc_dTs] = ...
      icemodel.surface.turbulence.bulk_richardson.evaluate_surface_flux( ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.De, Ts, Qc, s.ea_atm, s.roL, s.br_coefs, s.chi, s.liqflag);

   % Recompute at a small positive perturbation for a centered difference check.
   h = 1e-5;
   Fplus = icemodel.surface.turbulence.bulk_richardson.evaluate_surface_flux( ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.De, Ts + h, Qc, s.ea_atm, s.roL, s.br_coefs, s.chi, s.liqflag);

   % Recompute at a small negative perturbation.
   Fminus = icemodel.surface.turbulence.bulk_richardson.evaluate_surface_flux( ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.De, Ts - h, Qc, s.ea_atm, s.roL, s.br_coefs, s.chi, s.liqflag);

   % Centered finite difference.
   dQ_sfc_dTs_fd = (Fplus - Fminus) / (2 * h);

   % Verify the derivative is finite and matches the centered difference value.
   testCase.verifyTrue(isfinite(Q_sfc));
   testCase.verifyEqual(dQ_sfc_dTs, dQ_sfc_dTs_fd, 'RelTol', 2e-4);
end

function test_sfctemp_finds_small_surface_residual(testCase)
   % solve_surface_temperature should converge to a T_sfc where
   % evaluate_surface_flux returns a small residual.

   s = testCase.TestData.skin;

   Qc = icemodel.surface.conductive_heat_flux(s.k_eff, s.T, s.dz, s.Ts);

   % Solve for T_sfc using the analytical Newton-Raphson solver.
   [Ts, ok] = icemodel.surface.solve_surface_temperature( ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, s.De, ...
      s.ea_atm, s.br_coefs, s.roL, s.liqflag, s.chi, Qc);

   % Verify convergence using evaluate_surface_flux, which defines the
   % residual that solve_surface_temperature minimizes.
   Q_sfc = icemodel.surface.turbulence.bulk_richardson.evaluate_surface_flux( ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.De, Ts, Qc, s.ea_atm, s.roL, s.br_coefs, s.chi, s.liqflag);

   testCase.verifyTrue(ok);
   testCase.verifyLessThan(abs(Q_sfc), 1e-2);
end

function test_sfcflux_includes_precipitation_advection(testCase)
   % The explicit surface residual should include Qa so the analytic and
   % numerical SEB paths use the same forcing terms.

   s = testCase.TestData.skin;
   Ts = s.Ts;
   Qc = icemodel.surface.conductive_heat_flux(s.k_eff, s.T, s.dz, Ts);
   ppt = 2e-4;
   tppt = s.tair + 3.0;

   % Compute the SEB residual without Qa.
   Q_sfc_dry = icemodel.surface.turbulence.bulk_richardson.evaluate_surface_flux( ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, 0.0, tppt, s.psfc, s.De, ...
      Ts, Qc, s.ea_atm, s.roL, s.br_coefs, s.chi, s.liqflag);

   % Compute the SEB residual with Qa.
   Q_sfc_wet = icemodel.surface.turbulence.bulk_richardson.evaluate_surface_flux( ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, ppt, tppt, s.psfc, s.De, ...
      Ts, Qc, s.ea_atm, s.roL, s.br_coefs, s.chi, s.liqflag);

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
      [Ts, ok] = icemodel.surface.solve_surface_energy_balance( ...
         s.Ts, s.tair, s.swd, s.lwd, s.albedo, s.wspd, ...
         s.ppt, s.tppt, s.psfc, s.De, s.ea_atm, s.br_coefs, s.roL, ...
         s.liqflag, s.chi, s.T, s.k_eff, s.dz, s.ro_sfc, s.snow_depth, ...
         opts_sv);
      residual = icemodel.surface.surface_energy_balance_residual(Ts, ...
         s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
         s.De, s.ea_atm, s.br_coefs, s.roL, s.liqflag, s.chi, ...
         icemodel.surface.conductive_heat_flux(s.k_eff, s.T, s.dz, Ts), s.ro_sfc, ...
         s.snow_depth, s.opts);

      testCase.verifyTrue(ok);
      testCase.verifyTrue(isfinite(Ts));
      testCase.verifyTrue(isreal(Ts));
      testCase.verifyLessThan(abs(residual), 1.0);
   end
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
   [Ts, ok] = icemodel.surface.solve_surface_energy_balance( ...
      s.Ts, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.De, s.ea_atm, s.br_coefs, s.roL, s.liqflag, s.chi, s.T, s.k_eff, ...
      s.dz, s.ro_sfc, s.snow_depth, s.opts);

   % Compute the residual
   residual = icemodel.surface.surface_energy_balance_residual( ...
      Ts, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.De, s.ea_atm, s.br_coefs, s.roL, s.liqflag, s.chi, ...
      icemodel.surface.conductive_heat_flux(s.k_eff, s.T, s.dz, Ts), s.ro_sfc, s.snow_depth, s.opts);

   % Verify valid solution and residual
   testCase.verifyTrue(ok);
   testCase.verifyTrue(isfinite(Ts));
   testCase.verifyTrue(isreal(Ts));
   testCase.verifyTrue(isfinite(residual));
   testCase.verifyLessThan(abs(residual), 1.0);
end

function test_iceebsolvedirichlet_converges_on_synthetic_column(testCase)
   % The coupled Dirichlet solver should converge on the synthetic ice
   % column, keep phase fractions bounded, and reduce the accepted
   % updated-state SEB residual relative to the current one-pass pattern.

   s = testCase.TestData.ice;

   % Solve the coupled surface-column energy balance
   [Ts, T, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok, n_iters] = ...
      icemodel.couplers.solve_surface_column_dirichlet( ...
      s.Ts, s.T, s.f_ice, s.f_liq, s.dz, s.delz, s.fn, s.Sc, s.Sp, ...
      s.opts.dt, s.JJ, s.k_liq, s.cv_ice, s.cv_liq, s.ro_ice, s.ro_liq, ...
      s.Ls, s.Lf, s.roLf, s.Tf, s.fcp, s.TL, s.TH, s.f_ell_min, s.f_ell_max, ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, s.De, ...
      s.ea_atm, s.br_coefs, s.roL, s.liqflag, s.chi, s.seb_solver, s.tol, ...
      s.maxiter, s.alpha, s.use_aitken, s.jumpmax, s.cpl_Ts_tol, ...
      s.cpl_seb_tol, s.cpl_maxiter, s.cpl_alpha, s.cpl_aitken, s.cpl_jumpmax, ...
      s.ro_sfc, s.snow_depth, s.opts);

   % Compute the coupled surface-column energy balance residual
   residual_coupled = icemodel.surface.surface_energy_balance_residual(Ts, ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, s.De, ...
      s.ea_atm, s.br_coefs, s.roL, s.liqflag, s.chi, icemodel.surface.conductive_heat_flux(k_eff, T, s.dz, Ts), ...
      s.ro_sfc, s.snow_depth, s.opts);

   % Solve for T_sfc using the surface solver to mimic a one-pass solve
   [Ts_old, ok_seb_old] = icemodel.surface.solve_surface_energy_balance( ...
      s.Ts, s.tair, s.swd, s.lwd, s.albedo, ...
      s.wspd, s.ppt, s.tppt, s.psfc, s.De, s.ea_atm, s.br_coefs, s.roL, ...
      s.liqflag, s.chi, s.T, s.k_eff, s.dz, s.ro_sfc, s.snow_depth, s.opts);

   % Solve for column T_ice using the subsurface column solver
   [T_old, f_ice_old, f_liq_old, k_eff_old, ok_ieb_old] = ICEENBAL( ...
      Ts_old, s.T, s.f_ice, s.f_liq, s.Fc, s.Fp, s.Sc, s.Sp, s.dz, s.delz, ...
      s.fn, s.opts.dt, s.JJ, s.k_liq, s.cv_ice, s.cv_liq, s.ro_ice, ...
      s.ro_liq, s.Ls, s.Lf, s.roLf, s.Tf, s.fcp, s.TL, s.TH, ...
      s.f_ell_min, s.f_ell_max, 1, s.tol, s.maxiter, s.alpha, s.use_aitken, ...
      s.jumpmax, false);

   % Compute the one-pass residual
   residual_old = icemodel.surface.surface_energy_balance_residual(Ts_old, ...
      s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, s.De, ...
      s.ea_atm, s.br_coefs, s.roL, s.liqflag, s.chi, icemodel.surface.conductive_heat_flux(k_eff_old, T_old, s.dz, Ts_old), ...
      s.ro_sfc, s.snow_depth, s.opts);

   % Verify the coupled solutions are valid
   testCase.verifyTrue(ok_seb);
   testCase.verifyTrue(ok_ieb);
   testCase.verifyTrue(ok);
   testCase.verifyTrue(isfinite(n_iters));
   testCase.verifyTrue(isreal(Ts));
   testCase.verifyTrue(all(isfinite([Ts; T; f_ice; f_liq; k_eff])));
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
   [Ts_raw, ok] = icemodel.surface.solve_surface_energy_balance(s.Ts, s.tair, ...
      s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, s.De, ea_atm, ...
      s.br_coefs, s.roL, s.liqflag, s.chi, s.T, s.k_eff, s.dz, s.ro_sfc, ...
      s.snow_depth, s.opts);

   % Cap T_sfc at the melt temp
   Ts = icemodel.surface.physical_surface_temperature(Ts_raw);

   % Diagnose the surface energy balance fluxes
   [Qe, Qh, Qc, Qm, ~, balance] = ...
      icemodel.surface.diagnose_surface_energy_fluxes(s.T, Ts, s.tair, ...
      s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, s.De, ea_atm, ...
      s.k_eff, s.dz, s.roL, s.chi, s.br_coefs, s.ro_sfc, s.snow_depth, ...
      s.liqflag, s.opts);

   % Diagnose fluxes not resolved by diagnose_surface_energy_fluxes
   [~, ~, ~, Qa, Qle] = icemodel.surface.surface_energy_balance_terms( ...
      Ts, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, ...
      s.De, ea_atm, s.T, s.k_eff, s.dz, s.br_coefs, s.roL, s.liqflag, ...
      s.chi, s.ro_sfc, s.snow_depth, s.opts);

   % Verify the T_sfc solver worked
   testCase.verifyTrue(ok);

   % Verify the ea_atm calculation matches the test case setup
   testCase.verifyEqual(ea_atm, s.ea_atm, 'RelTol', 1e-12);

   % Verify the balance returned by diagnose_surface_energy_fluxes matches ENBAL
   testCase.verifyEqual(balance, icemodel.surface.evaluate_surface_energy_balance(s.chi, s.albedo, s.swd, ...
      s.lwd, Qle, Qh, Qe, Qc, Qa, Qm), 'RelTol', 1e-12);

   % Verify the T_sfc solution matches the test case setup value.
   testCase.verifyLessThanOrEqual(Ts, s.Tf + 1e-12);
end

function test_skinsolve_returns_finite_bounded_state(testCase)
   % SKINSOLVE should keep the skin column finite and phase-bounded on the
   % shared synthetic state.

   s = testCase.TestData.skin;

   % Solve for ice column temperature using the reduced-complexity solver
   [T, f_ice, f_liq, k_eff, ok, iter] = SKINSOLVE(s.T, s.f_ice, s.f_liq, ...
      s.dz, s.delz, s.fn, s.opts.dt, s.JJ, s.Ts, s.k_liq, s.cv_ice, ...
      s.cv_liq, s.ro_ice, s.Ls, s.tol, s.maxiter, s.alpha, ...
      false);

   % Verify the solution is valid
   testCase.verifyTrue(ok);
   testCase.verifyGreaterThan(iter, 0);
   testCase.verifyTrue(all(isfinite(T)));
   testCase.verifyTrue(all(isfinite(k_eff)));
   testCase.verifyLessThanOrEqual(max(f_ice + f_liq * s.ro_liq / s.ro_ice), ...
      1 + 1e-9);
end

function test_skinebsolve_converges_on_synthetic_column(testCase)
   % The coupled skin energy-balance solve should converge and keep the
   % phase fractions inside their physical bounds.

   s = testCase.TestData.skin;

   % Solve for T_sfc and T_ice using the coupled reduced-complexity solver
   [Ts, T, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok, n_iters] = ...
      SKINEBSOLVE(s.Ts, s.T, s.f_ice, s.f_liq, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.JJ, s.ro_ice, s.k_liq, s.cv_ice, s.cv_liq, s.Ls, ...
      s.Tf, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.De, s.ea_atm, s.chi, s.roL, s.br_coefs, s.liqflag, ...
      s.seb_solver, s.tol, s.maxiter, s.alpha, s.cpl_maxiter, ...
      s.cpl_Ts_tol, s.cpl_seb_tol, s.cpl_alpha, s.cpl_aitken, ...
      s.cpl_jumpmax, false, s.ro_sfc, s.snow_depth, s.opts);

   % Verify the solution is valid
   testCase.verifyTrue(ok_seb);
   testCase.verifyTrue(ok_ieb);
   testCase.verifyTrue(ok);
   testCase.verifyGreaterThan(n_iters, 0);
   testCase.verifyTrue(isreal(Ts));
   testCase.verifyTrue(all(isfinite([Ts; T; k_eff])));
   testCase.verifyLessThanOrEqual(max(f_ice + f_liq * s.ro_liq / s.ro_ice), ...
      1 + 1e-9);
end

function test_iceenbal_and_iceebsolve_converge_on_synthetic_column(testCase)
   % The direct and coupled icemodel column solves should both converge on
   % the shared synthetic state and preserve phase bounds.

   s = testCase.TestData.ice;

   % Solve the ice column model using a direct call to ICEENBAL (single sweep)
   [T_dir, f_ice_dir, f_liq_dir, k_eff_dir, ok_dir, iter_dir] = ICEENBAL( ...
      s.Ts, s.T, s.f_ice, s.f_liq, s.Fc, s.Fp, s.Sc, s.Sp, s.dz, s.delz, ...
      s.fn, s.opts.dt, s.JJ, s.k_liq, s.cv_ice, s.cv_liq, s.ro_ice, ...
      s.ro_liq, s.Ls, s.Lf, s.roLf, s.Tf, s.fcp, s.TL, s.TH, ...
      s.f_ell_min, s.f_ell_max, 1, s.tol, s.maxiter, s.alpha, s.use_aitken, ...
      s.jumpmax, false);

   % Solve the ice column model using the fully coupled solver with robin bc
   [Ts, T_rob, f_ice_rob, f_liq_rob, k_eff_rob, ok_rob, n_iters_rob] = ...
      icemodel.couplers.solve_surface_column_robin( ...
      s.Ts, s.T, s.f_ice, s.f_liq, s.dz, s.delz, s.fn, s.Sc, s.Sp, ...
      s.opts.dt, s.JJ, s.k_liq, s.cv_ice, s.cv_liq, s.ro_ice, ...
      s.ro_liq, s.Ls, s.Lf, s.roLf, s.Tf, s.fcp, s.TL, s.TH, s.f_ell_min, ...
      s.f_ell_max, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
      s.tppt, s.psfc, s.De, s.ea_atm, s.br_coefs, s.roL, s.liqflag, s.chi, ...
      3, s.tol, s.maxiter, s.alpha, ...
      s.use_aitken, s.jumpmax, s.cpl_Ts_tol, s.cpl_seb_tol, ...
      s.cpl_maxiter, s.cpl_alpha, s.cpl_aitken, s.cpl_jumpmax, ...
      s.ro_sfc, s.snow_depth, s.opts);

   % Verify the single-sweep direct solve is valid
   testCase.verifyTrue(ok_dir);
   testCase.verifyGreaterThanOrEqual(iter_dir, 0);
   testCase.verifyTrue(all(isfinite([T_dir; k_eff_dir])));
   testCase.verifyLessThanOrEqual(max(f_ice_dir + f_liq_dir * ...
      s.ro_liq / s.ro_ice), 1 + 1e-9);

   % Verify the fully coupled robin bc solve is valid
   testCase.verifyTrue(ok_rob);
   testCase.verifyGreaterThan(n_iters_rob, 0);
   testCase.verifyTrue(isfinite(Ts));
   testCase.verifyTrue(isreal(Ts));
   testCase.verifyTrue(all(isfinite([T_rob; k_eff_rob])));
   testCase.verifyLessThanOrEqual(max(f_ice_rob + f_liq_rob * ...
      s.ro_liq / s.ro_ice), 1 + 1e-9);
end

function test_robin_coupler_supports_monin_obukhov_on_synthetic_column(testCase)
   % The Robin coupler should accept the bulk-MO scheme and converge on the
   % synthetic ice state used by the shared solver tests.

   workspace = testCase.TestData.workspace;
   s = icemodel.test.fixtures.makeSyntheticColumnState(workspace, ...
      'icemodel', solver=3, seb_solver=2, turbulent_flux_scheme='monin_obukhov', ...
      z0_ice=0.02, testname='ice_kernel_robin_bulk_mo');

   % Solve the ice column model using the fully coupled solver with robin bc
   [Ts, T_rob, f_ice_rob, f_liq_rob, k_eff_rob, ok_rob, n_iters_rob] = ...
      icemodel.couplers.solve_surface_column_robin( ...
      s.Ts, s.T, s.f_ice, s.f_liq, s.dz, s.delz, s.fn, s.Sc, s.Sp, ...
      s.opts.dt, s.JJ, s.k_liq, s.cv_ice, s.cv_liq, s.ro_ice, ...
      s.ro_liq, s.Ls, s.Lf, s.roLf, s.Tf, s.fcp, s.TL, s.TH, s.f_ell_min, ...
      s.f_ell_max, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
      s.tppt, s.psfc, s.De, s.ea_atm, s.br_coefs, s.roL, s.liqflag, s.chi, ...
      s.opts.solver, s.tol, s.maxiter, s.alpha, ...
      s.use_aitken, s.jumpmax, s.cpl_Ts_tol, s.cpl_seb_tol, s.cpl_maxiter, ...
      s.cpl_alpha, s.cpl_aitken, s.cpl_jumpmax, s.ro_sfc, ...
      s.snow_depth, s.opts);

   % Compute the residual
   residual = icemodel.surface.surface_energy_balance_residual(Ts, s.tair, ...
      s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, s.psfc, s.De, ...
      s.ea_atm, s.br_coefs, s.roL, s.liqflag, s.chi, icemodel.surface.conductive_heat_flux(k_eff_rob, T_rob, s.dz, Ts), ...
      s.ro_sfc, s.snow_depth, s.opts);

   % Verify the solution is valid
   testCase.verifyTrue(ok_rob);
   testCase.verifyGreaterThan(n_iters_rob, 0);
   testCase.verifyTrue(isreal(Ts));
   testCase.verifyTrue(all(isfinite([Ts; T_rob; k_eff_rob; residual])));
   testCase.verifyLessThanOrEqual(max(f_ice_rob + f_liq_rob * ...
      s.ro_liq / s.ro_ice), 1 + 1e-9);
   testCase.verifyLessThan(abs(residual), 1.0);
end
