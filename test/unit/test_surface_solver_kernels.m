function tests = test_surface_solver_kernels
   %TEST_SURFACE_SOLVER_KERNELS Verify local surface and column solver kernels.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Build paired skinmodel and icemodel synthetic columns so the surface
   % and coupled column solvers can be exercised on matched states.

   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(2016, ...
      configure=true, nsteps=24, dt_seconds=900);
   testCase.TestData.workspace = workspace;
   testCase.TestData.skin = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'skinmodel', solver=1, testname='skin_kernel');
   testCase.TestData.ice = icemodel.test.fixtures.makeSyntheticColumnState( ...
      workspace, 'icemodel', solver=3, testname='ice_kernel');
end

function teardown(testCase)
   % Remove the shared synthetic columns after the file-level tests end.

   icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      testCase.TestData.workspace);
end

function test_sfcflux_derivative_matches_finite_difference(testCase)
   % SFCFLUX should return a derivative consistent with a centered finite
   % difference about the same surface state.

   s = testCase.TestData.skin;
   Ts = s.Ts;
   Qc = CONDUCT(s.k_eff, s.T, s.dz, Ts);
   [Fsfc, Fdot] = SFCFLUX(s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.psfc, ...
      s.De, Ts, Qc, s.ea, s.cv_air, s.emiss, s.SB, s.roL, s.Tf, s.scoef, ...
      s.chi, s.liqflag);

   h = 1e-5;
   Fplus = SFCFLUX(s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.psfc, s.De, ...
      Ts + h, Qc, s.ea, s.cv_air, s.emiss, s.SB, s.roL, s.Tf, s.scoef, ...
      s.chi, s.liqflag);
   Fminus = SFCFLUX(s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.psfc, s.De, ...
      Ts - h, Qc, s.ea, s.cv_air, s.emiss, s.SB, s.roL, s.Tf, s.scoef, ...
      s.chi, s.liqflag);
   Fdot_fd = (Fplus - Fminus) / (2 * h);

   testCase.verifyTrue(isfinite(Fsfc));
   testCase.verifyEqual(Fdot, Fdot_fd, 'RelTol', 2e-4);
end

function test_sfctemp_finds_small_surface_residual(testCase)
   % SFCTEMP should find a temperature that leaves only a small residual in
   % the explicit surface-flux balance.

   s = testCase.TestData.skin;
   Qc = CONDUCT(s.k_eff, s.T, s.dz, s.Ts);
   [Ts, ok] = SFCTEMP(s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.psfc, ...
      s.De, s.ea, s.cv_air, s.emiss, s.SB, s.Tf, s.chi, s.roL, s.scoef, ...
      s.liqflag, Qc);
   residual = SFCFLUX(s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.psfc, ...
      s.De, Ts, Qc, s.ea, s.cv_air, s.emiss, s.SB, s.roL, s.Tf, s.scoef, ...
      s.chi, s.liqflag);

   testCase.verifyTrue(ok);
   testCase.verifyLessThan(abs(residual), 1e-2);
end

function test_sebsolve_converges_across_root_finders(testCase)
   % All standalone SEBSOLVE root-finder modes should converge on the same
   % synthetic forcing state.

   s = testCase.TestData.skin;
   for seb_solver = 0:2
      [Ts, ok] = SEBSOLVE(s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
         s.tppt, s.psfc, s.De, s.ea, s.cv_air, s.cv_liq, s.emiss, s.SB, ...
         s.Tf, s.chi, s.roL, s.scoef, s.liqflag, s.Ts, s.T, s.k_eff, ...
         s.dz, seb_solver, false);
      residual = fSEB(Ts, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, ...
         s.tppt, s.psfc, s.De, s.ea, s.cv_air, s.cv_liq, s.emiss, s.SB, ...
         s.Tf, s.chi, s.roL, s.scoef, CONDUCT(s.k_eff, s.T, s.dz, Ts), ...
         s.liqflag);

      testCase.verifyTrue(ok);
      testCase.verifyTrue(isfinite(Ts));
      testCase.verifyLessThan(abs(residual), 1.0);
   end
end

function test_enbalance_matches_energy_balance_residual(testCase)
   % ENBALANCE should report the same residual as the explicit ENBAL helper
   % built from its returned component fluxes.

   s = testCase.TestData.skin;
   rh = s.met.rh(s.metstep);
   [Qm, ~, Qh, Qe, Qc, Qle, balance, Ts, ea] = ENBALANCE(s.tair, s.swd, ...
      s.lwd, s.albedo, s.wspd, rh, s.ppt, s.tppt, s.psfc, s.De, s.T, ...
      s.k_eff, s.Tf, s.dz, s.chi, s.Ts, s.cv_air, s.cv_liq, s.emiss, ...
      s.SB, s.roL, s.scoef, s.epsilon, s.liqflag, true, s.seb_solver);
   Qa = QADVECT(s.ppt, s.tppt, s.cv_liq);

   testCase.verifyEqual(ea, s.ea, 'RelTol', 1e-12);
   testCase.verifyEqual(balance, ENBAL(s.albedo, s.emiss, s.chi, s.swd, ...
      s.lwd, Qle, Qh, Qe, Qc, Qa, Qm), 'RelTol', 1e-12);
   testCase.verifyLessThanOrEqual(Ts, s.Tf + 1e-12);
end

function test_skinsolve_returns_finite_bounded_state(testCase)
   % SKINSOLVE should keep the skin column finite and phase-bounded on the
   % shared synthetic state.

   s = testCase.TestData.skin;
   [T, f_ice, f_liq, k_eff, ok, iter] = SKINSOLVE(s.T, s.f_ice, s.f_liq, ...
      s.dz, s.delz, s.fn, s.opts.dt, s.JJ, s.Ts, s.k_liq, s.cv_ice, ...
      s.cv_liq, s.ro_ice, s.Ls, s.Rv, s.Tf, s.tol, s.maxiter, s.alpha, ...
      false);

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
   [Ts, T, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok, n_iters] = ...
      SKINEBSOLVE(s.Ts, s.T, s.f_ice, s.f_liq, s.dz, s.delz, s.fn, ...
      s.opts.dt, s.JJ, s.ro_ice, s.k_liq, s.cv_ice, s.cv_liq, s.Ls, ...
      s.Rv, s.Tf, s.tair, s.swd, s.lwd, s.albedo, s.wspd, s.ppt, s.tppt, ...
      s.psfc, s.De, s.ea, s.cv_air, s.emiss, s.SB, s.chi, s.roL, s.scoef, ...
      s.liqflag, s.seb_solver, s.tol, s.maxiter, s.alpha, s.cpl_maxiter, ...
      s.cpl_Ts_tol, s.cpl_seb_tol, s.cpl_alpha, s.cpl_aitken, ...
      s.cpl_jumpmax, false);

   testCase.verifyTrue(ok_seb);
   testCase.verifyTrue(ok_ieb);
   testCase.verifyTrue(ok);
   testCase.verifyGreaterThan(n_iters, 0);
   testCase.verifyTrue(all(isfinite([Ts; T; k_eff])));
   testCase.verifyLessThanOrEqual(max(f_ice + f_liq * s.ro_liq / s.ro_ice), ...
      1 + 1e-9);
end

function test_iceenbal_and_iceebsolve_converge_on_synthetic_column(testCase)
   % The direct and coupled icemodel column solves should both converge on
   % the shared synthetic state and preserve phase bounds.

   s = testCase.TestData.ice;
   [T_dir, f_ice_dir, f_liq_dir, k_eff_dir, ok_dir, iter_dir] = ICEENBAL( ...
      s.T, s.f_ice, s.f_liq, s.dz, s.delz, s.fn, s.Sc, s.opts.dt, s.JJ, ...
      s.Ts, s.k_liq, s.cv_ice, s.cv_liq, s.ro_ice, s.ro_liq, s.Ls, s.Lf, ...
      s.roLf, s.Rv, s.Tf, s.fcp, s.TL, s.TH, s.f_ell_min, s.f_ell_max, ...
      s.Fc, s.Fp, 1, s.tol, s.maxiter, s.alpha, s.use_aitken, s.jumpmax, ...
      false);

   [Ts, T_rob, f_ice_rob, f_liq_rob, k_eff_rob, ok_rob, n_iters_rob] = ...
      ICEEBSOLVE(s.T, s.f_ice, s.f_liq, s.dz, s.delz, s.fn, s.Sc, ...
      s.opts.dt, s.JJ, s.Ts, s.k_liq, s.cv_ice, s.cv_liq, s.ro_ice, ...
      s.ro_liq, s.Ls, s.Lf, s.roLf, s.Rv, s.Tf, s.fcp, s.TL, s.TH, ...
      s.f_ell_min, s.f_ell_max, s.tair, s.swd, s.lwd, s.albedo, s.wspd, ...
      s.ppt, s.tppt, s.psfc, s.De, s.ea, s.cv_air, s.emiss, s.SB, s.roL, ...
      s.scoef, s.chi, s.liqflag, 3, s.tol, s.maxiter, s.alpha, ...
      s.use_aitken, s.jumpmax, s.cpl_Ts_tol, s.cpl_seb_tol, ...
      s.cpl_maxiter, s.cpl_alpha, s.cpl_aitken, s.cpl_jumpmax, false);

   testCase.verifyTrue(ok_dir);
   testCase.verifyGreaterThanOrEqual(iter_dir, 0);
   testCase.verifyTrue(all(isfinite([T_dir; k_eff_dir])));
   testCase.verifyLessThanOrEqual(max(f_ice_dir + f_liq_dir * ...
      s.ro_liq / s.ro_ice), 1 + 1e-9);

   testCase.verifyTrue(ok_rob);
   testCase.verifyGreaterThan(n_iters_rob, 0);
   testCase.verifyTrue(isfinite(Ts));
   testCase.verifyTrue(all(isfinite([T_rob; k_eff_rob])));
   testCase.verifyLessThanOrEqual(max(f_ice_rob + f_liq_rob * ...
      s.ro_liq / s.ro_ice), 1 + 1e-9);
end
