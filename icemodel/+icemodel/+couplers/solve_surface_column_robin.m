function [T_sfc, T_ice, f_ice, f_liq, k_eff, ok_ieb, n_iters] = ...
      solve_surface_column_robin(T_sfc, xT_ice, xf_ice, xf_liq, dz, delz, fn, ...
      Sc, Sp, dt, JJ, k_liq, cv_ice, cv_liq, ro_ice, ro_liq, Ls, Lf, roLf, Tf, ...
      fcp, TL, TH, f_ell_min, f_ell_max, tair, swd, lwd, albedo, wspd, ppt, tppt, ...
      psfc, De, ea_atm, br_coefs, roL, liqflag, chi, solver, tol, maxiter, alpha, ...
      use_aitken, jumpmax, cpl_Ts_tol, cpl_seb_tol, cpl_maxiter, cpl_alpha, ...
      cpl_aitken, cpl_jumpmax, ro_sfc, snow_depth, opts)
   %SOLVE_SURFACE_COLUMN_ROBIN Coupled icemodel Robin SEB solve.
   %
   % solver = 2 is the single-sweep special case of this Robin coupler
   % (cpl_maxiter = 1). solver = 3 runs the full outer Ts-T iterations.
   %
   %#codegen

   debug = opts.debug;

   % Initial values for SEB linearization coefficients Fc, Fp.
   [Fc, Fp] = icemodel.surface.surface_flux_linearization(T_sfc, tair, swd, ...
      lwd, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, roL, liqflag, ...
      chi, ro_sfc, snow_depth, opts);

   % Initial past Picard iterates for Aitken-acceleration.
   Ts_1 = nan;
   Ts_2 = nan;

   % Initial values for convergence checks.
   ok_cpl = false;
   Ts_old = T_sfc;
   Ts_diag = T_sfc;
   seb_res = nan;

   % Run outer Ts-T convergence loop (iterative block/Picard coupling).
   for cpliter = 1:cpl_maxiter
      Ts_old = T_sfc;

      % Inner subsurface solve with updated Ts, Fc, Fp and checkpoint state.
      [T_ice, f_ice, f_liq, k_eff, ok_ieb, n_iters, a1] = ICEENBAL(T_sfc, ...
         xT_ice, xf_ice, xf_liq, Fc, Fp, Sc, Sp, dz, delz, fn, dt, JJ, k_liq, ...
         cv_ice, cv_liq, ro_ice, ro_liq, Ls, Lf, roLf, Tf, fcp, TL, TH, ...
         f_ell_min, f_ell_max, solver, tol, maxiter, alpha, use_aitken, ...
         jumpmax, debug);

      if ~ok_ieb
         if debug
            dumpIceEbSolveFailure("iceenbal_failed", T_sfc, Ts_diag, Ts_old, ...
               T_ice, f_ice, f_liq, k_eff, Sc, dt, Fc, Fp, cpliter, cpl_maxiter, ...
               cpl_Ts_tol, cpl_seb_tol, seb_res, ok_cpl, n_iters);
         end
         break
      end

      % Diagnose Ts from frozen coefficients and updated (a1, T1).
      T_sfc = (Fc + a1 * T_ice(1)) / (a1 - Fp);
      Ts_diag = T_sfc;
      ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));

      % Update SEB linearization coefficients Fc, Fp.
      [Fc, Fp] = icemodel.surface.surface_flux_linearization(T_sfc, tair, swd, ...
         lwd, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, roL, liqflag, ...
         chi, ro_sfc, snow_depth, opts);

      % Diagnose SEB residual.
      seb_res ...
         = icemodel.surface.surface_energy_balance_residual(T_sfc, tair, swd, ...
         lwd, albedo, wspd, ppt, tppt, psfc, De, ea_atm, br_coefs, roL, liqflag, ...
         chi, icemodel.surface.conductive_heat_flux(k_eff, T_ice, dz, T_sfc), ...
         ro_sfc, snow_depth, opts);

      % Check convergence (bypass coupler if cpl_maxiter == 1).
      if (cpl_maxiter == 1) || ...
            (abs(T_sfc - Ts_old) < cpl_Ts_tol && abs(seb_res) < cpl_seb_tol)
         ok_cpl = true;
         break
      end

      % Aitken acceleration w/ relaxation-fallback (on failure or ~cpl_aitken).
      Ts_0 = T_sfc;
      T_sfc = icemodel.numerics.aitkenscalar(Ts_2, Ts_1, Ts_0, ...
         (1.0 - cpl_alpha) * Ts_old + cpl_alpha * T_sfc, ... % relaxation
         cpl_jumpmax, cpl_aitken);
      Ts_2 = Ts_1;
      Ts_1 = Ts_0;

      % The accelerated iterate is the actual state carried into the next
      % coupling sweep. Accept convergence here as well so a stationary
      % Aitken iterate cannot fall through the loop and fail spuriously at
      % cpl_maxiter.
      if abs(T_sfc - Ts_old) < cpl_Ts_tol && abs(seb_res) < cpl_seb_tol
         ok_cpl = true;
         break
      end
   end

   % Hitting max coupling iterations without ok_cpl is a substep fail
   ok_ieb = ok_ieb && ok_cpl;

   if ~ok_ieb && debug
      dumpIceEbSolveFailure("coupler_nonconvergence", T_sfc, Ts_diag, Ts_old, ...
         T_ice, f_ice, f_liq, k_eff, Sc, dt, Fc, Fp, cpliter, cpl_maxiter, ...
         cpl_Ts_tol, cpl_seb_tol, seb_res, ok_cpl, n_iters);
   end
end

function dumpIceEbSolveFailure(reason, Ts, Ts_diag, Ts_old, T, f_ice, ...
      f_liq, k_eff, Sc, dt, Fc, Fp, cpliter, cpl_maxiter, cpl_Ts_tol, ...
      cpl_seb_tol, seb_res, ok_cpl, n_iters)
   %DUMPICEEBSOLVEFAILURE Save coupled ice-SEB solver diagnostics on demand.

   debug_file = getenv('ICEMODEL_DEBUG_ICEEBSOLVE_FILE');
   if isempty(debug_file)
      return
   end

   debug_state = struct();
   debug_state.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   debug_state.reason = reason;
   debug_state.Ts = Ts;
   debug_state.Ts_diag = Ts_diag;
   debug_state.Ts_old = Ts_old;
   debug_state.T = T;
   debug_state.f_ice = f_ice;
   debug_state.f_liq = f_liq;
   debug_state.k_eff = k_eff;
   debug_state.Sc = Sc;
   debug_state.dt = dt;
   debug_state.Fc = Fc;
   debug_state.Fp = Fp;
   debug_state.cpliter = cpliter;
   debug_state.cpl_maxiter = cpl_maxiter;
   debug_state.cpl_Ts_tol = cpl_Ts_tol;
   debug_state.cpl_seb_tol = cpl_seb_tol;
   debug_state.seb_res = seb_res;
   debug_state.ok_cpl = ok_cpl;
   debug_state.n_iters = n_iters;

   save(debug_file, 'debug_state');
end
