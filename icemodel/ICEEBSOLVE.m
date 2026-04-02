function [Ts, T, f_ice, f_liq, k_eff, ok, n_iters] = ICEEBSOLVE( ...
      xT, xf_ice, xf_liq, dz, delz, fn, Sc, dt, JJ, Ts, k_liq, cv_ice, ...
      cv_liq, ro_ice, ro_liq, Ls, Lf, roLf, Rv, Tf, fcp, TL, TH, f_ell_min, ...
      f_ell_max, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, De, ea, ...
      cv_air, emiss, SB, roL, scoef, chi, liqflag, solver, tol, maxiter, ...
      alpha, use_aitken, jumpmax, cpl_Ts_tol, cpl_seb_tol, cpl_maxiter, ...
      cpl_alpha, cpl_aitken, cpl_jumpmax, debug, ro_sfc, snow_depth, opts)
   %ICEEBSOLVE Coupled surface/subsurface solve for Robin coupling modes.
   %
   % solver = 2 is the single-sweep special case of this Robin coupler
   % (cpl_maxiter = 1). solver = 3 runs the full outer Ts-T iterations.
   %
   %#codegen

   % Initial values for SEB linearization coefficients Fc, Fp.
   [Fc, Fp] = SFCFLIN(tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, De, ...
      ea, cv_air, cv_liq, emiss, SB, roL, scoef, chi, Tf, Ts, liqflag);

   % Initial past Picard iterates for Aitken-acceleration.
   Ts_1 = nan;
   Ts_2 = nan;

   % Initial values for convergence checks.
   ok_cpl = false;
   Ts_old = Ts;
   Ts_diag = Ts;
   seb_res = nan;

   % Run outer Ts-T convergence loop (iterative block/Picard coupling).
   for cpliter = 1:cpl_maxiter
      Ts_old = Ts;

      % Inner subsurface solve with updated Ts, Fc, Fp and checkpoint state.
      [T, f_ice, f_liq, k_eff, ok, n_iters, a1] = ICEENBAL(Ts, xT, xf_ice, ...
         xf_liq, Fc, Fp, Sc, [], dz, delz, fn, dt, JJ, k_liq, cv_ice, ...
         cv_liq, ro_ice, ro_liq, Ls, Lf, roLf, Tf, fcp, TL, TH, ...
         f_ell_min, f_ell_max, solver, tol, maxiter, use_aitken, jumpmax, ...
         debug);

      if ~ok
         if debug
            dumpIceEbSolveFailure("iceenbal_failed", Ts, Ts_diag, Ts_old, T, ...
               f_ice, f_liq, k_eff, Sc, dt, Fc, Fp, cpliter, cpl_maxiter, ...
               cpl_Ts_tol, cpl_seb_tol, seb_res, ok_cpl, n_iters);
         end
         break
      end

      % Diagnose Ts from frozen coefficients and updated (a1, T1).
      Ts = (Fc + a1 * T(1)) / (a1 - Fp);
      Ts_diag = Ts;
      ro_sfc = icemodel.surface.surface_bulk_density(f_ice(1), f_liq(1));

      % Update SEB linearization coefficients Fc, Fp.
      [Fc, Fp] = SFCFLIN(tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, ...
         De, ea, cv_air, cv_liq, emiss, SB, roL, scoef, chi, Tf, Ts, liqflag);

      % Diagnose SEB residual.
      seb_res = fSEB(Ts, tair, swd, lwd, albedo, wspd, ppt, tppt, ...
         psfc, De, ea, chi, roL, scoef, CONDUCT(k_eff, T, dz, Ts), ...
         liqflag, ro_sfc, snow_depth, opts);

      % Check convergence (bypass coupler if cpl_maxiter == 1).
      if (cpl_maxiter == 1) || ...
            (abs(Ts - Ts_old) < cpl_Ts_tol && abs(seb_res) < cpl_seb_tol)
         ok_cpl = true;
         break
      end

      % Aitken acceleration w/ relaxation-fallback (on failure or ~cpl_aitken).
      Ts_0 = Ts;
      Ts = aitkenscalar(Ts_2, Ts_1, Ts_0, ...
         (1.0 - cpl_alpha) * Ts_old + cpl_alpha * Ts, ... % relaxation
         cpl_jumpmax, cpl_aitken);
      Ts_2 = Ts_1;
      Ts_1 = Ts_0;

      % The accelerated iterate is the actual state carried into the next
      % coupling sweep. Accept convergence here as well so a stationary
      % Aitken iterate cannot fall through the loop and fail spuriously at
      % cpl_maxiter.
      if abs(Ts - Ts_old) < cpl_Ts_tol && abs(seb_res) < cpl_seb_tol
         ok_cpl = true;
         break
      end
   end

   % Hitting max coupling iterations without ok_cpl is a substep fail
   ok = ok && ok_cpl;
   
   if ~ok && debug
      dumpIceEbSolveFailure("coupler_nonconvergence", Ts, Ts_diag, Ts_old, T, ...
         f_ice, f_liq, k_eff, Sc, dt, Fc, Fp, cpliter, cpl_maxiter, ...
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
