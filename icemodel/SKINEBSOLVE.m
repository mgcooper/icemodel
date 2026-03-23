function [Ts, T, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok, n_iters] = ...
      SKINEBSOLVE(xTs, xT, xf_ice, xf_liq, dz, delz, fn, dt, JJ, ro_ice, ...
      k_liq, cv_ice, cv_liq, Ls, Rv, Tf, tair, swd, lwd, albedo, wspd, ...
      ppt, tppt, psfc, De, ea, cv_air, emiss, SB, chi, roL, scoef, ...
      liqflag, seb_solver, tol, maxiter, alpha, cpl_maxiter, cpl_Ts_tol, ...
      cpl_seb_tol, cpl_alpha, cpl_aitken, cpl_jumpmax, debug)
   %SKINEBSOLVE Coupled skin-subsurface Ts-T solve.
   %
   %#codegen

   % Pre-coupler Ts predictor using checkpoint state
   k_eff = GETGAMMA(xT, xf_ice, xf_liq, ro_ice, k_liq, Ls, Rv, Tf);
   [Ts, ok_seb] = SEBSOLVE(tair, swd, lwd, albedo, wspd, ppt, tppt, ...
      psfc, De, ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, ...
      liqflag, xTs, xT, k_eff, dz, seb_solver, debug);
   Ts = MELTTEMP(Ts, Tf);

   % Initialize past Picard iterates for Aitken-acceleration
   Ts_1 = nan;
   Ts_2 = nan;

   % Initial values for convergence checks
   ok_cpl = false;
   Ts_old = Ts;
   Ts_diag = Ts;
   seb_res = nan;

   % Run outer Ts-T convergence loop (iterative block/Picard coupling).
   for cpliter = 1:cpl_maxiter

      % Inner subsurface solve from checkpoint state w/o physical advancement
      [T, f_ice, f_liq, k_eff, ok_ieb, n_iters] = SKINSOLVE(xT, ...
         xf_ice, xf_liq, dz, delz, fn, dt, JJ, Ts, k_liq, cv_ice, ...
         cv_liq, ro_ice, Ls, Rv, Tf, tol, maxiter, alpha, debug);

      if ~ok_ieb
         if debug
            dumpSkinEbSolveFailure("skinsolve_failed", Ts, Ts_diag, Ts_old, ...
               T, f_ice, f_liq, k_eff, dt, cpliter, cpl_maxiter, cpl_Ts_tol, ...
               cpl_seb_tol, seb_res, n_iters, ok_seb, ok_ieb, ok_cpl);
         end
         break
      end

      % Inner surface solve (in-loop corrector using updated trial state)
      Ts_old = Ts;
      [Ts, ok_seb] = SEBSOLVE(tair, swd, lwd, albedo, wspd, ppt, tppt, ...
         psfc, De, ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, ...
         liqflag, Ts, T, k_eff, dz, seb_solver, debug);
      Ts = MELTTEMP(Ts, Tf);
      Ts_diag = Ts;

      if not(ok_seb)
         if debug
            dumpSkinEbSolveFailure("sebsolve_failed", Ts, Ts_diag, Ts_old, T, ...
               f_ice, f_liq, k_eff, dt, cpliter, cpl_maxiter, cpl_Ts_tol, ...
               cpl_seb_tol, seb_res, n_iters, ok_seb, ok_ieb, ok_cpl);
         end
         break
      end

      % SEB residual
      seb_res = abs(fSEB(Ts, tair, swd, lwd, albedo, wspd, ppt, tppt, ...
         psfc, De, ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, ...
         CONDUCT(k_eff, T, dz, Ts), liqflag));

      % At melt cap (Ts ~= Tf), positive residual is melt energy (Qm).
      if Ts >= Tf
         seb_res = 0.0;
      end

      % Check convergence (bypass coupler if cpl_maxiter == 1)
      if (cpl_maxiter == 1) || ...
            abs(Ts - Ts_old) < cpl_Ts_tol && seb_res < cpl_seb_tol
         ok_cpl = true;
         break
      end

      % Aitken acceleration with relaxation-fallback
      Ts_0 = Ts;
      Ts = MELTTEMP(aitkenscalar(Ts_2, Ts_1, Ts_0, ...
         (1.0 - cpl_alpha) * Ts_old + cpl_alpha * Ts, ... % relaxation
         cpl_jumpmax, cpl_aitken), Tf);
      Ts_2 = Ts_1;
      Ts_1 = Ts_0;

      if abs(Ts - Ts_old) < cpl_Ts_tol && seb_res < cpl_seb_tol
         ok_cpl = true;
         break
      end
   end

   % Hitting max coupling iterations without ok_cpl is a substep fail
   ok = ok_seb && ok_ieb && ok_cpl;

   if ~ok && debug
      dumpSkinEbSolveFailure("coupler_nonconvergence", Ts, Ts_diag, Ts_old, ...
         T, f_ice, f_liq, k_eff, dt, cpliter, cpl_maxiter, cpl_Ts_tol, ...
         cpl_seb_tol, seb_res, n_iters, ok_seb, ok_ieb, ok_cpl);
   end
end

function dumpSkinEbSolveFailure(reason, Ts, Ts_diag, Ts_old, T, f_ice, ...
      f_liq, k_eff, dt, cpliter, cpl_maxiter, cpl_Ts_tol, cpl_seb_tol, ...
      seb_res, n_iters, ok_seb, ok_ieb, ok_cpl)
   %DUMPSKINEBSOLVEFAILURE Save coupled skin-model solver diagnostics.

   debug_file = getenv('ICEMODEL_DEBUG_SKINEBSOLVE_FILE');
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
   debug_state.dt = dt;
   debug_state.cpliter = cpliter;
   debug_state.cpl_maxiter = cpl_maxiter;
   debug_state.cpl_Ts_tol = cpl_Ts_tol;
   debug_state.cpl_seb_tol = cpl_seb_tol;
   debug_state.seb_res = seb_res;
   debug_state.n_iters = n_iters;
   debug_state.ok_seb = ok_seb;
   debug_state.ok_ieb = ok_ieb;
   debug_state.ok_cpl = ok_cpl;

   save(debug_file, 'debug_state');
end
