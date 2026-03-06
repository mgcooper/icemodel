function [Ts, T, f_ice, f_liq, k_eff, ok_seb, ok_ieb, ok, n_iters] = ...
      SKINEBSOLVE(xTs, xT, xf_ice, xf_liq, dz, delz, fn, dt, JJ, ro_ice, ...
      k_liq, cv_ice, cv_liq, Ls, Rv, Tf, tair, swd, lwd, albedo, wspd, ...
      ppt, tppt, psfc, De, ea, cv_air, emiss, SB, chi, roL, scoef, ...
      liqflag, seb_solver, tol, maxiter, alpha, cpl_maxiter, cpl_Ts_tol, ...
      cpl_seb_tol, cpl_alpha, cpl_aitken, cpl_jumpmax)
   %SKINEBSOLVE Coupled skin-subsurface Ts-T solve.
   %
   %#codegen

   % Pre-coupler Ts predictor using checkpoint state
   k_eff = GETGAMMA(xT, xf_ice, xf_liq, ro_ice, k_liq, Ls, Rv, Tf);
   [Ts, ok_seb] = SEBSOLVE(tair, swd, lwd, albedo, wspd, ppt, tppt, ...
      psfc, De, ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, ...
      liqflag, xTs, xT, k_eff, dz, seb_solver);
   Ts = MELTTEMP(Ts, Tf);

   % Initialize past Picard iterates for Aitken-acceleration
   Ts_1 = nan;
   Ts_2 = nan;

   % Run the coupler
   ok_cpl = false;
   for cpliter = 1:cpl_maxiter

      % Inner subsurface solve from checkpoint state w/o physical advancement
      [T, f_ice, f_liq, k_eff, ok_ieb, n_iters] = SKINSOLVE(xT, ...
         xf_ice, xf_liq, dz, delz, fn, dt, JJ, Ts, k_liq, cv_ice, ...
         cv_liq, ro_ice, Ls, Rv, Tf, tol, maxiter, alpha);
      if not(ok_ieb)
         break
      end

      % Inner surface solve (in-loop corrector using updated trial state)
      old = Ts;
      [Ts, ok_seb] = SEBSOLVE(tair, swd, lwd, albedo, wspd, ppt, tppt, ...
         psfc, De, ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, ...
         liqflag, Ts, T, k_eff, dz, seb_solver);
      Ts = MELTTEMP(Ts, Tf);

      if not(ok_seb)
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
            abs(Ts - old) < cpl_Ts_tol && seb_res < cpl_seb_tol
         ok_cpl = true;
         break
      end

      % Aitken acceleration with relaxation-fallback
      Ts_0 = Ts;
      Ts = MELTTEMP(aitkenscalar(Ts_2, Ts_1, Ts_0, ...
         (1.0 - cpl_alpha) * old + cpl_alpha * Ts, ... % relaxation
         cpl_jumpmax, cpl_aitken), Tf);
      Ts_2 = Ts_1;
      Ts_1 = Ts_0;
   end

   ok = ok_seb && ok_ieb && ok_cpl;
end
