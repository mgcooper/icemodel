function [Ts, Fc, Fp, T, f_ice, f_liq, k_eff, ok, n_iters, a1] = ICEEBSOLVE( ...
      xT, xf_ice, xf_liq, dz, delz, fn, Sc, dt, JJ, Ts, k_liq, cv_ice, ...
      cv_liq, ro_ice, ro_liq, Ls, Lf, roLf, Rv, Tf, fcp, TL, TH, f_ell_min, ...
      f_ell_max, tair, swd, lwd, albedo, wspd, ppt, tppt, psfc, De, ea, ...
      cv_air, emiss, SB, roL, scoef, chi, liqflag, bc_type, tol, maxiter, ...
      alpha, use_aitken, jumpmax, cpl_Ts_tol, cpl_seb_tol, cpl_maxiter, ...
      cpl_alpha, cpl_aitken, cpl_jumpmax)
   %ICEEBSOLVE Coupled surface/subsurface solve for inner Robin coupling.
   %
   %#codegen

   % Initial values for SEB linearization coefficients Fc, Fp
   [Fc, Fp] = SFCFLIN(tair, swd, lwd, albedo, wspd, psfc, De, ea, cv_air, ...
      emiss, SB, roL, scoef, chi, Tf, Ts, liqflag);

   % Initial past Picard iterates for Aitken-acceleration
   Ts_1 = nan;
   Ts_2 = nan;

   % Run outer Ts-T convergence loop (iterative block/Picard coupling)
   ok_cpl = false;
   for cpliter = 1:cpl_maxiter
      old = Ts;

      % Inner subsurface solve with updated Ts, Fc, Fp
      [T, f_ice, f_liq, k_eff, ok, n_iters, a1] = ICEENBAL(xT, xf_ice, xf_liq, ...
         dz, delz, fn, Sc, dt, JJ, Ts, k_liq, cv_ice, cv_liq, ro_ice, ...
         ro_liq, Ls, Lf, roLf, Rv, Tf, fcp, TL, TH, f_ell_min, f_ell_max, ...
         Fc, Fp, bc_type, tol, maxiter, alpha, use_aitken, jumpmax);

      if not(ok)
         break
      end

      % Diagnose Ts from frozen coefficients and updated (a1,T1)
      Ts = (Fc + a1 * T(1)) / (a1 - Fp);

      % Update SEB linearization coefficients Fc, Fp
      [Fc, Fp] = SFCFLIN(tair, swd, lwd, albedo, wspd, psfc, De, ea, ...
         cv_air, emiss, SB, roL, scoef, chi, Tf, Ts, liqflag);

      % Diagnose SEB residual
      seb_res = fSEB(Ts, tair, swd, lwd, albedo, wspd, ppt, tppt, ...
         psfc, De, ea, cv_air, cv_liq, emiss, SB, Tf, chi, roL, scoef, ...
         CONDUCT(k_eff, T, dz, Ts), liqflag);

      % Check convergence
      if abs(Ts - old) < cpl_Ts_tol && abs(seb_res) < cpl_seb_tol
         ok_cpl = true;
         break
      end

      % Aitken acceleration with relaxation-fallback (on failure or ~cpl_aitken)
      Ts_0 = Ts;
      Ts = aitkenscalar(Ts_2, Ts_1, Ts_0, ...
         (1.0 - cpl_alpha) * old + cpl_alpha * Ts, ... % relaxation
         cpl_jumpmax, cpl_aitken);
      Ts_2 = Ts_1;
      Ts_1 = Ts_0;
   end

   % Hitting max coupling iterations without ok_cpl is a substep fail
   ok = ok && ok_cpl;
end
