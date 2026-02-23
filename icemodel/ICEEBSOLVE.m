function [Ts, Fc, Fp, T, f_ice, f_liq, k_eff, ok, N, a1] = ICEEBSOLVE( ...
      xT, xf_ice, xf_liq, dz, delz, fn, Sc, dt, JJ, Ts, k_liq, cv_ice, ...
      cv_liq, ro_ice, ro_liq, Ls, Lf, roLf, Rv, Tf, fcp, TL, TH, f_ell_min, ...
      f_ell_max, tair, swd, lwd, albedo, wspd, psfc, De, ea, cv_air, emiss, ...
      SB, roL, scoef, chi, liqflag, bc, omega, maxcpliter, cpltol)
   %ICEEBSOLVE Coupled surface/subsurface solve for inner Robin coupling.
   %
   %#codegen

   % Initial values for Fc, Fp
   [Fc, Fp] = SFCFLIN(tair, swd, lwd, albedo, wspd, psfc, De, ea, cv_air, ...
      emiss, SB, roL, scoef, chi, Tf, Ts, liqflag);

   % Outer Ts convergence
   ok_cpl = false;
   for cpliter = 1:maxcpliter
      old = Ts;

      % Inner subsurface solve with updated Ts, Fc, Fp
      [T, f_ice, f_liq, k_eff, ok, N, a1] = ICEENBAL(xT, xf_ice, xf_liq, ...
         dz, delz, fn, Sc, dt, JJ, Ts, k_liq, ...
         cv_ice, cv_liq, ro_ice, ...
         ro_liq, Ls, Lf, roLf, Rv, Tf, fcp, TL, TH, f_ell_min, ...
         f_ell_max, Fc, Fp, bc);

      if not(ok)
         break
      end

      % Diagnose Ts from frozen coefficients and updated (a1,T1)
      Ts = (Fc + a1 * T(1)) / (a1 - Fp);

      % Update coefficients
      [Fc, Fp] = SFCFLIN(tair, swd, lwd, albedo, wspd, psfc, De, ea, ...
         cv_air, emiss, SB, roL, scoef, chi, Tf, Ts, liqflag);

      % Check convergence
      if abs(Ts - old) < cpltol
         ok_cpl = true;
         break
      end

      % Apply relaxation
      Ts = (1-omega) * old + omega * Ts;
   end

   % Hitting max coupling iterations without cplOK is a substep fail
   ok = ok && ok_cpl;
end
