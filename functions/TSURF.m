function Ts = TSURF(T, ro_sno, cp_sno, f_liq, f_ice, Ls, Lf, dz, dt, dFdT, ...
      ro_liq, drovdT, H, H_old, T_old, Sc, k_eff, fn, delz, JJ)
   %TSURF Compute surface temperature from coefficients for the top node
   %
   % Note: This assumes the LHS of the equation is aP0(T1-T1o), where T1o is the
   % initial T1 at the start of iterations. If T1o on the last iteration is
   % used, then dH_1 should be H(n) - H(n-1) where n = the number of iterations.
   %
   % aP0 (T1 - T1o) = a2(T2-T1) - a1(T1-Ts) + Scdz - dH

   % coefficients for nodes below the melt zone: (W/m2/K)
   f_air = (1.0 - f_ice(1) - f_liq(1));
   aP0 = (ro_sno(1) * cp_sno(1) + Lf * ro_liq * dFdT(1) + ...
      Ls * f_air * drovdT(1)) * dz(1) / dt;

   % compute gamma at the control volume interfaces (eq. 4.9, p. 45) (JJ+1)
   g_b_ns = 1 ./ ((1 - fn) ./ [k_eff(1); k_eff] + fn ./ [k_eff; k_eff(JJ)]);

   % compute the aN and aS coefficients
   aN = g_b_ns(1:JJ)   ./ delz(1:JJ);
   aS = g_b_ns(2:JJ+1) ./ delz(2:JJ+1);
   % note that dely_p(1) and dely_p(end) are 1/2 CVs, which is correct

   aN = aN(1);
   aS = aS(1);

   % compute the aP coefficient and solution vector b
   aP = aN + aS + aP0; % -Sp.*dz;
   b = aP0 * T_old(1) + Sc(1) * dz(1) - (H(1) - H_old(1));

   Ts = (aP * T(1) - aS * T(2) - b) / aN;
end

