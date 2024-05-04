function [T, f_ice, f_liq, k_eff, OK, iter, a1, err] = ICEENBAL(T, f_ice, ...
      f_liq, dz, delz, fn, Sc, dt, JJ, Ts, k_liq, cv_ice, cv_liq, ro_ice, ...
      ro_liq, Ls, Lf, roLf, Rv, Tf, fcp, TL, TH, f_ell_min, f_ell_max, ...
      Fc, Fp, bc)
   %ICEENBAL Solve the ice energy balance.
   %
   %#codegen

   % Solver options
   tol = 1e-3;
   maxiter = 100;
   % alpha = 0.8;

   % Update the water fraction
   f_wat = f_liq + f_ice * ro_ice / ro_liq;

   % Update the melt-zone boundaries. These are the volumetric liquid 
   % fractions at T=TL and T=TH, given the current total water fraction.
   f_liq_min = f_wat .* f_ell_min;
   f_liq_max = f_wat .* f_ell_max;

   % Compute vapor density [kg m-3]
   ro_vap = VAPORHEAT(T, f_ice, f_liq, Tf, Rv, Ls);

   % Compute enthalpy [J m-3]
   H_old = TOTALHEAT(T, f_ice, f_liq, cv_ice, cv_liq, roLf, Ls * ro_vap, Tf);

   % Store past values
   T_old = T;
   f_liq_old = f_liq;

   % Initialize current values
   T_iter = T_old + 2 * tol;

   % Iterate to solve the nonlinear heat equation
   OK = true;
   for iter = 0:maxiter-1

      % Update vapor heat
      [ro_vap, drovdT, k_vap] = VAPORHEAT(T, f_ice, f_liq, Tf, Rv, Ls);

      % Update total enthalpy
      H = TOTALHEAT(T, f_ice, f_liq, cv_ice, cv_liq, roLf, Ls * ro_vap, Tf);

      % Update thermal conductivity
      k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, k_vap);

      % Update the derivative of enthalpy wrt temperature
      dHdT = cv_ice * f_ice + cv_liq * f_liq;
      dLdT = 2.0 * fcp ^ 2.0 * (Tf - min(T, Tf)) .* f_wat ...
         ./ (1.0 + fcp ^ 2.0 * (Tf - min(T, Tf)) .^ 2.0) .^ 2.0;

      % Update the general equation coefficients
      [aN, aP, aS, b, iM, a1, a2] = GECOEFS(T, f_ice, f_liq, dHdT, dLdT, ...
         drovdT, H-H_old, Sc, k_eff, delz, fn, dz, dt, Ts, Ls, Lf, ro_liq, ...
         TL, JJ, Fc, Fp, bc);

      % % Check diagonal dominance and condition number
      % icemodel.checkdiags(aP, aN, aS)

      % Exit here so the state variables are updated on the final iteration
      if all(abs(T - T_iter) < tol)
         break
      end

      % Capture past values
      T_iter = T;

      % Solve the equation (predictor step)
      T = TRISOLVE(-aN, aP, -aS, b);

      % Update the temperature-enthalpy relationship (corrector step)
      [T, f_ice, f_liq, OK] = MZTRANSFORM(T, T_iter, f_ice, f_liq, f_wat, ...
         ro_ice, ro_liq, Tf, TL, TH, fcp, f_liq_min, f_liq_max, iM, OK);

      % If failure, return to the main program and shorten the timestep
      if OK == false; return; end

      % Mass conservation check
      assertF(@() all(f_ice + f_liq * ro_liq / ro_ice <= 1 + eps))

      % Relaxation. Doesn't improve convergence but keep for reference.
      % T = alpha * T + (1 - alpha) * T_iter;
   end

   OK = iter < maxiter;

   % Surface energy balance linearization error [K]
   % err0 = -(Fc + Fp * Ts) / a1 - (T(1) - Ts);

   % Subsurface energy balance linearization error [K]
   err = (dt/dz(1) ...
      * (a2 * (T(2) - T(1)) ...
      + Fc + Fp * (Fc + a1 * T(1)) / (a1 - Fp) ...
      + Sc(1) * dz(1)) ...
      - ro_liq * Lf * (f_liq(1) - f_liq_old(1)) ) ...
      / (dHdT(1) + Ls * drovdT(1) * (1 - f_ice(1) - f_liq(1)))...
      - (T(1) - T_old(1));
end
