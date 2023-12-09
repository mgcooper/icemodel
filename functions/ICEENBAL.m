function [T, f_ice, f_liq, OK, errT, errH, iter] = ICEENBAL(T, f_ice, f_liq, ...
      dz, delz, fn, Sc, dt, JJ, Ts, k_liq, cv_ice, cv_liq, ro_ice, ro_liq, ...
      ro_air, Ls, Lf, roLf, Rv, Tf, fcp, TL, TH, flmin, flmax, ro_iwe, ro_wie)
   %ICEENBAL Solve the ice energy balance.
   %
   %#codegen

   % Solver options
   persistent tol maxiter % alpha
   if isempty(tol)
      tol = 1e-3;
      maxiter = 100;
      % alpha = 0.8;
   end

   % Update the water fraction
   f_wat = f_liq + f_ice * ro_iwe;

   % Update the melt-zone boundaries. These are the liquid fractions at T=TL
   % and T=TH, given the current total water fraction.
   fliqmin = f_wat .* flmin;
   fliqmax = f_wat .* flmax;

   % Compute vapor density [kg m-3]
   ro_vap = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls);

   % Compute enthalpy [J m-3]
   H_old = TOTALHEAT(T, f_ice, f_liq, cv_ice, cv_liq, roLf, Ls * ro_vap, Tf);

   % Iterate to solve the nonlinear heat equation
   errT = NaN; errH = NaN; OK = true;

   for iter = 1:maxiter

      T_iter = T;

      % Update vapor heat
      [ro_vap, drovdT, k_vap] = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls);

      % Update total enthalpy
      H = TOTALHEAT(T, f_ice, f_liq, cv_ice, cv_liq, roLf, Ls * ro_vap, Tf);
      
      % Update thermal conductivity
      k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, k_vap);

      % Update density, heat capacity, and d(f_liq)/dT
      ro_sno = ro_ice * f_ice + ro_liq * f_liq + ro_air * (1.0 - f_liq - f_ice);
      cp_sno = (cv_ice * f_ice + cv_liq * f_liq) ./ ro_sno;
      dFdT = 2.0 * fcp ^ 2.0 * (Tf - min(T, Tf)) .* f_wat ...
         ./ (1.0 + fcp ^ 2.0 * (Tf - min(T, Tf)) .^ 2.0) .^ 2.0;

      % Update the general equation coefficients
      [aN, aP, aS, b, iM] = GECOEFS(T, ro_sno, cp_sno, f_liq, f_ice, Ls, ...
         Lf, ro_liq, dz, dt, dFdT, drovdT, TL, H, H_old, Sc, k_eff, fn, ...
         delz, Ts, JJ);

      % % Check diagonal dominance and condition number
      % icemodel.checkdiags(aP, aN, aS)

      % Solve the equation
      T = TRISOLVE(-aN, aP, -aS, b);

      % Back-transform the meltzone transformation
      [T, f_liq, f_ice, OK] = MZTRANSFORM(T, T_iter, f_liq, f_wat, ro_liq, ...
         ro_wie, Tf, TL, TH, fcp, fliqmin, fliqmax, iM, OK);

      assertF(@() all(f_ice + f_liq * ro_wie <= 1))

      % If failure, return to the main program and adjust the timestep
      if OK == false; return; end

      % When a node first enters the melt zone, it will not be tagged by iM at
      % the start of this loop, so this checks for that and updates f_ice/liq/T
      [f_liq, f_ice, ~, T] = MELTCURVE(T, f_liq, f_ice, ro_wie, ro_iwe, Tf, fcp);

      assertF(@() all(f_ice + f_liq * ro_wie <= 1))

      % Relaxation does not tend to improve the solution but keep for reference.
      % T = alpha * T + (1 - alpha) * T_iter;

      % Prep for next iteration
      errT = T - T_iter;
      if all(abs(errT) < tol)
         break
      end
   end
   OK = iter < maxiter;

   % Compute the enthalpy residual
   errH = TOTALHEAT(T, f_ice, f_liq, cv_ice, cv_liq, roLf, ...
      Ls * VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls), Tf) - H;
end
