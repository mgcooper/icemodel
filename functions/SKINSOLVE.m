function [T, OK] = SKINSOLVE(T, f_ice, f_liq, ro_sno, cp_sno, k_eff, fn, ...
      delz, dz, dt, JJ, Ts, Tf, Rv, Ls)
    
   %SKINSOLVE Solve the 1-dimensional heat conduction equation

   % Solver options
   persistent tol maxiter alpha
   if isempty(tol); tol = 1e-2; end
   if isempty(maxiter); maxiter = 1000; end
   if isempty(alpha); alpha = 1.8; end
  
   % Solve the nonlinear heat equation by iteration (p. 47)
   for iter = 1:maxiter

      T_iter = T;
      
      % Compute gamma at the control volume interfaces (eq. 4.9, p. 45) (JJ+1)
      g_ns = [k_eff(1); k_eff(1:JJ); k_eff(JJ)];
      gb_ns = 1.0 ./ ( (1.0 - fn) ./ g_ns(1:JJ+1) + fn ./ g_ns(2:JJ+2));

      % Compute the enthalpy coefficient for each c.v. for the current timestep
      [~, drovdT] = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls);
      aP0 = (ro_sno .* cp_sno + Ls * (1 - f_liq - f_ice) .* drovdT) .* dz / dt;

      % Compute the aN and aS coefficients
      aN = gb_ns(1:JJ)   ./ delz(1:JJ);
      aS = gb_ns(2:JJ+1) ./ delz(2:JJ+1);

      % Account for the boundary conditions.
      bc_N = aN(1) * Ts;
      bc_S = 0.0;
      aS(JJ) = 0.0;

      % Compute the aP coefficient and solution vector b
      aP = aN(1:JJ) + aS(1:JJ) + aP0(1:JJ); % - Sp(1:N) .* dz(1:N);
      b = aP0(1:JJ) .* T(1:JJ); % + Sc(1:N) .* dz(1:N);

      % Account for Dirichlet upper and Neumann lower boundary conditions
      b(1) = b(1) + bc_N;
      b(JJ) = b(JJ) + bc_S;

      % Solve the equation
      T = TRISOLVE(-aN, aP, -aS, b);

      % Apply over-relaxation
      T = alpha * T + (1 - alpha) * T;

      % Prep for next iteration
      if all(abs(T - T_iter) < tol); break; end
   end
   OK = iter < maxiter;
end
