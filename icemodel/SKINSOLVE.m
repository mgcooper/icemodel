function [T, OK, iter] = SKINSOLVE(T, f_ice, f_liq, dz, delz, fn, dt, JJ, Ts, ...
      k_liq, cv_ice, cv_liq, ro_ice, Ls, Rv, Tf, varargin)
   %SKINSOLVE Solve the 1-dimensional heat conduction equation
   %
   %#codegen

   % Solver options
   tol = 1e-2;
   maxiter = 100;
   alpha = 1.8;
   if ~isempty(varargin)
      switch numel(varargin)
         case 1
            tol = varargin{1};
         case 2
            tol = varargin{1};
            maxiter = varargin{2};
         case 3
            tol = varargin{1};
            maxiter = varargin{2};
            alpha = varargin{3};
         otherwise
      end
   end
   if maxiter == 1
      alpha = 1;
   end
   drovdT = 0 * T;
   k_vap = 0 * T;

   % Solve the nonlinear heat equation by iteration (p. 47)
   for iter = 1:maxiter

      T_iter = T;

      % Update thermal conductivity
      % [~, drovdT, k_vap] = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls);
      % k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, k_vap);

      % Update thermal conductivity
      k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);

      % Compute gamma at the control volume interfaces (eq. 4.9, p. 45) (JJ+1)
      g_ns = [k_eff(1); k_eff(1:JJ); k_eff(JJ)];
      gb_ns = 1.0 ./ ( (1.0 - fn) ./ g_ns(1:JJ+1) + fn ./ g_ns(2:JJ+2));

      % Compute the enthalpy coefficient for each c.v. for the current timestep
      aP0 = (cv_ice * f_ice + cv_liq * f_liq ...
         + Ls * (1 - f_liq - f_ice) .* drovdT) .* dz / dt;

      % Compute the aN and aS coefficients
      aN = gb_ns(1:JJ)   ./ delz(1:JJ);
      aS = gb_ns(2:JJ+1) ./ delz(2:JJ+1);

      % Account for the boundary conditions.
      bc_N = aN(1) * Ts;
      bc_S = 0.0;
      aS(JJ) = 0.0;

      % Compute the aP coefficient and solution vector b
      aP = aN(1:JJ) + aS(1:JJ) + aP0(1:JJ);
      b = aP0(1:JJ) .* T(1:JJ);

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
   OK = iter <= maxiter;
end
