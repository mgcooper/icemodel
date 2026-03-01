function [T, f_ice, f_liq, k_eff, ok, iter] = SKINSOLVE(T, f_ice, f_liq, dz, ...
      delz, fn, dt, JJ, Ts, k_liq, cv_ice, cv_liq, ro_ice, Ls, Rv, Tf, tol, ...
      maxiter, alpha)
   %SKINSOLVE Solve the 1-dimensional heat conduction equation
   %
   %#codegen

   debug = false;

   % Solver options
   if maxiter == 1
      alpha = 1;
   end

   % Update thermal conductivity (vapor heat currently held constant)
   % [~, drovdT, k_vap] = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls);
   % k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, k_vap);
   drovdT = 0;

   % Update thermal conductivity
   k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);

   % Initial past Picard iterates for Aitken-acceleration
   % T_1 = nan(size(T));
   % T_2 = nan(size(T));

   % Iterate to solve the nonlinear heat equation (p. 47)
   ok = false;
   for iter = 1:maxiter

      % Capture current T iterate
      T_iter = T;

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

      if debug == true
         plot_temp(T, T_iter, Ts, dz)
      end

      % Prep for next iteration
      if all(abs(T - T_iter) < tol)
         ok = true;
         break
      end

      % Apply relaxation
      T = alpha * T + (1 - alpha) * T_iter;

      % Aitken acceleration (node-by-node) with relaxed value as fallback.
      % if use_aitken
      %    T_0 = T;
      %    for mm = 1:numel(T)
      %       T(mm) = aitkenscalar(T_2(mm), T_1(mm), T_0(mm), T(mm), ...
      %          jumpmax);
      %    end
      %    T_2 = T_1;
      %    T_1 = T_0;
      % end

      % Update thermal conductivity (T-k_eff consistency on final iteration)
      k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);
   end
end

function plot_temp(T, T_iter, Ts, dz)
   Z = cumsum(dz) - dz / 2;
   figure; hold on
   plot(T, Z)
   plot(T_iter, Z, '--')
   scatter(Ts, 0, 'filled')
   set(gca, 'YDir', 'reverse')
   legend('T', 'T iter', 'Ts')
end
