function [T, errH, errT, f_ice, f_liq, OK] = ICEENBAL(T, f_ice, f_liq, ...
      k_liq, cv_ice, cv_liq, ro_ice, ro_liq, ro_sno, cp_sno, Ls, Lf, roLf, ...
      Rv, Tf, dz, delz, fn, dt, JJ, Tsfc, Sc, fcp, TL, TH, flmin, flmax, ...
      ro_iwe, ro_wie) %#codegen
   %ICEENBAL Solve the ice energy balance.
   
   % Compute water fraction and the derivative of f_liq wrt to temperature
   f_wat = f_liq + f_ice * ro_iwe;
   dFdT = 2.0 * fcp ^ 2.0 * (Tf - min(T, Tf)) .* f_wat ...
      ./ (1.0 + fcp ^ 2.0 * (Tf - min(T, Tf)) .^ 2.0) .^ 2.0;
   
   % Update the melt-zone boundaries
   fliqmin = f_wat .* flmin;
   fliqmax = f_wat .* flmax;

   % Compute vapor diffusion
   ro_vap = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls);

   % Compute total enthalpy in J/m3
   H_old = TOTALHEAT(T, f_liq, f_ice, cv_liq, cv_ice, roLf, Ls * ro_vap, Tf);

   % Iterate to solve the nonlinear heat equation
   errT = 2e-2; iter = 0; errH = 0; OK = true; alpha = 0.8;

   while any(errT > 1e-2) && iter < 50

      T_old = T;
      iter = iter + 1;

      % Update thermal conductivity
      k_eff = GETGAMMA(T, f_liq, f_ice, ro_ice, k_liq, Ls, Rv, Tf);

      % Update vapor heat
      [ro_vap, drovdT] = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls);

      % Update total enthalpy
      H = TOTALHEAT(T, f_liq, f_ice, cv_liq, cv_ice, roLf, Ls * ro_vap, Tf);
      
      % Update the general equation coefficients
      [aN, aP, aS, b, iM] = GECOEFS(T, ro_sno, cp_sno, f_liq, f_ice, Ls, ...
         Lf, ro_liq, dz, dt, dFdT, drovdT, TL, H, H_old, Sc, k_eff, fn, ...
         delz, Tsfc, JJ);

      % % Check diagonal dominance and condition number
      % if ~checkdiags(aP, aN, aS)
      %    warning(['The system matrix is not diagonally dominant. ' ...
      %       'Convergence might be slow or not guaranteed.']);
      % end

      % Solve the equation
      T = TRISOLVE(-aN, aP, -aS, b);

      % Back-transform the meltzone transformation
      [T, f_liq, f_ice, OK] = MZTRANSFORM(T, T_old, f_liq, f_wat, ro_liq, ...
         ro_wie, Tf, TL, TH, fcp, fliqmin, fliqmax, iM, OK);

      % assertF(@() all(f_ice + f_liq - 1.0 <= 0))

      % If failure, return to the main program and adjust the timestep
      if OK == false; return; end

      % When a node first enters the melt zone, it will not be tagged by iM at
      % the start of this loop, so this checks for that and updates f_ice/liq/T
      [f_liq, f_ice, ~, T] = MELTCURVE(T, f_liq, f_ice, ro_wie, ro_iwe, Tf, fcp);

      % assertF(@() all(f_ice + f_liq - 1.0 <= 0))

      % Apply relaxation
      T = alpha * T + (1 - alpha) * T_old;

      % prep for next iteration
      errT = abs(T - T_old);

      % compute residual (need to add Sp)
      % resN = aN .* ([TN T(1:end-1)] - T);
      % resS = aS .* ([T(2:end) 0] - T);
      % res = (resN + resS + Sc .* dz) ./ aP0;
   end

   % assertF(@() all(f_ice + f_liq - 1.0 <= 0))

   % keep track of the residual. To accurately compute this, it is necessary to
   % update all terms (need to add that):
   errH = H - H_old;

   % % check Tsfc
   % Tsfc = TSURF(T, ro_sno, cp_sno, f_liq, f_ice, Ls, Lf, dz, dt, dFdT, ...
   %    drovdT, H, H_old, Sc, k_eff, fn, delz, JJ);
end

% Various debugging snippets
% find(iM)
% find(f_ice + f_liq - 1.0 > 0)
% find( f_wat - (f_liq + f_ice * ro_iwe))
%
% max(f_ice + f_liq - 1.0)
%
%

function isDominant = checkdiags(aP, aN, aS)
   % % Check diagonal dominance for internal nodes

   % Need to confirm if aN and aS need a neg sign in the isDominant checks
   % isDominant = all(aP(2:end-1) >= abs(aN(2:end-1)) + abs(aS(2:end-1)));
   %
   % % Check diagonal dominance for boundary nodes
   % isDominant = isDominant && (aP(1) >= abs(aN(1)));
   % isDominant = isDominant && (aP(end) >= abs(aS(end)));
   %
   % % Construct the tridiagonal matrix A and compute the condition number
   % A = diag(aP) + diag(-aN(2:end), 1) + diag(-aS(1:end-1), -1);

   % C = 1 / rcond(A);
   % % C = cond(A);
   %
   % % Compare with alternative solvers
   % A = diag(aP) + diag(-aN(2:end), 1) + diag(-aS(1:end-1), -1);
   % T1 = A \ b;
   % T2 = bicgstab(A, b);
   % T3 = gmres(A, b);
   %
   % % Solve the preconditioned system
   % M = diag(diag(A));
   % Minv = diag(1 ./ diag(M));
   % T4 = (Minv * A) \ (Minv * b);
end
