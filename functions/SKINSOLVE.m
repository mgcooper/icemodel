function [T, OK] = SKINSOLVE(Tsfc, T, k_eff, ro_sno, cp_sno, dz, dt, N, ...
      fn, delz, f_liq, f_ice, Tf, Rv, Ls)
   %SKINSOLVE Solve the 1-dimensional heat conduction equation
   %--------------------------------------------------------------------------
   %
   % [  1    0                                      ] * [ Tsfc ] = [ Tsfc ]
   % |-c1   a1   -b1                                |   | T1   |   | d1   |
   % |     -c2    a2   -b2                          |   | T2   |   | d2   |
   % |           -c3    a3   -b3                    |   | T3   |   | d3   |
   % |                  .     .     .               |   | .    |   |   .  |
   % |                        .     .     .         |   | .    |   |   .  |
   % |                             -cK-1  aK-1 -bK-1|   | TK-1 |   | dK-1 |
   % [                                   -cK    aK  ]   [ TK   ]   [ dK   ]
   %
   %      0        +   a1*T1       +   -b1*T2                    = d1 + c1*Tsfc
   %   -c2*T1      +   a2*T2       +   -b2*T3                    = d2
   %   -c3*T2      +   a3*T3       +   -b3*T4                    = d3
   %     ...           ...             ...                        ...
   %     ...           ...             ...                        ...
   %   -cK-1*TK-2  +   aK-1*TK-1   +   -bK-1*TK                  = d(K-1)
   %   -cK*TK-1    +   aK*TK       +   0*TK+1                    = d(K)
   %
   %--------------------------------------------------------------------------

   % Solver options
   persistent tol maxiter alpha
   if isempty(tol); tol = 1e-2; end
   if isempty(maxiter); maxiter = 1000; end
   if isempty(alpha); alpha = 1.8; end
  
   % Solve the nonlinear heat equation by iteration (p. 47)
   dif = 2*tol; iter = 0; OK = true;

   while any(dif > tol) && iter < maxiter

      % Compute gamma at the control volume interfaces (eq. 4.9, p. 45) (JJ+1)
      g_ns = [k_eff(1); k_eff(1:N); k_eff(N)];
      gb_ns = 1.0 ./ ( (1.0 - fn) ./ g_ns(1:N+1) + fn ./ g_ns(2:N+2));

      % Compute the enthalpy coefficient for each c.v. for the current timestep
      [~, drovdT] = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls);
      aP0 = (ro_sno .* cp_sno + Ls * (1 - f_liq - f_ice) .* drovdT) .* dz / dt;

      % Compute the aN and aS coefficients
      aN = gb_ns(1:N)   ./ delz(1:N);
      aS = gb_ns(2:N+1) ./ delz(2:N+1);

      % Account for the boundary conditions.
      bc_N = aN(1) * Tsfc;
      bc_S = 0.0;
      aS(N) = 0.0;

      % Compute the aP coefficient and solution vector b
      aP = aN(1:N) + aS(1:N) + aP0(1:N); % - Sp(1:N) .* dz(1:N);

      % % Check diagonal dominance and condition number
      % checkdiags(aP, aN, aS)

      b = aP0(1:N) .* T(1:N); % + Sc(1:N) .* dz(1:N);

      % Account for Dirichlet upper and Neumann lower boundary conditions
      b(1) = b(1) + bc_N;
      b(N) = b(N) + bc_S;

      % Solve the equation
      x = TRISOLVE(-aN, aP, -aS, b);

      % Apply over-relaxation
      x = alpha * x + (1 - alpha) * T;

      % Prep for next iteration
      dif = abs(T - x);
      T = x;
      iter = iter + 1;
   end
   if any(dif > tol)
      OK = false;
   end
end

function isDominant = checkdiags(aP, aN, aS)
   % Check diagonal dominance for internal nodes
   isDominant = all(aP(2:end-1) >= abs(-aN(2:end-1)) + abs(-aS(2:end-1)));

   % Check diagonal dominance for boundary nodes
   isDominant = isDominant && (aP(1) >= abs(aN(1)));
   isDominant = isDominant && (aP(end) >= abs(aS(end)));

   % Construct the tridiagonal matrix A and compute the condition number
   A = diag(aP) + diag(-aN(2:end), 1) + diag(-aS(1:end-1), -1);
   C = 1 / rcond(A);
   % C = cond(A);
   
   if ~isDominant
      warning( ...
         ['The system matrix is not diagonally dominant. ' ...
         'Convergence might be slow or not guaranteed.']);
   end
end

% function omega = optimalOmega(A)
%    % Given A, the coefficient matrix
%    D = diag(diag(A));
%    L = tril(A, -1);
%    U = triu(A, 1);
%
%    % Calculate the Gauss-Seidel iteration matrix G
%    G = -(D + L) \ U;
%
%    % Calculate the spectral radius of G
%    rho = max(abs(eig(G)));
%
%    % Calculate the optimal omega for SOR
%    omega = 4 / (2 - sqrt(4 - rho^2));
% end
