function [T] = SKINSOLVE(Tsfc,T,k_eff,ro_sno,cp_sno,dz,dt,N,fn,delz)
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
%    0*To       +   1*Tsfc      +     0*T1                    = Tsfc
%   -c1*Tsfc    +   a1*T1       +   -b1*T2                    = d1
%   -c2*T1      +   a2*T2       +   -b2*T3                    = d2
%   -c3*T2      +   a3*T3       +   -b3*T4                    = d3
%     ...           ...             ...                        ...
%     ...           ...             ...                        ...
%   -cK-1*TK-2  +   aK-1*TK-1   +   -bK-1*TK                  = d(K-1)
%   -cK*TK-1    +   aK*TK       +   0*TK+1                    = d(K)
% 
%--------------------------------------------------------------------------

% save a copy of the surface temperature
xT = Tsfc;

% set Sc/Sp to zero for pure 'skin' melt model
Sc = 0.0.*dz;
Sp = 0.0.*dz;

% solve the nonlinear heat equation by iteration (p. 47)
tol = 1e-1; dif = 2*tol; iter = 0; maxiter = 100;

while any(dif > tol) && iter < maxiter

   % compute gamma at the control volume interfaces (eq. 4.9, p. 45) (JJ+1)
   g_ns = [k_eff(1); k_eff(1:N); k_eff(N)];
   gb_ns = 1.0./( (1.0-fn)./g_ns(1:N+1)+fn./g_ns(2:N+2));

   % compute the enthalpy coefficient for each c.v. for the current timestep
   aP0 = ro_sno .* cp_sno .* dz ./ dt;

   % if including vapor heat:
   %[ ~,                                                      ...
   %  drovdT] = VAPORHEAT(T,f_liq,f_ice,Tf,Rv,Ls)
   %aP0 = (ro_sno.*cp_sno+Ls.*(1-f_liq-f_ice).*drovdT).*dz./dt;

   % compute the aN and aS coefficients
   aN = gb_ns(1:N)   ./ delz(1:N);
   aS = gb_ns(2:N+1) ./ delz(2:N+1);

   % Account for the boundary conditions.
   T_N = xT;
   bc_N = aN(1) * T_N;
   bc_S = 0.0;
   aS(N) = 0.0;

   % compute the aP coefficient and solution vector b
   aP = aN(1:N) + aS(1:N) + aP0(1:N) - Sp(1:N) .* dz(1:N);
   b = Sc(1:N) .* dz(1:N) + aP0(1:N) .* T(1:N);

   % modify b to account for Dirichlet boundary conditions
   b(1) = b(1) + bc_N;
   b(N) = b(N) + bc_S;

   % solve the equation
   x = TRISOLVE(-aN,aP,-aS,b);
   
   % prep for next iteration    
   dif = abs(T-x);
   T = x;
   iter = iter+1;
end
