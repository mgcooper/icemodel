function [aN, aP, aS, b, iM] = GECOEFS(T, ro_sno, cp_sno, f_liq, f_ice, Ls, Lf, ...
      ro_liq, dz, dt, dFdT, drovdT, TL, H, H_old, Sc, k_eff, fn, delz, Ts, JJ)
   %GECOEFS Compute the general equation coefficients
   %
   %  This function constructs the lower, middle, and upper diagonals of the 
   %  A matrix in a form compatible with TRISOLVE.
   % 
   %  Note: ro_sno * cp_sno = (cv_ice * f_ice + cv_liq * f_liq)
   %  See UPDATESTATE (or UPDATESUBSTEP) for how ro_sno and cp_sno are computed
   % 
   %  Subtle point: Pmelt here is identical to SNTHRM:
   %     P = g_liq - g_liq_o
   %       = ro_liq * (f_liq - f_liq_o)
   %       = g_wat * (f_ell - f_ell_o)
   % 
   %  Jordan uses P = g_wat * (f_ell - f_ell_o) 
   %  I use P = ro_liq * (f_liq - f_liq_o)
   %  This leads to different gv/gk switches, since my model is defined in terms
   %  of volumetric fraction f_liq and Jordan's is bulk density g_liq:
   %  T = dT + To
   %  T = gv * P + gk
   %    = 1 / (ro_liq * df_liq_dT) * (ro_liq * (f_liq - f_liq_o)) + To
   %    = 1 / (f_liq - f_liq_o) / dT * (f_liq - f_liq_o) + To
   %  T = dT + To
   % 
   % The same result is found using Jordan's definition of Pmelt and gv/gk.
   % 
   % See also: ICEENBAL, TRISOLVE

   % Commented statements are kept for reference. In some cases they are needed
   % but are set when initialized e.g. the gv/gk/dFdT BCs, in other cases they
   % are not used but would be for a frozen soil model.
   
   % Melt zone indices
   iM = TL <= T;

   % For a soil model, would need indices above the melt zone
   % iM = TL <= T & T <= TH;
   % iL = T > TH;
   
   % Compute gamma at the control volume interfaces (eq. 4.9, p. 45) (JJ+1)
   g_b_ns = 1 ./ ((1 - fn) ./ [k_eff(1); k_eff] + fn ./ [k_eff; k_eff(JJ)]);
   
   % Coefficients for nodes below the melt zone [W m-2 K-1]
   f_air = (1.0 - f_ice - f_liq);
   aP0 = ro_sno .* cp_sno + Lf * ro_liq * dFdT + Ls * f_air .* drovdT;
   gv = ones(JJ, 1);     % Eq 123
   gk = zeros(JJ, 1);    % Eq 123
   LfMZ = zeros(JJ, 1);  % Eq 123, melt-zone latent heat switch

   % % If using g_liq instead of f_liq in the definition of dFdT as in SNTHERM:
   % aP0 = (ro_sno .* cp_sno + Lf * ro_sno .* dFdT + Ls * f_air .* drovdT)

   if sum(iM) > 0 % nodes inside the melt zone:
      aP0(iM) = ro_sno(iM) .* cp_sno(iM) + Ls * f_air(iM) .* drovdT(iM);
      gv(iM) = 1 ./ (ro_liq * dFdT(iM)); % Eq 122b
      gk(iM) = T(iM);
      LfMZ(iM) = Lf * dz(iM) / dt;

      % % If using g_liq instead of f_liq as in SNTHERM:
      % gv(iM) = 1 ./ (ro_sno(iM) .* dFdT(iM)); % Eq 122b
      % gv(iM) = 1 ./ (ro_liq * f_wat(iM) .* dFdT(iM)); % Eq 122b
   end

   % % For a soil model, would need indices above the melt zone
   % if sum(iL)>0   % nodes above the melt zone:
   %    aP0(iL) = cv_liq;
   %  % gv(iM) = 1.0;
   %  % gk(iM) = 0.0;
   %  % dFdT(iL) = 0.0;
   % end

   % Convert from [J m-3 K-1] to [W m-2 K-1]
   aP0 = aP0 .* dz / dt;

   % Adjust gv/gk in terms of N/P/S
   gvN = vertcat(1, gv(1:JJ-1));
   gvS = vertcat(gv(2:JJ), 0);
   gkN = vertcat(0, gk(1:JJ-1)); % zero accounts for upper dirichlet BC
   gkS = vertcat(gk(2:JJ), 0);
   % gkP = gk;
   % gvP = gv;
   % Todo: Confirm if gkN(1) is zero if a Neumann condition is used up top.

   % Compute the aN and aS coefficients [W m-2 K-1]
   aN = g_b_ns(1:JJ)   ./ delz(1:JJ);
   aS = g_b_ns(2:JJ+1) ./ delz(2:JJ+1);
   % note that delz(1) and delz(end) are 1/2 CVs

   % Account for Neumann boundary conditions before constructing A
   bc_N = Ts;               % Dirichlet: TN = known
   bc_S = 0.0;                % Neumann: dT/dz = 0.0
   aS(JJ) = 0.0;              % Neumann: dT/dz = 0.0
   % aN(1) = 0.0;             % Neumann: qB = known

   % Compute the aP coefficient and solution vector b
   aP = aN + aS + aP0; % -Sp.*dz;
   b = aP0 .* T + Sc .* dz - (H - H_old) .* dz / dt; % [W m-2]

   % Modify b to account for boundary conditions
   b(1) = b(1) + aN(1) * bc_N;
   b(JJ) = b(JJ) + aS(JJ) * bc_S;

   % Apply the melt zone (enthalpy) transformation (Eq. 128/29)
   % b = b + aN .* gkN - (aN + aS + aP0) .* gk + aS .* gkS ;
   b = b + aN .* gkN - aP .* gk + aS .* gkS ;
   aN = aN .* gvN;
   aS = aS .* gvS;
   aP = aP .* gv + LfMZ;
   
   % Account for the upper boundary condition. Not strictly necessary b/c aN(1)
   % is not involved in the TRISOLVE solution, but technically correct. Must be
   % done here, after all other terms involving aN are computed, unlike aS(JJ),
   % which must be done prior to constructing A.
   % Update: Unless the original value is needed e.g., to update the wall temp
   % aN(1) = 0.0;
end

% NOTES:

% n = surface interface
% j+1 = N   (n-1)
% j+1/2 = N-P interface
% j = P
% j-1/2 - P-S interface
% j-1 = S

% A3 = aN
% A2 = aP
% A1 = aS
% Qs = aP0

% switches, interior nodes
% aN = aN*gvN;                       (inside and outside the meltzone)
% aS = aS*gvS;                       (inside and outside the meltzone)
% aP = (aP0+aN+aS)*gvP + Lf*dz/dt;   (inside)
% aP = (aP0+aN+aS);                  (outside)
% b = aP0*T_old + aS*gkS - aP0*gkP - aN*gkP - aS*gkP + aN*gkN + Sc*dx;
%                                       (inside and out)

% aP0 = (ro_sno*cp_sno + Lf*ro_liq*Fbar + Lv*f_air*CkT)*dz/dt (outside)
% aP0 = (ro_sno*cp_sno +      0         + Lv*f_air*CkT)*dz/dt (inside)

% gv = 1./ro_liq./Fbar;  (inside)
% gk = T_old             (inside)
% gv = 1;                (outside)
% gv = 0;                (outside)

% T_old = T_old;            (outside)
% T_old = P_melt;           (inside) (convert to get T_old)

% top node:
% aN(1) = 0.0;                                        (outside and inside)
% aP(1) = aP0(1) + g_b_ns(1) - Sp_1                   (outside)
% aP(1) = (aP0(1) + g_b_ns(1) - Sp_1)*gvP + Lf*dz/dt  (inside)
% aS(1) = -g_b_ns(1)*gvS                              (outside and inside)

% note: for these, I am using aS(1) as in aS(1) = -g_b_ns(1) i.e. b4 *gvS
% b(1) = aP0(1)*T_old(1) + Sc(1)*dz + aS(1)*gkS - aP0(1)*gkP - aS(1)*gkP

% note: n-1 = j+1 in the sense that n is the surface layer, which she shows
% as a layer, and n-1 is the next layer down, which is j+1 in the sense
% that j has an upper and lower neighbor for defining the equations,
