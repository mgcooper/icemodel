function [aN, aP, aS, b, iM, a1, a2, aP01] = GECOEFS(T, f_ice, f_liq, dHdT, ...
      dFdT, drovdT, dH, Sc, k_eff, delz, fn, dz, dt, Ts, Ls, Lf, ro_liq, TL, ...
      JJ, Fc, Fp, bc)
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
   aP0 = (dHdT + Lf * ro_liq * dFdT + Ls * f_air .* drovdT) .* dz / dt;
   gv = ones(JJ, 1);     % Eq 123
   gk = zeros(JJ, 1);    % Eq 123
   LfMZ = zeros(JJ, 1);  % Eq 123, melt-zone latent heat switch

   aP01 = aP0(1);

   % % If using g_liq instead of f_liq in the definition of dFdT as in SNTHERM:
   % aP0 = (dHdT + Lf * ro_sno .* dFdT + Ls * f_air .* drovdT)

   % Cofficients for nodes inside the melt zone [W m-2 K-1]
   if sum(iM) > 0
      aP0(iM) = (dHdT(iM) + Ls * f_air(iM) .* drovdT(iM)) .* dz(iM) / dt;
      gv(iM) = 1 ./ (ro_liq * dFdT(iM)); % Eq 122b
      gk(iM) = T(iM);
      LfMZ(iM) = Lf * dz(iM) / dt;

      % % If using g_liq instead of f_liq as in SNTHERM:
      % gv(iM) = 1 ./ (ro_sno(iM) .* dFdT(iM)); % Eq 122b
      % gv(iM) = 1 ./ (ro_liq * f_wat(iM) .* dFdT(iM)); % Eq 122b
   end

   % % For a soil model, would need indices above the melt zone
   % if sum(iL)>0 % nodes above the melt zone:
   %    aP0(iL) = cv_liq * dz / dt;
   %  % gv(iM) = 1.0;
   %  % gk(iM) = 0.0;
   %  % dFdT(iL) = 0.0;
   % end

   % Adjust gv/gk in terms of N/P/S. Set gkN(1) = 0 and gkS(JJ+1) = 0 to account
   % for the upper left and lower right off-diagonal BCs. Set gvS(JJ+1) = 0 to
   % account for the zero-flux lower BC. Set gvN(1) = 1 to account for the
   % non-zero upper BC (for both Dirichlet and Robin).
   gvN = vertcat(1, gv(1:JJ-1));
   gvS = vertcat(gv(2:JJ), 0);
   gkN = vertcat(0, gk(1:JJ-1));
   gkS = vertcat(gk(2:JJ), 0);
   % gkP = gk;
   % gvP = gv;

   % Compute the aN and aS conductances [W m-2 K-1]. Note that delz(1) &
   % delz(end) are 1/2 CVs.
   aN = g_b_ns(1:JJ)   ./ delz(1:JJ);
   aS = g_b_ns(2:JJ+1) ./ delz(2:JJ+1);

   % Keep the upper left and right conductances
   a1 = aN(1); % Note: astar = a1 / (a1 - Fp);
   a2 = aS(1);

   % Account for the upper boundary condition
   switch bc
      case 1
         % Dirichlet: Ts = known
         bc_N = a1 * Ts;
      case 2
         % Robin: qB = f(Ts)
         bc_N = a1 * Fc / (a1 - Fp);
         aN(1) = 0.0;
      case 3
         % Neumann: qB = known
         % bc_N = qN;
         % aN(1) = 0.0;
   end

   % Account for the lower boundary condition
   aS(JJ) = 0.0;                       % Neumann: dT/dz = 0.0
   bc_S = aS(JJ) * 0.0;                % Neumann: dT/dz = 0.0

   % Compute the aP coefficient and solution vector b
   aP = aN + aS + aP0; % -Sp.*dz;
   b = aP0 .* T + Sc .* dz - dH .* dz / dt; % [W m-2]

   % Modify b to account for boundary conditions
   b(1) = b(1) + bc_N;
   b(JJ) = b(JJ) + bc_S;

   % Robin:
   if bc == 2 || bc == 3
      aP(1) = aP(1) - Fp * a1 / (a1 - Fp);
   end

   % Apply the melt zone (enthalpy) transformation (Eq. 128/29)
   % b = b + aN .* gkN - (aN + aS + aP0) .* gk + aS .* gkS ;
   b = b + aN .* gkN - aP .* gk + aS .* gkS ;
   aN = aN .* gvN;
   aS = aS .* gvS;
   aP = aP .* gv + LfMZ;
end

% NOTES:

% n = surface layer
% j+1    = N   (n-1)
% j+1/2  = N-P interface
% j      = P
% j-1/2  = P-S interface
% j-1    = S

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

% Qnet are the past net heat fluxes which include the conductive flux into the
% surface layer

% Qnet_old(1) = chi * Qsi * (1 - albedo) + cp_liq * U_liq * T_old(1)  ...
%    - aS(1) * (T_old(1) - T_old(2));

% Her equation:
% Qnet_old(1) = Itop_old + cliq * Uliq * T_old(1) - ...
%    - ke_nm_1o2
