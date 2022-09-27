%------------------------------------------------------------------------------
% get the general equation coefficients
%------------------------------------------------------------------------------
function [aN,aP,aS,b,iM] = GECOEFS(T,ro_sno,cp_sno,f_liq,f_ice,Ls,Lf,   ...
                           dz,dt,dFdT,drovdT,TL,H,H_old,Sc,k_eff,fn,    ...
                           delz,Tsfc,JJ)
% was H1d_GETCOEFS_MZ

% melt zone indices
   iM       = TL <= T;
   
% for a soil model, would need indices above the melt zone   
 % iM       = TL <= T & T <= TH;
 % iL       = T>TH;

% coefficients for nodes below the melt zone: (W/m2/K)
   f_air    = (1.0-f_ice-f_liq);
   aP0      = (ro_sno.*cp_sno+Lf.*ro_sno.*dFdT+Ls.*f_air.*drovdT).*dz./dt;
   gv       = ones(JJ,1);                                % Eq 123
   gk       = zeros(JJ,1);
   LfMZ     = zeros(JJ,1);       % for melt-zone latent heat switch

% commented statements are for reference, they are set when initialized

if sum(iM)>0   % nodes inside the melt zone:
   aP0(iM)  = (ro_sno(iM).*cp_sno(iM)+Ls.*f_air(iM).*drovdT(iM)).*dz(iM)./dt;
   gv(iM)   = 1./dFdT(iM)./ro_sno(iM);                   % Eq 122b
 % gv(iM)   = 1./dFdT(iM)./f_wat(iM)./ro_liq;            % Eq 122b
   gk(iM)   = T(iM);
   LfMZ(iM) = Lf.*dz(iM)./dt;
end

% % for a soil model, would need indices above the melt zone
% if sum(iL)>0   % nodes above the melt zone:
%    aP0(iL)  = ro_sno(iL).*cp_liq.*dz(iL)./dt;
%  % gv(iM)   = 1.0;
%  % gk(iM)   = 0.0;
%  % dFdT(iL) = 0.0;
% end

% adjust gv/gk in terms of N/P/S
   gvN      =  vertcat(1,gv(1:JJ-1));
   gvS      =  vertcat(gv(2:JJ),0);
   gkN      =  vertcat(0,gk(1:JJ-1));
   gkS      =  vertcat(gk(2:JJ),0);
 % gkP      =  gk;
 % gvP      =  gv;
 
% compute gamma at the control volume interfaces (eq. 4.9, p. 45) (JJ+1)
   g_b_ns   =  1./((1-fn)./[k_eff(1);k_eff]+fn./[k_eff;k_eff(JJ)]); 

% compute the aN and aS coefficients
   aN       =  g_b_ns(1:JJ)   ./ delz(1:JJ);
   aS       =  g_b_ns(2:JJ+1) ./ delz(2:JJ+1);
% note that dely_p(1) and dely_p(end) are 1/2 CVs, which is correct

% Account for the boundary conditions.
   %bc_N    =  aN(1) * TN;      % Dirichlet: TN = known
   %bc_S    =  0.0;             % Neumann: dT/dz = 0.0
   aS(JJ)   =  0.0;             % Neumann: dT/dz = 0.0
   %aN(1)   =  0.0;             % Neumann: qB = known
   
% compute the aP coefficient and solution vector b
   aP       =  aN + aS + aP0; % -Sp.*dz;
   b        =  aP0.*T + Sc.*dz - (H - H_old);

% modify b to account for Dirichlet boundary conditions
   b(1)     =   b(1) + aN(1)*Tsfc;
    
% adjust coefficients for melt zone switches
   b        =   b - aP0.*gk+aS.*gkS-aS.*gk+aN.*gkN-aN.*gk;
   aN       =   aN.*gvN;
   aS       =   aS.*gvS;
   aP       =   aP.*gv + LfMZ;
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %   reference:

% n     = surface interface
% j+1   = N   (n-1)
% j+1/2 = N-P interface
% j     = P
% j-1/2 - P-S interface
% j-1   = S

% A3    = aN
% A2    = aP
% A1    = aS
% Qs    = aP0

% switches, interior nodes
% aN    = aN*gvN;                       (inside and outside the meltzone)
% aS    = aS*gvS;                       (inside and outside the meltzone)
% aP    = (aP0+aN+aS)*gvP + Lf*dz/dt;   (inside)
% aP    = (aP0+aN+aS);                  (outside)
% b     = aP0*T_old + aS*gkS - aP0*gkP - aN*gkP - aS*gkP + aN*gkN + Sc*dx;
%                                       (inside and out)

% aP0   = (ro_sno*cp_sno + Lf*ro_sno*Fbar + Lv*f_air*CkT)*dz/dt (outside)
% aP0   = (ro_sno*cp_sno +      0         + Lv*f_air*CkT)*dz/dt (inside)

% gv    = 1./ro_sno./Fbar;  (inside)
% gk    = T_old             (inside)
% gv    = 1;                (outside)
% gv    = 0;                (outside)

% T_old = T_old;            (outside)
% T_old = P_melt;           (inside) (convert to get T_old)

% top node:
% aN(1) = 0.0;                                        (outside and inside)
% aP(1) = aP0(1) + g_b_ns(1) - Sp_1                   (outside)
% aP(1) = (aP0(1) + g_b_ns(1) - Sp_1)*gvP + Lf*dz/dt  (inside)
% aS(1) = -g_b_ns(1)*gvS                              (outside and inside)

% note: for these, I am using aS(1) as in aS(1)=-g_b_ns(1) i.e. b4 *gvS
% b(1)  = aP0(1)*T_old(1) + Sc(1)*dz + aS(1)*gkS - aP0(1)*gkP - aS(1)*gkP 

% note: n-1 = j+1 in the sense that n is the surface layer, which she shows
% as a layer, and n-1 is the next layer down, which is j+1 in the sense
% that j has an upper and lower neighbor for defining the equations, 



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %   older notes
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
% below was my intuition, but it's wrong, because gvN(1) is the layer above
% the top layer, which is always outside the melt zone, so gvN=1 for the
% top layer, and gvS for the top layer is gv of layer 2, so gvS(1) = gv(2),
% similarly for gkN/gkS, except gvS = gv(j+1) is zero for the top layer
% (j=1). These are set up this way so they can be multiplied by aN/aS/aP in
% heatsolve just like g_b_ns is done
%     gvN         = [gv(1:end-1) 0];
%     gvP         = gv;
%     gvS         = [0 gv(2:end)];
%     
%     gkN         = [gk(1:end-1) 0];
%     gkP         = gk;
%     gkS         = [0 gk(2:end)];

    
% the net change in stored heat is explained by three processes: 1) change
% in the ice+liq+air mixture temperature, 2) change in the liquid mass
% fraction, and 3) changes in the water vapor mass fraction. If the mixture
% temperature is below the freezing depression point, then only 1) and 3)
% can occur, so the enthalpy coefficients only contain terms for change in
% specific heat capacity and change in water vapor heat capacity. If the
% mixture temperature is within the "melt zone", there is an extra term for
% change in latent heat stored in liquid water. One question that requires
% answering is "how much of the net change in stored heat goes into the
% latent heat term vs the other two terms?". The answer is given by the
% "freezing curve" which is the derivative of liquid water mass fraction
% wrt temperature. Similarly, the amount of stored heat explained by water
% vapor is provided by the clausius clapeyron relationship, which describes
% the change in water vapor mass fraction wrt temperature.

% NOTE: I confirmed that in the melt zone, aP0 does not have the Lf*dbldT
% term, which can be confirmed in thparam_ice. This seems strange, but I
% think the idea is that once we get to the melt zone, the latent heat is
% used up, and we use the 'extra' energy to change phase just like in my
% original method, and maybe the whole idea is to have a constraint on the
% temperarture zone so dbldT works, since in my method, the temperature can
% go anywhere, and in capp, we know the temperature changed according to
% dbldT, and so once in the melt zone, we can convert direclty to latent
% heat via fliq. 

% This is how Jordan does it, Fbar is in terms of gamma:
% P_melt = gam_wat*Fbar*delT        = gam_liq-gam_liq_old;
% This is in terms of theta, Fbar is in terms of theta
% P_melt = f_wat*Fbar_liq*delT        = gam_liq-gam_liq_old;
% P_melt = ro_liq*Fbar*delT   = f_liq-f_liq_old;
% P_melt = -ro_liq*d_f_ice_dT*delT  = -(f_ice-f_ice_old);

% where Fbar = d_f_liq_dT or d_f_ice_dT for the liq/ice version s