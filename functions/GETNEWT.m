%--------------------------------------------------------------------------
%   Repeat HEATSOLVE, dealing with latent heat
%--------------------------------------------------------------------------
function [  T_old,                                                      ...
            aP0_old ]   =   GETNEWT(Tsfc,T_old,T_tmp,ro_sno,Cp_sno,     ...
                            aP0_old,Sc,Sp,frac_ice,frac_liq,frac_air,   ...
                            ro_ice,xk_liq,xLs,Rv,Tfp,dt,dy_p,dely_p,    ...
                            f_n,JJ,opts)
%--------------------------------------------------------------------------

% copied this from v9

% recompute the temperature profile if the ice is melting
    xTsfc       =   Tsfc;
   %xTsfc       =   T_old(1);

% keep a copy of the T profile, use it to force Tf where water exists
    T_old0      =   T_old;

% set the T profile back to its value before solar radiation penetration
    T_old       =   T_tmp;
    
% update the vapor heat enthalpy term
    d_ro_vap_dT =   VAPOR_HEAT(T_old,Tfp,Rv,xLs,frac_liq,frac_air);
% NOTE: based on Eq. 131, this shouldn't be updated each iteration, and
% neither should frac_air, Cp_sno, or Fhatbar (dFliq/dT), which means i can
% move the enthalpy coefficient out of the while loop. i tested and it very
% slightly decreases runoff when its moved out like this, but there's also
% some increased instability in 2016 during the first few days 

% compute the enthalpy coefficient for each c.v. for the current timestep
    aP0         =   (ro_sno.*Cp_sno+xLs.*frac_air.*d_ro_vap_dT).*dy_p./dt;

% force the source terms to produce Tf at the positions with water,
%   including ice that was warmed to Tf by refreezing. Let the solar
%   radiation exist in other regions.
    Sc(T_old0==Tfp)  =   1e30*Tfp;
    Sp(T_old0==Tfp)  =   -1e30;    

% upper boundary condition
    T_N         =   xTsfc;

% the heat equation is nonlinear, iterate to solve it (p. 47)
    tol         =   1e-1; T_dif = 2*tol; iter = 0; maxiter = 10;
    
while any(T_dif > tol) && iter < maxiter    
    
% update the heat transfer coefficients.
    xk_eff      =   GETGAMMA(T_old,frac_ice,frac_liq,frac_air,xk_liq,   ...
                    ro_ice,ro_sno,xLs,Rv,Tfp,opts);

% update the vapor heat enthalpy term
%     d_ro_vap_dT =   VAPOR_HEAT(T_old,Tfp,Rv,xLs,frac_liq,frac_air);
    
% % compute the enthalpy coefficient for each c.v. for the current timestep
%     aP0         =   (ro_sno.*Cp_sno+xLs.*frac_air.*d_ro_vap_dT).*dy_p./dt;

% compute gamma at the c.v. interfaces (eq. 4.9, p. 45) (JJ+1)
    g_ns        =   [xk_eff(1) xk_eff(1:JJ) xk_eff(JJ)];
    g_b_ns      =   1.0./( (1.0-f_n)./g_ns(1:JJ+1)+f_n./g_ns(2:JJ+2));
    
% compute the lower (aN) and upper (aS) diagonal coefficients
    aN          =   g_b_ns(1:JJ)   ./ dely_p(1:JJ);
	aS          =   g_b_ns(2:JJ+1) ./ dely_p(2:JJ+1);
    
% Account for the boundary conditions.    
   %T_N         =   xTsfc;
    bc_N        =   aN(1) * T_N;
   %bc_S        =   0.0;
    aS(JJ)      =   0.0;  
    
% compute the middle (aP) coefficients and solution vector b
    aP          =   aN(1:JJ) + aS(1:JJ) + aP0(1:JJ) - Sp(1:JJ) .* dy_p(1:JJ);
    b           =   Sc(1:JJ) .* dy_p(1:JJ) + aP0_old(1:JJ) .* T_old(1:JJ);

% modify b to account for Dirichlet boundary conditions
    b(1)        =   b(1) + bc_N;
   %b(JJ)       =   b(JJ) + bc_S;
    
% solve the equation
    x           =   TRISOLVE(-aN,aP,-aS,b);

% prep for next iteration
    T_dif       =   abs(T_old-x);
    T_old       =   x;
    iter        =   iter+1;
end

% keep a copy of the enthalpy term for the next timestep
    aP0_old     =   (ro_sno.*Cp_sno+xLs.*frac_air.*d_ro_vap_dT).*dy_p./dt;
    
% % for water infiltration, would be:
%     Cp          =   ro_sno.*Cp_sno;
%     Cv          =   xLs.*d_ro_vap_dT.*frac_air;
%     Cu          =   xLf.*ro_sno.*dFdT;
%     aP0_old     =   (Cp+Cv+Cu).*dy_p./dt;
% % but this only accounts for aP0, b also needs to be modified
    
% note that T_tmp has changes in cold content in layers where all water
% froze and the temperature was reduced to below freezing. Otherwise,
% layers that experienced melt/freeze are forced to Tfp by the Sc/Sp=1e30,
% and all other layers are set to the temperature before solar radiation
% and melt/freeze was computed. Therefore, the iterations below, performed
% on T_old=T_tmp, account for the transfer of all heat produced and
% consumed during this timestep, and hte

% the key thing here is that the coefficients for the areas where T!=Tfp
% are for the T profile before solar radiation, since we are starting over
% with that profile, but we are also using T_tmp, which has changes in cold
% content due to solar radiation ... since we starting over, and allowing
% the original solar radiation to exist except where we already know T went
% to Tfp, do we need the cold content? I think yes, since what we're doing
% is putting the melt/freeze heat into the cv's, and now we're letting the
% heat from solar radiation and that melt/freeze heat transfer within the
% ice until convergence is achieved, so that means we can use the
% T_old=T_tmp to compute the coefficients on the first iteration, and then
% let them get updated, but do we want the new frac's, and cp_snow, and
% ro_snow? those will only change where melt/freeze occured, and if those
% areas go to Tfp, those coefficients won't affect the solution, so only
% areas where melt/freeze occured but T did not go to Tfp will have new
% coefficients, and those areas are also changed in T_tmp, so I think yes,
% we want the new coefficients, but its unclear what ap0_old should be    
    