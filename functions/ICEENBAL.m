function [  T,                                                          ...
            errH,                                                       ...
            errT,                                                       ...
            f_ice,                                                      ...
            f_liq,                                                      ...
            OK    ]  =  ICEENBAL(T,f_ice,f_liq,k_liq,cv_ice,cv_liq,     ...
                        ro_ice,ro_liq,ro_sno,cp_sno,Ls,Lf,roLf,Rv,      ...
                        Tf,dz,delz,fn,dt,JJ,Tsfc,Sc,fcp,TL,TH,flmin,    ...
                        flmax,ro_iwe,ro_wie) %#codegen

% % liquid water fraction and the derivative of f_liq wrt to temperature
   f_wat    =  f_liq+f_ice.*ro_iwe;           % frac_wat_old
   dFdT     =  2*fcp^2.*f_wat.*(Tf-min(T,Tf))./(1+fcp^2.*(Tf-min(T,Tf)).^2).^2;
    
% update the melt-zone boundaries
   fliqmin  =  f_wat.*flmin;
   fliqmax  =  f_wat.*flmax;
   
% vapor diffusion heat
   H_vap    =  VAPORHEAT(T,f_liq,f_ice,Tf,Rv,Ls);

% total enthalpy (in W/m2 for comparison with fluxes)
   H_old    = ((cv_ice.*f_ice+cv_liq.*f_liq).*(T-Tf)+roLf.*f_liq+H_vap).*dz./dt;
                     
% iterate to solve the nonlinear heat equation (p. 47)
   errT = 2e-2; iter = 0; errH = 0; OK = true;

while any(errT > 1e-2) && iter < 50
   
   T_old    =  T;
   iter     =  iter+1;
    
% update thermal conductivity
   k_eff    =  GETGAMMA(T,f_liq,f_ice,ro_ice,k_liq,Ls,Rv,Tf);
   
% update vapor heat   
[  H_vap,                                                                ...
   drovdT ] =  VAPORHEAT(T,f_liq,f_ice,Tf,Rv,Ls);

% update total enthalpy
   H        =  ((cv_ice.*f_ice+cv_liq.*f_liq).*(T-Tf)+roLf.*f_liq+H_vap).*dz./dt;
                 
% update the general model coefficients
[  aN,aP,                                                               ...
   aS,b,                                                                ...
   iM  ]    =  GECOEFS(T,ro_sno,cp_sno,f_liq,f_ice,Ls,Lf,dz,dt,         ...
               dFdT,drovdT,TL,H,H_old,Sc,k_eff,fn,delz,Tsfc,JJ);

% underrelaxation (doesn't work with melt-zone transformation)
% if iter>1
%     aP          =   aP./alpha;
%     b           =   b + beta.*aP.*T;
% end

% solve the equation
   T        =  TRISOLVE(-aN,aP,-aS,b);

% back-transform the meltzone preconditioner
[  T,                                                                   ...
   f_liq,                                                               ...
   f_ice,                                                               ...
   OK ]     =  MZTRANSFORM(T,T_old,f_liq,f_wat,ro_liq,ro_wie,Tf,TL,TH,  ...
               fcp,fliqmin,fliqmax,iM,OK);

% if failure, return to the main program and adjust the timestep
   if OK == false; return; end
               
% when a node first enters the melt zone, it will not be tagged by iM at
% the start of this loop, so this checks for that and updates f_ice/liq/T
[  f_liq,                                                               ...
   f_ice,                                                               ...
   ~,                                                                   ...
   T  ]     =  MELTCURVE(T,f_liq,f_ice,ro_wie,ro_iwe,Tf,fcp);
 
% prep for next iteration
   errT     =  abs(T_old-T);

% compute residual (need to add Sp)
%    resN   =  aN.*([TN T(1:end-1)]-T);
%    resS   =  aS.*([T(2:end) 0]-T);
%    res    =  (resN+resS+Sc.*dz)./aP0;

end

% keep track of the residual
   errH  =  H-H_old;
   
% % check Tsfc
%    Tsfc  =  TSURF(T,ro_sno,cp_sno,f_liq,f_ice,Ls,Lf,dz,dt,dFdT,   ...
%             drovdT,H,H_old,Sc,k_eff,fn,delz,JJ);
   
   
   
   
   