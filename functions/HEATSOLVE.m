function [  H,                                                             ...
            T,                                                             ...
            errH,                                                          ...
            errT,                                                          ...
            OK    ]  =  HEATSOLVE(T_old,H_old,f_wat,f_ice,f_liq,k_liq,     ...
                        cv_ice,cv_liq,ro_ice,ro_liq,ro_sno,cp_liq,cp_sno,  ...
                        dFdT,Ls,Lf,roLf,Rv,Tf,dz,delz,fn,dzdt,JJ_therm,    ...
                        Tsfc,Sc,fcp,TL,TH,flmin,flmax,fcpsq,ro_iwe,ro_wie)

%------------------------------------------------------------------------------

% iterate to solve the nonlinear heat equation (p. 47)
%     tol = 1e-2; T_dif = 2*tol; iter = 0; maxiter = 50;

% just for comparison with pre-combined speed test
   tol = 5e-2; errT = 2*tol; iter = 0; maxiter = 5; errH = 0; OK = true;

while any(errT > tol) && iter < maxiter
    
% update thermal conductivity
   k_eff       =  GETGAMMA(T_old,f_liq,f_ice,ro_ice,k_liq,Ls,Rv,Tf);
   
% update vapor heat   
[  H_vap,                                                                ...
   drovdT ]    =  VAPORHEAT(T_old,f_liq,f_ice,Tf,Rv,Ls);

% update total enthalpy
   H           =  TOTALHEAT(T_old,f_liq,f_ice,cv_liq,cv_ice,roLf,       ...
                  Tf,H_vap,dzdt);
                 
% update the general model coefficients
[  aN,aP,                                                               ...
   aS,b,                                                                ...
   iM  ]       =  GECOEFS(T_old,ro_sno,cp_sno,f_liq,f_ice,cp_liq,       ...
                  Ls,Lf,dz,dzdt,dFdT,drovdT,TL,TH,H,H_old,Sc,k_eff,     ...
                  fn,delz,Tsfc,JJ_therm);

   T           =  TRISOLVE(-aN,aP,-aS,b);

% back-transform the meltzone preconditioner
[  T,f_liq,                                                             ...
   OK ]        =  MZTRANSFORM(T,T_old,f_liq,f_wat,ro_liq,Tf,            ...
                  TL,TH,fcp,flmin,flmax,iM);

% if failure, return to the main program and adjust the timestep
   if OK == false; return; end
               
% when a node first enters the melt zone, it will not be tagged by iM at
% the start of this loop, so this checks for that and updates everything
[  f_liq,                                                               ...
   f_ice,                                                               ...
   ~,                                                                   ...
   T     ]     =  MELTCURVE(T,f_liq,f_ice,ro_wie,ro_iwe,Tf,fcp,fcpsq);
 
% prep for next iteration
   errT         =  abs(T_old-T);
   T_old       =  T;
   iter        =  iter+1;

% compute residual (need to add Sp)
%     resN        = aN.*([TN T(1:end-1)]-T);
%     resS        = aS.*([T(2:end) 0]-T);
%     res         = (resN+resS+Sc.*dz)./aP0;

end

% keep track of the residual
   errH = H-H_old;
   
% % sparse testing:
%    A           =  sparse(diag(-aN(2:end),-1)+diag(-aS(1:end-1),1)+diag(aP,0));
% 
% % [  A,b,                                                                 ...
% %    iM  ]       =  GECOEFS2(T_old,ro_sno,cp_sno,f_liq,f_ice,cp_liq,      ...
% %                   Ls,Lf,dz,dzdt,dFdT,drovdT,TL,TH,H,H_old,Sc,k_eff,     ...
% %                   fn,delz,Tsfc,JJ,A);
% 
% % solve the system
%    T           =  A \ b;
%    tic; T = sparse(diag(-aN(2:end),-1)+diag(-aS(1:end-1),1)+diag(aP,0))\b; toc;
%    tic; T = spdiags([-aS,aP,-aN],-1:1,JJ,JJ) \ b; toc;
%    tic; T = TRISOLVE(-aN,aP,-aS,b); toc;
%    tic; T = (diag(-aN(2:end),-1)+diag(-aS(1:end-1),1)+diag(aP,0))\b;toc;
%    
% %    T           =  spdiags([-aS,aP,-aN],-1:1,JJ,JJ) \ b;

% % total enthalpy
%     H0          =   ((ro_ice.*cp_ice.*f_ice + ro_liq.*cp_liq.*f_liq).*  ...
%                     (T_old-Tf) +  ro_ice.*Lf.*f_liq).*dz./dt;
% just keeping this version for now 
%     H0          =   ((ro_ice.*cp_ice.*f_ice + ro_liq.*cp_liq.*f_liq).*  ...
%                     (T_old-Tf)+ro_ice.*Lf.*f_liq+H_vap).*dz./dt;

% TN    = T_old(1);

%     % source term relaxation
%     Sc_prime    =   alpha.*Sc;
%     Sc_new      =   Sc_old + Sc_prime;
%     
%     Sc_use      =   alpha.*Sc_new + (1-alpha).*Sc_old;
%     
%     Sc_new      =   Sc_old + alpha.*Sc;
%     φnew,used = φold +λ(φnew,predicted −φold)


% % for water infiltration, would be:
%     Cp          =   ro_sno.*Cp_sno;
%     Cv          =   xLs.*d_ro_vap_dT.*frac_air;
%     Cu          =   xLf.*ro_sno.*dFdT;
%     aP0_old     =   (Cp+Cv+Cu).*dy_p./dt;
% % but this only accounts for aP0, b also needs to be modified
