function [  H,                                                          ...
            T,                                                          ...
            errH,                                                       ...
            errT,                                                       ...
            OK    ]  =  HEATSOLVE(T_old,H_old,f_wat,f_ice,f_liq,k_liq,  ...
                        cv_ice,cv_liq,ro_ice,ro_liq,ro_sno,cp_liq,cp_sno,...
                        dFdT,Ls,Lf,roLf,Rv,Tf,dz,delz,fn,dzdt,JJ_therm, ...
                        Tsfc,Sc,fcp,TL,TH,flmin,flmax,fcpsq,ro_iwe,ro_wie)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
   
% NOTE: still need to figrue out why i had the version w/o H_vap just in
% case but must have just been a test
% not sure why Hvap isn't in here ... actuallly i think it may be because
% my porosity doesn't change so there's always .1 air which is unrealistic  
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

% key notes: flmin/max must be at least 0.95/1.05, tol must be >=5e-2 (not
% 1e-2), maxiter might need to be capped at 5, or some combination of the
% aforementioned, otherwise Qc blows up. Might need to figure out if some
% combination, such as 1d-2 and maxiter>10, works such that Qc doesn't need
% to be included in Tsfc calculation

% these are the min/max controls in case i want to go back
% % f_liq(iM)   =   max(0.0,min(1.0,T(iM)./ro_liq + f_liq(iM)));
%   f_liq(iM)   =   min(1.0,T(iM)./ro_liq + f_liq(iM));
%   f_ice(iM)   =   max(0.0,(f_wat(iM)-f_liq(iM).*ro_liq./ro_ice));
  
% this is how she says she updates it (eq. 133a)
%   T(iM)       =   Tf - 1/fcp.*sqrt(ro_sno(iM)./(f_liq(iM).*ro_liq)-1);

% but this is the direct translation of eq 122a
%   T(iM)       =   T(iM)./ro_sno(iM)./dFdT(iM) + gkP(iM);
%   T(iM)       =   T(iM).*gvP(iM) + gkP(iM);

% and this is how it's coded in ftemp.f (line 79 and 117)
%   g_liq(iM)   =   T(iM)+g_liq_old(iM)*dz(iM);
%   T(iM)       =   Tf-sqrt(ratio-1.0)/a1;
% which is identical to what i have
    
% check if the phase boundary is overshot by more than 5%    
   %if any( f_liq(iM) < 0.950*flmin) || any(f_liq(iM) > 1.050*flmax)
%     if  ( T(iM)<T_old(iM) & any( f_liq(iM) < 1e-3) )  |      ...
%         ( T(iM)>T_old(iM) & any( f_liq(iM) > 1.001) ) 
%         asflag = true;
%         return;
%     end
% the reason this doesn't work is b/c if the prior step was f_liq==1e-3 or
% whatever is used, then the next step gets stuck trying to get f_liq
% smaller but that's impossible unless dt=0.0, so i'm not sure how she got
% around that, I think by just setting f_liq to the small value and
% continuing




% THIS IS HOW I HAVE IT IN v9e/HEATSOLVE_MZ
% x=Pmelt in the melt zone, so those layers need to be transformed to T 
% this method first updates f_liq, then uses f_liq to update x(iM), but
% note that f_ice must also be updated, but i can probably get around it
% with some algebra involving f_wat
%     f_wat       =   f_liq(iM)+f_ice(iM).*ro_ice./ro_liq;
%     f_liq(iM)   =   min(1.0,x(iM)./ro_liq + f_liq(iM));
%     f_ice(iM)   =   max(0.0,(f_wat-f_liq(iM)).*ro_liq./ro_ice);
%     x(iM)        =  Tf-sqrt(f_ice(iM).*ro_ice./f_liq(iM)./ro_liq)./fcp;
    
% note: i think the method above is correct but not right b/c it masks the
% situation where a phase boundary is overshot, instead, i think i need to
% use eq. 122a, which will return T>Tf so i can check overshoot correctly
% T(iM) = T(iM).*gvP(iM) + gkP(iM);

% this is how she says she updates it (eq. 133a)
%   x(iM)       =   Tf - 1/fcp.*sqrt(ro_sno(iM)./(f_liq(iM).*ro_liq)-1);

% but this is the direct translation of 
%   x(iM)       =   x(iM)./ro_sno(iM)./dFdT(iM) + gkP(iM);


% after tridiag is finished, this is what she does:
%     flmin     =   fliquid(f_wat(iM),tdl(m,1),tdl13(m,1),flglim(m,1))
%     flmax     =   fliquid(f_wat(iM),tdl(m,2),tdl13(m,2),flglim(m,2))
%     [f_liq,T] =   ftemp(T(iM),dt,f_liq(iM),f_wat(iM),dz(iM),TL(iM),TH(iM),...
%                   f_liq_old(i),bmelt(i),flgo(i),flmin,flmax,dtmin,dFdT);
            
% since her model has water fluxes, flimin/max need to be updated b/c they
% are functions of total water content, so i can use the precomputed
% values. ftemp returns liquid water content and temperature BUT ALSO the
% overshoot condition, so 
     
% flglim are the unfrozen water fraction values 'flg' at phase boundaries
% MELTZONE determines temperature limits th and tl for the meltzone,
% and computes the corresponding fractional water content limits
% flgliml and flglimh for the capillary portion of the freezing curv

% This was in HEATSOLVE_CAPP, which I deleted
% % for water infiltration, would be:
%     Cp          =   ro_sno.*Cp_sno;
%     Cv          =   xLs.*d_ro_vap_dT.*frac_air;
%     Cu          =   xLf.*ro_sno.*dFdT;
%     aP0_old     =   (Cp+Cv+Cu).*dy_p./dt;
% % but this only accounts for aP0, b also needs to be modified
