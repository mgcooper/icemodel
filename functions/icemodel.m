
function [enbal,ice2,met,opts] = icemodel(opts)
   
%------------------------------------------------------------------------------
%   INITIALIZE THE MODEL
%------------------------------------------------------------------------------
   
% load the physical constants to be used.
   load( 'PHYSCONS', 'cp_ice','cv_air','cv_liq','cv_ice','emiss',       ...
                     'SB','epsilon','fcp','k_liq','Lf','Ls','Lv',       ...
                     'ro_air','ro_ice','ro_liq','roLf','roLs','roLv',   ...
                     'Rv','Tf','TINY');
                  
% open the input meteorological data file
   [met,opts]  = METINIT(opts);
   
   if strcmp(opts.error,'metfile does not exist')
      error(opts.error);
   end


%  INITIALIZE THE ICE COLUMN
[  f_ice,                                                               ...
   f_liq,                                                               ...
   T,                                                                   ...
   TL,                                                                  ...
   TH,                                                                  ...
   flmin,                                                               ...
   flmax,                                                               ...
   cp_sno,                                                              ...
   k_eff,                                                               ...
   dz,                                                                  ...
   fn,                                                                  ...
   delz,                                                                ...
   grid_therm,                                                          ...
   dz_therm,                                                            ...
   dz_spect,                                                            ...
   JJ_therm,                                                            ...
   JJ_spect,                                                            ...
   Sc,                                                                  ...
   Sp,                                                                  ...
   wcoef,                                                               ...
   scoef,                                                               ...
   ro_sno,                                                              ...
   ro_iwe,                                                              ...
   ro_wie,                                                              ...
   enbal,                                                               ...
   ice2  ]     =   ICEINIT(opts,met);

%------------------------------------------------------------------------------
%   INITIALIZE THE EXTINCTION COEFFICIENTS
%------------------------------------------------------------------------------

% load the spectral mie values, solar spectrum, and abs. coefficients
[  radii,                                                            ...
   scattercoefs,                                                     ...
   solar,                                                            ...
   kabs,                                                             ...
   kice    ]         =     SPECINIT(opts);
   
% Initialize the extinction coefficients that control the solar radiation
%   source term within the snow/ice matrix.
[  total_solar,                                                         ...
   grid_spect,                                                          ...
   z_walls,                                                             ...
   spect_lower,                                                         ...
   spect_upper,                                                         ...
   solardwavl  ]     =     EXTCOEFSINIT(opts,radii,scattercoefs,        ...
                           solar,kabs,kice,dz_spect,JJ_spect,ro_ice);
%------------------------------------------------------------------------------
%   START THE MODEL
%------------------------------------------------------------------------------
   
   % extract struct values that are used frequently
   Tsfc        =  T(1);
   itime       =  met.Time(1);             % initialize model time
   
   nloops      =  opts.annual_loops;
   dt          =  opts.dt;
   fopts       =  opts.fzero;
   f_min       =  opts.f_ice_min;
%  qsfactor    =  opts.qsfactor;
   
%------------------------------------------------------------------------------
% substepping settings
%------------------------------------------------------------------------------
   maxiter     =  opts.maxiter;
   maxsubiter  =  200;
   minsubiter  =  1;
   dt_min      =  dt/maxsubiter;  % or: dt_min = 5; maxsubiter  = dt/dt_min;
   dt_max      =  dt/minsubiter;  % or: dt_min = 5; maxsubiter  = dt/dt_min;
   dt_new      =  dt_max;
   iter        =  1;
   subiter     =  1;

% we need these on the first iteration, they are updated after each substep
   xf_liq      =  f_liq;
   xTsfc       =  Tsfc;
   roL         =  roLs;
   liqflag     =  false;
   Qc          =  CONDUCT(k_eff,T,dz,xTsfc);
   
%  zD          =  sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1)));

for MM = 1:nloops
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
   
while iter <= maxiter
   
% reset the change in liq water content
   d_liq       =  0.0.*f_liq;
   d_drn       =  0.0.*f_liq;

%  load the atmospheric forcing data for this time step.
   Tair        =  met.tair(iter);
   rh          =  met.rh(iter);
   wspd        =  met.wspd(iter);
   Qsi         =  met.swd(iter);
   Qli         =  met.lwd(iter);
   Pa          =  met.psfc(iter);
   albedo      =  met.albedo(iter);
%    snowd       =  met.snowd(iter);
   De          =  wspd*wcoef;
   ea          =  VAPPRESS(rh,Tair,liqflag);
   
% update the upper layer diffusion length scale
%  zD          =  sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1)));
%    if zD > dz(1)
%       pause;
%    end   
   
%  Update the subsurface absorption profile and solar radiation source-term
if Qsi>0
[  Sc,                                                                  ...
   chi   ]     =  UPDATEEXTCOEFS(grid_spect,JJ_spect,dz_spect,          ...
                  grid_therm,JJ_therm,dz_therm,z_walls,dz,ro_sno,       ...
                  Qsi,total_solar,spect_lower,spect_upper,albedo,       ...
                  solardwavl);
else
   Sc          =  0.0.*Sc;
   chi         =  0.0;
end

% [Qsi*(1-albedo) sum(Sc.*dz)]                        % total
% [chi*Qsi*(1-albedo) Sc(1).*dz(1)]                   % top node
% [(1-chi)*Qsi*(1-albedo) sum(Sc(2:end).*dz(2:end))] 	% interior nodes

% if the Sc(1)=0 method is used, then this should hold:
% Qsi*(1-albedo) - ( chi*Qsi*(1-albedo) + sum(Sc.*dz) ) = 0
% (1-chi)*Qsi*(1.0-albedo) - Sc(1)*dz(1) = 0

% if using the Qseb method, the total absorbed radiation should equal the
% portion allocated to the surface plus the sum of the subsurface absorption
% [Qsi*(1-albedo) sum(Sc.*dz)+Qseb]                   % total
% [chi*Qsi*(1.0-albedo) Qseb]                         % SEB 
% [(1-chi)*Qsi*(1.0-albedo) sum(Sc.*dz)]              % interior

% so we can either pass Qseb or chi out of extcoefs 

%------------------------------------------------------------------------------
%   Initialize dynamic timestepping
%------------------------------------------------------------------------------
   dt_sum      =  0.0;
   subfail     =  0;                      % keep track of failed substeps
   dt_flag     =  false;                  %
%    LCflag      =  false(size(f_ice));     % layer combination flag
   OK          =  false;                  % assume sub-stepping is needed

   % could use this to get a better estimate of Tsfc if convergence fails
%    xTsfc       =  T(1) + SEB/(k_eff(1)/(dz/2)); % T1 + SEB/a1
   
   while OK == false || dt_sum < dt

% SURFACE TEMPERATURE
   [  Tsfc  ]  =  SFCTEMP(Tair,Qsi,Qli,ea,albedo,De,Pa,wspd,cv_air,     ...
                  emiss,SB,Tf,Qc,xTsfc,chi,roL,scoef,fopts,liqflag);

%       if Tflag == false
%          fprintf('iter = %d (%.2f%%), dt = %.0f, Tsfc not converged\n',            ...
%             iter,100*iter/maxiter,dt_new)
%       end

% SUBSURFACE ENERGY BALANCE
   [  T,                                                                ...
      errH,                                                             ...
      errT,                                                             ...
      f_ice,                                                            ...
      f_liq,                                                            ...
      OK    ]  =  ICEENBAL(T,f_ice,f_liq,k_liq,cv_ice,cv_liq,           ...
                  ro_ice,ro_liq,ro_sno,cp_sno,Ls,Lf,roLf,Rv,Tf,         ...
                  dz,delz,fn,dt_new,JJ_therm,Tsfc,Sc,fcp,TL,TH,         ...
                  flmin,flmax,ro_iwe,ro_wie);

% this will slow down the code a lot but is useful for debugging               
%     fprintf('iter = %d (%.2f%%), dt = %.0f, success = %s\n',            ...
%             iter,100*iter/maxiter,dt_new,mat2str(OK))

% if a phase boundary was overshot, decrease the time step and start over
      if OK == false && subfail < maxsubiter
         subfail     =  subfail+1;
         subiter     =  min(subiter+1,maxsubiter);
         dt_new      =  dt_max/subiter;
         T           =  xT;
         Tsfc        =  xTsfc;
         f_ice       =  xf_ice;
         f_liq       =  xf_liq;
         continue;
      else
         
      % need to update this for Qc
         k_eff       =  GETGAMMA(T,f_liq,f_ice,ro_ice,k_liq,Ls,Rv,Tf);
         
      % Substep was successful, recompute top flux and linearization error
      [  Qe,                                                            ...
         Qh,                                                            ...
         Qc,                                                            ...
         Qm,                                                            ...
         Qf,                                                            ...
         ea,                                                            ...
         balance ]   =  SEBFLUX(T,Tair,MELTTEMP(Tsfc,Tf),Qsi,Qli,albedo,...
                        wspd,rh,Pa,Tf,k_eff,cv_air,roL,emiss,SB,epsilon,...
                        scoef,De,dz,liqflag,chi);

% %    % faster option to send an updated Qe to ICEMF to compute subl:
%          ea       =  VAPPRESS(rh,Tair,Tsfc,Tf,liqflag);
%          S        =  STABLEFN(Tair,Tsfc,wspd,scoef);
%          es0      =  VAPOR(Tsfc,Tf,liqflag);
%          Qe       =  LATENT(De,S,ea,es0,roL,epsilon,Pa);

      % Substep was successful, compute fluxes and update state
      [  T,                                                             ...
         f_ice,                                                         ...
         f_liq,                                                         ...
         d_liq,                                                         ...
         d_drn ]  =  ICEMF(T,f_ice,f_liq,ro_ice,ro_liq,cp_ice,Lf,Ls,    ...
                     Lv,Tf,TL,fcp,xf_liq,Sc,Sp,JJ_therm,f_min,fopts,    ...
                     dz_therm,dt_new,Qe,liqflag,ro_iwe,d_liq,d_drn);
                        
      % update density (kg/m3), heat capacity (J/kg/K), thermal K (W/m/K)
       % k_eff    =  GETGAMMA(T,f_liq,f_ice,ro_ice,k_liq,Ls,Rv,Tf);
         ro_sno   =  f_ice.*ro_ice+f_liq.*ro_liq+(1.0-f_liq-f_ice).*ro_air;
         cp_sno   =  (cv_ice.*f_ice+cv_liq.*f_liq)./ro_sno;
         
%       % update the upper layer diffusion length scale
%          zD       =  sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1)));
%          if zD > dz(1)
%          end
         
      % allocate this step to the running total
         dt_sum   =  dt_sum + dt_new;
         itime    =  itime + seconds(dt_new);

      % this gets the last successful substep AND the last full-step
         xT       =  T;
         xTsfc    =  Tsfc;
         xf_liq   =  f_liq;
         xf_ice   =  f_ice;
         
      % if the top node contains >2% liquid water, use Lv for latent heat flux
         if f_liq(1)>0.02
            roL      = roLv;  % ro_air*Lv
            liqflag  = true; 
         else
            roL      = roLs;  % ro_air*Ls
            liqflag  = false;
         end
         % first is true if step incomplete, second if overage will occur
         if dt-dt_sum > TINY && (dt_sum+dt_new)-dt > TINY
            dt_new   =  max(dt-dt_sum,dt_min);
            dt_flag  =  true;
         end

      end
      
% NOTE: this can go outside of substepping unless we want it during subs
%       LCflag  =   logical(LCflag+lcflag);    % track layer combinations

   end % substepping
   
   % i think this check can be deleted b/c next check gets it
   if dt_flag == true 
      dt_new  =   dt_max/subiter;
   end

   % If convergence is accepted, try increasing the timestep if very small
   if OK == true
      subiter =   max(1,subiter-1);   % reverse the subiter decrement
      dt_new  =   dt_max/subiter;
   end

%------------------------------------------------------------------------------
%   Save the output
%------------------------------------------------------------------------------

% although i could compute the enbal after the model is finished, we need
% all components to compute the balance if used for error control, so we
% might as well save the components here rather than compute in POSTPROC.
% but should test to see if memory overhead exceeds compute.

   if MM==nloops
      
%       for a stripped-down run this is all that's needed:
%       df_liq(:,iter)    =  d_liq;
%       frac_ice(:,iter)  =  f_ice;         % fraction ice
%       T_ice(:,iter)     =  T;
%       T_srf(iter,1)     =  Tsfc;
      
      % save the surface energy balance
      enbal.Tsfc(iter,1)      =  Tsfc;           % surface temp
      enbal.Qh(iter,1)        =  Qh;             % sensible
      enbal.Qe(iter,1)        =  Qe;             % latent
      enbal.Qc(iter,1)        =  Qc;             % conductive into surface
      enbal.Qm(iter,1)        =  Qm;             % melt energy
      enbal.Qf(iter,1)        =  Qf;             % freeze deficit
      enbal.chi(iter,1)       =  chi;            % skin parameter
      enbal.balance(iter,1)   =  balance;        % SEB residual
      enbal.dt(iter,1)        =  dt_sum;        % check
      %enbal.zD(iter,1)        =  zD;
      
      % Save the ice column data
      ice2.Tice(:,iter)       =  T;             % ice temperature
      ice2.f_ice(:,iter)      =  f_ice;         % fraction ice
      ice2.f_liq(:,iter)      =  f_liq;         % fraction liq
      ice2.Sc(:,iter)         =  Sc;            % source term
      ice2.errH(:,iter)       =  errH;          % enthalpy error
      ice2.errT(:,iter)       =  errT;          % temperature error
      ice2.df_liq(:,iter)     =  d_liq;
      ice2.df_drn(:,iter)     =  d_drn;

      % save the diagnostic data
      % diags.Tflag(iter,1)     =  Tflag; 
      % diags.LCflag(iter,1)    =  LCflag(1);

   end
   
   % move to next timestep
   iter     =  iter + 1;

end
   % restart the model
   iter     =  1;
   subiter  =  1;
   itime    =  met.Time(1);             % initialize model time
end