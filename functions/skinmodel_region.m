
function [ice1,ice2,met,opts] = skinmodel_region(opts)
   
%-------------------------------------------------------------------------------
%   INITIALIZE THE MODEL
%-------------------------------------------------------------------------------
   
% load the physical constants to be used.
load( 'PHYSCONS', 'cp_ice','cv_air','cv_liq','cv_ice','emiss',       ...
                  'SB','epsilon','fcp','k_liq','Lf','Ls','Lv',       ...
                  'ro_air','ro_ice','ro_liq','roLf','roLs','roLv',   ...
                  'Rv','Tf','TINY','cp_liq');
                  

%  LOAD THE FORCING DATA
[  met,opts ] = METINIT(opts);
                  
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
   ice1,                                                                ...
   ice2  ]    =  ICEINIT(opts,met);


% % init output arrays - for region
%    df_liq      =  nan(JJ_therm,maxiter);
%    frac_ice    =  nan(JJ_therm,maxiter);
%    frac_liq    =  nan(JJ_therm,maxiter);
%    T_ice       =  nan(JJ_therm,maxiter);
%    T_sfc       =  nan(maxiter,1);
   

%    surf_melt      = 0.0;
%    surf_freeze    = 0.0;
%    surf_subl      = 0.0;
%    surf_cond      = 0.0;
%    surf_runoff    = 0.0;
%                         
%    ice1.surf_melt    = zeros(opts.maxiter,1);
%    ice1.surf_freeze  = zeros(opts.maxiter,1);
%    ice1.surf_subl    = zeros(opts.maxiter,1);
%    ice1.surf_cond    = zeros(opts.maxiter,1);
%    ice1.surf_runoff  = zeros(opts.maxiter,1);
      
%-------------------------------------------------------------------------------
%   INITIALIZE THE EXTINCTION COEFFICIENTS
%-------------------------------------------------------------------------------

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

%    clear radii scattercoefs solar
%-------------------------------------------------------------------------------
%  initialize timestepping and initial state values
%-------------------------------------------------------------------------------
   
% timestepping settings
   dt          =  opts.dt;
   maxiter     =  opts.maxiter;
   maxsubiter  =  200;
   minsubiter  =  1;
   dt_min      =  dt/maxsubiter;  % or: dt_min = 5; maxsubiter  = dt/dt_min;
   dt_max      =  dt/minsubiter;  % or: dt_min = 5; maxsubiter  = dt/dt_min;
   dt_new      =  dt_max;
   iter        =  1;
   subiter     =  1;

% extract parameter values that are used frequently
   nloops      =  opts.annual_loops;
   fopts       =  opts.fzero;
   f_min       =  opts.f_ice_min;

   simyears    =  opts.simyears;

% we need these on the first iteration, they are updated after each substep
   Tsfc        =  T(1);
   xf_liq      =  f_liq;
   xTsfc       =  Tsfc;
   roL         =  roLs;
   liqflag     =  false;
   Qc          =  CONDUCT(k_eff,T,dz,xTsfc);

   
if nloops > 1
   simyears    =  [repmat(simyears(1),nloops,1) simyears(2:end)];
end

   numyears    =  numel(simyears);
   
for MM = 1:numyears
   
   % set the output folder and load the met data
   simyear     = int2str(simyears(MM));
   pathoutput  = [opts.pathoutput simyear '/'];
   
   t1          = datetime(simyears(MM),1,1,0,0,0);
   t2          = datetime(simyears(MM),12,31,23,45,0);
   t1.TimeZone = 'UTC';
   t2.TimeZone = 'UTC';
   Time        = met.Time(isbetween(met.Time,t1,t2));
   iter0       = find(met.Time==t1);

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
   
   while iter <= maxiter
   
   metiter     =  iter+iter0-1;

% reset the change in liq water content
   d_liq       =  0.0.*f_liq;
   d_drn       =  0.0.*f_liq;

%  load the atmospheric forcing data for this time step.
   Tair        =  met.tair(metiter);
   rh          =  met.rh(metiter);
   wspd        =  met.wspd(metiter);
   Qsi         =  met.swd(metiter);
   Qli         =  met.lwd(metiter);
   Pa          =  met.psfc(metiter);
   albedo      =  met.albedo(metiter);
   De          =  wspd*wcoef;
   
%  Update the subsurface absorption profile and solar radiation source-term
   if opts.meltmodel == "skinmodel"
      chi      =  1.0;
   elseif Qsi>0
   [  Sc,                                                               ...
      chi   ]  =  UPDATEEXTCOEFS(grid_spect,JJ_spect,dz_spect,          ...
                  grid_therm,JJ_therm,dz_therm,z_walls,dz,ro_sno,       ...
                  Qsi,total_solar,spect_lower,spect_upper,albedo,       ...
                  solardwavl);
   else
      Sc       =  0.0.*Sc;
      chi      =  0.0;
   end
   
%-------------------------------------------------------------------------------
%   Initialize dynamic timestepping
%-------------------------------------------------------------------------------
   dt_sum      =  0.0;
   subfail     =  0;                      % keep track of failed substeps
   dt_flag     =  false;                  %
   %OK         =  false;                  % assume sub-stepping is needed
   OK          =  true;                   % assume sub-stepping is needed
   
   % while OK == false || dt_sum < dt
   while dt_sum < dt

% SURFACE ENERGY BALANCE
   [  Qm,                                                               ...
      Qf,                                                               ...
      Qh,                                                               ...
      Qe,                                                               ...
      Qc,                                                               ...
      ~,                                                                ...
      balance,                                                          ...
      Tsfc  ]  =  ENBALANCE(Tair,wspd,rh,Qsi,Qli,albedo,k_eff,Pa,       ...
                  T,Tf,dz,chi,xTsfc,cv_air,emiss,SB,roL,De,scoef,       ...
                  epsilon,fopts,liqflag,opts);

% SUBSURFACE ENERGY BALANCE
      T        =  SKINSOLVE(Tsfc,T,k_eff,ro_sno,cp_sno,dz,dt_new,       ...
                  JJ_therm,fn,delz);

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

%          [ surf_melt,                                                   ...
%            surf_freeze,                                                 ...
%            surf_subl,                                                   ...
%            surf_cond,                                                   ...
%            surf_runoff]  =  ICEABLATIONtmp(Qm,Qe,Qf,surf_melt,          ...
%                             surf_freeze,surf_runoff,surf_subl,          ...
%                             surf_cond,ro_ice,ro_liq,Lf,Ls,dt,opts);
                         
      % recompute top fluxes
%          k_eff    =  GETGAMMA(T,f_liq,f_ice,ro_ice,k_liq,Ls,Rv,Tf);
%          Qc       =  CONDUCT(k_eff,T,dz,Tsfc);
%          ea       =  VAPPRESS(rh,Tair,liqflag);
%          S        =  STABLEFN(Tair,Tsfc,wspd,scoef);
%          es0      =  VAPOR(Tsfc,Tf,liqflag);
%          Qe       =  LATENT(De,S,ea,es0,roL,epsilon,Pa);
%          Qh       =  SENSIBLE(De,S,Tair,Tsfc,cv_air);
%          Qle      =  LONGOUT(Tsfc,emiss,SB);
% 
%       % Compute the energy flux available for melting or freezing.
%          [Qm,Qf]  =  MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Tsfc,Tf,      ...
%                      Tair,wspd,De,ea,roL,Pa,cv_air,emiss,SB,k_eff,      ...
%                      T,dz,epsilon,scoef,chi);
%          balance  =  ENBAL(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,emiss,chi);

%       % Substep was successful, compute fluxes and update state
%       [  T,                                                             ...
%          f_ice,                                                         ...
%          f_liq,                                                         ...
%          d_liq,                                                         ...
%          d_drn ]  =  ICEMF(T,f_ice,f_liq,ro_ice,ro_liq,cp_ice,Lf,Ls,    ...
%                      Lv,Tf,TL,fcp,xf_liq,Sc,Sp,JJ_therm,f_min,fopts,    ...
%                      dz_therm,dt_new,Qe,liqflag,ro_iwe,d_liq,d_drn);
                        
      % update density (kg/m3), heat capacity (J/kg/K), thermal K (W/m/K)
      [  k_eff,                                                         ...
         ro_sno,                                                        ...
         cp_sno ] =  UPDATESTATE(T,f_ice,f_liq,ro_ice,ro_liq,ro_air,    ...
                     cp_ice,cp_liq,k_liq,Ls,Rv,Tf);
         
      % allocate this step to the running total
         dt_sum   =  dt_sum + dt_new;

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

%-------------------------------------------------------------------------------
%   Save the output
%-------------------------------------------------------------------------------

   % NEED TO GO BACK TO ARRAYS SO THEY CAN BE WRITTEN OVER EACH YEAR or USE
   % THE KEEP VERSIONS BELOW ... BUT ALSO NEED TO REVISIT SURFRUNOFF TO SEE
   % HOW I DEALT WITH FREEZE, IT CAN ONLY EXIST IF WATER IS AVAILABLE SO AS
   % IT IS IT IS WRONG

   if MM-nloops >= 0

      % save the surface energy balance
      ice1.Tsfc(iter,1)       =  Tsfc;          % surface temp
      ice1.Qm(iter,1)         =  Qm;            % melt energy
      ice1.Qf(iter,1)         =  Qf;            % freeze deficit
      ice1.Qe(iter,1)         =  Qe;            % latent
      ice1.Qh(iter,1)         =  Qh;            % sensible
      ice1.Qc(iter,1)         =  Qc;            % conductive into surface
      ice1.chi(iter,1)        =  chi;           % skin parameter
      ice1.balance(iter,1)    =  balance;       % SEB residual
      ice1.dt(iter,1)         =  dt_sum;        % check
      %ice1.zD(iter,1)        =  zD;
      
      % Save the ice column data
      ice2.Tice(:,iter)       =  T;             % ice temperature
    % ice2.f_ice(:,iter)   =  f_ice;         % fraction ice
    % ice2.f_liq(:,iter)   =  f_liq;         % fraction liq
    % ice2.df_liq(:,iter)  =  d_liq;
    % ice2.df_drn(:,iter)  =  d_drn; 
    % ice2.Sc(:,iter)      =  Sc;            % source term
    % ice2.errH(:,iter)    =  errH;          % enthalpy error
    % ice2.errT(:,iter)    =  errT;          % temperature error
      
      % save the surf data
%       ice1.surf_melt(iter,1)     =  surf_melt;
%       ice1.surf_freeze(iter,1)   =  surf_freeze;    % skin freeze
%       ice1.surf_subl(iter,1)     =  surf_subl;      % sublimation
%       ice1.surf_cond(iter,1)     =  surf_cond;      % condensation
%       ice1.surf_runoff(iter,1)   =  surf_runoff;    % skin runoff
   
   end

   % move to next timestep
   iter     =  iter + 1;
   
   end % timesteps (one year)
   
   % option to spin up year 1 by looping over it 
   if MM-nloops<0
      % restart the counters
      iter     =  1;
      subiter  =  1;
      continue
   end
   
   ice1keep = ice1;
   ice2keep = ice2;
   
   % post process
   % [ice1,ice2,met] = POSTPROC(ice1,ice2,met,opts);
   
   [ice1,ice2] = POSTPROC2(ice1,ice2,Time,opts);
   
%-------------------------------------------------------------------------------
%  save the data
%-------------------------------------------------------------------------------
   if opts.savedata == true
      if ~exist(pathoutput,'dir'); mkdir(pathoutput); end
      save([pathoutput 'ice1_' int2str(opts.ipt)],'ice1','-v7.3');
      save([pathoutput 'ice2_' int2str(opts.ipt)],'ice2','-v7.3');
   end
   
   % restart the counters
   iter     =  1;
   subiter  =  1;
   
   ice1     =  ice1keep;
   ice2     =  ice2keep;
end % numyears
