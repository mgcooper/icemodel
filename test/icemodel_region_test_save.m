
function [ice1,ice2,times] = icemodel_region_test_save(opts)
   
%--------------------------------------------------------------------------
%   INITIALIZE THE MODEL
%--------------------------------------------------------------------------
   
% load the physical constants to be used.
   load( 'PHYSCONS', 'cp_ice','cv_air','cv_liq','cv_ice','emiss',       ...
                     'SB','epsilon','fcp','k_liq','Lf','Ls','Lv',       ...
                     'ro_air','ro_ice','ro_liq','roLf','roLs','roLv',   ...
                     'Rv','Tf','TINY');
                  
% load the first met file
   simyear  = int2str(opts.simyears(1));
   load([opts.pathmet simyear '/met_' int2str(opts.ipt) '.mat'],'met');

   if isleap(opts.simyears(1))
      met = rmleapinds(met);
   end
   maxiter = size(met,1);
   
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
   ro_wie ]    =  ICEINIT(opts,met);
      

% init output arrays
   df_liq      =  nan(JJ_therm,maxiter);
   frac_ice    =  nan(JJ_therm,maxiter);
   frac_liq    =  nan(JJ_therm,maxiter);
   T_ice       =  nan(JJ_therm,maxiter);
   T_sfc       =  nan(maxiter,1);
   
%--------------------------------------------------------------------------
%   INITIALIZE THE EXTINCTION COEFFICIENTS
%--------------------------------------------------------------------------

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
%--------------------------------------------------------------------------
%   START THE MODEL
%--------------------------------------------------------------------------
   
   % extract struct values that are used frequently
   dt          =  opts.dt;
   fopts       =  opts.fzero;
   f_min       =  opts.f_ice_min;
   
%--------------------------------------------------------------------------
% substepping settings
%--------------------------------------------------------------------------
   maxsubiter  =  200;
   minsubiter  =  1;
   dt_min      =  dt/maxsubiter;  % or: dt_min = 5; maxsubiter  = dt/dt_min;
   dt_max      =  dt/minsubiter;  % or: dt_min = 5; maxsubiter  = dt/dt_min;
   dt_new      =  dt_max;
   iter        =  1;
   subiter     =  1;

% we need these on the first iteration, they are updated after each substep
   Tsfc        =  T(1);
   xf_liq      =  f_liq;
   xTsfc       =  Tsfc;
   roL         =  roLs;
   liqflag     =  false;
   Qc          =  CONDUCT(k_eff,T,dz,xTsfc);

for MM = 1:opts.numyears
   
   % set the output folder and load the met data
   simyear     = int2str(opts.simyears(MM));
   pathoutput  = [opts.pathsave simyear '/'];
   
   load([opts.pathmet simyear '/met_' int2str(opts.ipt) '.mat'],'met');
   if isleap(opts.simyears(MM))
      met = rmleapinds(met);
   end
   
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
   
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
   De          =  wspd*wcoef;
   ea          =  VAPPRESS(rh,Tair,liqflag);
   
   
%  Update the subsurface absorption profile and solar radiation source-term
   if Qsi>0
   [  Sc,                                                               ...
      chi   ]  =  UPDATEEXTCOEFS(grid_spect,JJ_spect,dz_spect,          ...
                  grid_therm,JJ_therm,dz_therm,z_walls,dz,ro_sno,       ...
                  Qsi,total_solar,spect_lower,spect_upper,albedo,       ...
                  solardwavl);
   else
      Sc       =  0.0.*Sc;
      chi      =  0.0;
   end

%--------------------------------------------------------------------------
%   Initialize dynamic timestepping
%--------------------------------------------------------------------------
   dt_sum      =  0.0;
   subfail     =  0;                      % keep track of failed substeps
   dt_flag     =  false;                  %
   OK          =  false;                  % assume sub-stepping is needed
   
   while OK == false || dt_sum < dt

% SURFACE TEMPERATURE
   [  Tsfc  ]  =  SFCTEMP(Tair,Qsi,Qli,ea,albedo,De,Pa,wspd,cv_air,     ...
                  emiss,SB,Tf,Qc,xTsfc,chi,roL,scoef,fopts,liqflag);

% SUBSURFACE ENERGY BALANCE
   [  T,                                                                ...
      ~,                                                                ...
      ~,                                                                ...
      f_ice,                                                            ...
      f_liq,                                                            ...
      OK    ]  =  ICEENBAL(T,f_ice,f_liq,k_liq,cv_ice,cv_liq,           ...
                  ro_ice,ro_liq,ro_sno,cp_sno,Ls,Lf,roLf,Rv,Tf,         ...
                  dz,delz,fn,dt_new,JJ_therm,Tsfc,Sc,fcp,TL,TH,         ...
                  flmin,flmax,ro_iwe,ro_wie);

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

      % faster option to send an updated Qe to ICEMF to compute subl:
         k_eff    =  GETGAMMA(T,f_liq,f_ice,ro_ice,k_liq,Ls,Rv,Tf);
         Qc       =  CONDUCT(k_eff,T,dz,MELTTEMP(Tsfc,Tf));
         ea       =  VAPPRESS(rh,Tair,liqflag);
         S        =  STABLEFN(Tair,MELTTEMP(Tsfc,Tf),wspd,scoef);
         es0      =  VAPOR(MELTTEMP(Tsfc,Tf),Tf,liqflag);
         Qe       =  LATENT(De,S,ea,es0,roL,epsilon,Pa);

      % Substep was successful, compute fluxes and update state
      [  T,                                                             ...
         f_ice,                                                         ...
         f_liq,                                                         ...
         d_liq,                                                         ...
         d_drn ]  =  ICEMF(T,f_ice,f_liq,ro_ice,ro_liq,cp_ice,Lf,Ls,    ...
                     Lv,Tf,TL,fcp,xf_liq,Sc,Sp,JJ_therm,f_min,fopts,    ...
                     dz_therm,dt_new,Qe,liqflag,ro_iwe,d_liq,d_drn);
                        
      % update density (kg/m3), heat capacity (J/kg/K), thermal K (W/m/K)
         ro_sno   =  f_ice.*ro_ice+f_liq.*ro_liq+(1.0-f_liq-f_ice).*ro_air;
         cp_sno   =  (cv_ice.*f_ice+cv_liq.*f_liq)./ro_sno;
         
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

%--------------------------------------------------------------------------
%   Save the output
%--------------------------------------------------------------------------

      % for a stripped-down run this is all that's needed:
      T_sfc    (iter,1)    =  Tsfc;
      T_ice    (:,iter)    =  T;
      frac_ice (:,iter)    =  f_ice;
      frac_liq (:,iter)    =  f_liq;
      df_liq   (:,iter)    =  d_liq;

      % move to next timestep
      iter     =  iter + 1;

   end % timesteps (one year)

%--------------------------------------------------------------------------
%%  save the data
%--------------------------------------------------------------------------
   if opts.savedata == true
      if ~exist(pathoutput,'dir'); mkdir(pathoutput); end
      
      versions = {'-v4','-v6','-v7','-v7.3'};
      times    = nan(3,numel(versions)+2);
      
      for n = 1:numel(versions)
         
         v     = versions{n};
         fsfx  = strrep(strrep(v,'-',''),'.','');

         % tt
         [ice1,ice2] = POSTPROC2(T_sfc,T_ice,frac_ice,frac_liq,...
                                       df_liq,met,opts);
         ice1 = removevars(ice1,{'surf_runoff','column_runoff'});
         tic;  save([pathoutput 'tt_' fsfx],'ice1','ice2',v);
         times(1,n) = toc;
         
         if n == numel(versions)
            tic;  save([pathoutput 'tt_v7_no'],'ice1','ice2','-v7','-nocompression');
            times(1,n+1) = toc;
            tic;  save([pathoutput 'tt_v73_no'],'ice1','ice2','-v7.3','-nocompression');
            times(1,n+2) = toc;
         end
         
         % struct
         [ice1,ice2] = POSTPROC2_test(T_sfc,T_ice,frac_ice,frac_liq,...
                                       df_liq,met,opts);
         tic;  save([pathoutput 'struct_' fsfx],'ice1','ice2',v);
         times(2,n) = toc;
         
         if n == numel(versions)
            tic;  save([pathoutput 'struct_v7_no'],'ice1','ice2','-v7','-nocompression');
            times(2,n+1) = toc;
            tic;  save([pathoutput 'struct_v73_no'],'ice1','ice2','-v7.3','-nocompression');
            times(2,n+2) = toc;
         end
         
         % arrays
         [Tsfc,runoff,melt,freeze,                                   ...
         Tice,f_ice,f_liq,d_liq] = POSTPROC2_test2(T_sfc,T_ice,    ...
                                 frac_ice,frac_liq,df_liq,met,opts);
      
         tic; save([pathoutput 'vec_' fsfx], ...
         'Tsfc','runoff','melt','freeze', ...
         'Tice','f_ice','f_liq','d_liq',v);
         times(3,n) = toc;
         
         if n == numel(versions)
            tic;  save([pathoutput 'vec_v7_no'], ...
            'Tsfc','runoff','melt','freeze', ...
            'Tice','f_ice','f_liq','d_liq','-v7','-nocompression');
            times(3,n+1) = toc;
            tic;  save([pathoutput 'vec_v73_no'], ...
            'Tsfc','runoff','melt','freeze', ...
            'Tice','f_ice','f_liq','d_liq','-v7.3','-nocompression');
            times(3,n+2) = toc;
         end
         
      end
      times = array2table(times,'VariableNames',{'v4','v6','v7','v73','v7no','v73no'}, ...
      'RowNames',{'timetable','struct','arrays'});
      
   end
   
   % restart the counters
   iter     =  1;
   subiter  =  1;
   
end % numyears