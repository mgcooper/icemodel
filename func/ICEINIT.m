%--------------------------------------------------------------------------
%   ICEINIT - initialize the model
%--------------------------------------------------------------------------

function [  f_ice,                                                      ...
            f_liq,                                                      ...
            T,                                                          ...
            TL,                                                         ...
            TH,                                                         ...
            flmin,                                                      ...
            flmax,                                                      ...
            cp_sno,                                                     ...
            k_eff,                                                      ...
            dz,                                                         ...
            fn,                                                         ...
            delz,                                                       ...
            grid_therm,                                                 ...
            dz_therm,                                                   ...
            dz_spect,                                                   ...
            JJ_therm,                                                   ...
            JJ_spect,                                                   ...
            Sc,                                                         ...
            Sp,                                                         ...
            wcoef,                                                      ...
            scoef,                                                      ...
            ro_sno,                                                     ...
            ro_iwe,                                                     ...
            ro_wie ]   =   ICEINIT(opts,met)
%--------------------------------------------------------------------------

% load the physical constants to be used.
   load('PHYSCONS','cp_ice','cp_liq','fcp','kappa','Lf','ro_ice',       ...
      'ro_air','ro_liq','k_liq','Tf','Ls','Rv');

   % pull out anything in opts used below

   % # of nodes in the thermal and spectral grid
   JJ_therm       =  opts.z0_thermal/opts.dz_thermal;
   JJ_spect       =  opts.z0_spectral/opts.dz_spectral;
   dz_therm       =  opts.dz_thermal;
   dz_spect       =  opts.dz_spectral;
   z0_therm       =  opts.z0_thermal;
   ro_i           =  opts.ro_snow_i;
%    maxiter        =  opts.maxiter;
%    use_ro_glc     =  opts.ro_glc;
   
%    % test
%    A1 = ones(JJ_therm,1);
%    A = spdiags([A1,A1,A1],-1:1,JJ_therm,JJ_therm);
   
   % Compute c.v. size information.
[  dz,                                                               ...
   delz,                                                             ...
   fn,                                                               ...
   grid_therm ]   =  CVTHERMAL(z0_therm,dz_therm);
   
   % Density for melt/freeze
%    ro_glc      =  (917+1000)/2;                    % [kg/m3]
%    if use_ro_glc == true
%       ro_ice   =  ro_glc;
%       ro_liq   =  ro_glc;
%    end
   ro_wie      =  ro_liq./ro_ice;
   ro_iwe      =  ro_ice./ro_liq;

%    % Init ice temperature
%    if opts.use_init == true
%       load(opts.f_init);
%       T        =  init.T+met.tair(1)-init.T(1);
%    else
%       T        =  (met.tair(1)-1).*ones(JJ_therm,1);
%    end

   T           =  (met.tair(1)-1).*ones(JJ_therm,1);

   % Init liquid/ice water fraction, bulk densities, and aP coeff
   Tdep        =  Tf-T(:,1);                       % [K]
   fliq        =  1./(1+(fcp.*Tdep).^2);           % [-]
   bd_ice      =  ro_i;                            % at t1, bd_ice=ro_sno
   bd_liq      =  bd_ice.*(fliq./(1-fliq));
   f_liq       =  bd_liq./ro_liq;
   f_ice       =  bd_ice./ro_ice.*ones(JJ_therm,1);


[  k_eff,                                                            ...
   ro_sno,                                                           ...
   cp_sno   ]  =  UPDATESTATE(T,f_ice,f_liq,ro_ice,ro_liq,ro_air,    ...
                  cp_ice,cp_liq,k_liq,Ls,Rv,Tf);

   % lower and upper melt zone temperature and water fractions
   TL          =  Tf-(2.0*Lf/(fcp*fcp*cp_ice))^(1.0/3.0);   % Eq 120
   TH          =  Tf-cp_liq/(Lf*2.0*fcp*fcp);
   flmin       =  1/(1+(fcp*(Tf-TL))^2);
   flmax       =  1/(1+(fcp*(Tf-TH))^2);
   
   % precomputed stability coefficients to speed up the code
   z_obs       =  opts.z_obs;
   z_0         =  opts.z_0;
   wcoef       =  (kappa/log(z_obs/z_0))^2;        % wind transfer coef
   scoef(1)    =  5.3*9.4*wcoef*sqrt(z_obs/z_0);   % gamma Eq. A15
   scoef(2)    =  9.4*9.81*z_obs;                  % 9.81=gravity
   scoef(3)    =  scoef(1)*sqrt(9.81*z_obs);
   clear kappa z_obs z_0
   
   % source term linearization vectors
   Sc          =  zeros(JJ_therm,1);
   Sp          =  zeros(JJ_therm,1);
   
   % Initialize the output structures
%    enbal.Tsfc        = nan(maxiter,1);
%    enbal.Qle         = nan(maxiter,1);
%    enbal.Qh          = nan(maxiter,1);
%    enbal.Qe          = nan(maxiter,1);
%    enbal.Qc          = nan(maxiter,1);
%    enbal.Qm          = nan(maxiter,1);
%    enbal.Qf          = nan(maxiter,1);
%    enbal.chi         = nan(maxiter,1);
%    enbal.balance     = nan(maxiter,1);
%    enbal.dt          = nan(maxiter,1);
%    enbal.zD          = nan(maxiter,1);
%    
%    % Save the vertical ice column data   
%    ice2.Tice         = nan(JJ_therm,maxiter);
%    ice2.f_ice        = nan(JJ_therm,maxiter);
%    ice2.f_liq        = nan(JJ_therm,maxiter);
%    ice2.Sc           = nan(JJ_therm,maxiter);
%    ice2.errH         = nan(JJ_therm,maxiter);
%    ice2.errT         = nan(JJ_therm,maxiter);
%    ice2.df_liq       = nan(JJ_therm,maxiter);
%    ice2.df_drn       = nan(JJ_therm,maxiter);
   
%    diags.Tflag       = false(maxiter,1);
%    diags.LCflag      = false(maxiter,1);
   

%       dz       =  opts.dz_thermal;
%       Z        =  opts.z0_thermal;
%       z        =  (dz/2:dz:Z-dz/2)';      
%       figure; plot([init.T init.T+dT],z); set(gca,'YDir','reverse'); 
      
%       cp_sno(:,1) =  init.cp_sno;
%       k_eff(:,1)  =  init.k_eff;
%       ro_sno(:,1) =  init.ro_sno;
%       f_liq(:,1)  =  init.f_liq;
%       f_ice(:,1)  =  init.f_ice;
   
% normally f_liq = f_wat/(1+(fcp*(Tf-Td))^2) 
% so I use f_liq_min = f_wat/(1+(fcp*(Tf-TL))^2) = flmin
% and similar for f_liq_max
% but f_liq_min/max should change as f_wat changes. 
% so instead I should define flmin/max as 
%    flmin    = 1/(1+(fcp*(Tf-TL))^2)
%    flmax    = 1/(1+(fcp*(Tf-TH))^2)
% where flmin/max are now the same quantities defined in sntherm as
% flglim(1) and flglim(2) (i.e., f_sub_ell in eq. 67), rather than
% frac_liq_min/max as I have them above returned from MELTCURVE. I would
% then update frac_liq_min/max at the top of ICEENBAL as:
%    fliqmin = f_wat./flmin;
%    fliqmax = f_wat./flmax;
% and send these into MZTRANSFORM
   
% this is how jordan does it:
% first compute TL,TH (as above)
% then compute flmin/flmax as:
%    flmin    = 1/(1+(fcp*(Tf-TL))^2)
%    flmax    = 1/(1+(fcp*(Tf-TH))^2)
% note that f_liq/f_wat = g_liq/g_wat 


% % confirm that MELTCURVE recovers the same frac's
%    f_wat    = f_liq+f_ice.*ro_ice./ro_liq;              % frac_wat_old
%    f_liq2   = f_wat./(1.0+(fcp.*Tdep).^2.0);            % eq 67, Jordan
%    f_ice2   = (f_wat - f_liq) .* ro_liq ./ ro_ice;
%    f_air2   = 1.0-f_liq-f_ice;
   
%    % save the ice surface data
%    surf_mlt    =  0.0;
%    surf_frz    =  0.0;
%    surf_rof    =  0.0;
%    surf_sub    =  0.0;
%    surf_con    =  0.0;

%    ice1.surf_melt    = zeros(maxiter,1);
%    ice1.surf_freeze  = zeros(maxiter,1);
%    ice1.surf_subl    = zeros(maxiter,1);
%    ice1.surf_cond    = zeros(maxiter,1);
%    ice1.surf_runoff  = zeros(maxiter,1);
   
      
   % put the data into an output structure (unimplemented)
   % init.dy_p                       =   dy_p;
   % init.f_n                        =   f_n;
   % init.dely_p                     =   dely_p;
   % init.T_old                      =   T_old;
   % init.melt_frac                  =   melt_frac;
   % init.freeze                     =   freeze;
   % init.gamma                      =   gamma;
   % init.ro_snow_z                  =   ro_snow_z;
   % init.xk_snow                    =   xk_snow;
   % init.freeze_frac                =   freeze_frac;
   % init.water_depth_old            =   water_depth_old;
   % init.freeze_depth_old           =   freeze_depth_old;