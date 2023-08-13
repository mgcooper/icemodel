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
            scoef,                                                      ...
            ro_sno,                                                     ...
            ro_iwe,                                                     ...
            ro_wie,                                                     ...
            xTsfc,                                                      ...
            xf_liq,                                                     ...
            roL,                                                        ...
            Qc,                                                         ...
            f_min,                                                      ...
            fopts,                                                      ...
            liqflag,                                                    ... 
            ice1,                                                       ...
            ice2  ]   =   ICEINIT(opts, tair) %#codegen
%--------------------------------------------------------------------------

   % load the physical constants.
   [cp_ice,cp_liq,fcp,kappa,Lf,ro_ice,ro_air,ro_liq,k_liq,Tf,Ls,Rv,roLs] = ...
      icemodel.physicalConstant('cp_ice','cp_liq','fcp','kappa','Lf', ...
      'ro_ice','ro_air','ro_liq','k_liq','Tf','Ls','Rv','roLs');

   % pull out anything in opts used below

   % # of nodes in the thermal and spectral grid
   JJ_therm       =  opts.z0_thermal/opts.dz_thermal;
   JJ_spect       =  opts.z0_spectral/opts.dz_spectral;
   dz_therm       =  opts.dz_thermal;
   dz_spect       =  opts.dz_spectral;
   z0_therm       =  opts.z0_thermal;
   ro_i           =  opts.ro_snow_i;
   maxiter        =  opts.maxiter;
   use_ro_glc     =  opts.use_ro_glc;
   
   % Compute c.v. size information.
[  dz,                                                               ...
   delz,                                                             ...
   fn,                                                               ...
   grid_therm ]   =  CVTHERMAL(z0_therm,dz_therm);
   
   % Density for melt/freeze
   ro_glc      =  (917+1000)/2;                    % [kg/m3]
   if use_ro_glc == true
      ro_ice   =  ro_glc;
      ro_liq   =  ro_glc;
   end
   ro_wie      =  ro_liq./ro_ice;
   ro_iwe      =  ro_ice./ro_liq;

%    % Init ice temperature
%    if opts.use_init == true
%       load(opts.f_init);
%       T        =  init.T+met.tair(1)-init.T(1);
%    else
%       T        =  (met.tair(1)-1).*ones(JJ_therm,1);
%    end

   T        =  (tair(1)-1).*ones(JJ_therm,1);

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
   scoef       =  nan(1, 3);
   scoef(1)    =  5.3*9.4*wcoef*sqrt(z_obs/z_0);   % gamma Eq. A15
   scoef(2)    =  9.4*9.81*z_obs;                  % 9.81=gravity
   scoef(3)    =  scoef(1)*sqrt(9.81*z_obs);
   clear kappa z_obs z_0
   
   % source term linearization vectors
   Sc          =  zeros(JJ_therm,1);
   Sp          =  zeros(JJ_therm,1);

   % initialize state variables and parameters needed on the first iteration
   xTsfc       =  T(1);
   xf_liq      =  f_liq;
   roL         =  roLs;
   liqflag     =  false;
   Qc          =  CONDUCT(k_eff,T,dz,xTsfc);
   % zD        =  sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1)));
   fopts       =  opts.fzero;
   f_min       =  opts.f_ice_min;

% Initialize the output structures
if opts.sitename == "sector"
   ice1.Tsfc         = nan(maxiter,1);
   ice2.Tice         = nan(JJ_therm,maxiter);
   ice2.f_ice        = nan(JJ_therm,maxiter);
   ice2.f_liq        = nan(JJ_therm,maxiter);
   ice2.df_liq       = nan(JJ_therm,maxiter);
else
   ice1.Tsfc        = nan(maxiter,1);
   ice1.Qm          = nan(maxiter,1);
   ice1.Qf          = nan(maxiter,1);
   ice1.Qe          = nan(maxiter,1);
   ice1.Qh          = nan(maxiter,1);
   ice1.Qc          = nan(maxiter,1);
   ice1.chi         = nan(maxiter,1);
   ice1.balance     = nan(maxiter,1);
   ice1.dt          = nan(maxiter,1);
   ice1.zD          = nan(maxiter,1);

   % Save the vertical ice column data   
   ice2.Tice         = nan(JJ_therm,maxiter);
   ice2.f_ice        = nan(JJ_therm,maxiter);
   ice2.f_liq        = nan(JJ_therm,maxiter);
   ice2.Sc           = nan(JJ_therm,maxiter);
   ice2.errH         = nan(JJ_therm,maxiter);
   ice2.errT         = nan(JJ_therm,maxiter);
   ice2.df_liq       = nan(JJ_therm,maxiter);
   ice2.df_drn       = nan(JJ_therm,maxiter);
   ice2.df_evp       = nan(JJ_therm,maxiter);
   % diags.Tflag       = false(maxiter,1);
   % diags.LCflag      = false(maxiter,1);
end
   
%    % test
%    A1 = ones(JJ_therm,1);
%    A = spdiags([A1,A1,A1],-1:1,JJ_therm,JJ_therm);
