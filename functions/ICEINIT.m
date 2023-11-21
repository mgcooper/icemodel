function [f_ice, f_liq, T, TL, TH, flmin, flmax, cp_sno, k_eff, dz, fn, ...
      delz, z_therm, dz_therm, dz_spect, JJ_therm, JJ_spect, Sc, Sp, scoef, ...
      ro_sno, ro_iwe, ro_wie, xTsfc, xf_liq, roL, Qc, f_min, fopts, ...
      liqresid, liqflag, ice1, ice2] = ICEINIT(opts, tair)
   %ICEINIT initialize the 1-d ice column
   %
   %#codegen

   % load the physical constants.
   [cp_ice, cp_liq, fcp, kappa, Lf, ro_ice, ro_air, ro_liq, k_liq, Tf, Ls, ...
      Rv, roLs, gravity] ...
      = icemodel.physicalConstant('cp_ice','cp_liq','fcp','kappa','Lf', ...
      'ro_ice','ro_air','ro_liq','k_liq','Tf','Ls','Rv','roLs','gravity');

   % COMPUTE THE # OF NODES IN THE THERMAL AND SPECTRAL GRIDS
   JJ_therm = opts.z0_thermal / opts.dz_thermal;
   JJ_spect = opts.z0_spectral / opts.dz_spectral;
   dz_therm = opts.dz_thermal;
   dz_spect = opts.dz_spectral;
   z0_therm = opts.z0_thermal;

   % NUMBER OF TIMESTEPS TO INITIALIZE OUTPUTS
   maxiter = numel(tair) / opts.numyears;

   % GENERATE A THERMAL MESH
   [dz, delz, z_therm, ~, fn] = CVMESH(z0_therm, dz_therm);

   % INITIALIZE DENSITIES
   ro_glc = (ro_ice + ro_liq) / 2; % [kg/m3]
   if opts.use_ro_glc == true
      ro_ice = ro_glc;
      ro_liq = ro_glc;
   end
   ro_wie = ro_liq / ro_ice;
   ro_iwe = ro_ice / ro_liq;

   % INITIALIZE ICE TEMPERATURE TO AIR TEMPERATURE
   T = (tair(1) - 1) * ones(JJ_therm, 1);

   % INITIALIZE LIQUID/ICE WATER FRACTION AND BULK DENSITIES
   T_dep = Tf - T;                           % [K]
   fmliq = 1 ./ (1 + (fcp * T_dep) .^ 2);    % [-]
   g_ice = opts.ro_snow_i;
   g_liq = g_ice .* (fmliq ./ (1 - fmliq));
   f_liq = g_liq ./ ro_liq;
   f_ice = g_ice ./ ro_ice .* ones(JJ_therm, 1);

   % INITIAL THERMAL CONDUCTIVITY, BULK DENSITY, AND HEAT CAPACITY
   [k_eff, ro_sno, cp_sno] = UPDATESTATE(T, f_ice, f_liq, ro_ice, ro_liq, ...
      ro_air, cp_ice, cp_liq, k_liq, Ls, Rv, Tf);

   % LOWER AND UPPER MELT ZONE TEMPERATURE AND WATER FRACTIONS
   TL = Tf - (2.0 * Lf / (fcp ^ 2.0 * cp_ice)) ^ (1.0 / 3.0); % Eq 120
   TH = Tf - cp_liq / (Lf * 2.0 * fcp ^ 2.0);
   flmin = 1 / (1 + (fcp * (Tf - TL)) ^ 2.0);
   flmax = 1 / (1 + (fcp * (Tf - TH)) ^ 2.0);

   % PRECOMPUTE STABILITY COEFFICIENTS
   z_0 = opts.z_0;
   z_obs = opts.z_obs;
   wcoef = (kappa / log(z_obs / z_0)) ^ 2; % wind transfer coef
   scoef = nan(1, 3);
   scoef(1) = 5.3 * 9.4 * wcoef * sqrt(z_obs / z_0); % gamma Eq. A15
   scoef(2) = 9.4 * gravity * z_obs;
   scoef(3) = scoef(1) * sqrt(gravity * z_obs);

   % SOURCE TERM LINEARIZATION VECTORS
   Sc = zeros(JJ_therm, 1);
   Sp = zeros(JJ_therm, 1);

   % STATE VARIABLES AND PARAMETERS NEEDED ON THE FIRST ITERATION
   fopts = opts.fzero;
   f_min = opts.f_ice_min;
   liqresid = opts.liqresid;

   roL = roLs;
   xTsfc = T(1);
   xf_liq = f_liq;
   liqflag = false;
   Qc =  CONDUCT(k_eff, T, dz, xTsfc);
   % zD = sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1)));

   % INITIALIZE THE OUTPUT STRUCTURES
   for n = 1:numel(opts.vars1) % ice1 = 1-d data
      ice1.(opts.vars1{n}) = nan(maxiter, 1);
   end
   for n = 1:numel(opts.vars2) % ice2 = 2-d data
      ice2.(opts.vars2{n}) = nan(JJ_therm, maxiter);
   end

   % diags.Tflag = false(maxiter,1);
   % diags.LCflag = false(maxiter,1);
   % A1 = ones(JJ_therm,1);
   % A = spdiags([A1,A1,A1],-1:1,JJ_therm,JJ_therm);
end
