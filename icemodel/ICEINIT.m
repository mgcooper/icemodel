function [ice1, ice2, T, f_ice, f_liq, k_eff, fn, dz, delz, roL, liqflag, ...
      Ts, JJ, Sc, Sp, Fc, Fp, TL, TH, f_ell_min, f_ell_max, f_ice_min, ...
      f_liq_res, ro_iwe, ro_wie] = ICEINIT(opts, tair)
   %ICEINIT initialize the 1-d ice column
   %
   %#codegen

   debug = false;

   % LOAD PHYSICAL CONSTANTS AND PARAMETERS
   [cp_ice, cp_liq, fcp, Lf, ro_ice, ro_liq, k_ice, k_liq, Tf, Ls, Rv, roLs] ...
      = icemodel.physicalConstant( ...
      'cp_ice','cp_liq', 'fcp','Lf', 'ro_ice','ro_liq', 'k_ice', 'k_liq', ...
      'Tf', 'Ls','Rv','roLs');

   % GENERATE A THERMAL MESH
   dz_therm = opts.dz_thermal;
   z0_therm = opts.z0_thermal;
   [dz, delz, ~, ~, fn] = CVMESH(z0_therm, dz_therm);

   % NUMBER OF TIMESTEPS TO INITIALIZE OUTPUTS
   maxiter = numel(tair) / opts.numyears;

   % NUMBER OF NODES IN THE THERMAL MESH
   JJ = numel(dz);

   % INITIALIZE DENSITIES
   ro_glc = (ro_ice + ro_liq) / 2; % [kg/m3]
   if opts.use_ro_glc == true
      ro_ice = ro_glc;
      ro_liq = ro_glc;
   end
   ro_iwe = ro_ice / ro_liq;
   ro_wie = ro_liq / ro_ice;

   % LOWER AND UPPER MELT ZONE TEMPERATURE AND WATER (MASS) FRACTION
   TL = Tf - (2.0 * Lf / (fcp ^ 2.0 * cp_ice)) ^ (1.0 / 3.0); % Eq 120
   TH = Tf - cp_liq / (Lf * 2.0 * fcp ^ 2.0);
   f_ell_min = 1 / (1 + (fcp * (Tf - TL)) ^ 2.0);
   f_ell_max = 1 / (1 + (fcp * (Tf - TH)) ^ 2.0);

   % INITIALIZE ICE TEMPERATURE

   % Initialize with a physically scaled exponential profile anchored to the
   % top node so T(1) = T_ref exactly. Use a slight cold offset from air to
   % avoid a zero-gradient startup.
   T_deep = min(TL - 1, Tf + opts.T_ice_init);
   T_ref = min(TL - 1, min(tair(1), Tf) - 1.0);
   Z = cumsum(dz) - dz / 2;
   omega_yr = 2.0 * pi / (365.0 * 86400.0);
   kappa_ice = k_ice / (ro_ice * cp_ice);
   z_scale = sqrt(2.0 * kappa_ice / omega_yr);
   T = T_deep + (T_ref - T_deep) .* exp(-(Z - Z(1)) ./ z_scale);

   % INITIALIZE LIQUID/ICE WATER FRACTION (f) AND BULK DENSITIES (g)
   T_dep = Tf - T;                           % [K]
   f_ell = 1 ./ (1 + (fcp * T_dep) .^ 2);    % [-], f_ell = liquid mass fraction
   g_ice = opts.ro_ice_init;
   g_liq = g_ice .* (f_ell ./ (1 - f_ell));
   f_liq = g_liq ./ ro_liq;
   f_ice = g_ice ./ ro_ice .* ones(JJ, 1);

   % INITIAL THERMAL CONDUCTIVITY
   k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);

   % SOURCE TERM LINEARIZATION VECTORS
   Sc = zeros(JJ, 1);
   Sp = zeros(JJ, 1);

   % BOUNDARY FLUX LINEARIZATION SCALARS
   Fc = 1;
   Fp = 1;

   % STATE VARIABLES AND PARAMETERS NEEDED ON THE FIRST ITERATION
   f_ice_min = opts.f_ice_min;
   f_liq_res = opts.f_liq_resid;

   Ts = (min(tair(1), Tf) + T(1)) / 2;
   roL = roLs;
   liqflag = false;
   % zD = sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1)));

   % INITIALIZE THE OUTPUT STRUCTURES
   for n = 1:numel(opts.vars1) % ice1 = 1-d data
      ice1.(opts.vars1{n}) = nan(maxiter, 1);
   end
   for n = 1:numel(opts.vars2) % ice2 = 2-d data
      if strcmp(opts.vars2{n}, 'lcflag')
         ice2.(opts.vars2{n}) = false(JJ, maxiter);
      else
         ice2.(opts.vars2{n}) = nan(JJ, maxiter);
      end
   end

   % diags.Tflag = false(maxiter,1);
   % diags.LCflag = false(maxiter,1);
   % A1 = ones(JJ_therm,1);
   % A = spdiags([A1,A1,A1],-1:1,JJ_therm,JJ_therm);

   function plot_T_init(T, Z, Ta, Ts)
      figure; hold on
      plot(T, Z)
      scatter(Ts, 0, 'filled')
      scatter(Ta, 0, 'filled')
      set(gca, 'YDir', 'reverse')
      legend('T', 'Ts', 'Ta')
   end
   if debug == true
      plot_T_init(T, Z, tair(1), Ts);
   end
end
