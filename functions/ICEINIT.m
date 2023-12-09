function [ice1, ice2, T, f_ice, f_liq, k_eff, fn, dz, delz, roL, liqflag, ...
      Ts, JJ_therm, JJ_spect, z_therm, dz_therm, dz_spect, Sc, Sp, Fc, Fp, ...
      TL, TH, flmin, flmax, f_min, liqresid, ro_iwe, ro_wie] = ICEINIT(opts, tair)
   %ICEINIT initialize the 1-d ice column
   %
   %#codegen

   % LOAD PHYSICAL CONSTANTS AND PARAMETERS
   [cp_ice, cp_liq, fcp, Lf, ro_ice, ro_liq, k_liq, Tf, Ls, Rv, roLs] ...
      = icemodel.physicalConstant( ...
      'cp_ice','cp_liq', 'fcp','Lf', 'ro_ice','ro_liq', 'k_liq','Tf', ...
      'Ls','Rv','roLs');

   % COMPUTE THE # OF NODES IN THE THERMAL AND SPECTRAL MESHES
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

   % LOWER AND UPPER MELT ZONE TEMPERATURE AND WATER FRACTIONS
   TL = Tf - (2.0 * Lf / (fcp ^ 2.0 * cp_ice)) ^ (1.0 / 3.0); % Eq 120
   TH = Tf - cp_liq / (Lf * 2.0 * fcp ^ 2.0);
   flmin = 1 / (1 + (fcp * (Tf - TL)) ^ 2.0);
   flmax = 1 / (1 + (fcp * (Tf - TH)) ^ 2.0);

   % INITIALIZE ICE TEMPERATURE TO AIR TEMPERATURE
   T = min(TL - 1, (tair(1) + 1) * ones(JJ_therm, 1));

   % INITIALIZE LIQUID/ICE WATER FRACTION AND BULK DENSITIES
   T_dep = Tf - T;                           % [K]
   fmliq = 1 ./ (1 + (fcp * T_dep) .^ 2);    % [-]
   g_ice = opts.ro_snow_i;
   g_liq = g_ice .* (fmliq ./ (1 - fmliq));
   f_liq = g_liq ./ ro_liq;
   f_ice = g_ice ./ ro_ice .* ones(JJ_therm, 1);

   % INITIAL THERMAL CONDUCTIVITY
   k_eff = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, Ls, Rv, Tf);

   % SOURCE TERM LINEARIZATION VECTORS
   Sc = zeros(JJ_therm, 1);
   Sp = zeros(JJ_therm, 1);

   % BOUNDARY FLUX LINEARIZATION SCALARS
   Fc = 1;
   Fp = 1;

   % STATE VARIABLES AND PARAMETERS NEEDED ON THE FIRST ITERATION
   f_min = opts.f_ice_min;
   liqresid = opts.liqresid;

   Ts = T(1);
   roL = roLs;
   liqflag = false;
   % zD = sqrt(k_eff(1)*dt/(ro_sno(1)*cp_sno(1)));

   % INITIALIZE THE OUTPUT STRUCTURES
   for n = 1:numel(opts.vars1) % ice1 = 1-d data
      ice1.(opts.vars1{n}) = nan(maxiter, 1);
   end
   for n = 1:numel(opts.vars2) % ice2 = 2-d data
      if strcmp(opts.vars2{n}, 'lcflag')
         ice2.(opts.vars2{n}) = false(JJ_therm, maxiter);
      else
         ice2.(opts.vars2{n}) = nan(JJ_therm, maxiter);
      end
   end

   % diags.Tflag = false(maxiter,1);
   % diags.LCflag = false(maxiter,1);
   % A1 = ones(JJ_therm,1);
   % A = spdiags([A1,A1,A1],-1:1,JJ_therm,JJ_therm);

   %% CV balance notes

   % Note the volumetric water fraction is ~0.9:
   %
   % f_wat_1 = f_liq(1) + f_ice(1) * ro_ice / ro_liq
   % printf(f_wat_1, 10)
   %
   % And the maximum mass fraction (flmax) is ~1 but <1:
   % printf(flmax, 10)
   %
   % Thus the max volumetric water fraction (the water fraction at T = TH):
   % fliqmax = f_wat_1 * flmax
   % printf(fliqmax, 10)
   %
   % Cannot exceed the initial fraction (f_wat_1):
   % printf(f_wat_1 - fliqmax , 10)
   % assert(fliqmax < f_wat_1)
   %
   % These are not needed but demonstrate cvconvert:
   %
   % [fm_ice, fm_liq] = icemodel.cvconvert("volumefraction", "massfraction", ...
   %    dz(1), [ro_ice, ro_liq], f_ice(1), f_liq(1));
   % printf(fm_ice, 10)
   % printf(fm_liq, 10)
   % printf(fm_ice + fm_liq, 10)
end
