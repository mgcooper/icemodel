function k_ice = THERMALK(T, f_ice, ro_ice)
   %THERMALK Compute porous ice thermal conductivity (Calonne 2019 Eq. 5).
   %
   %  k_ice = THERMALK(T, f_ice, ro_ice)
   %
   %  Computes the temperature-dependent thermal conductivity of snow/firn/ice
   %  using Calonne et al. (2019) Eq. 5, which blends a low-density (snow) and
   %  high-density (firn) parameterization via a logistic transition function.
   %
   %  The parameterization is:
   %
   %     k = (kiceT / ki_ref) * [(1-theta) * ksnow + theta * kfirn]
   %
   %  where:
   %     ro     = ro_ice * f_ice                bulk ice density [kg m-3]
   %     theta  = 1 / (1 + exp(-a*(ro-ro_mid))) logistic snow/firn transition [1]
   %     kfirn  = ki_ref + c_firn*(ro-ro_ref)   firn conductivity [W m-1 K-1]
   %     ksnow  = ka_ref + c1*ro + c2*ro^2      snow conductivity [W m-1 K-1]
   %     kiceT  = c_kT * exp(-c_kT_exp * T)     temperature-dependent ice k [W m-1 K-1]
   %
   %  All parameters are loaded from icemodel.parameterLookup (cal_* namespace).
   %  The (kiceT / ki_ref) factor rescales the density-dependent blended
   %  conductivity to the actual ice thermal conductivity at the local
   %  temperature, since kfirn and ksnow were regressed at a reference
   %  temperature -3°C where kiceT = ki_ref = 2.107 W m-1 K-1.
   %
   %  Inputs:
   %     T      - Temperature [K]
   %     f_ice  - Ice volume fraction [1]
   %     ro_ice - Intrinsic density of ice [kg m-3]
   %
   %  Output:
   %     k_ice  - Snow/ice thermal conductivity [W m-1 K-1]
   %
   %  Reference:
   %     Calonne et al. (2019), "Thermal Conductivity of Snow, Firn, and
   %     Porous Ice from 3-D Image-Based Computations." Geophysical Research
   %     Letters, 46(22), 13079-13089.
   %
   % See also: BULKTHERMALK, icemodel.parameterLookup
   %
   %#codegen

   persistent ki_ref ka_ref ro_ref ro_mid a c_firn c_sno1 c_sno2 c_kT c_kT_exp
   if isempty(ki_ref)
      [ki_ref, ka_ref, ro_ref, ro_mid, a, c_firn, c_sno1, c_sno2, ...
         c_kT, c_kT_exp] = icemodel.parameterLookup( ...
         'cal_ki_ref', 'cal_ka_ref', 'cal_rho_ref', 'cal_rho_mid', ...
         'cal_a_theta', 'cal_c_firn', 'cal_c_sno1', 'cal_c_sno2', ...
         'cal_c_kT', 'cal_c_kT_exp');
   end

   % Bulk density of ice (dry snow density) [kg m-3]
   g_ice = ro_ice * f_ice;

   % Logistic transition between snow and firn regimes
   theta = 1 ./ (1 + exp(-a * (g_ice - ro_mid)));

   % Density-dependent reference conductivities [W m-1 K-1]
   kfirn = ki_ref + c_firn * (g_ice - ro_ref);
   ksnow = ka_ref + c_sno1 * g_ice + c_sno2 * g_ice .^ 2;

   % Temperature-dependent ice thermal conductivity [W m-1 K-1]
   kiceT = c_kT * exp(-c_kT_exp * T);

   % Blended thermal conductivity, rescaled by temperature [W m-1 K-1]
   k_ice = (kiceT / ki_ref) .* ((1 - theta) .* ksnow + theta .* kfirn);
end
