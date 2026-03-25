function k_sno = GETKTHERMAL(T, f_ice, ro_ice)
   %GETKTHERMAL Compute snow/ice thermal conductivity (Calonne 2019 Eq. 5).
   %
   %  k_sno = GETKTHERMAL(T, f_ice, ro_ice)
   %
   %  Computes the temperature-dependent thermal conductivity of snow/ice
   %  using Calonne et al. (2019) Eq. 5, which blends a low-density (snow)
   %  and high-density (firn) parameterization via a logistic transition.
   %
   %  Inputs:
   %     T      - Temperature [K]
   %     f_ice  - Ice volume fraction [1]
   %     ro_ice - Intrinsic density of ice [kg m-3]
   %
   %  Output:
   %     k_sno  - Snow/ice thermal conductivity [W m-1 K-1]
   %
   % See also: GETGAMMA, icemodel.kernels.snowThermalConductivity
   %
   %#codegen

   % Bulk density of ice (dry snow density) [kg m-3]
   g_ice = ro_ice * f_ice;

   % Logistic transition between snow and firn regimes
   theta = 1 ./ (1 + exp(-0.04 * (g_ice - 450.0)));

   % Reference conductivities [W m-1 K-1]
   kfirn = 2.107 + 0.003618 * (g_ice - 917.0);
   ksnow = 0.024 - 0.000123 * g_ice + 2.5e-6 * g_ice .^ 2;

   % Temperature-dependent ice thermal conductivity [W m-1 K-1]
   kiceT = 9.828 * exp(-0.0057 * T);

   % Blended thermal conductivity [W m-1 K-1]
   k_sno = 0.47461 * (1 - theta) .* kiceT .* ksnow ...
      + theta .* kiceT / 2.107 .* kfirn;
end
