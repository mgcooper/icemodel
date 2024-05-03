function [dFdT, f_wat] = FREEZECURVE(T, f_ice, f_liq, ro_ice, ro_liq, fcp, Tf)
   %FREEZECURVE Derivative of water fraction with respect to temperature
   %
   %  [DFDT, F_WAT] = FREEZECURVE(T, F_ICE, F_LIQ, RO_ICE, RO_LIQ, FCP, TF)
   %
   % See also: MELTCURVE
   %
   %#codegen

   T_dep = Tf - min(T, Tf);

   % In terms of volumetric liquid fraction:
   f_wat = f_liq + f_ice * ro_ice / ro_liq;
   dFdT = 2.0 * fcp ^ 2 * T_dep .* f_wat ...
      ./ (1.0 + fcp ^ 2 * T_dep .^ 2) .^ 2;

   % In terms of the volumetric liquid mass fraction f_ell = f_liq / f_wat:
   % dFdT = 2.0 * fcp ^ 2.0 * Tdep ./ (1.0 + (fcp * Tdep) .^ 2.0) .^ 2.0;
end
