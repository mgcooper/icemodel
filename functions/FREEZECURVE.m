function [dFdT, f_wat] = FREEZECURVE(T, f_liq, f_ice, ro_iwe, Tf, fcp) %#codegen
   %FREEZECURVE Derivative of water fraction with respect to temperature
   %
   %  [dFdT, f_wat] = FREEZECURVE(T, f_liq, f_ice, ro_iwe, Tf, fcp)
   %
   % See also: MELTCURVE

   T_dep = Tf - min(T, Tf);

   % In terms of volumetric liquid fraction:
   f_wat = f_liq + ro_iwe * f_ice;
   dFdT = 2.0 * fcp ^ 2 * T_dep .* f_wat ...
      ./ (1.0 + fcp ^ 2 * T_dep .^ 2) .^ 2;
   
   % In terms of the volumetric liquid mass fraction f_ell = f_liq / f_wat:
   % dFdT = 2.0 * fcp ^ 2.0 * Tdep ./ (1.0 + (fcp * Tdep) .^ 2.0) .^ 2.0;
end
