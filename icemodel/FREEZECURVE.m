function [dFdT, f_wat] = FREEZECURVE(T, ro_ice, ro_liq, fcp, Tf, ...
      f_ice, f_liq, f_wat)
   %FREEZECURVE Derivative of water fraction with respect to temperature
   %
   %  [DFDT, F_WAT] = FREEZECURVE(T, RO_ICE, RO_LIQ, FCP, TF, F_ICE, F_LIQ)
   %  [DFDT, F_WAT] = FREEZECURVE(..., F_WAT)
   %
   % If F_WAT is provided, the derivative is evaluated against that total water
   % fraction instead of recomputing it from F_ICE and F_LIQ. This is useful
   % inside melt-zone corrector iterations where F_WAT is held fixed.
   %
   % See also: MELTCURVE
   %
   %#codegen

   T_dep = Tf - min(T, Tf);

   if nargin < 8 || isempty(f_wat)
      % In terms of volumetric liquid fraction:
      f_wat = f_liq + f_ice * ro_ice / ro_liq;
   end
   dFdT = 2.0 * fcp ^ 2 * T_dep .* f_wat ...
      ./ (1.0 + fcp ^ 2 * T_dep .^ 2) .^ 2;
   dFdT = max(dFdT, sqrt(eps));

   % In terms of the volumetric liquid mass fraction f_ell = f_liq / f_wat:
   % dFdT = 2.0 * fcp ^ 2.0 * Tdep ./ (1.0 + (fcp * Tdep) .^ 2.0) .^ 2.0;
end
