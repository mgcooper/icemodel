function [dFdT, f_wat] = liquid_fraction_derivative(T_ice, f_ice, f_liq, f_wat)
   %LIQUID_FRACTION_DERIVATIVE Liquid fraction derivative wrt temperature.
   %
   %  [DFDT, F_WAT] = liquid_fraction_derivative(T_ICE, F_ICE, F_LIQ)
   %  [DFDT, F_WAT] = liquid_fraction_derivative(..., F_WAT)
   %
   % If you supply F_WAT, the function evaluates the derivative against that
   % total water fraction and does not recompute it from F_ICE and F_LIQ.
   % Melt-zone corrector iterations use this form because they hold F_WAT fixed.
   %
   % See also: icemodel.column.liquid_fraction_function
   %
   %#codegen

   persistent Tf fcp
   if isempty(Tf)
      Tf = icemodel.physicalConstant('Tf');
      fcp = icemodel.parameterLookup('fcp');
   end

   T_dep = Tf - min(T_ice, Tf);

   if nargin < 4 || isempty(f_wat)
      % In terms of volumetric liquid fraction:
      f_wat = icemodel.column.water_fraction(f_ice, f_liq);
   end

   % Differentiate the liquid fraction-temperature relationship [K^-1].
   dFdT = 2.0 * fcp ^ 2 * T_dep .* f_wat ./ (1.0 + fcp ^ 2 * T_dep .^ 2) .^ 2;
   dFdT = max(dFdT, sqrt(eps));

   % In terms of the volumetric liquid mass fraction f_ell = f_liq / f_wat:
   % dFdT = 2.0 * fcp ^ 2.0 * Tdep ./ (1.0 + (fcp * Tdep) .^ 2.0) .^ 2.0;
end
