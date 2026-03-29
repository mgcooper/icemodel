function [H, dHdT, dLdT, dVdT, k_vap] = TOTALHEAT2(T, f_ice, f_liq, cv_ice, ...
      cv_liq, ro_liq, Lf, Ls, Tf)
   %TOTALHEAT Compute total enthalpy (J/m3) and it's temperature derivative
   %
   %#codegen

   % Experimental helper. Keep the calculation shape aligned with ICEENBAL.
   fcp = icemodel.parameterLookup('fcp');
   ro_ice = icemodel.physicalConstant('ro_ice');
   f_wat = f_liq + f_ice * ro_ice / ro_liq;

   % Compute the derivative of enthalpy wrt temperature
   dHdT = cv_ice * f_ice + cv_liq * f_liq;
   dLdT = 2.0 * fcp ^ 2.0 * (Tf - min(T, Tf)) .* f_wat ...
      ./ (1.0 + fcp ^ 2.0 * (Tf - min(T, Tf)) .^ 2.0) .^ 2.0;

   % Compute vapor density, its derivative wrt temperature, and k_vap
   [ro_vap, dVdT] = VAPORDENSITY(T, f_liq);
   k_vap = VAPORK(T, f_liq, dVdT);

   % Compute total enthalpy wrt the reference temperature.
   H = ( ...
      dHdT .* (T - Tf) ...                               % specific heat
      + Lf * ro_liq .* f_liq ...                         % latent heat
      + Ls * ro_vap .* (1.0 - f_liq - f_ice) ...         % vapor heat
      );
end
