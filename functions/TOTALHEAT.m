function H = TOTALHEAT(T, f_liq, f_ice, cv_liq, cv_ice, roLf, roLs, Tf)
   %TOTALHEAT Compute total enthalpy (J/m3)   
   H = ( ...
      (cv_ice * f_ice + cv_liq * f_liq) .* (T - Tf) ...  % specific heat
      + roLf .* f_liq ...                                % latent heat
      + roLs .* (1.0 - f_liq - f_ice) ...                % vapor heat
      );
end
