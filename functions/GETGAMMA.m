function [k_eff, k_vap] = GETGAMMA(T, f_liq, f_ice, ro_ice, k_liq, Ls, Rv, Tf)
   %GETGAMMA Compute thermal conductivity (gamma)
   %
   % gamma = ki + kv, where ki is thermal conductivity of ice and kv is vapor
   % diffusivity. gamma is used as ke in Patankar (e.g. eq. 4.9 pg. 45)
   %
   %  A = 611.15; % reference vapor pressure [Pa]
   %  B = 22.452; % unitless coefficient
   %  C = 272.55; % unitless coefficient
   %
   %  kfirn is a reference conductivity for firn
   %  ksnow is a reference conductivity for snow
   %  kiceT is a temperature-dependent ice thermal conductivity
   %  ki_ref (2.107) is ice thermal conductivity
   %  ka_ref (0.024) is air thermal conductivity
   %  k_airT (0.024) is temperature-dependent air thermal conductivity, taken
   %  equal to the reference air thermal conductivity
   %
   %  k_sno = (1-th).*(k_iceT.*k_airT)./(ki_ref*ka_ref).*k_snow ...
   %     + th.*(k_iceT./ki_ref).*k_firn;
   %
   % See also: GETKTHERMAL, VAPORHEAT

   % Compute snow thermal k
   g_ice = ro_ice * f_ice;
   theta = 1 ./ (1 + exp(-0.04 * (g_ice - 450.0)));
   kfirn = 2.107 + 0.003618 * (g_ice - 917.0);
   ksnow = 0.024 - 1.23e-4 * g_ice + 2.5e-6 * g_ice .^ 2;
   kiceT = 9.828 * exp(-5.7e-3 * T);
   k_sno = 0.47461 * (1 - theta) .* kiceT .* ksnow ...
      + theta .* kiceT / 2.107 .* kfirn;

   % Compute snow vapor k
   es = 611.15 * exp((22.452 * (T - Tf)) ./ (272.55 + T - Tf)); % [Pa]
   k_vap = Ls * 9e-5 * (T / Tf) .^ 14 ./ (Rv * T) ...
     * 22.452 * 272.55 .* es ./ ((272.55 + T - Tf) .^ 2); % des/dT [Pa K-1]
   
   % Combine them into a bulk value
   k_sno = f_liq .* k_liq + f_ice .* k_sno;
   k_eff = (f_liq + f_ice) .* k_sno + (1.0 - f_liq - f_ice) .* k_vap;
end
