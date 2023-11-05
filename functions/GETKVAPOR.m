function k_vap = GETKVAPOR(T, Ls, Rv, Tf)
   %GETKVAPOR Compute the vapor heat diffusion coefficient
   %
   %  k_vap = GETKVAPOR(T, Ls, Rv, Tf)
   %
   % This follows Anderson (1976)
   %  k = Ls * De / (Rv * T) * de / dT
   %
   % where:
   %  De = 9e-5 * (T / Tf) ^ nd;
   %
   % with nd the 'temperature exponent', Anderson (1976) eq. 3.13
   %
   % See also:

   % Coefficients over ice.
   A = 611.15; % reference vapor pressure [Pa]
   B = 22.452; % unitless coefficient
   C = 272.55; % unitless coefficient
   
   esi = A .* exp((B * (T - Tf)) ./ (C + T - Tf)); % [Pa]
   desi_dT = B * C * esi ./ ((C + T - Tf) .^ 2);   % [Pa K-1]
   De = 9e-5 * (T / Tf) .^ 14;                     % [m2 s-1]
   k_vap = Ls * De ./ (Rv * T) .* desi_dT;         % [W m-1 K-1]
   
   %k_vap = Ls * 9e-5 * (T / Tf) .^ 6 ./ (Rv * T) ...
   %   * B * C * esi ./ ((C + T - Tf) .^ 2);
   
   % % snthrm version:
   % Lv = 2.500e6;
   % iliq = T >= TL;
   % 
   % c1_ice = 8.047e9;                                  % [kg/m3 K]
   % c1_liq = 5.726e8;                                  % [kg/m3 K]
   % 
   % Tice = Ls ./ (Rv * T);                             % [-]
   % Tliq = Lv ./ (Rv * T(iliq));                       % [-]
   % 
   % CkT = c1_ice / T .^2 .* (Tice - 1) .* exp(-Tice);  % d_ro_vap_dT
   % 
   % CkT(iliq) = c1_liq / T(iliq) .^ 2 .* (Tliq - 1) .* exp(-Tliq);
   % 
   % k_vap = Ls * De .* CkT;
end
