function k_vap = GETKVAPOR(T, Ls, Rv, Tf)
   %GETKVAPOR Compute the vapor heat diffusion coefficient
   %
   %  k_vap = GETKVAPOR(T, Ls, Rv, Tf)
   %
   % This follows Anderson (1976)
   %  k = Ls * De / (Rv * T) * de_s / dT
   %
   % where:
   %  De = 9e-5 * (T / Tf) ^ nd;
   %
   % with nd the 'temperature exponent', Anderson (1976) Eq. 3.13 (& Fig. 4.3)
   %
   % See also:

   % Coefficients over ice.
   A = 611.15; % reference vapor pressure [Pa]
   B = 22.452; % unitless coefficient
   C = 272.55; % unitless coefficient
   
   % Compute saturation vapor pressure and its derivative wrt temperature
   es = A * exp(B * (T - Tf) ./ (C + T - Tf));     % [Pa]
   des_dT = B * C * es ./ (C + T - Tf) .^ 2;       % [Pa K-1]
   
   % Compute vapor diffusivity and the vapor thermal diffusion coefficient
   De = 9e-5 * (T / Tf) .^ 14;                     % [m2 s-1]
   k_vap = Ls * De ./ (Rv * T) .* des_dT;          % [W m-1 K-1]
end
