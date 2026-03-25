function k_vap = GETKVAPOR(T, Ls, Rv, Tf)
   %GETKVAPOR Compute the vapor heat diffusion coefficient.
   %
   %  k_vap = GETKVAPOR(T, Ls, Rv, Tf)
   %
   % This follows Anderson (1976):
   %  k = Ls * De / (Rv * T) * de_s / dT
   %
   % where:
   %  De = De0 * (T / Tf) ^ nd;
   %
   % with nd the 'temperature exponent', Anderson (1976) Eq. 3.13 (& Fig. 4.3)
   %
   % See also: VAPORHEAT, GETGAMMA, icemodel.kernels.buckVaporPressure
   %
   %#codegen

   % Ambaum (2020) Rankine-Kirchhoff coefficients over ice
   persistent ai bi ci nd De0
   if isempty(ai)
      [ai, bi, ci, nd, De0] = icemodel.parameterLookup( ...
         'ai', 'bi', 'ci', 'nd', 'De0');
   end

   % Saturation vapor pressure over ice [Pa]: es = a * exp(b/T) * T^c
   es = ai * exp(bi ./ T) .* T .^ ci;

   % Derivative of es wrt temperature [Pa K-1]
   des_dT = es ./ T .* (ci - bi ./ T);

   % Vapor diffusivity [m2 s-1]
   De = De0 * (T / Tf) .^ nd;

   % Vapor thermal diffusion coefficient [W m-1 K-1]
   k_vap = Ls * De ./ (Rv * T) .* (des_dT - es ./ T);

   % % From GETGAMMA before I added k_vap to VAPORHEAT:
   %
   % % Compute snow vapor k: (Ls·De)/(Rv·T)(∂es/∂T-es/T)
   % es = 611.15 * exp((22.452 * (T - Tf)) ./ (272.55 + T - Tf)); % [Pa]
   % k_vap = Ls * 9e-5 * (T / Tf) .^ 6 ./ (Rv * T) ...
   %    .* (22.452 * 272.55 .* es ./ ((272.55 + T - Tf) .^ 2) ... % ∂es/∂T
   %    - es ./ T); % es/T
end
