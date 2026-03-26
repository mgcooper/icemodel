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
   % See also: VAPORHEAT, GETGAMMA, icemodel.kernels.buckVaporModel
   %
   %#codegen

   % Ambaum (2020) / Romps (2021) Rankine-Kirchhoff coefficients over ice
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
end
