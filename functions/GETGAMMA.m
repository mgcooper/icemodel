function [k_eff, k_vap] = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, varargin)
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
   %  ksnow / (ki_ref * ka_ref)
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

   switch nargin
      case 6
         % k_vap provided by an external model.
         k_vap = varargin{1};

      case 8
         % Compute dry snow vapor k: (Ls·De)/(Rv·T)(∂es/∂T-es/T)
         Ls = varargin{1};
         Rv = varargin{2};
         Tf = varargin{3};
         es = 611.15 * exp((22.452 * (T - Tf)) ./ (272.55 + T - Tf)); % [Pa]
         k_vap = Ls * 9e-5 * (T / Tf) .^ 6 ./ (Rv * T) ...
            .* (22.452 * 272.55 .* es ./ ((272.55 + T - Tf) .^ 2) ... % ∂es/∂T
            - es ./ T); % es/T
      otherwise
         % do not compute k_vap
         k_vap = 0;
   end

   % Combine the dry snow and liquid water values into a bulk wet-snow value:
   k_sno = f_liq .* k_liq + f_ice .* k_sno;

   % Combine the wet snow and vapor values:
   k_eff = (f_liq + f_ice) .* k_sno + (1.0 - f_liq - f_ice) .* k_vap;

   % Combine the dry snow, liquid water, and vapor values into a bulk value:
   % k_eff = f_liq .* k_liq + f_ice .* k_sno + (1.0 - f_liq - f_ice) .* k_vap;
end
