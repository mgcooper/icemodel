function [k_eff, k_vap] = GETGAMMA(T, f_ice, f_liq, ro_ice, k_liq, varargin)
   %GETGAMMA Compute effective thermal conductivity (gamma).
   %
   % gamma = ki + kv, where ki is thermal conductivity of ice and kv is vapor
   % diffusivity. gamma is used as ke in Patankar (e.g. eq. 4.9 pg. 45)
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
   %
   %#codegen

   % Compute snow thermal conductivity (Calonne 2019 Eq. 5)
   k_sno = GETKTHERMAL(T, f_ice, ro_ice);

   switch nargin
      case 6
         % k_vap provided by an external model.
         k_vap = varargin{1};

      case 8
         % Compute dry snow vapor k via GETKVAPOR
         Ls = varargin{1};
         Rv = varargin{2};
         Tf = varargin{3};
         k_vap = GETKVAPOR(T, Ls, Rv, Tf);

      otherwise
         % do not compute k_vap
         k_vap = 0;
   end

   % Combine the dry snow and liquid water values into a bulk wet-snow value:
   k_sno = f_liq .* k_liq + f_ice .* k_sno;

   % % Combine the wet snow and vapor values:
   k_eff = (f_liq + f_ice) .* k_sno + (1.0 - f_liq - f_ice) .* k_vap;

   % Combine the dry snow, liquid water, and vapor values into a bulk value:
   % k_eff = f_liq .* k_liq + f_ice .* k_sno + (1.0 - f_liq - f_ice) .* k_vap;
end
