function [H_vap, dro_vapdT, ro_vap, bd_vap] = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls)
   %VAPORHEAT Compute saturation vapor density within porous ice
   %
   % [H_vap, dro_vapdT, ro_vap, bd_vap] = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls)
   % computes the saturation vapor density within the ice and the latent heat
   % due to water vapor diffusion within the ice. 
   %
   % ro_vap    - Vapor mass per vapor volume, a function of temperature.
   % dro_vapdT - The derivative of ro_vap with respect to temperature.
   % H_vap     - Vapor enthalpy per unit volume [J / m3]
   % k=i       - wrt ice 
   % k=l       - wrt liq
   %
   % ro_vap is analogous to an intrinsic density. The calculations assume the
   % air voids are saturated with respect to water vapor.
   % 
   % See also:

   % Locate the indices with and without water
   iM = f_liq > 0.02;

   % Saturation vapor pressure over ice (Pvk_sat in Jordan, but not used) [Pa]
   es = 611.15 * exp((22.452 * (T - Tf)) ./ (272.55 + (T - Tf)));

   % Saturation vapor pressure over water (Pvk_sat in Jordan) [Pa]
   if sum(iM) > 0
      es(iM) = 611.21 * exp(17.502 * (T(iM) - Tf) ./ (240.97 + (T(iM) - Tf)));
   end

   % Equilibrium water vapor density wrt phase k: [kg m-3] (ro_vk_sat)
   ro_vap = es ./ (Rv * T);

   % Diffusion of water vapor [kg/m^2/s]
   % De = 9.0e-5 * (T / Tf) .^ 14; % [m2 s-1]

   % Vapor diffusivity [W m-1 K-1]
   % k_vap = Ls * 9.0e-5 * (T / Tf) .^ 14  ...
   %    * 22.452 * 272.55 .* ro_vap ./ (272.55 + T - Tf) .^ 2; % des/dT [Pa K-1]
   % TODO: Add option to return one argument k_vap 
   
   % Bulk vapor density (Eq. 18) (assume f_rh = 1.0) % [kg m-3]
   bd_vap = (1.0 - f_liq - f_ice) .* ro_vap;

   % Derivative of vapor density wrt to temperature over ice (Eq. 20)
   dro_vapdT = ro_vap.*(22.452 * 272.55 ./ ((272.55 + (T - Tf)) .^ 2) - 1 ./ T);

   % dbd_vapdT = dro_vapdT * (1.0 - f_liq - f_ice);

   % Derivative of vapor density wrt to temperature over water
   if sum(iM) > 0
      dro_vapdT(iM) = ro_vap(iM) .* (17.502 * 240.97 ./ ...
         ((240.97 + (T(iM) - Tf)) .^ 2) - 1 ./ T(iM));
   end

   % Vapor heat (enthalpy) per unit volume:
   H_vap = Ls * bd_vap;

   % The change in stored heat due to water vapor diffusion (Eq. 74, term 2)
   % H_vap = Ls * f_air .* f_rh .* d_ro_vap_dT;           % [J/K/m3]

   % Vapor thermal diffusion coefficient % [W m-1 K-1]
   % k_vap = ro_vap .* Ls .* De .* 22.452 * 272.55 ./ (272.55 + Td) .^ 2;

   % equivalently:
   % k_vap = Ls * De .* (d_ro_vap_dT + ro_vap ./ T);
end
