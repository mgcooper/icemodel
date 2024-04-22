function [ro_vap, dro_vapdT, k_vap] = VAPORHEAT(T, f_liq, f_ice, Tf, Rv, Ls)
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
   % De0s = 9.2e-5 [m2 s-1] is the reference effective diffusion coefficient for
   % water vapor in snow at 1000 mb and 0oC. (See Jordan pg v. Nomenclature).
   %
   % De0g = 1.61e-5 * porosity is the reference effective diffusion coefficient
   % for water vapor in soil at 1000 mb and 0oC.
   %
   % See also:

   % Define coefficients over water and ice
   persistent aw bw cw ai bi ci nd
   if isempty(aw)
      aw = 611.21;
      bw = 17.502;
      cw = 240.97;
      ai = 611.15;
      bi = 22.452;
      ci = 272.55;
      nd = 14;
   end

   % Locate the indices with and without water
   iM = f_liq > 0.02;

   % Saturation vapor pressure over ice [Pa]
   es = ai * exp(bi * (T - Tf) ./ (ci + T - Tf));

   % Saturation vapor pressure over water [Pa]
   if sum(iM) > 0
      es(iM) = aw * exp(bw * (T(iM) - Tf) ./ (cw + T(iM) - Tf));
   end

   % Equilibrium water vapor density wrt phase k: [kg m-3]
   ro_vap = es ./ (Rv * T);

   % Derivative of vapor density wrt to temperature over ice [kg m-3 K-1]
   dro_vapdT = ro_vap .* (bi * ci ./ (ci + T - Tf) .^ 2 - 1 ./ T);

   % Derivative of vapor density wrt to temperature over water
   if sum(iM) > 0
      dro_vapdT(iM) = ro_vap(iM) .* (bw * cw ./ ...
         (cw + T(iM) - Tf) .^ 2 - 1 ./ T(iM));
   end

   % Vapor diffusivity [m2 s-1]
   % De = 9.0e-5 * (T / Tf) .^ 6;

   % Vapor thermal diffusion coefficient [W m-1 K-1]: Ls * De * dro_vapdT
   k_vap = Ls * 9.0e-5 * (T / Tf) .^ nd .* dro_vapdT;

   % k_vap = Ls * 2.1664062e-19 * T .^ 14 .* dro_vapdT;

   % Below here not currently implemented, but would be used for grain growth

   % Bulk vapor density (Eq. 18) (assume f_rh = 1.0) [kg m-3]
   % bd_vap = ro_vap .* (1.0 - f_liq - f_ice);

   % Derivative of bulk vapor density wrt temperature [kg m-3 K-1]
   % dbd_vapdT = dro_vapdT .* (1.0 - f_liq - f_ice);

   % Vapor heat (enthalpy) per unit volume:
   % H_vap = Ls * bd_vap;

   % The change in stored heat due to water vapor diffusion (Eq. 74, term 2)
   % dH_vap_dT = Ls * f_air .* f_rh .* d_ro_vap_dT; % [J/K/m3]

   % Diffusion of water vapor [kg/m^2/s]
   % From here ... need dT/dz
   % U_vap = -De .* dro_vapdT .* dTdz;
   % U_vap = -k_vap ./ Ls .* dTdz; % equivalent to above
   % U_vap = -De .* dro_vapdz; % if dro_vapdz = dro_vapdT * dTdz is valid

   % Time derivative of grain growth, but only for dry snow:
   % dgdt = g1 ./ d .* De .* dro_vapdT .* dTdz;
   % dgdt = g1 ./ d(n) .* Uv(n);
   % for n = 1:numel(f_liq)
   %
   %    % For wet snow:
   %    if 0.02 <= f_liq(n) && f_liq(n) < 0.09
   %
   %       dgdt(n) = g2 / d(n) * (f_liq(n) + 0.05);
   %
   %    elseif 0.09 <= f_liq(n)
   %
   %       dgdt(n) = g2 / d(n) * 0.14;
   %    else
   %       % For dry snow:
   %       dgdt(n) = g1 ./ d(n) .* Uv(n);
   %    end
   % end
end
