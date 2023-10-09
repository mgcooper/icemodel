function [H_vap, dro_vapdT, ro_vap, bd_vap] = VAPORHEAT(T,f_liq,f_ice,Tf,Rv,Ls)
   %VAPORHEAT compute the saturation vapor density within the ice
   %
   % [H_vap,d_ro_vap_dT,ro_vap_sat,gamma_vap] =
   % VAPORHEAT(T,f_liq,f_ice,Tf,Rv,Ls) computes the saturation vapor density
   % within the ice and the release of latent heat due to water vapor diffusion
   % within the ice.
   %
   % k=i wrt ice, and l wrt liq, that's what the i_liq / i_ice do
   %
   % See also:

   % Locate the indices with and without water
   iM = f_liq > 0.02;

   % Saturation vapor pressure over ice (Pvk_sat in Jordan, but not used) [Pa]
   esi = 611.15.*exp((22.452.*(T-Tf))./(272.55+(T-Tf)));

   % Saturation vapor pressure over water (Pvk_sat in Jordan) [Pa]
   if sum(iM) > 0
      esi(iM) = 611.21.*exp(17.502.*(T(iM)-Tf)./(240.97+(T(iM)-Tf)));
   end
   % note: es0 is the surface value computed in VAPOR

   % equilibrium water vapor density wrt phase k: [kg m-3] (ro_vk_sat)
   ro_vap = esi ./ (Rv .* T);

   % bulk vapor density (Eq. 18) (assume f_rh = 1.0) % [kg m-3]
   bd_vap = (1.0-f_liq-f_ice).*ro_vap;

   % Compute the diffusion of water vapor [kg/m^2/s]
   %De = 9.0e-5.*(T./Tf).^14; % [m2 s-1]

   % Derivative of vapor density wrt to temperature over ice (Eq. 20)
   dro_vapdT = ro_vap.*(22.452*272.55./((272.55+(T-Tf)).^2)-1./T);

   % dbd_vapdT = dro_vapdT*(1.0-f_liq-f_ice);

   % Derivative of vapor density wrt to temperature over water
   if sum(iM) > 0
      dro_vapdT(iM) = ro_vap(iM).*(17.502*240.97./ ...
         ((240.97+(T(iM)-Tf)).^2)-1./T(iM));
   end

   % cast in terms of the actual heat (enthalpy) per unit volume:
   H_vap = Ls.*bd_vap;

   % the change in stored heat due to water vapor diffusion (Eq. 74, term 2)
   % H_vap = Ls.*f_air.*f_rh.*d_ro_vap_dT;           % [J/K/m3]

   % vapor thermal diffusion coefficient % [W m-1 K-1]
   % k_vap = ro_vap_sat.*Ls.*De.*22.452*272.55./(272.55+Td).^2;

   % equivalently:
   % k_vap = Ls.*De.*(d_ro_vap_dT + ro_vap_sat./T);
end

% NOTE: might change to dbvdT, for consistency with snthrm, i.e., convert
% ro_vap to bulk density, so the thermal coefficient aP0 is ro_sno*c_sno +
% Lf*dbldT + Ls*dbvdT, which is intuitive.

% ro_vap_sat is the 'intrinsic' density of water vapor, meaning the mass of
% water vapor per volume of water vapor, a function of temperature. This
% means ro_vap_sat is analogous to ro_ice/ro_water, not gamma_ice/water
% The derivative of ro_vap_sat wrt temperature dro_vap_sat_dT

% note: in Jordan, Pvk_sat is defined, but not used. In icemodel, this is
% es0 for the surface, and esi within the ice matrix over ice, and would be
% esl over liquid but I just call it esi. Similarly, Jordan defines and
% uses ro_vk_sat (ro_vap_sat wrt to constituent k), but I don't use it
% except here. It's easy to confuse ro_vk_sat and Pvk_sat in Jordan b/c the
% rho and P are similar.

% Pvk_sat       = es0 / esi         (could call this Pvap_sat)
% ro_vap_sat    = ro_vap_sat        (previously undefined in Glen's model)
% ro_vap_sat    = Pvk_sat/(Rv*T)    (conversion b/w esi and ro_vap_sat)
% CkT           = d_ro_vap_sat_dT   (previously undefined in Glen's model)
% desi_dT       = dPvk_sat_dT       (undefined in Jordan)
%                                   (could call this dPvap_sat_dT)
% Uvap          = -Deos*(T/Tfp)^n*CkT*dT/dz (see xk_vapor)
% for reference, 'kvap' meaning vapor conductivity from GETGAMMA:
% kvap          = De0*(T/Tfp)^n*Ls/(Rv*T)*desi/dT

% es0 = VAPOR = sat. water vapor pressure at the surface (wrt ice or water)
% ea  = VAPPRESS = same as above, but for the atmosphere above the surface
% esi = same as above, but over ice/water within the ice matrix, used in
% GETGAMMA

% to pick up from here ... need to figure out exactly how to replace
% GETGAMMA with this function, i.e., how to map Uvap above to what happens
% in GETGAMMA and/or vice versa,

% Coeffs for saturation vapor pressure over water (Buck 1981).
%   Note: temperatures for Buck's equations are in deg C, and
%   vapor pressures are in mb.  Do the adjustments so that the
%   calculations are done with temperatures in K, and vapor
%   pressures in Pa.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % In Jordan, it would be:
%     c1_ice      =   8.047e9;                                    % [kg/m3 K]
%     c1_liq      =   5.726e8;                                    % [kg/m3 K]
%     Tice        =   xLs./(Rv.*T_old);                           % [-]
%     Tliq        =   xLv./(Rv.*T_old);                           % [-]
% % over ice: (Eq. 20)
%     e           =   c1_ice./(T_old.^2).*(Tice-1).*exp(-Tice);   % [kg/m3/K]
% % over water:
%     e(i_liq)    =   c1_liq./(T_old.^2).*(Tliq-1).*exp(-Tliq);   % [kg/m3/K]
%     d_ro_vap_dT =   e;                                          % [kg/m3/K]
%     H_vap       =   xLs.*frac_air.*d_ro_vap_dT;                 % [J/K/m3]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Jordan calls d_ro_vap_dT 'CkT', where k refers to constituent ice, liq,
% Eq. 21: (I commented this out, I think drovap/dT is all I need)
%     dT          =   T_old(2:JJ) - T_old(1:JJ-1);
%     dT(JJ)      =   0.0;
%     dT_dz       =   dT ./ dy_p;
%     Uvapor      =   -De .* d_ro_vap_dT .* dT_dz;

% NEED TO PICK UP on equation 21, the remaining expressions for mass
% fluxes, involving melt/evap I think, and the condition that frac_liq be
% below a threshold for vapor diffusion to occur (no diffusion if
% saturated)

% A few notes though while this is fresh ... The quantity Uvapor is the
% diffusive vapor flux term in Eq. 17. I also have all the ingredients for
% the expulsion vapor flux, see equation 25:

%     Evapor      =   (porosity_sd - theta_liq) .* f_rh .* d_ro_vap_dT;

% This would become the coefficient in the tridiagonal matrix on T, whereas
% Uvapor would be differentiated i.e. Uv(j+1) - Uv(j-1) and put directly
% into the matrix (I think, probably a deltaz term needed). The remaining
% ingredient is the Mvi, Mvl terms, which are sublimation/condensation
% within the ice matrix.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% these notes were further up

% below I have a differet expression for Uvap but I think this is correct:
%   Uvap        =   -De0.*(1000/Pa).*(T_old./Tfp).^nd.*.*dT_dz;

% this is from GETGAMMA:
%    gamma_vap  =   De0.*(1000/Pa).*(T_old./Tfp).^nd.*xLv./Rv./T_old.*desi_dT;
% NOTE: this is a different gamma than the one above, which is

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Jordan notation:
% gamma_v = bulk vapor density = mass vapor / total volume (kg/m3)
% gamma_v = theta_v * ro_v = theta_v * frh * ro_vk_sat
% ro_vk_sat = equilibrium water vapor density wrt phase k (kg/m3)
% theta_v = fractional volume (m3/m3) (my frac_ice/frac_air/frac_liq)
