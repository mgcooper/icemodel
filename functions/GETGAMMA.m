function [k_eff,k_vap] = GETGAMMA(T,f_liq,f_ice,ro_ice,k_liq,Ls,Rv,Tf)
   %GETGAMMA Compute thermal conductivity (gamma)
   %
   % gamma = ki + kv, where ki is thermal conductivity of ice and kv is vapor
   % diffusivity. gamma is used as ke in Patankar (e.g. eq. 4.9 pg. 45)
   %
   % See also: GETKTHERMAL, VAPORHEAT

   % compute snow thermal k
   gam_ice  =   ro_ice.*f_ice;
   th       =   1./(1+exp(-0.04.*(gam_ice-450)));
   kfirn    =   2.107 + 0.003618.*(gam_ice-917);
   ksnow    =   0.024 - 1.23e-4.*gam_ice + 2.5e-6.*gam_ice.^2;
   kiceT    =   9.828 .* exp(-5.7e-3.*T);
   k_sno    =  (1-th).*0.47461.*kiceT.*ksnow+th.*kiceT./2.107.*kfirn;

   % compute snow vapor k
   esi      =   611.15.*exp((22.452.*(T-Tf))./(272.55+T-Tf));   % [Pa]
   desi_dT  =   esi.*22.452.*272.55./((272.55+T-Tf).^2);        % [Pa K-1]
   De       =   9e-5 .* (T./Tf).^14;                              % [m2 s-1]
   k_vap    =   Ls.*De./(Rv.*T).*desi_dT;                         % [W m-1 K-1]

   % combine them into a bulk value
   k_sno    =   (f_liq.*k_liq) + (f_ice.*k_sno);
   k_eff    =   (f_liq+f_ice) .* k_sno + (1.0-f_liq-f_ice).* k_vap;
end

% EXPLICIT VERSION:

% function k_eff = GETGAMMA(T,f_ice,f_liq,k_liq,ro_ice,Ls,Rv,Tf)
% % gamma = ki + kv, where ki is thermal conductivity of ice and kv is vapor
% % diffusivity. gamma is used as ke in Patankar (e.g. eq. 4.9 pg. 45)
%
% % Over ice.
%     A           =   611.15;         % reference vapor pressure [Pa]
%     B           =   22.452;         % unitless coefficient
%     C           =   272.55;         % unitless coefficient
%
% % compute snow thermal k
%     gam_ice     =   ro_ice.*f_ice;
%     th          =   1./(1+exp(-0.04.*(gam_ice-450)));
%     k_firn      =   2.107 + 0.003618.*(gam_ice-917);               %ref
%     k_snow      =   0.024 - 1.23e-4.*gam_ice + 2.5e-6.*gam_ice.^2; %ref
%     k_iceT      =   9.828 .* exp(-5.7e-3.*T);
%     ki_ref      =   2.107;
%     ka_ref      =   0.024;
%     k_airT      =   0.024;
%     k_sno       =   (1-th).*(k_iceT.*k_airT)./(ki_ref*ka_ref).*k_snow ...
%                         + th.*(k_iceT./ki_ref).*k_firn;
% % compute snow vapor k
%     nd          =   14; % 'temperature exponent', Anderson (1976) eq. 3.13
%     deltaT      =   T-Tf;
%     esi         =   A .* exp((B.*deltaT)./(C+deltaT));      % [Pa]
%     desi_dT     =   esi .* B.*C./((C+deltaT).^2);           % [Pa K-1]
%     De          =   9e-5 .* (T./Tf).^nd;               % [m2 s-1]
%     k_vap       =   Ls.*De./(Rv.*T).*desi_dT;          % [W m-1 K-1]
%
% % combine them into a bulk value
%     k_sno       =   (f_liq .* k_liq) + (f_ice .* k_sno);
%     k_eff       =   (f_liq+f_ice) .* k_sno + (1.0-f_liq-f_ice).* k_vap;
