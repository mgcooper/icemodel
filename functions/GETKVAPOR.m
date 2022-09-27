%--------------------------------------------------------------------------
%   COMPUTE THE VAPOR HEAT DIFFUSION COEFFICIENT
%--------------------------------------------------------------------------

function xk_vap = GETKVAPOR(T,Ls,Rv,Tf)

% Over ice.
    A           =   6.1115 * 100.0;     % reference vapor pressure [Pa]
    B           =   22.452;             % unitless coefficient
    C           =   272.55;             % unitless coefficient
% compute gamma
    nd          =   14; % 'temperature exponent', Anderson (1976) eq. 3.13
    deltaT      =   T-Tf;
    esi         =   A .* exp((B.*deltaT)./(C+deltaT));      % [Pa]
    desi_dT     =   esi .* B.*C./((C+deltaT).^2);           % [Pa K-1]
    De          =   9e-5 .* (T./Tf).^nd;               % [m2 s-1]
    xk_vap      =   Ls.*De./(Rv.*T).*desi_dT;      % [W m-1 K-1]

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Jordan's version:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     Lv          =   2.500e6;
%     i_liq       =   T_old==Tfp;
%     c1_ice      =   8.047e9;                                    % [kg/m3 K]
%     c1_liq      =   5.726e8;                                    % [kg/m3 K]
%     Tice        =   xLs./(Rv.*T_old);                           % [-]
%     Tliq        =   Lv./(Rv.*T_old(i_liq));                    % [-]
%     CkT         =   c1_ice./T_old.^2.*(Tice-1).*exp(-Tice);
%     CkT(i_liq)  =   c1_liq./T_old(i_liq).^2.*(Tliq-1).*exp(-Tliq);
%     gamma0      =   xLs.*De.*CkT;
%     gamma       =   10.*xLs.*De.*CkT;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Yen's diffusivity enhancement
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % increase the vapor diffusivity to account for ventilation in the very
% % low-density near-surface ice following Yen et al. 1981
%     Dep         =   10; % 6.74
%     xk_vap      =   Dep.*xLs.*De./(Rv.*T_old).*desi_dT;    % [W m-1 K-1]

% % equivalent to this:
% %   G           =   40e-3;
% %   De          =   De.*33.511.*sqrt(G+4.56e-4);
    
% for testing, use the density-dependent equation
%     G           =   40e-3;
%     De          =   9e-5 .* (T_old./Tfp).^nd;
%     ro_s        =   ro_snow_dry./1000;
%     idx         =   ro_s < 0.9;
%     De(idx)     =   De(idx)+0.0513.*(ro_s(idx)).^3.2.*G^0.615;
%     xk_vapor    =   xLs.*De./(Rv.*T_old).*desi_dT;          % [W m-1 K-1]
%     gamma       =   xk_vapor;
    
% Yosida, 1950: De = 0.7 to 1.0e-4 m2/s, about 4-5x larger than D0
% where D0 is water vapor diffusivity in air. Forced convection near
% the surface increases De substantially for unconsolidated snow, which
% is analogous to weatherign crust, which is totally unconsolidated.
% For this case, 