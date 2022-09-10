clean
savedata    =  true;

%--------------------------------------------------------------------------
%   physical constants
%--------------------------------------------------------------------------
Lv          =  2.500e6;       % latent heat vapor,             [J kg-1]
Lf          =  3.34e5;        % latent heat fusion,            [J kg-1]
Ls          =  Lv+Lf;         % latent heat subl. (2.83e6)     [J kg-1]
Tf          =  273.15;        % freezing point,                [K]
ro_air      =  1.275;         % density air,                   [kg m-3]
ro_ice      =  917.0;         % density pure ice,              [kg m-3]
ro_liq      =  1000.0;        % density water,                 [kg m-3]
cp_air      =  1004;          % sp. heat cap. air,             [J kg-1 K-1]
cp_liq      =  4180.0;        % sp. heat cap. water,           [J kg-1 K-1]
cp_ice      =  2102.0;        % sp. heat cap. ice,             [J kg-1 K-1]
k_liq       =  0.552;         % therm. cond. water at 0oC,     [W m-1 K-1]
k_ice       =  2.10;          % therm. cond. ice at 0oC,       [W m-1 K-1]
Rv          =  461.0;         % gas constant water vapor       [J kg-1 K-1]
emiss       =  0.98;          % surface emissivity,            [-]
SB          =  5.6696e-8;     % stefan boltzman constant,      [W m-2 K-4]
gravity     =  9.81;          % gravitational accel.,          [m s-2]
kappa       =  0.4;           % von Karman constant,           [-]
one_atmos   =  101300.0;      % one atm,                       [Pa = J m-3 = kg m-1 s-2]
epsilon     =  0.622;         % Rd/Rv                          [-]       
es0         =  611;           % Ref sat vap press in C-C Eqn   [Pa]
S0          =  1367;          % Solar constant in              [W m-2]
N0          =  0.08;          % Marshall-Palmer parameter      [cm-4]
psychro     =  64.6;          % Psychrometric constant         [Pa K-1]
dalr        =  9.800;         % dry adiabatic lapse rate       [K/km]
malr        =  5.1;           % moist adiabatice lapse rate    [K/km]
cv_air      =  ro_air*cp_air; % volumetric heat cap dry air    [J/m3/K]
cv_ice      =  ro_ice*cp_ice; % volumetric heat cap ice        [J/m3/K]
cv_liq      =  ro_liq*cp_liq; % volumetric heat cap liquid h20 [J/m3/K]
fcp         =  100;           % freezing curve parameter
scale_ht    =  8500.0;
hrsperday   =  24;

% combined values
roLv        =  ro_air*Lv;
roLs        =  ro_air*Ls;
roLf        =  ro_ice*Lf;
ro_iwe      =  ro_ice./ro_liq;                   % ice water equivalent
ro_wie      =  ro_liq./ro_ice;                   % liq ice equivalent
emissSB     =  emiss*SB;
fcpsq       =  fcp*fcp;
TINY        =  1e-8;

if savedata == true
   clear savedata 
   save('PHYSCONS');
end

% function consts = PHYSICAL_CONSTS()
%
% consts.xLv         =   2.500e6;     % latent heat vapor,            [J kg-1]
% consts.xLf         =   3.34e5;      % latent heat fusion, 			[J kg-1]
% consts.xLs         =   xLv+xLf;     % latent heat subl. (2.83e6) 	[J kg-1]
% consts.Tfp         =   273.15;      % freezing point,               [K]
% consts.ro_air      =   1.275;       % density air,                  [kg m-3]
% consts.ro_ice      =   917.0; 		% density pure ice,             [kg m-3]
% consts.ro_liq      =   1000.0;      % density water,                [kg m-3]
% consts.Cp_air      =   1004;        % sp. heat cap. air,            [J kg-1 K-1]
% consts.Cp_liq      =   4180.0;      % sp. heat cap. water,          [J kg-1 K-1]
% consts.Cp_ice      =   2102.0;      % sp. heat cap. ice,            [J kg-1 K-1]
% consts.xk_liq      =   0.552;       % therm. cond. water at 0oC,    [W m-1 K-1]
% consts.xk_ice      =   2.10;        % therm. cond. ice at 0oC,      [W m-1 K-1]
% consts.Rv          =   461.0; 		% gas constant water vapor      [J kg-1 K-1]
% consts.emiss_sfc   =   0.98;        % surface emissivity,           [-]
% consts.Stef_Boltz  =   5.6696e-8; 	% stefan boltzman constant, 	[W m-2 K-4]
% consts.gravity     =   9.81; 		% gravitational accel.,         [m s-2]
% consts.xkappa      =   0.4; 		% von Karman,                   [-]
% consts.one_atmos   =   101300.0; 	% one atm,                      [Pa = J m-3 = kg m-1 s-2]
% consts.scale_ht    =   8500.0;
% consts.ihrs_day    =   24;
