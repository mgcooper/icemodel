
% 

T15m     = ice1.Time(1):minutes(15):ice1.Time(end)+minutes(45);
T1hr     = ice1.Time;


keff     = interp1(T15m,ice2.k_eff(1,:),ice1.Time);
cp       = interp1(T15m,ice2.cp_sno(1,:),ice1.Time);
rho      = interp1(T15m,ice2.ro_sno(1,:),ice1.Time);

dt       = 15*60;                               % sec
deltaz   = sqrt(keff./(rho.*cp).*dt);

figure; plot(T1hr,deltaz)




% the 'radiation to conduction' parameter:
load('physicalconstants.mat','SB_const')
L        = 0.06;                             % layer thickness, m
T        = 273;                              % temperature, K
keff     = 2;                                % eff. thermal k, W/m/K
K        = 5;                                % abs. coeff. 1/m
tau      = K*L;                              % opacity [-]
N        = keff*K/(4*SB_const*T^3);          

% the Rosseland conductivity:
eta      = 1.3;                              % refractive index [-]
betaR    = 5;                                % Rosseland ext. coeff. eq 11.188
kR       = (16*eta^2*SB_const*T^3)/(3*betaR)

% kR








