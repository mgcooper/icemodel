function [AAA,CCC,DDD,EEE,FFF,B1,B2] = PREP_SEBSOLVE(ro_air,Cp_air,De_h,...
                                            Pa,emiss_sfc,Stef_Boltz,chi,...    
                                            albedo,Qsi,Qli,Qc,xLs,z_obs,...
                                            z_0,xkappa,gravity,Tair,wspd)
    
    AAA     =   ro_air * Cp_air * De_h;                         % [W m-2 K-1]
    CCC     =   0.622 / Pa;                                     % [Pa-1] = [m3 J-1]
    DDD     =   emiss_sfc * Stef_Boltz;                         % [W m-2 K-4]
    EEE     =   chi*(1.0-albedo) * Qsi + Qli + Qc;              % [W m-2]
    FFF     =   ro_air * xLs * De_h;                            % [W m-2]

    % Compute the constants used in the stability coefficient computations
    a1      =   5.3*9.4;                                        % [-]
    z1      =   z_obs/z_0;                                      % [-]
    C1      =   a1 * (xkappa/(log(z1)))^2 * sqrt(z1);           % [-]
    C2      =   gravity * z_obs/(Tair*wspd^2);                  % [K-1]
    B1      =   9.4 * C2;                                       % [K-1]
    B2      =   C1 * sqrt(C2);