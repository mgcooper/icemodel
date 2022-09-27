%--------------------------------------------------------------------------
%   SOLVE THE ENERGY BALANCE FOR SURFACE TEMPERATURE
%--------------------------------------------------------------------------

function [Tsfc,OK] = SFCTEMP(Tair,Qsi,Qli,ea,albedo,De_h,Pa,wspd,cv_air,...
                     emissSB,Tf,Qc,xTsfc,chi,roL,scoef,fopts,liqflag)
%--------------------------------------------------------------------------

% define these in PHYSICAL_CONSTS and/or put in an array S since they won't
% be accessed often and obscurity is probably acceptable
% rocp_air        = ro_air*cp_air;
% roLs_air        = ro_air*Ls;
% sigma_epsilon   = emiss_sfc*stef_boltz;
% p622_Pa         = 0.622/Pa;
% S_a1            = 5.3*9.4;
% S_z1            = z_obs/z_0;
% S_C1            = S_a1*(kapp/(log(S_z1)))^2*sqrt(S_z1);    % [-]

% gather terms in the SEB equation. Note that if icond_flag == 0, Qc = 0
    AAA     =   cv_air * De_h;                        % [W m-2 K-1]
    CCC     =   0.622 / Pa;                           % [Pa-1] = [m3 J-1]
    DDD     =   emissSB;                              % [W m-2 K-4]
    EEE     =   chi*(1.0-albedo) * Qsi + Qli + Qc;    % [W m-2]
    FFF     =   roL * De_h;                           % [W m-2]

% % Compute the constants used in the stability coefficient computations
% % (following are needed for SOLVE but not for fzero)
%   a1      =   5.3*9.4;                                % [-]
%   z1      =   z_obs/z_0;                              % [-]
%   C1      =   a1 * (kapp/(log(z1)))^2 * sqrt(z1);     % [-]
%   C2      =   grav * z_obs/(Tair*wspd^2);             % [K-1]
%   B1      =   9.4 * C2;                               % [K-1]
%   B2      =   C1 * sqrt(C2);                          % [K-1]

  % these should replace the block above
   B1       =  scoef(2)/(Tair*wspd^2);
   B2       =  scoef(3)/(sqrt(Tair)*wspd);
   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% fzero
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     Sfnc    =   @STABLEFN;
%     Vfnc    =   @VAPOR;
%     fSEB    =   @(Tsfc) EEE - DDD.*Tsfc.^4 +                            ...
%                 AAA.*Sfnc(Tair,Tsfc,wspd,scoef).*(Tair-Tsfc) +  ...
%                 FFF.*CCC.*Sfnc(Tair,Tsfc,wspd,scoef) .* ...
%                 (ea-Vfnc(Tsfc,Tf,liqflag));
%                 % + cp_liq*ppt*Tppt; % ppt in kg/m2/s
% 
% % solve the equation, using Tair as an intial guess for Tsfc
%   % [Tsfc,~,~]      =   fzero(fSEB,Tair,fopts);
%     [Tsfc,~,~]      =   fzero(fSEB,xTsfc,fopts);
%   % [Tsfc,~,~]      =   fzero(fSEB,[Tair-5 Tair+5],fopts);
%   
% %   SEB     =   fSEB(250:1:300);
% %   SEB     =   derivative(fSEB(250:1:300));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% hand made
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  [ Tsfc,                                                               ...
    OK   ]  =   SOLVE(Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf,xTsfc,liqflag); % [K]

%%%%%%%%%%%%%%%%%%%%%%%
% use this for the standard method above
   if OK == false || isnan(Tsfc)
      Tsfc  = Tair;
      OK    = false;
%    else
%       OK    = true;
   end
%%%%%%%%%%%%%%%%%%%%%%%
