%--------------------------------------------------------------------------
%   SOLVE THE ENERGY BALANCE FOR SURFACE TEMPERATURE
%--------------------------------------------------------------------------

function [Tsfc,OK] = SFCTEMP(Tair,Qsi,Qli,ea,albedo,De_h,Pa,wspd,cv_air,...
                     emissSB,Tf,Qc,xTsfc,chi,roL,scoef,fopts,liqflag)
%--------------------------------------------------------------------------

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
% hand made
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   if xTsfc>=Tf && liqflag == true
      % Over water.
      A     =  611.21;
      B     =  17.502;
      C     =  240.97;
   else
      % Over ice.
      A     =  611.15;
      B     =  22.452;
      C     =  272.55;
   end

   old   = Tair;
   OK    = false;

   for i=1:2000

      % commented vars can be activated for clarity / debugging
      % This accounts for an increase in turbulent fluxes under unstable conditions.
      % other1    =  AAA * (Tair-old);
      es0         =  A * exp((B*(old-Tf))/(C+(old-Tf)));
      % other2    =  FFF*CCC*(ea-es0);
      % dother1   =  -AAA;
      % dother2   =  -FFF*CCC*es0*B*C/((C+old-Tf)^2);

      if (old>Tair)                          % Unstable case.
         B3    =  1.0+B2*sqrt(old-Tair);
         S     =  1.0+B1*(old-Tair)/B3;
         dS    =  B1/B3-(B1*B2*(old-Tair))/(2.0*B3*B3*sqrt(old-Tair));
         df1   =  -4.0*DDD*old^3;
         df2   =  S*-AAA+AAA*(Tair-old)*dS;
         df3   =  S*-FFF*CCC*es0*B*C/((C+old-Tf)^2)+FFF*CCC*(ea-es0)*dS;
         df4   =  -0.0;
      elseif (old<Tair)                      % Stable case.
         B8    =  B1/2.0;
         S     =  1.0/((1.0+B8*(Tair-old))^2);
         dS    =  2.0*B8/((1.0+B8*(Tair-old))^3);
         df1   =  -4.0*DDD*old^3;
         df2   =  S*-AAA+AAA*(Tair-old)*dS;
         df3   =  S*-FFF*CCC*es0*B*C/((C+old-Tf)^2)+FFF*CCC*(ea-es0)*dS;
         df4   =  -0.0;
      else                                   % Neutrally stable case.
         S     =  1.0;
         df1   =  -4.0*DDD*old^3;
         df2   =  -AAA;
         df3   =  -FFF*CCC*es0*B*C/((C+old-Tf)^2);
         df4   =  -0.0;
      end

      % f1-5 can be activated for clarity / debugging
      % f1  =   EEE - DDD*old^4;    % net radiation + conduction (Qc is in EEE)
      % f2  =   AAA*(Tair-old)*S;   % sensible flux
      % f3  =   FFF*CCC*(ea-es0)*S; % latent flux
      % f4  =   0.0;                % for the case where Qc is taken out from EEE

      f     = EEE-DDD*old^4+AAA*(Tair-old)*S+FFF*CCC*(ea-es0)*S; % + f4
      Tsfc  = old-f/(df1 + df2 + df3 + df4);

      % if xnew converges pass it to SFCTEMP
      if (abs(Tsfc - old)<1e-4)
         OK = true;
         return
      end
      old   = Tsfc;
   end
   
   % if it doesn't converge it goes to tair
   if OK == false || isnan(Tsfc)
      Tsfc  = Tair;
      OK    = false;
   end
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
