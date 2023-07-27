function [Tsfc,OK] = SFCTEMP(Tair,Qsi,Qli,ea,albedo,De,Pa,wspd,cv_air,...
                     emiss,SB,Tf,Qc,xTsfc,chi,roL,scoef,fopts,liqflag)
%SFCTEMP solve the energy balance for surface temperature

% debug =  false;
% dflag =  true;

% gather terms in the SEB equation. Note that if icond_flag == 0, Qc = 0
   AAA   =  cv_air * De;                           % [W m-2 K-1]
   CCC   =  0.622 / Pa;                            % [Pa-1] = [m3 J-1]
   EEE   =  chi*(1.0-albedo)*Qsi + emiss*Qli + Qc; % [W m-2]
   FFF   =  roL * De;                              % [W m-2]

% % Compute the constants used in the stability coefficient computations
% % (following are needed for SOLVE but not for fzero)
%   a1   =  5.3*9.4;                               % [-]
%   z1   =  z_obs/z_0;                             % [-]
%   C1   =  a1 * (kapp/(log(z1)))^2 * sqrt(z1);    % [-]
%   C2   =  grav * z_obs/(Tair*wspd^2);            % [K-1]
%   B1   =  9.4 * C2;                              % [K-1]
%   B2   =  C1 * sqrt(C2);                         % [K-1]

  % these replace the block above
   B1    =  scoef(2)/(Tair*wspd^2);
   B2    =  scoef(3)/(sqrt(Tair)*wspd);


   % SOLVE
   if liqflag == true
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

   old = Tair; errT = 2e-3; iter = 0; OK = false;   
   
   while errT>1e-3 && iter<200

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
         df1   =  -4.0*emiss*SB*old^3;
         df2   =  S*-AAA+AAA*(Tair-old)*dS;
         df3   =  S*-FFF*CCC*es0*B*C/((C+old-Tf)^2)+FFF*CCC*(ea-es0)*dS;
         df4   =  -0.0;
      elseif (old<Tair)                      % Stable case.
         B8    =  B1/2.0;
         S     =  1.0/((1.0+B8*(Tair-old))^2);
         dS    =  2.0*B8/((1.0+B8*(Tair-old))^3);
         df1   =  -4.0*emiss*SB*old^3;
         df2   =  S*-AAA+AAA*(Tair-old)*dS;
         df3   =  S*-FFF*CCC*es0*B*C/((C+old-Tf)^2)+FFF*CCC*(ea-es0)*dS;
         df4   =  -0.0;
      else                                   % Neutrally stable case.
         S     =  1.0;
         df1   =  -4.0*emiss*SB*old^3;
         df2   =  -AAA;
         df3   =  -FFF*CCC*es0*B*C/((C+old-Tf)^2);
         df4   =  -0.0;
      end

      % f1-5 can be activated for clarity / debugging
      % f1  =   EEE - emiss*SB*old^4;    % net radiation + conduction (Qc is in EEE)
      % f2  =   AAA*(Tair-old)*S;   % sensible flux
      % f3  =   FFF*CCC*(ea-es0)*S; % latent flux
      % f4  =   0.0;                % for the case where Qc is taken out from EEE

      f     =  EEE-emiss*SB*old^4+AAA*(Tair-old)*S+FFF*CCC*(ea-es0)*S; % + f4
      Tsfc  =  old-f/(df1 + df2 + df3 + df4);

%       if debug == true
%          if dflag == true
%             figure; 
%             s1 = subplot(1,2,1); plot(old,f,'ko'); hold on;
%             xylabel('T','F(T)');
%             legend('f(Tair)','AutoUpdate','off'); 
%             s2 = subplot(1,2,2); plot(old,(df1 + df2 + df3 + df4),'ko');
%             xylabel('T','dF/dT');
%             hold on; dflag = false;
%          else
%             plot(s1,old,f,'o'); plot(s2,old,(df1 + df2 + df3 + df4),'o');
%             pause;
%          end
%       end
      
      % prep for next iteration
      errT  =  abs(Tsfc-old);
      old   =  Tsfc;
      iter  =  iter+1;
   end
   
   % if converges pass it to SFCTEMP
   if errT<1e-3 && abs(xTsfc-Tsfc) < 10
      OK    =  true;  return
   else
      
      % cycle through various solver opts until we get a good solution
      dif      =  20;
      tryflag  =  0;
      
      while dif > 10 && tryflag < 3
         
         tryflag = tryflag + 1;
         
         % try fsolve
         try %#ok<TRYNC>
            Tsfc  = SFCTMPFSOLVE(EEE,AAA,FFF,CCC,Tair,ea,wspd,emiss,SB,Tf, ...
                    xTsfc,scoef,fopts,liqflag,tryflag);
%          catch ME
%             if strcmp(ME.identifier,'MATLAB:fzero:ValuesAtEndPtsSameSign')
%                msg      =  'Tsfc failed using endpoints, using midpoint instead';
%                causeE   =  MException('icemodel:SFCTEMP:rootfinding',msg);
%                ME       =  addCause(ME,causeE); % let it go
%             end

         end
         dif = abs(xTsfc-Tsfc);
      end
      
      if isnan(Tsfc) || dif>10
         Tsfc  =  Tair; % if it doesn't converge use tair
      else
         OK    = true;
      end
   end
   
end
   

% fzero
function Tsfc = SFCTMPFSOLVE(EEE,AAA,FFF,CCC,Tair,ea,wspd,emiss,SB,Tf,  ...
                  xTsfc,scoef,fopts,liqflag,tryflag)
                  
   Sfnc    =   @STABLEFN;
   Vfnc    =   @VAPOR;
   fSEB    =   @(Tsfc) EEE - emiss*SB.*Tsfc.^4 +                       ...
               AAA.*Sfnc(Tair,Tsfc,wspd,scoef).*(Tair-Tsfc) +  ...
               FFF.*CCC.*Sfnc(Tair,Tsfc,wspd,scoef) .* ...
               (ea-Vfnc(Tsfc,Tf,liqflag));
               % + cp_liq*ppt*Tppt; % ppt in kg/m2/s

% solve the equation, using Tair as an intial guess for Tsfc
      
if tryflag == 1         % first try bounding by +/- 5 deg past Tsfc
   [Tsfc,~,~]  = fzero(fSEB,[xTsfc-5 xTsfc+5],fopts);
elseif tryflag == 2     % next try Tair in case of bad past Tsfc
   [Tsfc,~,~]  = fzero(fSEB,[Tair-5 Tair+5],fopts);
elseif tryflag == 3     % next try Tair in case of bad past Tsfc
   [Tsfc,~,~]  = fzero(fSEB,Tair,fopts); 
end 

end