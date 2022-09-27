function [Tsfc,OK] = SFCTEMP_TEST(Tair,Qsi,Qli,ea,albedo,De,Pa,wspd,    ...
                     cv_air,emissSB,Tf,Qc,xTsfc,chi,roL,scoef,liqflag)

% % this tests defining fewer new variables up top, but doing so seems to
% slow down the code just slightly
                  
   EEE   =  chi*(1.0-albedo) * Qsi + Qli + Qc;     % [W m-2]
   B1    =  scoef(2)/(Tair*wspd^2);
   B2    =  scoef(3)/(sqrt(Tair)*wspd);
   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SOLVE
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

   old = Tair; errT = 2e-4; iter = 0; OK = false;

   while errT>1e-4 && iter<2000
      % This accounts for an increase in turbulent fluxes under unstable conditions.
      es0         =  A * exp((B*(old-Tf))/(C+(old-Tf)));

      if (old>Tair)                          % Unstable case.
         B3    =  1.0+B2*sqrt(old-Tair);
         S     =  1.0+B1*(old-Tair)/B3;
         dS    =  B1/B3-(B1*B2*(old-Tair))/(2.0*B3*B3*sqrt(old-Tair));
         df1   =  -4.0*emissSB*old^3;
         df2   =  cv_air*De*((Tair-old)*dS-S);
         df3   =  S*-roL*De*0.622/Pa*es0*B*C/((C+old-Tf)^2)+roL*De*0.622/Pa*(ea-es0)*dS;
         df4   =  -0.0;
      elseif (old<Tair)                      % Stable case.
         B8    =  B1/2.0;
         S     =  1.0/((1.0+B8*(Tair-old))^2);
         dS    =  2.0*B8/((1.0+B8*(Tair-old))^3);
         df1   =  -4.0*emissSB*old^3;
         df2   =  cv_air*De*((Tair-old)*dS-S);
         df3   =  S*-roL*De*0.622/Pa*es0*B*C/((C+old-Tf)^2)+roL*De*0.622/Pa*(ea-es0)*dS;
         df4   =  -0.0;
      else                                   % Neutrally stable case.
         S     =  1.0;
         df1   =  -4.0*emissSB*old^3;
         df2   =  -cv_air*De;
         df3   =  -roL*De*0.622/Pa*es0*B*C/((C+old-Tf)^2);
         df4   =  -0.0;
      end
      f     =  EEE-emissSB*old^4+cv_air*De*(Tair-old)*S+roL*De*0.622/Pa*(ea-es0)*S; % + f4
      Tsfc  =  old-f/(df1+df2+df3+df4);

      % prep for next iteration
      errT  =  abs(Tsfc-old);
      old   =  Tsfc;
      iter  =  iter+1;
   end
   
   % if xnew converges pass it to SFCTEMP
   if errT<1e-4
      OK    =  true;  return
   else
      Tsfc  =  Tair; % if it doesn't converge use tair
   end