function [Tsfc,OK] = SOLVE(Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf,xTsfc,liqflag)
   %SOLVE NEWTON-RHAPSON SOLVER for TSFC

   if xTsfc >= Tf && liqflag == true
      % Over water.
      A = 611.21;
      B = 17.502;
      C = 240.97;
   else
      % Over ice.
      A = 611.15;
      B = 22.452;
      C = 272.55;
   end

   old = Tair;
   OK = false;

   for i = 1:2000

      % commented vars can be activated for clarity / debugging
      % This accounts for an increase in turbulent fluxes under unstable
      % conditions. 
      % other1 = AAA * (Tair-old);
      es0 = A * exp((B*(old-Tf))/(C+(old-Tf)));
      % other2 = FFF*CCC*(ea-es0);
      % dother1 = -AAA;
      % dother2 = -FFF*CCC*es0*B*C/((C+old-Tf)^2);

      if (old>Tair)                          % Unstable case.
         B3 = 1.0+B2*sqrt(old-Tair);
         S = 1.0+B1*(old-Tair)/B3;
         dS = B1/B3-(B1*B2*(old-Tair))/(2.0*B3*B3*sqrt(old-Tair));
         df1 = -4.0*DDD*old^3;
         df2 = S*-AAA+AAA*(Tair-old)*dS;
         df3 = S*-FFF*CCC*es0*B*C/((C+old-Tf)^2)+FFF*CCC*(ea-es0)*dS;
         df4 = -0.0;
      elseif (old<Tair)                      % Stable case.
         B8 = B1/2.0;
         S = 1.0/((1.0+B8*(Tair-old))^2);
         dS = 2.0*B8/((1.0+B8*(Tair-old))^3);
         df1 = -4.0*DDD*old^3;
         df2 = S*-AAA+AAA*(Tair-old)*dS;
         df3 = S*-FFF*CCC*es0*B*C/((C+old-Tf)^2)+FFF*CCC*(ea-es0)*dS;
         df4 = -0.0;
      else                                   % Neutrally stable case.
         S = 1.0;
         df1 = -4.0*DDD*old^3;
         df2 = -AAA;
         df3 = -FFF*CCC*es0*B*C/((C+old-Tf)^2);
         df4 = -0.0;
      end

      % f1-5 can be activated for clarity / debugging
      % f1 = EEE - DDD*old^4;    % net radiation + conduction (Qc is in EEE)
      % f2 = AAA*(Tair-old)*S;   % sensible flux
      % f3 = FFF*CCC*(ea-es0)*S; % latent flux
      % f4 = 0.0;                % for the case where Qc is taken out from EEE

      f = EEE-DDD*old^4+AAA*(Tair-old)*S+FFF*CCC*(ea-es0)*S; % + f4
      Tsfc = old-f/(df1 + df2 + df3 + df4);

      % if xnew converges pass it to SFCTEMP
      if (abs(Tsfc - old)<1e-4)
         OK = true;
         return
      end
      old = Tsfc;  % Tsfc = Tair; % moved this out of the loop
   end
   Tsfc = Tair; % if it doesn't converge it goes to tair
end
