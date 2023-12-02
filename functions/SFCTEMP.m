function [Ts, ok] = SFCTEMP(Ta, Qsi, Qli, albedo, wspd, Pa, De, ...
      ea, cv_air, emiss, SB, Tf, chi, roL, scoef, Qc, liqflag)
   %SFCTEMP Solve the energy balance for surface temperature

   persistent tol maxiter
   if isempty(tol)
      tol = 1e-3;
      maxiter = 100;
   end

   % Gather terms in the SEB equation.
   AAA = cv_air * De;                                    % [W m-2 K-1]
   CCC = 0.622 / Pa;                                     % [Pa-1] = [m3 J-1]
   EEE = chi * (1.0 - albedo) * Qsi + emiss * Qli + Qc;  % [W m-2]
   FFF = roL * De;                                       % [W m-2]

   % Compute the constants used in the stability coefficient computations.
   B1 = scoef(2) / (Ta * wspd ^ 2);
   B2 = scoef(3) / (sqrt(Ta) * wspd);
         
   % Define the vapor pressure coefficients.
   if liqflag == true
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

   old = Ta;
   
   for iter = 1:maxiter

      % Account for an increase in turbulent fluxes under unstable conditions.
      es = A * exp(B * (old - Tf) / (C + old - Tf));
      
      if old > Ta 
         % Unstable case.
         S = 1.0 + B1 * (old - Ta) / (1.0 + B2 * sqrt(old - Ta));
         dS = B1 / (1.0 + B2 * sqrt(old - Ta)) - (B1 * B2 * (old - Ta)) ...
            / (2.0 * (1.0 + B2 * sqrt(old - Ta)) ^ 2 * sqrt(old - Ta));
         
      elseif old < Ta
         % Stable case.
         S = 1.0 / (1.0 + B1 / 2.0 * (Ta - old)) ^ 2;
         dS = 2.0 * B1 / 2.0 / (1.0 + B1 / 2.0 * (Ta - old)) ^ 3;
      else
         % Neutrally stable case.
         S = 1.0;
         dS = 0;
      end

      dfdT = -4.0 * emiss * SB * old ^ 3 ...
         + S * -AAA + AAA * (Ta - old) * dS ...
         + S * -FFF * CCC * es * B * C / (C + old - Tf) ^ 2 ...
         + FFF * CCC * (ea - es) * dS;

      f = EEE - emiss * SB * old ^ 4 + AAA * (Ta - old) * S ...
         + FFF * CCC * (ea - es) * S;
      
      Ts = old - f / dfdT;

      if abs(Ts - old) < tol
         ok = true;
         return
      elseif isnan(Ts) || iter == maxiter
         ok = false;
         Ts = Ta;
         return
      end
      old = Ts;
   end
   
   %    % Check convergence
   %    if err < tol && abs(xTs - Ts) < 10
   %       OK = true;
   %       return
   %    else
   %
   %       % Try brent's method
   %       fSEB = @(Ts) EEE - emiss * SB * Ts ^ 4 + ...
   %       AAA * STABLEFN(Ta, Ts, wspd, scoef) * (Ta - Ts) + ...
   %       FFF * CCC * STABLEFN(Ta, Ts, wspd, scoef) * ...
   %       (ea - VAPPRESS(Ts, Tf, liqflag));
   %       % + cp_liq*ppt*Tppt; % ppt in kg/m2/s
   %
   %       [Ts, ~, ok] = fsearchzero(fSEB, ...
   %          xTs, xTs-50, xTs+50, Ta, fopts.TolX);
   %
   %       if not(ok) || abs(xTs - Ts) > 10
   %          Ts = Ta;
   %       else
   %          OK = true;
   %       end
   %
   %       % If newton and brent fail, could try fzero bounded and unconstrained:
   %       % [Ts, ~, ~] = fzero(fSEB, [Ta-50 Ta+50], fopts);
   %       % [Ts, ~, ~] = fzero(fSEB, Ta, fopts);
   %
   %    end
end
