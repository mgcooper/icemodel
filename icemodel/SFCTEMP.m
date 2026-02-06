function [Ts, ok] = SFCTEMP(Ta, Qsi, Qli, albedo, wspd, Pa, De, ea, cv_air, ...
      emiss, SB, Tf, chi, roL, scoef, liqflag, varargin)
   %SFCTEMP Solve the energy balance for surface temperature
   %
   % This function uses a traditional Newton-Rhapson iteration to find Tsfc
   %
   %#codegen

   persistent tol maxiter
   if isempty(tol)
      tol = 1e-3;
      maxiter = 100;
   end
   
   Ts = nan;
   ok = false;

   switch numel(varargin)
      case 1
         Qc = varargin{1};
         a1 = 0.0;
      case 3
         % Qc = a1 * T(1) - a1 * Ts;
         [k_eff, T, dz] = deal(varargin{:});
         a1 = k_eff(1) / (dz(1) / 2);
         Qc = 0.0;
      case 5
         % use this to test w/wo the Qc derivative
         [Qc, k_eff, T, dz, flag] = deal(varargin{:});
         if flag
            a1 = k_eff(1) / (dz(1) / 2);
            Qc = 0.0;
         else
            a1 = 0.0;
         end
      otherwise
         error('unrecognized number of inputs')
   end

   % Gather terms in the SEB equation.
   AAA = cv_air * De;   % [W m-2 K-1]
   CCC = 0.622 / Pa;    % [Pa-1] = [m3 J-1]
   EEE = chi * (1.0 - albedo) * Qsi + emiss * Qli + Qc + a1 * T(1); % [W m-2]
   FFF = roL * De;      % [W m-2]

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

      % Update surface saturation vapor pressure
      es = A * exp(B * (old - Tf) / (C + old - Tf));

      % Account for an increase in turbulent fluxes under unstable conditions.
      if old < Ta
         % Stable case.
         S = 1.0 / (1.0 + B1 / 2.0 * (Ta - old)) ^ 2;
         dS = 2.0 * B1 / 2.0 / (1.0 + B1 / 2.0 * (Ta - old)) ^ 3;

      elseif old > Ta
         % Unstable case.
         S = 1.0 + B1 * (old - Ta) / (1.0 + B2 * sqrt(old - Ta));
         dS = B1 / (1.0 + B2 * sqrt(old - Ta)) - (B1 * B2 * (old - Ta)) ...
            / (2.0 * (1.0 + B2 * sqrt(old - Ta)) ^ 2 * sqrt(old - Ta));
      else
         % Neutrally stable case.
         S = 1.0;
         dS = 0;
      end

      f = EEE - emiss * SB * old ^ 4 + AAA * (Ta - old) * S ...
         + FFF * CCC * (ea - es) * S - a1 * old;

      dfdT = -4.0 * emiss * SB * old ^ 3 ...
         + S * -AAA + AAA * (Ta - old) * dS ...
         + S * -FFF * CCC * es * B * C / (C + old - Tf) ^ 2 ...
         + FFF * CCC * (ea - es) * dS ...
         - a1;

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
end
