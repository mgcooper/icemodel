function [Twb, flag] = SOLVEWB(Tair, rh, Ls, cp_air, P, liqflag)
   %SOLVEWB Compute wet bulb from air temperature and relative humidity
   %
   % This function computes the wet bulb temperature based on air temperature,
   % relative humidity, latent heat of sublimation, specific heat capacity of
   % air, and atmospheric pressure using an iterative method. If the solution
   % fails to converge, the wet bulb temperature is to the air temperature.
   %
   % Syntax:
   %   Twb = SOLVEWB(Tair, rh, Ls, cp_air, P)
   %
   % Inputs:
   %   Tair   - Air temperature in Kelvins.
   %   rh     - Relative humidity in percentage (0 - 100).
   %   Ls     - Latent heat of sublimation (J/kg).
   %   cp_air - Specific heat capacity of air at constant pressure (J/(kg*K)).
   %   P      - Atmospheric pressure in Pa.
   %
   % Output:
   %   Twb    - Wet bulb temperature in Kelvins. If the function fails to
   %            converge to a solution, Twb is set to the air temperature.
   %   flag   - Exit flag indicating if the solution converged (FLAG=1) or if
   %            the maximum iterations were exceeded when solving for Twb
   %            (FLAG=0).
   %
   % Example:
   %   Twb = SOLVEWB(293, 50, 2.83e6, 1005, 101325)
   %
   % See also VAPPRESS
   
   if nargin < 6
      liqflag = false;
   end
   
   % Compute atmospheric vapor pressure from relative humidity data.
   ea = VAPPRESS(Tair, 273.16, liqflag) * rh / 100.0;
   
   % Compute dew point as a fall back if the solution does not converge
   % Tdew = 273.16 + C * log(ea / A) / (B - log(ea / A));
   
   % Use this as a fall-back on Twb, but I kept it here to clarify the correct
   % Tdw formula for the A, B, C 
   % Tdw = Tf + C * log(ea / A) / (B - log(ea / A));
   % Twb = tair(metiter) - (tair(metiter) - Tdew) / 3

   % Solver options.
   maxiter = 1000;
   tol = 1.0e-3;
   old = Tair;
   
   % Default return value if solution fails to converge
   Twb = Tair;

   % Solve for the wet bulb temperature.
   for n = 1:maxiter

      f = old - Tair ...
         + Ls/cp_air * 0.622/P * (10.0 ^ (11.40 - 2353.0 / old) - ea);

      fprime = 1.0 ...
         + Ls/cp_air * 0.622/P * log(10.0) * 2353.0 ...
         * (10.0 ^ (11.40 - 2353.0 / old)) / old ^ 2;

      new = old - f/fprime;

      if (abs(new - old) < tol)
         Twb = new;
         break
      end
      old = new;
   end

   % Exit flag
   flag = n < maxiter;
end
