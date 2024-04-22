function [Tw, flag] = SOLVEWB(Ta, rh, Ls, cp, Pa, liqflag)
   %SOLVEWB Compute wet bulb from air temperature and relative humidity
   %
   % This function computes the wet bulb temperature based on air temperature,
   % relative humidity, latent heat of sublimation, specific heat capacity of
   % air, and atmospheric pressure using an iterative method. If the solution
   % fails to converge, the wet bulb temperature is to the air temperature.
   %
   % Syntax:
   %   Tw = SOLVEWB(Ta, rh, Ls, cp, Pa)
   %
   % Inputs:
   %   Ta     - Air temperature in Kelvins.
   %   rh     - Relative humidity in percentage (0 - 100).
   %   Ls     - Latent heat of sublimation (J/kg).
   %   cp     - Specific heat capacity of air at constant pressure (J/(kg*K)).
   %   Pa     - Atmospheric pressure in Pa.
   %
   % Output:
   %   Tw    - Wet bulb temperature in Kelvins. If the function fails to
   %            converge to a solution, Tw is set to the air temperature.
   %   flag   - Exit flag indicating if the solution converged (FLAG=1) or if
   %            the maximum iterations were exceeded when solving for Tw
   %            (FLAG=0).
   %
   % Example:
   %   Tw = SOLVEWB(293, 50, 2.83e6, 1005, 101325)
   %
   % See also VAPPRESS

   if nargin < 6
      liqflag = false;
   end

   % Compute atmospheric vapor pressure from relative humidity data.
   ea = VAPPRESS(Ta, 273.16, liqflag) * rh / 100.0;

   % Compute dew point as a fall back if the solution does not converge
   % Tdew = 273.16 + C * log(ea / A) / (B - log(ea / A));

   % Use this as a fall-back on Tw, but I kept it here to clarify the correct
   % Tdw formula for the A, B, C
   % Tdw = Tf + C * log(ea / A) / (B - log(ea / A));
   % Tw = tair(metiter) - (tair(metiter) - Tdew) / 3

   % Solver options.
   maxiter = 1000;
   tol = 1.0e-3;
   old = Ta;

   % Default return value if solution fails to converge
   Tw = Ta;

   % Solve for the wet bulb temperature.
   for n = 1:maxiter

      f = old - Ta ...
         + Ls/cp * 0.622/Pa * (10.0 ^ (11.40 - 2353.0 / old) - ea);

      fprime = 1.0 ...
         + Ls/cp * 0.622/Pa * log(10.0) * 2353.0 ...
         * (10.0 ^ (11.40 - 2353.0 / old)) / old ^ 2;

      new = old - f/fprime;

      if (abs(new - old) < tol)
         Tw = new;
         break
      end
      old = new;
   end

   % Exit flag
   flag = n < maxiter;
end
