function [Tw, ok] = TWETBULB(Ta, rh, Pa, liqflag)
   %TWETBULB Compute wet-bulb temperature from air temperature and humidity.
   %
   %  [Tw, flag] = TWETBULB(Ta, rh, Pa)
   %  [Tw, flag] = TWETBULB(Ta, rh, Pa, liqflag)
   %
   %  Solves the psychrometric equation for wet-bulb temperature using
   %  Newton iteration with analytic derivatives from VAPPRESS:
   %
   %     f(Tw) = Tw - Ta + (Ls/cp_air) * (epsilon/Pa) * (es(Tw) - ea) = 0
   %
   %  If Newton iteration fails to converge, falls back to the Stull (2011)
   %  heuristic: Tw = Ta - (Ta - Tdew) / 3.
   %
   %  Inputs:
   %     Ta      - Air temperature [K]
   %     rh      - Relative humidity [%] (0-100)
   %     Ls      - Latent heat of sublimation [J kg-1]
   %     cp_air  - Specific heat capacity of air [J kg-1 K-1]
   %     Pa      - Atmospheric pressure [Pa]
   %     liqflag - (optional) Phase flag for saturation: true = liquid,
   %               false = ice (default: false)
   %
   %  Outputs:
   %     Tw   - Wet-bulb temperature [K]
   %     flag - true if Newton converged, false if heuristic was used
   %
   %  Reference:
   %     Stull (2011), "Wet-Bulb Temperature from Relative Humidity and Air
   %     Temperature." Journal of Applied Meteorology and Climatology, 50(11).
   %
   % See also: VAPPRESS, TDEWPOINT
   %
   %#codegen

   if nargin < 4
      liqflag = false;
   end

   persistent Tf Ls cp_air epsilon
   if isempty(Tf)
      [Tf, Ls, cp_air, epsilon] = icemodel.physicalConstant(...
         'Tf', 'Ls', 'cp_air', 'epsilon');
   end

   % Actual vapor pressure from relative humidity [Pa]
   ea = VAPPRESS(Ta, liqflag) * rh / 100.0;

   % Solver options
   maxiter = 1000;
   tol = 1.0e-3;
   old = Ta;
   ok = false;

   % Solve for wet-bulb temperature using Newton iteration
   % Psychrometric equation: Tw - Ta + (Ls/cp_air)*(epsilon/Pa)*(es(Tw)-ea)=0
   coeff = Ls / cp_air * epsilon / Pa;
   for n = 1:maxiter

      [es_wb, des_dT_wb] = VAPPRESS(old, liqflag);

      f = old - Ta + coeff * (es_wb - ea);
      fprime = 1.0 + coeff * des_dT_wb;

      new = old - f / fprime;

      if abs(new - old) < tol
         ok = true;
         Tw = new;
         return
      end
      old = new;
   end

   % Heuristic fallback (Stull 2011): Tw = Ta - (Ta - Tdew) / 3
   Tdew = TDEWPOINT(Ta, rh, liqflag);
   Tw = Ta - (Ta - Tdew) / 3;
end
