function Tdew = TDEWPOINT(T, rh, liqflag)
   %TDEWPOINT Compute dew point temperature from air temperature and humidity.
   %
   %  Tdew = TDEWPOINT(T, rh)
   %  Tdew = TDEWPOINT(T, rh, liqflag)
   %
   %  Dew point is the temperature at which air reaches saturation. Uses
   %  Newton iteration on the Rankine-Kirchhoff saturation vapor pressure
   %  formula (Ambaum 2020 / Romps 2021) via VAPPRESS.
   %
   %  The actual vapor pressure is ea = es(T) * rh/100, and the dew point
   %  solves es(Tdew) = ea via Newton's method.
   %
   %  Inputs:
   %     T       - Temperature [K] (nominally air temperature)
   %     rh      - Relative humidity [%] (0-100)
   %     liqflag - (optional) Reference surface phase: true = liquid water,
   %               false = ice (default: true)
   %
   %  Output:
   %     Tdew - Dew point temperature [K]
   %
   % See also: VAPPRESS, TWETBULB
   %
   %#codegen

   if nargin < 3
      liqflag = true;
   end

   persistent al bl ai bi
   if isempty(al)
      [al, bl, ai, bi] = icemodel.parameterLookup('al', 'bl', 'ai', 'bi');
   end

   % Actual vapor pressure [Pa]
   ea = VAPPRESS(T, liqflag) .* rh / 100;

   % Initial guess: invert dominant exponential (ignoring T^c term)
   if liqflag
      Tdew = bl ./ log(ea ./ al);
   else
      Tdew = bi ./ log(ea ./ ai);
   end

   % Newton iteration (typically converges in 3-5 steps)
   for iter = 1:10
      [es_guess, des_dT_guess] = VAPPRESS(Tdew, liqflag);
      Tdew = Tdew - (es_guess - ea) ./ des_dT_guess;
   end
end
