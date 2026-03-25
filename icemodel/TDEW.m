function Tdew = TDEW(T, rh, liqflag)
   %TDEW Compute dew point temperature from air temperature and humidity.
   %
   % Dew point is the temperature at which the air reaches saturation. For
   % the Ambaum (2020) Rankine-Kirchhoff formula es = a * exp(b/T) * T^c,
   % the dew point is transcendental and solved by Newton iteration.
   %
   % Syntax:
   %   Tdew = TDEW(T, rh)
   %   Tdew = TDEW(T, rh, liqflag)
   %
   % Inputs:
   %   T       - Temperature in Kelvins (nominally air temperature).
   %   rh      - Relative humidity in percentage (0 - 100).
   %   liqflag - (optional) Flag indicating whether the reference surface is
   %              liquid water (liqflag = true) or ice (liqflag = false).
   %
   % Output:
   %   Tdew - Dew point temperature in Kelvins.
   %
   % See also: VAPPRESS
   %
   %#codegen

   if nargin < 3
      liqflag = true;
   end

   % Ambaum (2020) Rankine-Kirchhoff coefficients (i = ice, l = liquid)
   persistent al bl cl ai bi ci
   if isempty(al)
      [al, bl, cl, ai, bi, ci] = icemodel.parameterLookup( ...
         'al', 'bl', 'cl', 'ai', 'bi', 'ci');
   end

   persistent Tf
   if isempty(Tf)
      Tf = icemodel.physicalConstant('Tf');
   end

   % Select phase coefficients
   if liqflag == true || T >= Tf
      a = al; b = bl; c = cl;
   else
      a = ai; b = bi; c = ci;
   end

   % Actual vapor pressure [Pa]
   e_air = a * exp(b ./ T) .* T .^ c .* (rh / 100);

   % Newton iteration to solve: e_air = a * exp(b / Tdew) * Tdew ^ c
   % Initial guess: invert the dominant exponential term (ignoring T^c)
   Tdew = b ./ log(e_air ./ a);

   for iter = 1:10
      es_guess = a * exp(b ./ Tdew) .* Tdew .^ c;
      des_dT = es_guess ./ Tdew .* (c - b ./ Tdew);
      Tdew = Tdew - (es_guess - e_air) ./ des_dT;
   end
end
