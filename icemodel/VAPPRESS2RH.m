function rh = VAPPRESS2RH(ea, Ta, liqflag)
   %VAPPRESS2RH Relative humidity from atmospheric vapor pressure.
   %
   %  rh = VAPPRESS2RH(ea, Ta, liqflag) computes relative humidity [%] from
   %  atmospheric vapor pressure ea [Pa] and air temperature Ta [K] using the
   %  Ambaum (2020) / Romps (2021) Rankine-Kirchhoff saturation vapor pressure
   %  formula.
   %
   % See also: VAPPRESS
   %
   %#codegen

   % Ambaum (2020) / Romps (2021) Rankine-Kirchhoff coefficients (i=ice, l=liq)
   persistent al bl cl ai bi ci
   if isempty(al)
      [al, bl, cl, ai, bi, ci] = icemodel.parameterLookup( ...
         'al', 'bl', 'cl', 'ai', 'bi', 'ci');
   end

   if liqflag == true
      % Over liquid: es = al * exp(bl / T) * T ^ cl
      rh = 100 * ea ./ (al * exp(bl ./ Ta) .* Ta .^ cl);
   else
      % Over ice: es = ai * exp(bi / T) * T ^ ci
      rh = 100 * ea ./ (ai * exp(bi ./ Ta) .* Ta .^ ci);
   end
end
