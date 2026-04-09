function Ri = richardsonNumber(tsfc, tair, wspd, z_tair, z_wind)
   %RICHARDSON_NUMBER
   %
   % Note: The Louis parameterization is developed such that z_wind and z_tair
   % are identical (referred to as z_obs). The expression below was developed
   % for the case where z_wind and z_tair differ.
   % 
   % TODO: define a coefficient:
   %  coef = g / z_tair * (z_wind ./ wspd) .^ 2;
   % If nargin == 5, perform fully-explicit canonical calculation of Ri for
   % reference and for general use, then compute the coefficient and return it
   % from the function as the second output if nargin > 1. During model
   % initialization, retrieve it. During model runtime, to avoid extra
   % calculations on the hot path, do this:
   % if nargin == 3: 
   %  Ri = coef .* (1 - tsfc ./ tair);
   %  return
   % end
   % 
   % This way there's a canonical richardson_number function rather than one
   % function that creates the coefficients and another one that computes the
   % richardson number in a non-canonical way. The main body of the function
   % shows the canonical calculation but isn't touched when nargin == 3. Same
   % pattern could potentially be extended to the bulk Richardson
   % stability_factor function, but need to confirm.

   g = 9.81;
   Tf = 273.16;

   if nargin == 4
      z_wind = z_tair;
   end

   tair = tair(:);
   tsfc = tsfc(:);
   wspd = wspd(:);

   if min(tair) < 0
      tair = tair + Tf;
   end
   if min(tsfc) < 0
      tsfc = tsfc + Tf;
   end

   Ri = g / z_tair * (z_wind ./ wspd) .^ 2 .* (1 - tsfc ./ tair);

   % Note: if z_wind = z_tair:
   % Ri = g * z_tair / wspd ^ 2 * (1 - tsfc ./ tair)

   % took these out of stabelfn
   % if (scoef(2) / 9.4 / wspd ^ 2 * (1 - tsfc / tair)) > 0.1
   %    % this is the "very stable" regime
   % end

   % thought htis might be needed to avoide blowing up but i dont think theres
   % any chance of it except maybe in the derivative so i left it in sfctemp
   % if abs(tsfc - tair) < 1e-1 % Neutrally stable case.

end
