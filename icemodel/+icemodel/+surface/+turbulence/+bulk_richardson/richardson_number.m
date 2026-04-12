function [Ri, coef] = richardson_number(tsfc, tair, wspd_or_coef, z_tair, z_wind)
   %RICHARDSON_NUMBER Compute the bulk Richardson number.
   %
   %  Ri = icemodel.surface.turbulence.bulk_richardson.richardson_number(
   %     tsfc, tair, wspd, z_tair, z_wind)
   %
   %  Ri = icemodel.surface.turbulence.bulk_richardson.richardson_number(
   %     tsfc, tair, wspd, z_tair)
   %
   %  [Ri, coef] = icemodel.surface.turbulence.bulk_richardson.richardson_number(
   %     tsfc, tair, wspd, z_tair, z_wind)
   %
   %  Ri = icemodel.surface.turbulence.bulk_richardson.richardson_number(
   %     tsfc, tair, coef)
   %
   % The bulk Richardson number is:
   %
   %   Ri = g/z_tair * (z_wind/wspd)^2 * (1 - tsfc/tair)
   %      = coef * (1 - tsfc/tair)
   %
   % where coef = g/z_tair * (z_wind/wspd)^2 can be precomputed for
   % efficiency when wind speed is constant or varies slowly.
   %
   % Calling conventions:
   %
   %   Canonical/initialization mode (nargin == 4 or 5):
   %     Computes Ri from the fully-explicit formula using the observation
   %     heights z_tair [m] and z_wind [m] (defaults to z_wind = z_tair when
   %     nargin == 4). Returns coef as the second output for fast-path use.
   %
   %   Production/fast mode (nargin == 3):
   %     Accepts a precomputed coef (scalar or array) in place of wspd.
   %     Computes Ri = coef * (1 - tsfc/tair) efficiently.
   %
   % Note: The Louis parameterization assumes z_wind == z_tair (referred to
   % as z_obs). The 5-argument form supports differing measurement heights.
   %
   % Inputs are converted to Kelvin if passed in Celsius (tair or tsfc < 0).

   [tsfc, tair, wspd, coef, z_wind] = parseinputs( ...
      tsfc, tair, wspd_or_coef, z_tair, z_wind, nargin);

   persistent g
   if isempty(g)
      g = icemodel.physicalConstant('gravity');
   end

   % Fast mode: wspd_or_coef is a precomputed coefficient.
   if nargin == 3
      Ri = coef .* (1 - tsfc ./ tair);
      return
   end

   % Canonical mode (nargin == 4 or 5).
   coef = g / z_tair * (z_wind ./ wspd) .^ 2;
   Ri = coef .* (1 - tsfc ./ tair);
end

function [tsfc, tair, wspd, coef, z_wind] = parseinputs( ...
      tsfc, tair, wspd_or_coef, z_tair, z_wind, numargsin)

   persistent Tf
   if isempty(Tf)
      Tf = icemodel.physicalConstant('Tf');
   end

   tair = tair(:);
   tsfc = tsfc(:);
   if min(tair) < 0
      tair = tair + Tf;
   end
   if min(tsfc) < 0
      tsfc = tsfc + Tf;
   end

   if numargsin == 3
      % Fast mode: wspd_or_coef is a precomputed coefficient.
      coef = wspd_or_coef;
   else
      % Canonical mode (nargin == 4 or 5).
      wspd = wspd_or_coef;
      if numargsin < 5
         z_wind = z_tair;
      end
      wspd = wspd(:);
   end
end
