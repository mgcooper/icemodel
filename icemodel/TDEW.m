function Tdew = TDEW(T, rh, liqflag)
   %TDEW Compute dew point temperature from air temperature and humidity
   %
   % This function calculates dew point using the Magnus formula. Dew point is
   % the temperature at which the air reaches saturation. Thus to compute dew
   % point, the saturation vapor pressure equation is inverted to express
   % temperature as a function of saturation vapor pressure, and then solved by
   % substituting the actual vapor pressure for the saturation vapor pressure.
   % This can be described as finding the temperature at which the actual vapor
   % pressure would equal the saturation vapor pressure.
   %
   % Syntax:
   %   Tdew = TDEW(T, rh)
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
   % Example:
   %   Tdew = TDEW(20, 50);
   %
   % See also:
   %
   %#codegen

   if nargin < 3
      liqflag = true;
   end

   persistent Tf
   if isempty(Tf)
      Tf = 273.16;
   end

   % Direct calculation, adjusted to Celsius for ease, then back to Kelvins.
   %
   % T = T - Tf;
   % alpha = b * T / (c + T) + log(rh / 100)
   % Tdew = Tf + c * alpha / (b - alpha)
   %
   % Tdew = Tf + c * (b * T / (c + T) + log(rh / 100)) ...
   %    / (b - (b * T / (c + T) + log(rh / 100)));

   % Compute dew point temperature over water or ice.
   if liqflag == true || T >= Tf
      % Over water
      % b = 17.502
      % c = 240.97 (c - Tf = -32.19)
      Tdew = Tf + 240.97 * (17.502 * (T - Tf) / (T - 32.19) + log(rh / 100)) ...
         / (17.502 - (17.502 * (T - Tf) / (T - 32.19) + log(rh / 100)));
   else
      % Over ice
      % b = 22.452
      % c = 272.55 (c - Tf = -0.61)
      Tdew = Tf + 272.55 * (22.452 * (T - Tf) / (T - 0.61) + log(rh / 100)) ...
         / (22.452 - (22.452 * (T - Tf) / (T - 0.61) + log(rh / 100)));
   end

   % % For reference, using the vapor pressure:
   % ea = VAPPRESS(T, Tf, liqflag) * rh / 100;
   %
   % % Tdew = Tf + c * log(ea / es0) / (b - log(ea / es0));
   % if liqflag == true % || T >= Tf
   %    % Over water
   %    % a = 100 * 6.1121; % reference vapor pressure, es0 [Pa]
   %    % b = 17.502;
   %    % c = 240.97;
   %    Tdew = Tf +  (240.97 * log(ea / 611.21)) / (17.502 - log(ea / 611.21));
   % else
   %    % Over ice
   %    % A = 100 * 6.1115; % reference vapor pressure, es0 [Pa]
   %    % B = 22.452;
   %    % C = 272.55;
   %    Tdew = Tf +  (272.55 * log(ea / 611.15)) / (22.452 - log(ea / 611.15));
   % end
end
