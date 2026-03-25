function [es, des_dT, T_dew] = buckVaporPressure(T, Tf, liqflag, rh)
   %BUCKVAPORPRESSURE Buck (1981) saturation vapor pressure — reference archive.
   %
   %  ES = buckVaporPressure(T, TF, LIQFLAG) computes saturation vapor
   %  pressure using the Rankine-Kirchhoff formula with Buck (1981)
   %  coefficients:
   %
   %     ES = A * exp(B * (T - Tf) / (C + T - Tf))   [Pa]
   %
   %  where (al, bl, cl) are over liquid and (ai, bi, ci) are over ice, and
   %  Tf = 273.16 is the triple point temperature.
   %
   %  Note: temperatures for Buck's equations are in deg C, and vapor
   %  pressures are in mb. Here they have been converted so the calculations
   %  are done with temperatures in K, and vapor pressures in Pa. (1 mb =
   %  100 Pa). Curves ew1 and ei2 are adopted for water and ice, respectively
   %  (See Buck Eq. 3,8 and Table 2).
   %
   %  [ES, DES_DT] = buckVaporPressure(T, TF, LIQFLAG) Also computes the
   %  derivative of saturation vapor pressure with respect to temperature.
   %
   %  [ES, DES_DT, T_DEW] = buckVaporPressure(T, TF, LIQFLAG, RH) Also
   %  computes the dew point temperature in Kelvins.
   %
   %  This is the archived Buck (1981) implementation. Production code uses
   %  Ambaum (2020) via VAPPRESS / icemodel.parameterLookup.
   %
   % See also: VAPPRESS, icemodel.parameterLookup

   % Buck (1981) coefficients (i = ice, l = liquid)
   al = 611.21;
   bl = 17.502;  % 18.678
   cl = 240.97;
   ai = 611.15;
   bi = 22.452;
   ci = 272.55;

   if nargin < 2
      Tf = 273.16;
   end
   if nargin < 3
      liqflag = false;
   end

   if liqflag == true
      es = saturationVaporPressure(T, Tf, al, bl, cl);

      if nargout > 1
         des_dT = saturationVaporPressureDerivative(T, Tf, es, bl, cl);
      end
   else
      es = saturationVaporPressure(T, Tf, ai, bi, ci);

      if nargout > 1
         des_dT = saturationVaporPressureDerivative(T, Tf, es, bi, ci);
      end
   end

   if nargout == 3
      if liqflag == true
         T_dew = dewPointTemperature(es * rh / 100, Tf, al, bl, cl);
      else
         T_dew = dewPointTemperature(es * rh / 100, Tf, ai, bi, ci);
      end
   end
end

function es = saturationVaporPressure(T, Tf, a, b, c)
   %SATURATIONVAPORPRESSURE Buck (1981) saturation vapor pressure
   es = a * exp(b * (T - Tf) ./ (c + T - Tf)); % [Pa]
end

function des_dT = saturationVaporPressureDerivative(T, Tf, es, b, c)
   %SATURATIONVAPORPRESSUREDERIVATIVE Derivative of es wrt temperature
   des_dT = b * c * es ./ (c + T - Tf) .^ 2; % [Pa K-1]
end

function T_dew = dewPointTemperature(e_air, Tf, a, b, c)
   %DEWPOINTTEMPERATURE Buck (1981) dew point temperature
   T_dew = Tf + c * log(e_air / a) / (b - log(e_air / a));
end

function ro_vap = saturationVaporDensity(es, T, Rv) %#ok<*DEFNU>
   %SATURATIONVAPORDENSITY Saturation vapor density from ideal gas law
   ro_vap = es ./ (Rv .* T); % [kg m-3]
end

function dro_vapdT = saturationVaporDensityDerivative(ro_vap, T, Tf, b, c)
   %SATURATIONVAPORDENSITYDERIVATIVE Derivative of saturation vapor density
   dro_vapdT = ro_vap .* (b * c ./ (c + T - Tf) .^ 2 - 1 ./ T); % [kg m-3 K-1]

   % Equivalently:
   % dro_vapdT = es ./ (Rv * T) .* (b * c ./ (c + T - Tf) .^ 2 - 1 ./ T);
   % dro_vapdT = (des_dT - es ./ T) ./ (Rv * T);
end

%% Notes
%
% To convert from e_sat to T:
%
% T = Tf + c * log(e_sat / a) / (b - log(e_sat / a));
%
% Thus to compute Tdew, substitute e_air for e_sat:
%
% T = Tf + c * log(e_air / a) / (b - log(e_air / a));
%
% Buck suggests ew1, ei2, fw3, fi3. Below are expressions for f, not used
% here, but for reference:
% fw3 = 1.0007 + 3.46e-6 * P; % note: P is in mb, need to convert.
% fi3 = 1.0003 + 4.18e-6 * P;
%
% These enhancement factors are not implemented.
