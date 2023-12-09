function [es, des_dT, T_dew] = VAPPRESS(T, Tf, liqflag, rh)
   %VAPPRESS Compute saturation vapor pressure over water or ice.
   %
   %  E_SAT = VAPPRESS(T, TF, LIQFLAG) Computes saturation vapor pressure
   %  using the Magnus formula:
   %
   %     E_SAT = ES0 * exp(B * (T - Tf) / (C + T - Tf))
   %
   %  Where A, B, C are coefficients from Buck (1981), and Tf = 273.16 is
   %  the triple point temperature. Note that temperatures for Buck's
   %  equations are in deg C, and vapor pressures are in mb. Here they have
   %  been converted so the calculations are done with temperatures in K,
   %  and vapor pressures in Pa. (1 mb = 100 Pa). Curves ew1 and ei2 are
   %  adopted for water and ice, respectively (See Buck Eq. 3,8 and Table 2).
   %
   %  [E_SAT, E_AIR] = VAPPRESS(T, TF, LIQFLAG, RH) Also computes vapor
   %  pressure from relative humidity: E_AIR = E_SAT * RH / 100.
   %
   %  [E_SAT, E_AIR, T_DEW] = VAPPRESS(T, TF, LIQFLAG, RH) Also computes the
   %  dew point temperature in Kelvins.
   %
   %  Note that A is the reference vapor pressure, es0 [Pa]
   %
   % Units: [Pa = J m-3 = N m-2 = kg m-1 s-2]
   %
   % See also: VAPORHEAT

   % Not currently used, but could default to 100%
   % if nargin < 4
   %    rh = 100 * ones(size(T));
   % end

   % Define coefficients over water and ice
   persistent aw bw cw ai bi ci
   if isempty(aw)
      aw = 611.21;
      bw = 17.502;
      cw = 240.97;
      ai = 611.15;
      bi = 22.452;
      ci = 272.55;
   end

   if nargin < 2
      Tf = 273.16;
   end
   if nargin < 3
      liqflag = false;
   end

   if liqflag == true % || T >= Tf
      es = saturationVaporPressure(T, Tf, aw, bw, cw);

      % If requested, compute the derivative of es wrt temperature
      if nargout > 1
         des_dT = saturationVaporPressureDerivative(T, Tf, es, bw, cw);
      end
   else
      es = saturationVaporPressure(T, Tf, ai, bi, ci);

      % If requested, compute the derivative of es wrt temperature
      if nargout > 1
         des_dT = saturationVaporPressureDerivative(T, Tf, es, bi, ci);
      end
   end

   % If requested, compute dew point from relative humidity
   if nargout == 3
      if liqflag == true % || T >= Tf
         T_dew = dewPointTemperature(es * rh / 100, Tf, aw, bw, cw);
      else
         T_dew = dewPointTemperature(es * rh / 100, Tf, ai, bi, ci);
      end
   end
end

function es = saturationVaporPressure(T, Tf, a, b, c)
   %SATURATIONVAPORPRESSURE Compute saturation vapor pressure over water or ice
   es = a * exp(b * (T - Tf) ./ (c + T - Tf)); % [Pa]
end

function des_dT = saturationVaporPressureDerivative(T, Tf, es, b, c)
   %SATURATIONVAPORPRESSUREDERIVATIVE Derivative of es wrt temperature
   des_dT = b * c * es ./ (c + T - Tf) .^ 2; % [Pa K-1]
end

function T_dew = dewPointTemperature(e_air, Tf, a, b, c)
   %DEWPOINTTEMPERATURE Compute dew point temperature over water or ice
   T_dew = Tf + c * log(e_air / a) / (b - log(e_air / a));
end

function ro_vap = saturationVaporDensity(es, T, Rv) %#ok<*DEFNU>
   %SATURATIONVAPORDENSITY Compute saturation vapor density
   ro_vap = es ./ (Rv .* T); % [kg m-3]
end

function dro_vapdT = saturationVaporDensityDerivative(ro_vap, T, Tf, b, c)
   %SATURATIONVAPORDENSITYDERIVATIVE Derivative of saturation vapor density
   dro_vapdT = ro_vap .* (b * c ./ (c + T - Tf) .^ 2 - 1 ./ T); % [kg m-3 K-1]

   % % Equivalently:
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
