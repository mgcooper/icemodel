function [es, des_dT, T_dew] = VAPPRESS(T, Tf, liqflag, rh)
   %VAPPRESS Compute saturation vapor pressure over liquid or ice.
   %
   %  ES = VAPPRESS(T, TF, LIQFLAG) computes saturation vapor pressure
   %  using the Ambaum (2020) Rankine-Kirchhoff formula:
   %
   %     ES = A * exp(B / T) * T ^ C   [Pa]
   %
   %  where (al, bl, cl) are coefficients over liquid and (ai, bi, ci) are
   %  over ice, obtained from icemodel.parameterLookup.
   %
   %  [ES, DES_DT] = VAPPRESS(T, TF, LIQFLAG) Also computes the derivative
   %  of saturation vapor pressure with respect to temperature.
   %
   %  [ES, DES_DT, T_DEW] = VAPPRESS(T, TF, LIQFLAG, RH) Also computes the
   %  dew point temperature in Kelvins.
   %
   %  Note: the Tf argument is retained for signature compatibility but is
   %  not used by the Ambaum formula (which operates in absolute temperature).
   %
   % Units: [Pa = J m-3 = N m-2 = kg m-1 s-2]
   %
   % See also: VAPORHEAT, icemodel.kernels.buckVaporPressure
   %
   %#codegen

   % Ambaum (2020) Rankine-Kirchhoff coefficients (i = ice, l = liquid)
   persistent al bl cl ai bi ci
   if isempty(al)
      [al, bl, cl, ai, bi, ci] = icemodel.parameterLookup( ...
         'al', 'bl', 'cl', 'ai', 'bi', 'ci');
   end

   if nargin < 2
      Tf = 273.16; %#ok<NASGU>
   end
   if nargin < 3
      liqflag = false;
   end

   if liqflag == true
      es = saturationVaporPressure(T, al, bl, cl);

      % If requested, compute the derivative of es wrt temperature
      if nargout > 1
         des_dT = saturationVaporPressureDerivative(T, es, bl, cl);
      end
   else
      es = saturationVaporPressure(T, ai, bi, ci);

      % If requested, compute the derivative of es wrt temperature
      if nargout > 1
         des_dT = saturationVaporPressureDerivative(T, es, bi, ci);
      end
   end

   % If requested, compute dew point from relative humidity
   if nargout == 3
      if liqflag == true
         T_dew = dewPointTemperature(es * rh / 100, al, bl, cl);
      else
         T_dew = dewPointTemperature(es * rh / 100, ai, bi, ci);
      end
   end
end

function es = saturationVaporPressure(T, a, b, c)
   %SATURATIONVAPORPRESSURE Ambaum (2020) saturation vapor pressure
   es = a * exp(b ./ T) .* T .^ c; % [Pa]
end

function des_dT = saturationVaporPressureDerivative(T, es, b, c)
   %SATURATIONVAPORPRESSUREDERIVATIVE Derivative of es wrt temperature
   des_dT = es ./ T .* (c - b ./ T); % [Pa K-1]
end

function T_dew = dewPointTemperature(e_air, a, b, c)
   %DEWPOINTTEMPERATURE Dew point temperature from Ambaum (2020) inversion
   %
   %  For the Ambaum formula es = a * exp(b/T) * T^c, the dew point is the
   %  solution of: e_air = a * exp(b/T_dew) * T_dew^c
   %
   %  This is transcendental in T_dew. Use Newton iteration starting from
   %  an initial guess based on the ideal gas approximation.

   % Initial guess: invert the dominant exponential term
   T_dew = b ./ log(e_air ./ a); % ignoring T^c term

   % Newton iteration (typically converges in 3-5 steps)
   for iter = 1:10
      es_guess = a * exp(b ./ T_dew) .* T_dew .^ c;
      des_dT = es_guess ./ T_dew .* (c - b ./ T_dew);
      T_dew = T_dew - (es_guess - e_air) ./ des_dT;
   end
end

function ro_vap = saturationVaporDensity(es, T, Rv) %#ok<*DEFNU>
   %SATURATIONVAPORDENSITY Saturation vapor density from ideal gas law
   ro_vap = es ./ (Rv .* T); % [kg m-3]
end

function dro_vapdT = saturationVaporDensityDerivative(T, es, Rv, des_dT)
   %SATURATIONVAPORDENSITYDERIVATIVE Derivative of saturation vapor density
   dro_vapdT = (des_dT - es ./ T) ./ (Rv * T); % [kg m-3 K-1]
end
