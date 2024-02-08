function [es, des_dT, d2es_dT2, ro_vap, dro_vapdT, d2ro_vapdT2] = VAPPRESS2(T, liqflag)
   %VAPPRESS Compute saturation vapor pressure over liquid or solid water
   %
   %  ES = VAPPRESS(T, TF, LIQFLAG) Computes saturation vapor pressure
   %  using the Rankine–Kirchhoff formula, following Ambaum (2020):
   %
   %     ES = A * exp(B / T) * T ^ C
   %
   %  Where A, B, C are coefficients derived from physical constants given in
   %  Ambaum (2020).
   %
   %  [ES, dES_dT] = VAPPRESS(T, TF, LIQFLAG) Also computes the derivative of
   %  vapor with respect to temperature.
   %
   %  [ES, dES_dT, d2ES_dT2] = VAPPRESS(T, TF, LIQFLAG) Also computes the second
   %  derivative of vapor with respect to temperature.
   %
   %  Units: [Pa = J m-3 = N m-2 = kg m-1 s-2]
   %
   % See also: VAPORHEAT

   %  Tf = 273.16 is the triple point temperature.
   %  es = vapor pressure in Pa. (1 mb = 100 Pa).
   %
   %  https://romps.berkeley.edu/papers/pubdata/2021/ambaum/21ambaum.pdf
   %  https://romps.berkeley.edu/papers/pubdata/2020/dewpoint/20dewpoint.pdf

   % Define coefficients over water and ice
   if nargin < 2 || liqflag
      a = 9.7070e+24;
      b = -6716.9;
      c = -4.72885;
   else
      % Over ice
      a = 7.5948e+13;
      b = -6273.12;
      c = -0.45986;
   end
   
   % Compute saturation vapor pressure
   es = saturationVaporPressure(T, a, b, c);

   % Compute temperature derivatives of vapor pressure
   if nargout > 2
      [des_dT, d2es_dT2] = saturationVaporPressureDerivative(T, es, b, c);
   elseif nargout > 1
      des_dT = saturationVaporPressureDerivative(T, es, b, c);
   end
   
   % Compute vapor density and derivatives
   if nargout > 3
      Rv = 461;
      ro_vap = saturationVaporDensity(T, es, Rv);
      if nargout > 5
         [dro_vapdT, d2ro_vapdT2] = saturationVaporDensityDerivative( ...
            T, es, Rv, des_dT, b, c, ro_vap);
      elseif nargout > 4
         dro_vapdT = saturationVaporDensityDerivative( ...
            T, es, Rv, des_dT, b, c, ro_vap);
      end
   end
   
   % Actual expressions for latent enthalpies:
   % Lv = Lv0 + (cpv_l - cpl) * (T - T0);
   % Ls = Ls0 + (cpv_s - cps) * (T - T0);
end

%%
function es = saturationVaporPressure(T, a, b, c)
   %SATURATIONVAPORPRESSURE Compute saturation vapor pressure over water or ice
   es = a * exp(b ./ T) .* T .^ c; % [Pa]
end

function [des_dT, d2es_dT2] = saturationVaporPressureDerivative(T, es, b, c)
   %SATURATIONVAPORPRESSUREDERIVATIVE Temperature derivative of vapor pressure 
   
   % First derivative [Pa K-1]
   des_dT = es ./ T .* (c - b ./ T);

   if nargout == 2
      % Second derivative [Pa K-2]
      d2es_dT2 = c ./ T .* (des_dT + es ./ T) - des_dT ./ T .* (2 + b ./ T);
      
      % This form follows the form of des_dT somewhat but involves more powers
      % d2es_dT2 = es ./ T.^2 .* ((c-1) * (c - 2*b./T) + b^2 ./ T.^2);
   end
end

function ro_vap = saturationVaporDensity(T, es, Rv) %#ok<*DEFNU>
   %SATURATIONVAPORDENSITY Compute saturation vapor density
   ro_vap = es ./ (Rv .* T); % [kg m-3]
end

function [dro_vapdT, d2ro_vapdT2] = saturationVaporDensityDerivative(T, es, ...
      Rv, des_dT, b, c, ro_vap)
   %SATURATIONVAPORDENSITYDERIVATIVE Derivative of saturation vapor density
   dro_vapdT = (des_dT - es ./ T) ./ (Rv * T); % [kg m-3 K-1]
   
   % In terms of ro_vap (notice analogous form to des_dT):
   % dro_vapdT = ro_vap ./ T .* (c - b ./ T - 1);
   
   if nargout == 2
      % Second derivative % [kg m-3 K-2]
      
      % something is wrong with this one
      % d2ro_vapdT2 = (des_dT - es ./ T.^2 .* (2 * (c - b ./ T) - 1)) ./ (Rv * T);
      
      % In terms of ro_vap (notice analogous form to d2es_dT2):
      d2ro_vapdT2 = ro_vap ./ T.^2 .* ((c-2) * (c-1 - 2*b./T) + b^2 ./ T.^2);
   end
end
