function [es, des_dT, d2es_dT2, ro_vap, dro_vapdT, d2ro_vapdT2] = VAPPRESS2(T, liqflag)
   %VAPPRESS2 Compute saturation vapor pressure over liquid or ice.
   %
   %  ES = VAPPRESS2(T, LIQFLAG) computes saturation vapor pressure
   %  using the Ambaum (2020) Rankine-Kirchhoff formula:
   %
   %     ES = A * exp(B / T) * T ^ C   [Pa]
   %
   %  where (al, bl, cl) are coefficients over liquid and (ai, bi, ci) are
   %  over ice, obtained from icemodel.parameterLookup.
   %
   %  [ES, dES_dT] = VAPPRESS2(T, LIQFLAG) Also computes the derivative of
   %  saturation vapor pressure with respect to temperature.
   %
   %  [ES, dES_dT, d2ES_dT2] = VAPPRESS2(T, LIQFLAG) Also computes the
   %  second derivative.
   %
   %  [..., RO_VAP, DRO_VAPDT, D2RO_VAPDT2] = VAPPRESS2(T, LIQFLAG) Also
   %  computes saturation vapor density and its temperature derivatives.
   %
   %  Units: [Pa = J m-3 = N m-2 = kg m-1 s-2]
   %
   % See also: VAPPRESS, VAPORHEAT
   %
   %#codegen

   % Ambaum (2020) Rankine-Kirchhoff coefficients (i = ice, l = liquid)
   persistent al bl cl ai bi ci
   if isempty(al)
      [al, bl, cl, ai, bi, ci] = icemodel.parameterLookup( ...
         'al', 'bl', 'cl', 'ai', 'bi', 'ci');
   end

   persistent Rv
   if isempty(Rv)
      Rv = icemodel.physicalConstant('Rv');
   end

   % Select phase coefficients
   if nargin < 2 || liqflag
      a = al; b = bl; c = cl;
   else
      a = ai; b = bi; c = ci;
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
   % Lv = Lv0 + (cpv_l - cp_liq) * (T - T0);
   % Ls = Ls0 + (cpv_i - cp_ice) * (T - T0);
end

%%
function es = saturationVaporPressure(T, a, b, c)
   %SATURATIONVAPORPRESSURE Ambaum (2020) saturation vapor pressure
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
   %SATURATIONVAPORDENSITY Saturation vapor density from ideal gas law
   ro_vap = es ./ (Rv .* T); % [kg m-3]
end

function [dro_vapdT, d2ro_vapdT2] = saturationVaporDensityDerivative(T, es, ...
      Rv, des_dT, b, c, ro_vap)
   %SATURATIONVAPORDENSITYDERIVATIVE Derivative of saturation vapor density
   dro_vapdT = (des_dT - es ./ T) ./ (Rv * T); % [kg m-3 K-1]

   % In terms of ro_vap (notice analogous form to des_dT):
   % dro_vapdT = ro_vap ./ T .* (c - b ./ T - 1);

   if nargout == 2
      % Second derivative [kg m-3 K-2]
      d2ro_vapdT2 = ro_vap ./ T.^2 .* ((c-2) * (c-1 - 2*b./T) + b^2 ./ T.^2);
   end
end
