function [es, des_dT, d2es_dT2] = vappress(T, liqflag)
   %vappress Compute saturation vapor pressure over liquid or ice.
   %
   %  ES = icemodel.vapor.vappress(T, LIQFLAG) computes saturation vapor pressure
   %  using the Ambaum (2020) / Romps (2021) Rankine-Kirchhoff formula:
   %
   %     ES = A * exp(B / T) * T ^ C   [Pa]
   %
   %  where (al, bl, cl) are coefficients over liquid and (ai, bi, ci) are
   %  over ice, obtained from icemodel.parameterLookup.
   %
   %  [ES, DES_DT] = icemodel.vapor.vappress(T, LIQFLAG) also computes the first
   %  temperature derivative:
   %
   %     d(es)/dT = es / T * (c - b / T)   [Pa K-1]
   %
   %  [ES, DES_DT, D2ES_DT2] = icemodel.vapor.vappress(T, LIQFLAG) also computes the
   %  second temperature derivative:
   %
   %     d2(es)/dT2 = es / T^2 * ((c - 1) * (c - 2*b/T) + b^2/T^2)
   %
   %  This is equivalent to differentiating des_dT = es/T * (c - b/T)
   %  directly and is the pressure-side half of the derivative chain used
   %  by icemodel.vapor.vapordensity.
   %
   % Notes
   %
   %  des_dT      = es     / T   * (c - b/T)
   %  dro_vapdT   = ro_vap / T   * (c - b/T - 1)
   %
   %  d2es_dT2    = es     / T^2 * ((c-1) * (c - 2*b/T)     + b^2 / T^2)
   %  d2ro_vapdT2 = ro_vap / T^2 * ((c-2) * (c - 2*b/T - 1) + b^2 / T^2)
   %
   % Compute directly from es terms (forms above involve expensive ^'s):
   %
   %  dro_vapdT   = (des_dT - es / T) / (Rv * T)
   %
   %  d2es_dT2    = c/T     * (des_dT    + es/T)     - des_dT / T    * (2 + b/T)
   %  d2ro_vapdT2 = (c-1)/T * (dro_vapdT + ro_vap/T) - dro_vapdT / T * (2 + b/T)
   %              = (c-1) * des_dT / (Rv * T^2) - dro_vapdT / T * (2 + b/T)
   %              = (d2es_dT2 - 2*des_dT/T + 2*es/T^2) / (Rv*T)
   %
   % Units: [Pa = J m-3 = N m-2 = kg m-1 s-2]
   %
   % See also: icemodel.vapor.vapordensity, icemodel.vapor.tdewpoint,
   %  icemodel.kernels.buckVaporModel
   %
   %#codegen

   % Ambaum (2020) / Romps (2021) Rankine-Kirchhoff coefficients (i=ice, l=liq)
   persistent al bl cl ai bi ci
   if isempty(al)
      [al, bl, cl, ai, bi, ci] = icemodel.parameterLookup( ...
         'al', 'bl', 'cl', 'ai', 'bi', 'ci');
   end

   if nargin < 2
      liqflag = false;
   end

   % saturation vapor pressure over liquid
   if liqflag
      es = saturationVaporPressure(T, al, bl, cl);

      % derivative wrt temperature
      if nargout > 1
         des_dT = saturationVaporPressureDerivative(T, es, bl, cl);
      end
      if nargout > 2
         d2es_dT2 = saturationVaporPressureSecondDerivative(T, es, bl, cl);
      end
   else
      % saturation vapor pressure over ice
      es = saturationVaporPressure(T, ai, bi, ci);

      % derivative wrt temperature
      if nargout > 1
         des_dT = saturationVaporPressureDerivative(T, es, bi, ci);
      end
      if nargout > 2
         d2es_dT2 = saturationVaporPressureSecondDerivative(T, es, bi, ci);
      end
   end
end

function es = saturationVaporPressure(T, a, b, c)
   %SATURATIONVAPORPRESSURE Ambaum (2020) / Romps (2021) sat vapor pressure
   es = a * exp(b ./ T) .* T .^ c; % [Pa]
end

function des_dT = saturationVaporPressureDerivative(T, es, b, c)
   %SATURATIONVAPORPRESSUREDERIVATIVE Derivative of es wrt temperature
   des_dT = es ./ T .* (c - b ./ T); % [Pa K-1]
end

function d2es_dT2 = saturationVaporPressureSecondDerivative(T, es, b, c)
   %SATURATIONVAPORPRESSURESECONDDERIVATIVE Second derivative of es wrt T
   b_over_T = b ./ T;
   d2es_dT2 = es ./ T .^ 2 .* ((c - 1) .* (c - 2 * b_over_T) + b_over_T .^ 2);
end
