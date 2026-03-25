function es = saturationVaporPressure(T, liqflag)
   %SATURATIONVAPORPRESSURE Canonical Romps/Ambaum saturation vapor pressure.
   %
   %  ES = saturationVaporPressure(T) computes saturation vapor pressure over
   %  liquid using the Clausius-Clapeyron equation with temperature-dependent
   %  latent enthalpy following Romps (2017) / Ambaum (2020).
   %
   %  ES = saturationVaporPressure(T, LIQFLAG) computes over liquid (true) or
   %  ice (false).
   %
   %  This file is the canonical reference implementation for the Romps/Ambaum
   %  saturation vapor pressure formula. Paper-specific notation (cvl, cvi) is
   %  retained for traceability. Production code uses the Rankine-Kirchhoff
   %  coefficients derived in VAPORINIT from these same expressions.
   %
   %  References:
   %     Ambaum (2020): https://romps.berkeley.edu/papers/pubdata/2021/ambaum/21ambaum.pdf
   %     Romps (2017): https://romps.berkeley.edu/papers/pubdata/2020/dewpoint/20dewpoint.pdf
   %
   % See also: VAPORINIT, icemodel.kernels.latentEnthalpyWater

   % Regarding optimal values, Romps optimized cvl, cvi, and Rv simultaneously,
   % finding one set of values that minimize errors in sum over liquid and ice
   % whereas Ambaum optimized values separately for liquid and ice.

   % -----------------------------------------------------------------------
   % Constants following Romps (2017) — retain for reference
   % -----------------------------------------------------------------------
   T0 = 273.16;    % Triple point temperature [K]
   Rv = 461;       % Specific gas constant for water vapor [J/kg/K]
   es0 = 611.65;   % Saturation vapor pressure at the triple point [Pa]
   cvv = 1418;     % Specific heat of vapor at constant volume [J/kg/K]
   % cpl = 4119;   % Romps: specific heat of liquid at constant volume [J/kg/K]
   % cpi = 1861;   % Romps: specific heat of ice at constant volume [J/kg/K]
   % cpv_l = cvv + Rv = 1879;  % Romps: vapor cp (same for liquid and ice)
   % cpv_i = cvv + Rv = 1879;

   % -----------------------------------------------------------------------
   % Constants following Ambaum (2020) — adopted for production
   % -----------------------------------------------------------------------
   % These are the actual triple point values. cpv is optimized separately for
   % liquid and ice (unlike Romps who uses one cpv for both).
   cpl = 4220;     % Specific heat of liquid water [J/kg/K]
   cpi = 2097;     % Specific heat of ice [J/kg/K]
   cpv_l = 2040;   % Ambaum optimal vapor cp over liquid, yields cpl - cpv_l = 2180
   cpv_i = 1885;   % Ambaum optimal vapor cp over ice, yields cpi - cpv_i = 212

   % Note: strictly, cpl and cpi are cvl and cvi but they are taken as equal
   % so cp is used instead of cv as cp appears in the published works.
   %
   % In these expressions, cvl = cpl, but cpv != cvv. The Le expressions are
   % equivalent under a substitution of cvv = cpv - Rv.
   %
   % (T/T0) ^ ((cpv - cvl) / Rv)  % Romps
   % (T/T0) ^ ((cpv - cpl) / Rv)  % Ambaum (identical since cvl = cpl)
   %
   % Le = E0v + Rv * T + (cvv - cvl) * (T - T0);   % Romps
   % Le = L0 + (cpv_l - cpl) * (T - T0);            % Ambaum
   %
   % The two expressions are equivalent under a substitution of cvv = cpv-Rv:
   % Le = E0v + Rv * T + (cvv - cvl) * (T - T0);    % Romps
   % Le = E0v + Rv * T + (cpv - Rv - cvl) * (T - T0);
   %    = E0v + Rv * T + (cpv - cvl) * (T - T0) - Rv * (T - T0);
   %    = E0v + Rv * T + (cpv - cvl) * (T - T0) - Rv * T + Rv * T0;
   %    = E0v + Rv * T0 + (cpv - cvl) * (T - T0);
   %    = L0 + (cpv - cvl) * (T - T0);
   % Which equals Ambaum's expression given cvl = cpl.

   if nargin < 2 || liqflag
      es = saturationVaporPressureLiquid(T, T0, es0, Rv, cpl, cpv_l);
   else
      es = saturationVaporPressureIce(T, T0, es0, Rv, cpi, cpv_i);
   end
end

function es = saturationVaporPressureLiquid(T, T0, es0, Rv, cvl, cpv)
   % Saturation vapor pressure over liquid pv*,l
   % Uses temperature-dependent latent enthalpy via latentEnthalpyWater.

   Le_T0 = icemodel.kernels.latentEnthalpyWater(T0);
   Le_T = icemodel.kernels.latentEnthalpyWater(T);

   es = es0 * (T/T0).^((cpv - cvl)/Rv) .* exp(Le_T0/(Rv * T0) - Le_T./(Rv * T));
end

function es = saturationVaporPressureIce(T, T0, es0, Rv, cvi, cpv)
   % Saturation vapor pressure over ice pv*,i
   % Uses temperature-dependent latent enthalpy via latentEnthalpyWater.

   [~, Ls_T0] = icemodel.kernels.latentEnthalpyWater(T0);
   [~, Ls_T] = icemodel.kernels.latentEnthalpyWater(T);

   es = es0 * (T/T0).^((cpv - cvi)/Rv) .* exp(Ls_T0/(Rv * T0) - Ls_T./(Rv * T));
end
