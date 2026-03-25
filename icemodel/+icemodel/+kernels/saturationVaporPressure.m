function es = saturationVaporPressure(T, liqflag)
   %SATURATIONVAPORPRESSURE Compute saturation vapor pressure over water or ice
   %
   % https://romps.berkeley.edu/papers/pubdata/2021/ambaum/21ambaum.pdf
   % https://romps.berkeley.edu/papers/pubdata/2020/dewpoint/20dewpoint.pdf

   % Regarding optimal values, Romps optimized cvl, cvs, and Rv simultaneously,
   % finding one set of values that minimize errors in sum over liquid and solid
   % whereas Ambaum optimized values separately for liquid and solid.

   % Define the constants following Romps 2017
   T0 = 273.16; % Reference temperature (Triple point of water) [K]
   Rv = 461; % Specific gas constant for water vapor [J/kg/K]
   es0 = 611.65; % Reference saturation vapor pressure over liquid and solid ice [Pa]
   cvv = 1418; % Specific heat capacity of water vapor at constant volume [J/kg/K]
   cpl = 4119; % Specific heat capacity of liquid water at constant volume [J/kg/K]
   cps = 1861; % Specific heat capacity of ice (solid water) at constant volume [J/kg/K]
   cpv_l = 1879; % = cvv + Rv;
   cpv_s = 1879; % = cvv + Rv;

   % Define the constants following Ambaum 2020. These are the actual triple
   % point values, and cpv is optimized separately for liquid and solid (see
   % saturationVaporPressure subfunction for Ambaum cpv values).
   cpl = 4220;
   cps = 2097;
   cpv_l = 2040; % Ambaum, optimal value, yields cpl - cpv_l = 2180
   cpv_s = 1885; % Ambaum, optimal value, yields cps - cpv_s = 212


   % Note: strictly, cpl and cps are cvl and cvs but they are taken as equal so
   % cp is used instead of cv as cp appears in the published works.
   %
   % In these expressions, cvl = cpl, but cpv != cvv. The Le expressions are
   % equivelent under a substitution of cvv = cpv - Rv.
   %
   % (T/T0) ^ (cpv - cvl) / Rv % Romps
   % (T/T0) ^ (cpv - cpl) / Rv % Ambaum
   %
   % Le = E0v + Rv * T + (cvv - cvl) * (T - T0);   % Romps
   % Le = L0 + (cpv_l - cpl) * (T - T0);           % Ambaum
   %
   % The two expressions are equivelent under a substitution of cvv = cpv-Rv:
   % Le = E0v + Rv * T + (cvv - cvl) * (T - T0);   % Romps
   % Le = E0v + Rv * T + (cpv - Rv - cvl) * (T - T0);
   %    = E0v + Rv * T + (cpv - cvl) * (T - T0) - Rv * (T - T0);
   %    = E0v + Rv * T + (cpv - cvl) * (T - T0) - Rv * T + Rv * T0;
   %    = E0v + Rv * T0 + (cpv - cvl) * (T - T0);
   %    = L0 + (cpv - cvl) * (T - T0);
   % Which equals Ambaum's expression given cvl = cpl.

   if nargin < 2 || liqflag
      es = saturationVaporPressureLiquid(T, T0, es0, Rv, cpl, cpv_l);
   else
      es = saturationVaporPressureSolid(T, T0, es0, Rv, cps, cpv_s);
   end
end

function es = saturationVaporPressureLiquid(T, T0, es0, Rv, cvl, cpv)
   % This function computes the saturation vapor pressure over liquid pv∗,l
   % T - Temperature in K
   % T0, pv0_star_l, Rv, E0v, cvv, cvl are constants

   Le_T0 = snowphysics.latentEnthalpyWater(T0);
   Le_T = snowphysics.latentEnthalpyWater(T);

   es = es0 * (T/T0).^((cpv - cvl)/Rv) .* exp(Le_T0/(Rv * T0) - Le_T./(Rv * T));
end

function es = saturationVaporPressureSolid(T, T0, es0, Rv, cvs, cpv)
   % This function computes the saturation vapor pressure over solid pv∗,s
   % T - Temperature in K
   % T0, pv0_star_s, Rv, E0v, E0s, cvv, cvs are constants

   [~, Ls_T0] = snowphysics.latentEnthalpyWater(T0);
   [~, Ls_T] = snowphysics.latentEnthalpyWater(T);

   es = es0 * (T/T0)^((cpv - cvs)/Rv) * exp(Ls_T0/(Rv * T0) - Ls_T/(Rv * T));
end
