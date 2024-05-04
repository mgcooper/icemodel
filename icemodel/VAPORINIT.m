function [al, bl, cl, as, bs, cs] = VAPORINIT()
   %VAPORINIT Initialize the vapor model.
   %
   % T0 - Triple point temperature [K]
   % es0
   % Rv
   % cvv
   % cvl
   %
   % https://romps.berkeley.edu/papers/pubdata/2021/ambaum/21ambaum.pdf
   % https://romps.berkeley.edu/papers/pubdata/2020/dewpoint/20dewpoint.pdf
   %
   % See also:
   %
   %#codegen

   % Define the constants following Ambaum 2020.
   T0 = 273.16;   % Reference temperature (Triple point of water) [K]
   Rv = 461;      % Specific gas constant for water vapor [J/kg/K]
   Lv0 = 2.501e6; %
   Ls0 = 2.834e6; %
   es0 = 611.65;  % Reference saturation vapor pressure over water [Pa]
   cpl = 4220;    %
   cps = 2097;

   % Optimal values of cp_vap over liquid and solid water from Ambaum 2020
   cpv_l = 2040;
   cpv_s = 1885;

   % Combined values
   cpv_l_star = cpv_l - cpl;
   cpv_s_star = cpv_s - cps;
   Lv0_star = Lv0 - cpv_l_star * T0;
   Ls0_star = Ls0 - cpv_s_star * T0;

   % Coefficients over liquid
   bl = -Lv0_star / Rv;
   cl = cpv_l_star / Rv;
   al = es0 * exp(-bl / T0) / T0 ^ cl;

   % Coefficients over solid
   bs = -Ls0_star / Rv;
   cs = cpv_s_star / Rv;
   as = es0 * exp(-bs / T0) / T0 ^ cs;
end
