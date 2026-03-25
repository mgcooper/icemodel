function [Le, Ls] = latentEnthalpyWater(T)

   % Define the constants following Romps 2017
   cvv = 1418; % Specific heat capacity of water vapor at constant volume [J/kg/K]
   pv0 = 611.65; % Reference saturation vapor pressure over liquid and solid ice [Pa]
   T0 = 273.16; % Reference temperature (Triple point of water) [K]
   E0v = 2.3740e6; % Difference in specific internal energy between vapor and liquid at T0 [J/kg]
   E0s = 0.3337e6; % Difference in specific internal energy between liquid and solid at T0 [J/kg]
   Rv = 461; % Specific gas constant for water vapor [J/kg/K]
   cvl = 4119; % Specific heat capacity of liquid water at constant volume [J/kg/K]
   cvs = 1861; % Specific heat capacity of ice (solid water) at constant volume [J/kg/K]

   % Note: This is an exact implementation of Romps' expressions. However,
   % simpler expressions are given in Ambaum, and either can be made simpler by
   % combining constants.

   Le = latentHeatVaporization(T, E0v, Rv, cvv, cvl, T0);
   Ls = latentHeatSublimation(T, T0, E0v, E0s, Rv, cvv, cvs);
end

function Le = latentHeatVaporization(T, E0v, Rv, cvv, cvl, T0)
   %LATENTHEATVAPORIZATION Compute the specific latent enthalpy of vaporization
   %
   % This function computes the specific latent enthalpy of
   % evaporation/vaporization Le(T)
   % T - Temperature in K
   % T0, E0v, Rv, cvv, cvl are constants
   %
   % See also:

   % Romps expression:
   Le = E0v + Rv * T + (cvv - cvl) * (T - T0);

   % Ambaum:
   % L0 = E0v + Rv * T0;
   % Le = L0 + (cvv + Rv - cvl) * (T - T0);

   % Lv = Lv0 + (cpv_l - cpl) * (T - T0);
   % Ls = Ls0 + (cpv_s - cps) * (T - T0);

   % L0 = 2.501e6;
   % Le = L0 - (cpl - cpv) * (T - T0);
   % Le = 2.501e6 - 2180 * (T - T0);
end

function Ls = latentHeatSublimation(T, T0, E0v, E0s, Rv, cvv, cvs)
   %LATENTHEATSUBLIMATION Compute the specific latent enthalpy of sublimation
   %
   % This function computes the specific latent enthalpy of sublimation Ls(T)
   % T - Temperature in K
   % T0, E0v, E0s, Rv, cvv, cvs are constants
   %
   % See also:

   Ls = E0v + E0s + Rv * T + (cvv - cvs) * (T - T0);
end
