function [Le, Ls] = latentEnthalpyWater(T)
   %LATENTHALPYWATER Canonical Romps/Ambaum latent enthalpy reference.
   %
   %  [Le, Ls] = latentEnthalpyWater(T) computes the specific latent enthalpy
   %  of vaporization Le(T) and sublimation Ls(T) following Romps (2021).
   %
   %  This file is the canonical reference implementation for the Romps (2021)
   %  latent enthalpy expressions and their algebraic equivalence with the
   %  Ambaum (2020) formulation. Paper-specific notation (cvv, cvl, cvs) is
   %  retained here for direct traceability to the published equations.
   %  Production code uses icemodel-centric names (cp_ice, cp_liq, cpv_l,
   %  cpv_i) via VAPORINIT.
   %
   %  Note on temperature dependence: Le(T) and Ls(T) are linear in T,
   %  reflecting the Kirchhoff relation dL/dT = cpv - cp_phase. The
   %  Rankine-Kirchhoff saturation vapor pressure formula es = a*exp(b/T)*T^c
   %  absorbs this temperature dependence into the coefficients a, b, c. As a
   %  result, production code (VAPPRESS, VAPORHEAT, etc.) correctly treats L as
   %  constant at the reference value (Lv0 = 2.501e6, Ls0 = 2.834e6 J/kg) —
   %  the L(T) effect is already encoded in the Rankine-Kirchhoff exponents.
   %
   %  References:
   %     Ambaum (2020), "Accurate, simple equation for saturated vapour
   %        pressure over water and ice." QJRMS, 146, 4252-4258.
   %        DOI: 10.1002/qj.3899
   %     Romps (2021), "The Rankine-Kirchhoff approximations for moist
   %        thermodynamics." QJRMS, 147(741), 3493-3497.
   %        DOI: 10.1002/qj.4154
   %
   % See also: VAPORINIT, icemodel.kernels.saturationVaporPressure

   % -----------------------------------------------------------------------
   % Constants following Romps (2021)
   % -----------------------------------------------------------------------
   % These use constant-volume heat capacities (cv) as published in Romps.
   % Romps optimized cvl, cvs, and Rv simultaneously, finding one set of
   % values that minimizes errors in the sum over liquid and solid.

   cvv = 1418;     % Specific heat of vapor at constant volume [J/kg/K]
   pv0 = 611.65;   % Reference saturation vapor pressure at the triple point [Pa]
   T0 = 273.16;    % Triple point temperature [K]
   E0v = 2.3740e6; % Internal energy difference: vapor vs liquid at T0 [J/kg]
   E0s = 0.3337e6; % Internal energy difference: liquid vs solid at T0 [J/kg]
   Rv = 461;       % Specific gas constant for water vapor [J/kg/K]
   cvl = 4119;     % Specific heat of liquid water at constant volume [J/kg/K]
   cvs = 1861;     % Specific heat of ice at constant volume [J/kg/K]

   % -----------------------------------------------------------------------
   % Compute latent enthalpies (Romps exact expressions)
   % -----------------------------------------------------------------------
   Le = latentHeatVaporization(T, E0v, Rv, cvv, cvl, T0);
   Ls = latentHeatSublimation(T, T0, E0v, E0s, Rv, cvv, cvs);
end

function Le = latentHeatVaporization(T, E0v, Rv, cvv, cvl, T0)
   %LATENTHEATVAPORIZATION Specific latent enthalpy of vaporization Le(T).
   %
   % T - Temperature [K]
   % E0v, Rv, cvv, cvl, T0 - Constants (Romps 2021)

   % Romps (2021) expression:
   Le = E0v + Rv * T + (cvv - cvl) * (T - T0);

   % -----------------------------------------------------------------------
   % Ambaum (2020) equivalent expressions — retain for reference.
   %
   % Ambaum uses constant-pressure heat capacities (cp) and optimizes cpv
   % separately for liquid and ice (solid), yielding two optimal values:
   %   cpv_l = 2040 [J/kg/K] — optimal vapor cp over liquid
   %   cpv_i = 1885 [J/kg/K] — optimal vapor cp over ice
   %
   % In icemodel production code (VAPORINIT), cp_liq and cp_ice refer to the
   % measurement-consensus values at 0C (4218, 2093), not the Ambaum fitting
   % values (4220, 2097). The difference is negligible.
   %
   % The Romps and Ambaum expressions are algebraically equivalent under the
   % substitution cvv = cpv - Rv (since cv + R = cp for an ideal gas):
   %
   %   Romps:
   %     Le = E0v + Rv * T + (cvv - cvl) * (T - T0)
   %        = E0v + Rv * T + (cpv - Rv - cvl) * (T - T0)
   %        = E0v + Rv * T + (cpv - cvl) * (T - T0) - Rv * (T - T0)
   %        = E0v + Rv * T0 + (cpv - cvl) * (T - T0)
   %        = L0 + (cpv - cvl) * (T - T0)
   %
   %   which equals Ambaum's expression given cvl = cpl:
   %
   %   Ambaum: Le = Lv0 + (cpv_l - cp_liq) * (T - T0)
   %
   %   where Lv0 = E0v + Rv*T0 = 2.501e6 [J/kg]
   %
   %   Evaluating at T0: cpv_l - cp_liq = 2040 - 4220 = -2180
   %   So: Le = 2.501e6 - 2180*(T - T0)
   %
   % Romps:  (T/T0)^((cpv - cvl)/Rv) in saturation vapor pressure formula
   % Ambaum: (T/T0)^((cpv - cpl)/Rv) — identical since cvl = cpl
   %
   % The latent heat of sublimation follows the same pattern:
   %   Romps:  Ls = E0v + E0s + Rv*T + (cvv - cvs)*(T - T0)
   %   Ambaum: Ls = Ls0 + (cpv_i - cp_ice)*(T - T0)
   %           where Ls0 = E0v + E0s + Rv*T0 = 2.834e6 [J/kg]
   % -----------------------------------------------------------------------
end

function Ls = latentHeatSublimation(T, T0, E0v, E0s, Rv, cvv, cvs)
   %LATENTHEATSUBLIMATION Specific latent enthalpy of sublimation Ls(T).
   %
   % T - Temperature [K]
   % T0, E0v, E0s, Rv, cvv, cvs - Constants (Romps 2021)

   Ls = E0v + E0s + Rv * T + (cvv - cvs) * (T - T0);
end

% Older notes, may be clearer than above.
%
% Ambaum:
% L0 = E0v + Rv * T0;
% Le = L0 + (cvv + Rv - cvl) * (T - T0);

% Lv = Lv0 + (cpv_l - cpl) * (T - T0);
% Ls = Ls0 + (cpv_s - cps) * (T - T0);

% L0 = 2.501e6;
% Le = L0 - (cpl - cpv) * (T - T0);
% Le = 2.501e6 - 2180 * (T - T0);
