function [al, bl, cl, ai, bi, ci] = initialize_vapor_model()
   %INITIALIZE_VAPOR_MODEL Initialize Ambaum (2020) Rankine-Kirchhoff vapor
   % coefficients.
   %
   %  [AL, BL, CL, AI, BI, CI] = icemodel.vapor.initialize_vapor_model() computes
   %  the Rankine-Kirchhoff
   %  coefficients for the Ambaum (2020) saturation vapor pressure formula:
   %
   %     es = a * exp(b / T) * T ^ c   [Pa]
   %
   %  where (al, bl, cl) are over liquid and (ai, bi, ci) are over ice.
   %
   %  The physical constants come from icemodel.physicalConstant. The optimal
   %  vapor heat capacities (cpv_l, cpv_i) come from Ambaum (2020). This
   %  function uses the measurement-consensus values of cp_liq and cp_ice
   %  from physicalConstant (4218, 2093). They differ a little from the
   %  Ambaum fitting values (4220, 2097), and the difference in the derived
   %  coefficients is negligible.
   %
   %  For the algebraic derivation and equivalence with Romps (2021), see:
   %     icemodel.kernels.latentEnthalpyWater
   %     icemodel.kernels.saturationVaporPressure
   %
   %  References:
   %     Ambaum (2020), "Accurate, simple equation for saturated vapour
   %        pressure over water and ice." QJRMS, 146, 4252-4258.
   %        DOI: 10.1002/qj.3899
   %     Romps (2021), "The Rankine-Kirchhoff approximations for moist
   %        thermodynamics." QJRMS, 147(741), 3493-3497.
   %        DOI: 10.1002/qj.4154
   %
   % See also: icemodel.parameterLookup, icemodel.physicalConstant,
   %  icemodel.kernels.latentEnthalpyWater,
   %  icemodel.kernels.saturationVaporPressure
   %
   %#codegen

   % Obtain physical constants from the canonical source.
   [T0, Rv, Lv0, Ls0, es0, cp_liq, cp_ice, cpv_l, cpv_i] ...
      = icemodel.physicalConstant( ...
      'Tf', 'Rv', 'Lv', 'Ls', 'es0', 'cp_liq', 'cp_ice', 'cpv_l', 'cpv_i');

   % Combined heat capacity differences [J/kg/K]
   cpv_l_star = cpv_l - cp_liq;   % = 2040 - 4218 = -2178
   cpv_i_star = cpv_i - cp_ice;   % = 1885 - 2093 = -208

   % Adjusted reference latent heats [J/kg]
   Lv0_star = Lv0 - cpv_l_star * T0;
   Ls0_star = Ls0 - cpv_i_star * T0;

   % Rankine-Kirchhoff coefficients over liquid
   bl = -Lv0_star / Rv;
   cl = cpv_l_star / Rv;
   al = es0 * exp(-bl / T0) / T0 ^ cl;

   % Rankine-Kirchhoff coefficients over ice
   bi = -Ls0_star / Rv;
   ci = cpv_i_star / Rv;
   ai = es0 * exp(-bi / T0) / T0 ^ ci;
end
