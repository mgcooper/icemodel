function k_sat = saturated_hydraulic_conductivity(f_ice, f_liq, kwargs)
   %SATURATED_HYDRAULIC_CONDUCTIVITY Snow saturated hydraulic conductivity.
   %
   %  k_sat = saturated_hydraulic_conductivity(f_ice, f_liq)
   %  k_sat = saturated_hydraulic_conductivity(f_ice, f_liq, ...
   %     method="shimizu1970", grainsz=2e-3)
   %  k_sat = saturated_hydraulic_conductivity(f_ice, f_liq, ...
   %     method="darcy", permeability=2.967e-8)
   %
   %  Returns the snow saturated hydraulic conductivity [m s-1] under one of
   %  three closure schemes selected by the `method` keyword.
   %
   %  Methods:
   %
   %    "colbeck1972" (default)
   %
   %       k_sat = col_k0 * exp(col_a_phi * f_por)
   %
   %    where f_por = 1 - f_ice. The constants col_k0 and col_a_phi are
   %    Colbeck's (1972) empirical fit at 0 deg C; f_liq is unused on this path
   %    and is accepted only for parity with the Shimizu method so the dispatch
   %    surface is uniform.
   %
   %    "shimizu1970"
   %
   %       k_sat = (shi_k0 * ro_liq * g / mu_water_ref) * grainsz^2 ...
   %               * exp(-shi_a_wat * f_wat)
   %
   %    where f_wat is the total water fraction
   %    (icemodel.column.water_fraction). The Shimizu permeability times (ro_liq
   %    g / mu) gives a hydraulic conductivity at the reference dynamic
   %    viscosity of water at 0 deg C. Equivalent to the "darcy" form below with
   %    kappa = shi_k0 * grainsz^2 * exp(-shi_a_wat * f_wat) evaluated at the
   %    *current* state.
   %
   %    "darcy"
   %
   %       k_sat = kappa * ro_liq * g / mu_water_ref
   %
   %    where kappa [m^2] is an explicit intrinsic permeability supplied by the
   %    caller via the `permeability` keyword. This is the right path when the
   %    reference solver pins kappa to a precomputed value (e.g. the SUMMA
   %    Laugh-Tests Colbeck verification, where kappa is Shimizu evaluated at
   %    the initial state and held constant for the duration of the run). For
   %    the three Colbeck experiments at initial state f_wat = 0.300, the
   %    matching Shimizu permeabilities are 2.967e-8 m^2 (exp1, exp2: grainsz =
   %    2e-3 m) and 2.967e-10 m^2 (exp3: grainsz = 2e-4 m).
   %
   %  Inputs
   %    f_ice  Ice volume fraction [1]
   %    f_liq  Liquid water volume fraction [1] (unused for "colbeck1972"
   %           and "darcy")
   %
   %  Name-value
   %    method        "colbeck1972" | "shimizu1970" | "darcy"
   %                  (default "colbeck1972")
   %    grainsz       Grain diameter [m]; required for "shimizu1970"
   %    permeability  Intrinsic permeability [m^2]; required for "darcy"
   %
   %  Output
   %    k_sat  Saturated hydraulic conductivity [m s-1] (column-shaped)
   %
   %  Empirical coefficients (col_k0, col_a_phi, shi_k0, shi_a_wat,
   %  mu_water_ref) and physical constants (ro_liq, ro_ice via water_fraction,
   %  g) come from icemodel.parameterLookup and icemodel.physicalConstant.
   %
   %  References
   %    Colbeck, S. C. (1972). A theory of water percolation in snow.
   %      J. Glaciol., 11(63), 369-385.
   %    Shimizu, H. (1970). Air permeability of deposited snow.
   %      Contrib. Inst. Low Temp. Sci., Ser. A, 22, 1-32.
   %    Clark, M. P. et al. (2017). An analytical test case for snow
   %      models. Water Resour. Res., 53(1), 909-922.
   %
   % See also: icemodel.column.liquid_flux, icemodel.column.infiltration,
   %  icemodel.column.water_fraction, icemodel.parameterLookup
   %
   %#codegen

   arguments
      f_ice
      f_liq
      kwargs.method (1, 1) string {mustBeMember(kwargs.method, ...
         ["colbeck1972", "shimizu1970", "darcy"])} = "colbeck1972"
      kwargs.grainsz      (1, 1) double = NaN
      kwargs.permeability (1, 1) double = NaN
   end

   % Load empirical coefficients and physical constants once. col_* are Colbeck
   % 1972 prefactor / porosity exponent; shi_* are Shimizu 1970 prefactor /
   % water-fraction exponent; mu_water_ref is the reference water dynamic
   % viscosity at 0 deg C; ro_liq and grav are physical constants used to
   % convert intrinsic permeability to hydraulic conductivity (k_sat = kappa *
   % ro_liq * g / mu).
   persistent col_k0 col_a_phi shi_k0 shi_a_wat mu_water_ref ro_liq grav
   if isempty(col_k0)
      [col_k0, col_a_phi, shi_k0, shi_a_wat, mu_water_ref] = ...
         icemodel.parameterLookup( ...
         'col_k0', 'col_a_phi', 'shi_k0', 'shi_a_wat', 'mu_water_ref');
      [ro_liq, grav] = icemodel.physicalConstant('ro_liq', 'gravity');
   end

   switch kwargs.method

      case "colbeck1972"
         % Colbeck 1972 empirical fit. Density-dependent only; f_liq is ignored
         % on this path. f_por = 1 - f_ice is the snow pore fraction (volume of
         % voids per unit total volume).
         f_por = 1.0 - f_ice;
         k_sat = col_k0 * exp(col_a_phi * f_por);

      case "shimizu1970"
         % Shimizu 1970 grain-size and water-content dependent fit. The Shimizu
         % permeability is kappa = shi_k0 * grainsz^2 * exp(-shi_a_wat * f_wat);
         % the conversion to hydraulic conductivity uses k_sat = kappa * ro_liq
         % * g / mu.
         if isnan(kwargs.grainsz)
            error( ...
               'icemodel:column:saturated_hydraulic_conductivity:missingGrainsz', ...
               'shimizu1970 method requires the grainsz keyword input [m].');
         end
         f_wat = icemodel.column.water_fraction(f_ice, f_liq);
         k_sat = (shi_k0 * ro_liq * grav / mu_water_ref) ...
            * kwargs.grainsz ^ 2 * exp(-shi_a_wat * f_wat);

      case "darcy"
         % Darcy form with caller-supplied intrinsic permeability kappa [m^2].
         % Returns a constant k_sat shaped like the column so downstream code
         % can index by layer. Useful when the reference solver pins kappa to a
         % precomputed value (SUMMA Laugh-Tests parameter trials evaluate
         % Shimizu at the initial state and freeze kappa for the run).
         if isnan(kwargs.permeability)
            error( ...
               'icemodel:column:saturated_hydraulic_conductivity:missingPermeability', ...
               'darcy method requires the permeability keyword input [m^2].');
         end
         k_sat = (kwargs.permeability * ro_liq * grav / mu_water_ref) ...
            * ones(numel(f_ice), 1);
   end
end
