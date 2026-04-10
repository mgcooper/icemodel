function k = thermal_conductivity_firn(T, rho_firn, reference, depth)
   %THERMAL_CONDUCTIVITY_FIRN Firn conductivity from temperature and density.
   %
   %  K = icemodel.kernels.thermal_conductivity_firn(T, rho_firn)
   %  K = icemodel.kernels.thermal_conductivity_firn(T, rho_firn, reference)
   %  K = icemodel.kernels.thermal_conductivity_firn(T, rho_firn, reference, depth)
   %
   % Description:
   %  Returns firn thermal conductivity K [W m-1 K-1] using firn-specific
   %  or firn-applicable references. The primary interface is temperature
   %  plus firn density, so temperature-dependent porous-ice formulas can be
   %  compared directly. The Oster and Albert (2020) depth relation is kept
   %  as a special case via the optional depth input.
   %
   % Inputs:
   %  T         - Temperature [K]
   %  rho_firn  - Firn density [kg m-3]. Ignored when reference is
   %              "oster_2020_depth".
   %  reference - Conductivity reference selector
   %  depth     - Depth [m], required only when reference is
   %              "oster_2020_depth"
   %
   % Options:
   %  "calonne_2019_eq5"   - Calonne et al. (2019) Eq. 5 logistic blend
   %  "calonne_2019_eq1"   - Calonne et al. (2019) Eq. 1 high-density firn branch
   %  "schwerdtfeger_1963" - Schwerdtfeger (1963) mixture-theory form
   %  "yen_1981"           - Yen (1981) density-ratio fit
   %  "van_dusen_1929"     - van Dusen (1929) lower-limit density fit
   %  "oster_2020_density" - Oster and Albert (2020) density relation
   %  "oster_2020_depth"   - Oster and Albert (2020) depth relation
   %
   % References:
   %
   %   Oster, M. and Albert, M. (2020), The thermal conductivity of firn,
   %   Journal of Glaciology.
   %
   %   Calonne et al. (2019), Thermal Conductivity of Snow, Firn, and
   %   Porous Ice from 3-D Image-Based Computations, Geophysical Research
   %   Letters, 46(22), 13079-13089.
   %
   %   Schwerdtfeger, P. (1963), The thermal properties of snow and ice.
   %
   %   Yen, Y.C. (1981), Review of thermal properties of snow, ice, and sea ice.
   %
   %   van Dusen, M.S. (1929), Thermal conductivity of non-metallic solids.
   %
   % See also: icemodel.column.firn_thermal_conductivity,
   %  thermal_conductivity_snow

   arguments
      T {mustBeNumeric} = 263.15
      rho_firn {mustBeNumeric} = 600
      reference (1, 1) string {mustBeMember(reference, [ ...
         "calonne_2019_eq5", "calonne_2019_eq1", "schwerdtfeger_1963", ...
         "yen_1981", "van_dusen_1929", "oster_2020_density", ...
         "oster_2020_depth"])} = "calonne_2019_eq5"
      depth = []
   end

   switch reference
      case "oster_2020_density"
         % From Oster and Albert (2020), Journal of Glaciology.
         k = 0.144 .* exp(0.00308 .* rho_firn) + zeros(size(T));

      case "oster_2020_depth"
         % From Oster and Albert (2020), Journal of Glaciology.
         if isempty(depth)
            error('depth is required when reference is "oster_2020_depth".');
         end
         k = 0.536 .* exp(0.0144 .* depth);

      otherwise
         ro_ice = icemodel.physicalConstant('ro_ice');
         f_ice = rho_firn ./ ro_ice;
         k = icemodel.kernels.thermal_conductivity_snow( ...
            T, f_ice, ro_ice, rho_firn, reference);
   end
end
