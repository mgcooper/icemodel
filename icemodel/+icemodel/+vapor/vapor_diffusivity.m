function De = vapor_diffusivity(T, Pa)
   %VAPOR_DIFFUSIVITY Effective water-vapor diffusion coefficient in porous ice.
   %
   %  De = icemodel.vapor.vapor_diffusivity(T, Pa)
   %
   %  Computes the effective water vapor diffusion coefficient in snow/ice:
   %
   %     De = De0 * 100000 / Pa * (T / Tf) ^ nd   [m2 s-1]
   %
   %  where De0 = 9e-5 [m2 s-1] is the reference effective diffusivity for water
   %  vapor in snow and nd is the temperature exponent (Anderson 1976 Eq. 3.13,
   %  Fig. 4.3; Jordan 1991 pg v. Nomenclature).
   %
   %  Inputs:
   %     T  - Temperature [K]
   %     Pa - Atmospheric pressure [Pa]
   %
   %  Output:
   %     De - Effective vapor diffusivity [m2 s-1]
   %
   % See also: icemodel.vapor.vapor_thermal_conductivity,
   %  icemodel.column.vapor_mass_transfer
   %
   %#codegen

   persistent nd De0 Tf
   if isempty(nd)
      Tf = icemodel.physicalConstant('Tf');
      [nd, De0] = icemodel.parameterLookup('nd', 'De0');
   end

   % Every call site omits Pa, so this default of 100000 Pa makes the
   % pressure-correction term equal to one. The vapor diffusivity formulation
   % needs further review before any call site passes a real pressure.
   if nargin < 2
      Pa = 100000;
   end

   De = De0 * 100000 ./ Pa .* (T ./ Tf) .^ nd;
end
