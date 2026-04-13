function [Qe, dQe_dT_sfc] = latent_heat_flux(es_sfc, ea_atm, De, stability, ...
      psfc, ro_air_Lv, des_sfc_dT, dstability_dT_sfc)
   %LATENT_HEAT_FLUX Compute the turbulent latent heat flux.
   %
   %  Qe = icemodel.surface.turbulence.bulk_richardson.latent_heat_flux(...)
   %  [Qe, dQe_dT_sfc] = ...
   %     icemodel.surface.turbulence.bulk_richardson.latent_heat_flux(...)
   %
   %  Qe = ro_air * L * De_e * stability * (ea_atm - es_sfc) * epsilon/psfc;
   %  [W m-2] = [kg m-3] * [J kg-1] * [m s-1] * [-] * [Pa] * [Pa-1]
   %
   % epsilon = Rd/Rv where Rd = 287 J/K/kg is the gas constant for dry air and
   % Rv = 461.5 J/kg/K is the gas constant for water vapor. psfc is atmospheric
   % pressure, ea_atm is atmospheric (2-m) vapor pressure, and es_sfc is surface
   % saturation vapor pressure.
   %
   % In this implementation, ro_air_Lv is either ro_air * Ls or ro_air * Lv
   % depending on the surface liquid water state, updated each timestep in
   % icemodel.timestepping.updatesubstep. The density used is the fixed
   % reference air density ro_air (a constant). A more physics-consistent
   % alternative would be to use the local moist-air density, as the
   % monin_obukhov scheme already does via icemodel.vapor.moist_air_density;
   % that path divides the precomputed ro_air_Lv by ro_air_ref to recover the
   % specific latent heat and then multiplies by the local ro_atm. If a
   % moist-air correction is desired for the bulk-Richardson path, the same
   % approach could be applied here.
   %
   % The exchange coefficient De is identical for latent and sensible heat
   % fluxes because the scalar exchange roughness lengths for latent and
   % sensible heat are assumed equal to the momentum exchange roughness length.
   %
   % When the derivative is requested, provide the temperature derivatives of
   % both the surface saturation vapor pressure and the stability factor so
   % the returned derivative is the full dQe/dT_sfc.
   %
   % See also: icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux
   %
   %#codegen
   persistent epsilon
   if isempty(epsilon)
      epsilon = icemodel.physicalConstant('epsilon');
   end
   Qe = ro_air_Lv .* De .* stability .* (ea_atm - es_sfc) .* epsilon ./ psfc;

   if nargout > 1
      dQe_dT_sfc = ro_air_Lv .* De .* epsilon ./ psfc .* ...
         ((ea_atm - es_sfc) .* dstability_dT_sfc - stability .* des_sfc_dT);
   end
end
