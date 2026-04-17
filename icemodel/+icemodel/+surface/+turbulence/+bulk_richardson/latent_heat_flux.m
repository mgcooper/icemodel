function [Qe, dQe_dT_sfc] = latent_heat_flux(es_sfc, ea_atm, De, stability, ...
      psfc, ro_air_Lv, des_sfc_dT, dstability_dT_sfc)
   %LATENT_HEAT_FLUX Compute the turbulent latent heat flux.
   %
   %  Qe = icemodel.surface.turbulence.bulk_richardson.latent_heat_flux(...)
   %  [Qe, dQe_dT_sfc] = ...
   %     icemodel.surface.turbulence.bulk_richardson.latent_heat_flux(...)
   %
   %  Qe = H_e * stability * (ea_atm - es_sfc);
   %  [W m-2] = [W m-2 Pa-1] * [-] * [Pa]
   %
   % where H_e = hv_atm * De_e = ro_atm * L * De * epsilon / psfc is the
   % latent heat transport prefactor precomputed at each substep. This uses
   % the local moist-air density ro_atm rather than the dry-air reference
   % density, giving a physically consistent moist-air correction.
   %
   % When the derivative is requested, provide the temperature derivatives of
   % both the surface saturation vapor pressure and the stability factor so the
   % returned derivative is the full dQe/dT_sfc used in the newton solve rather
   % than only the fixed-stability partial used in the linearization.
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
