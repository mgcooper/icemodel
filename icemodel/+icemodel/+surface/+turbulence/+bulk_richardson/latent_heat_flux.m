function [Qe, dQe_dT_sfc] = latent_heat_flux(es_sfc, ea_atm, H_e, stability, ...
      des_sfc_dT, dstability_dT_sfc)
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
   % latent heat transport coefficient precomputed at each substep.
   %
   % Pass the temperature derivatives of the surface saturation vapor pressure
   % and stability factor when you request the derivative output. The function
   % then returns the full dQe/dT_sfc that the newton solve needs, not the
   % fixed-stability partial derivative that the linearization uses.
   %
   % See also: icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux
   %
   %#codegen

   Qe = H_e .* stability .* (ea_atm - es_sfc);

   if nargout > 1
      dQe_dT_sfc = H_e .* ...
         ((ea_atm - es_sfc) .* dstability_dT_sfc - stability .* des_sfc_dT);
   end
end
