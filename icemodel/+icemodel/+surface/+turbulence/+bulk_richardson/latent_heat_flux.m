function Qe = latent_heat_flux(es_sfc, ea_atm, De, stability, psfc, roL)
   %LATENT_HEAT_FLUX Compute the turbulent latent heat flux.
   %
   %  Qe = icemodel.surface.turbulence.bulk_richardson.latent_heat_flux(...)
   %
   %  Qe = ro_air * L * De_e * stability * (ea_atm - es_sfc) * epsilon/psfc;
   %  [W m-2] = [kg m-3] * [J kg-1] * [m s-1] * [-] * [Pa] * [Pa-1]
   %
   % epsilon = Rd/Rv where Rd = 287 J/K/kg is the gas constant for dry air and
   % Rv = 461.5 J/kg/K is the gas constant for water vapor. psfc is atmospheric
   % pressure, ea_atm is atmospheric (2-m) vapor pressure, and es_sfc is surface
   % saturation vapor pressure.
   %
   % In this implementation, roL is either ro_air * Ls or ro_air * Lv depending
   % on the surface liquid water state, updated each timestep in
   % icemodel.timestepping.updatesubstep.
   %
   % The exchange coefficient De is identical for latent and sensible heat
   % fluxes because the scalar exchange roughness lengths for latent and
   % sensible heat are assumed equal to the momentum exchange roughness length. 
   %
   % See also: icemodel.surface.turbulence.bulk_richardson.sensible_heat_flux
   %
   %#codegen
   persistent epsilon
   if isempty(epsilon)
      epsilon = icemodel.physicalConstant('epsilon');
   end
   Qe = roL * De * stability * (ea_atm - es_sfc) * epsilon / psfc;
end
