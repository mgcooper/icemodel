function U_top = surface_vapor_mass_flux(d_vap_sfc, dz_top, dt)
   %SURFACE_VAPOR_MASS_FLUX Convert a surface vapor fraction to a mass flux.
   %
   %  U_top = icemodel.surface.surface_vapor_mass_flux(d_vap_sfc, dz_top, dt)
   %
   % D_VAP_SFC is the realized liquid-water-equivalent volume fraction of the
   % mass moved through the surface. Demand partitioning and state limits
   % happen upstream, so this is a unit conversion and carries no latent heat:
   %
   %   U_top = d_vap_sfc * ro_liq * dz_top / dt   [kg m-2 s-1]
   %
   % Positive is downward, into the column.
   %
   % Splitting the conversion from the phase correction is what keeps the
   % driver's accumulator a fraction. The driver sums fractions across
   % substeps, and the kilogram basis appears only here, inside the function
   % that needs it.
   %
   % Inputs
   %   d_vap_sfc - Liquid-water volume fraction the surface exchanged [-].
   %   dz_top    - Thickness of the top control volume [m].
   %   dt        - Length of the step the fraction spans [s].
   %
   % Output
   %   U_top     - Vapor mass flux at the surface face [kg m-2 s-1].
   %
   % See also: icemodel.surface.potential_surface_vapor_exchange,
   %  icemodel.surface.apply_surface_vapor_exchange
   %
   %#codegen

   persistent ro_liq
   if isempty(ro_liq)
      ro_liq = icemodel.physicalConstant('ro_liq');
   end

   U_top = d_vap_sfc * ro_liq * dz_top / dt;
end
