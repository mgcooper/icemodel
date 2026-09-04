function [tf, f_res] = vapor_exchange_is_wet(f_ice, f_liq, f_res_por)
   %VAPOR_EXCHANGE_IS_WET Decide which phase a cell exchanges vapor with.
   %
   %  tf = icemodel.column.vapor_exchange_is_wet(f_ice, f_liq, f_res_por)
   %  [tf, f_res] = ...
   %     icemodel.column.vapor_exchange_is_wet(f_ice, f_liq, f_res_por)
   %
   % Returns true where a mobile liquid film exists, so vapor exchange goes
   % through the liquid phase at Lv. Elsewhere the exchange goes to ice at Ls.
   %
   % This function defines cell "wetness" for vapor exchange. A cell is "wet"
   % when its liquid fraction exceeds the residual floor set by the maximum of
   % capillary retention and the phase-fraction characteristic-curve
   % thermodynamic minimum. Liquid at or below that floor is immobile, so it
   % cannot supply evaporation.
   %
   % icemodel.surface.potential_surface_vapor_exchange calls this for the
   % surface demand partition, then applies the result through
   % icemodel.column.apply_vapor_transfer. icemodel.column.vapor_transport_terms
   % calls it during the enthalpy solve to pick each face's donor latent heat
   % (L_vap); icemodel.column.couple_vapor_step uses L_vap to ensure interior
   % transport uses the same phase for mass transfer as the donor latent heat.
   %
   % This is a different question from the one
   % icemodel.vapor.latent_enthalpy_switch answers. That function picks the
   % latent heat for vapor storage and conduction in the enthalpy solve, using a
   % fixed liquid-fraction threshold (f_liq_phase_switch_threshold). This
   % function picks the phase that supplies the mass exchange.
   %
   % icemodel.column.liquid_flux asks a third question and uses f_res_pore as
   % it's floor because capillarity holds water against flow. Its relative
   % saturation uses that same term. This function uses the maximum of the
   % capillary term and the phase fraction curve's thermodynamic minimum. That
   % value bounds what a cell can give up to phase change.
   %
   % Inputs
   %   f_ice     - Ice fraction [-].
   %   f_liq     - Liquid-water fraction [-].
   %   f_res_por - Residual liquid-water fraction per pore volume [-].
   %
   % Outputs
   %   tf        - True where the exchange goes through the liquid phase.
   %   f_res     - The volumetric residual floor used to decide 'tf' [-], for
   %               callers that need both.
   %
   % See also: icemodel.column.residual_water_fraction,
   %  icemodel.surface.potential_surface_vapor_exchange,
   %  icemodel.column.apply_vapor_transfer
   %
   %#codegen

   f_res = icemodel.column.residual_water_fraction(f_ice, f_liq, f_res_por);
   tf = f_liq > f_res;
end
