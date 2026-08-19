function [tf, f_res] = vapor_exchange_is_wet(f_ice, f_liq, f_res_por)
   %VAPOR_EXCHANGE_IS_WET Decide which phase a cell exchanges vapor with.
   %
   %  tf = icemodel.column.vapor_exchange_is_wet(f_ice, f_liq, f_res_por)
   %  [tf, f_res] = ...
   %     icemodel.column.vapor_exchange_is_wet(f_ice, f_liq, f_res_por)
   %
   % Returns true where a mobile liquid film exists, so vapor exchange goes
   % through the liquid phase at Lv. Elsewhere the exchange goes straight to
   % ice at Ls.
   %
   % A cell counts as wet when its liquid fraction exceeds the residual floor
   % that capillary retention and the Jordan thermodynamic minimum set
   % together. Liquid at or below that floor is held, not mobile, so it
   % cannot supply evaporation.
   %
   % One function owns this decision. The surface demand partition and
   % icemodel.column.couple_vapor_step consume it before they call the shared
   % icemodel.column.apply_vapor_transfer state mutator. The interior path
   % evaluates it on the state the solve converged on. Two different criteria
   % would make the selected
   % phase and face latent heat disagree. In the disagreement band, a cell
   % would gain or lose about twelve percent of its mass.
   %
   % icemodel.column.vapor_transport_faces also calls this predicate
   % to select the donor-cell latent heat. The solve energy and applied mass
   % then use the same latent heat in the disagreement band below.
   %
   % This is a different question from the one
   % icemodel.vapor.latent_enthalpy_switch answers. That function picks the
   % latent heat for vapor storage and conduction in the enthalpy solve, on a
   % fixed liquid-fraction threshold. This one picks the phase that supplies
   % a mass exchange, on the residual floor. They agree over most of the
   % state space and disagree in a band whose width depends on f_ice.
   %
   % icemodel.column.liquid_flux asks a third question and keeps its own
   % floor. Its hydraulic mask tests only the capillary term because
   % capillarity holds water against flow. Its relative saturation uses that
   % same term. This predicate uses the maximum of the capillary term and
   % Jordan's thermodynamic minimum. That value bounds what a cell can give
   % up to a phase change. The two floors represent different quantities,
   % not a duplicated rule. Closing the band between them would move
   % production results.
   %
   % Inputs
   %   f_ice     - Ice fraction [-].
   %   f_liq     - Liquid-water fraction [-].
   %   f_res_por - Residual liquid-water fraction per pore volume [-].
   %
   % Outputs
   %   tf        - True where the exchange goes through the liquid phase.
   %   f_res     - The volumetric residual floor this decision used [-]. A
   %               caller that needs the floor uses this output, so both use
   %               the same evaluation.
   %
   % See also: icemodel.column.residual_water_fraction,
   %  icemodel.surface.potential_surface_vapor_exchange,
   %  icemodel.column.apply_vapor_transfer
   %
   %#codegen

   f_res = icemodel.column.residual_water_fraction(f_ice, f_liq, f_res_por);
   tf = f_liq > f_res;
end
