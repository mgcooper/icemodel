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
   % One function owns this decision. Both appliers ask it.
   % icemodel.surface.apply_surface_vapor_exchange routes the surface
   % exchange with it. icemodel.column.apply_vapor_transport picks with it the
   % latent heat for the energy demand it records. That heat converts a
   % transported mass into the energy a cell could not supply. Two different
   % criteria would make those disagree. A cell in the disagreement band would
   % gain or lose about twelve percent of its mass.
   %
   % icemodel.column.couple_vapor_transport does not ask it. That function
   % moves mass and makes no phase decision.
   %
   % icemodel.column.vapor_face_conductance does not ask it either, and that
   % is a defect. It takes the donor phase from latent_enthalpy_switch, so in
   % the band below its energy and the applier's mass use different latent
   % heats. Bead icemodel-55x carries the fix.
   %
   % This is a different question from the one
   % icemodel.vapor.latent_enthalpy_switch answers. That function picks the
   % latent heat for vapor storage and conduction in the enthalpy solve, on a
   % fixed liquid-fraction threshold. This one picks the phase that supplies
   % a mass exchange, on the residual floor. They agree over most of the
   % state space and disagree in a band whose width depends on f_ice.
   %
   % icemodel.column.liquid_flux asks a third form of the question and keeps
   % its own floor. Its mask is hydraulic: it tests against the capillary term
   % alone, because capillarity is what holds water against flow, and its
   % relative saturation is built on that same term. This predicate takes the
   % maximum of the capillary term and Jordan's thermodynamic minimum, which
   % bounds what a cell can give up to a phase change. The two floors are
   % therefore different quantities, not a duplicated rule. Closing the band
   % between them would move default-mode results.
   %
   % Inputs
   %   f_ice     - Ice fraction [-].
   %   f_liq     - Liquid-water fraction [-].
   %   f_res_por - Residual liquid-water fraction per pore volume [-].
   %
   % Outputs
   %   tf        - True where the exchange goes through the liquid phase.
   %   f_res     - The volumetric residual floor this decision used [-]. A
   %               caller that also needs the floor takes it from here. That
   %               keeps one evaluation per cell, and the floor is then the
   %               one the decision used.
   %
   % See also: icemodel.column.residual_water_fraction,
   %  icemodel.surface.apply_surface_vapor_exchange,
   %  icemodel.column.apply_vapor_transport
   %
   %#codegen

   f_res = icemodel.column.residual_water_fraction(f_ice, f_liq, f_res_por);
   tf = f_liq > f_res;
end
