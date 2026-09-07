function budget = accumulate_vapor_transport( ...
      budget, d_vap_liq_applied, d_vap_ice_applied, dz)
   %ACCUMULATE_VAPOR_TRANSPORT Budget one substep's interior vapor transport.
   %
   %  budget = icemodel.column.accumulate_vapor_transport( ...
   %     budget, d_vap_liq_applied, d_vap_ice_applied, dz)
   %
   % The inputs are applied per-cell increments [-]. D_VAP_ICE_APPLIED uses
   % ice-fraction units and DZ is cell thickness [m]. A fully applied face
   % transfer has zero column total, but cv storage limits can make the total
   % nonzero if transfer is rejected.
   %
   % See also: icemodel.column.couple_vapor_step,
   %  icemodel.column.apply_vapor_transfer,
   %  icemodel.column.initialize_budget_state
   %
   %#codegen

   persistent ro_ice ro_liq
   if isempty(ro_ice)
      [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   end

   % Add this substep's liquid and ice transport changes in metres water
   % equivalent. Scale each cell by its dz before the sum because their
   % thicknesses can differ for a nonuniform mesh.
   budget.mass_budget_vapor_transport_liquid_mwe = ...
      budget.mass_budget_vapor_transport_liquid_mwe ...
      + sum(d_vap_liq_applied .* dz);
   budget.mass_budget_vapor_transport_solid_mwe = ...
      budget.mass_budget_vapor_transport_solid_mwe ...
      + sum(d_vap_ice_applied .* dz) * ro_ice / ro_liq;
end
