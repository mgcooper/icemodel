function budget = accumulate_vapor_exchange( ...
      budget, d_pevp, d_vap_liq, d_vap_ice, d_rof, dz)
   %ACCUMULATE_VAPOR_EXCHANGE Budget one substep's surface vapor exchange.
   %
   %  budget = icemodel.column.accumulate_vapor_exchange( ...
   %     budget, d_pevp, d_vap_liq, d_vap_ice, d_rof, dz)
   %
   % Accumulate applied liquid and ice increments [-], condensation overflow
   % [-], and potential demand on the Lv basis [-]. DZ is cell thickness [m].
   % The budget uses:
   %
   %   applied vapor energy = ro_liq * (Ls * vapor_solid
   %                          + Lv * vapor_liquid
   %                          + Lv * condensation_overflow)
   %
   % The difference between potential and applied exchange is the vapor-energy
   % residual (unapplied demand).
   %
   % See also: icemodel.column.budget_surface_mass_balance,
   %  icemodel.surface.apply_surface_vapor_exchange,
   %  icemodel.column.initialize_budget_state
   %
   %#codegen

   persistent ro_ice ro_liq Lv
   if isempty(ro_ice)
      [ro_ice, ro_liq, Lv] = ...
         icemodel.physicalConstant('ro_ice', 'ro_liq', 'Lv');
   end

   % Add the surface vapor demand for the top-cell volume.
   budget.mass_budget_vapor_potential_j_m2 = ...
      budget.mass_budget_vapor_potential_j_m2 ...
      + ro_liq * Lv * d_pevp * dz(1);

   % Budget the applied storage increments [mwe]. Scale each cell by its dz
   % before the sum because vapor exchange can reach more than one cell and
   % their thicknesses can differ for a nonuniform mesh.
   budget.mass_budget_vapor_liquid_mwe = budget.mass_budget_vapor_liquid_mwe ...
      + sum(d_vap_liq .* dz);
   budget.mass_budget_vapor_solid_mwe = budget.mass_budget_vapor_solid_mwe ...
      + sum(d_vap_ice .* dz) * ro_ice / ro_liq;

   % Assign the accumulated condensation overflow for this forcing step. Unlike
   % the increments budgeted above, d_rof is reset once per forcing step and
   % accumulated across substeps, so it arrives as the step total.
   budget.mass_budget_condensation_overflow_mwe = d_rof * dz(1);
end
