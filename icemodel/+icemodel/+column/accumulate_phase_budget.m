function budget = accumulate_phase_budget(budget, xT, xf_ice, xf_liq, ...
      T, f_ice, f_liq, dz)
   %ACCUMULATE_PHASE_BUDGET Add one substep's phase-change storage increments.
   %
   %  budget = icemodel.column.accumulate_phase_budget( ...
   %     budget, xT, xf_ice, xf_liq, T, f_ice, f_liq, dz)
   %
   % Integrates the prior checkpoint and updated column state to compute total
   % melt/freeze phase change over one substep. Call this once per accepted
   % substep, after the column enthalpy solve and before any surface vapor
   % exchange modifies the top layer, so the budget separates thermodynamic
   % phase change from vapor exchange. Vapor functions add their own increments
   % to separate budget fields.
   %
   % Inputs
   %   budget              - Forcing-step budget (see
   %                         icemodel.column.initialize_budget_state).
   %   xT, xf_ice, xf_liq  - Checkpointed column state entering the substep.
   %   T, f_ice, f_liq     - Column state after the accepted phase-change solve.
   %   dz                  - Control-volume thickness [m].
   %
   % Outputs
   %   budget              - Budget with the phase change increments added.
   %
   % See also: icemodel.column.initialize_budget_state,
   %  icemodel.column.integrate_column_budget,
   %  icemodel.column.accumulate_vapor_exchange
   %
   %#codegen

   % Integrate the prior checkpoint state and the updated accepted state.
   [solid_0, liquid_0] = icemodel.column.integrate_column_budget( ...
      xT, xf_ice, xf_liq, dz);
   [solid_p, liquid_p] = icemodel.column.integrate_column_budget( ...
      T, f_ice, f_liq, dz);

   % Difference them to get the column-integrated phase change increments for
   % this substep and accumulate the increments across substeps.
   budget.mass_budget_phase_solid_mwe = ...
      budget.mass_budget_phase_solid_mwe + (solid_p - solid_0);
   budget.mass_budget_phase_liquid_mwe = ...
      budget.mass_budget_phase_liquid_mwe + (liquid_p - liquid_0);
end
