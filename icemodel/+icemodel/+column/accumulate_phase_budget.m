function budget = accumulate_phase_budget(budget, xT, xf_ice, xf_liq, ...
      T, f_ice, f_liq, dz)
   %ACCUMULATE_PHASE_BUDGET Add one substep's phase-change storage increments.
   %
   %  budget = icemodel.column.accumulate_phase_budget( ...
   %     budget, xT, xf_ice, xf_liq, T, f_ice, f_liq, dz)
   %
   % Call this once per accepted substep, after the column enthalpy solve and
   % before any surface vapor exchange rewrites the top layer, so the budget
   % separates thermodynamic phase change from vapor exchange.
   %
   % The checkpoint baseline is integrated fresh here rather than carried
   % from the previous substep's post-remesh state, so the storage closure
   % keeps its power to detect any state mutation the budget did not see.
   %
   % Inputs
   %   budget              - Forcing-step budget (see
   %                         icemodel.column.initialize_budget_state).
   %   xT, xf_ice, xf_liq  - Checkpointed column state entering the substep.
   %   T, f_ice, f_liq     - Column state after the accepted phase-change solve.
   %   dz                  - Control-volume thickness [m].
   %
   % Outputs
   %   budget              - Budget with the accumulated phase increments and
   %                         the post-phase baselines in budget.substep
   %                         (solid_p, liquid_p), the reference state for
   %                         icemodel.column.accumulate_vapor_budget.
   %
   % See also: icemodel.column.initialize_budget_state
   %  icemodel.column.integrate_column_budget
   %  icemodel.column.accumulate_vapor_budget
   %
   %#codegen

   % Compute the checkpointed and solved total solid and liquid mass in mwe
   % on one fixed storage basis, so the increment cannot absorb a change of
   % reference density. Both calls must therefore use the same densities.
   [solid_0, liquid_0] = icemodel.column.integrate_column_budget( ...
      xT, xf_ice, xf_liq, dz);
   [solid_p, liquid_p] = icemodel.column.integrate_column_budget( ...
      T, f_ice, f_liq, dz);

   % Difference the checkpointed and solved states.
   budget.mass_budget_phase_solid_mwe = ...
      budget.mass_budget_phase_solid_mwe + (solid_p - solid_0);
   budget.mass_budget_phase_liquid_mwe = ...
      budget.mass_budget_phase_liquid_mwe + (liquid_p - liquid_0);

   % Hand the post-phase baselines to the surface vapor budget.
   budget.substep.solid_p = solid_p;
   budget.substep.liquid_p = liquid_p;
end
