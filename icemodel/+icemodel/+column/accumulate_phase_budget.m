function [ledger, solid_p, liquid_p, phase_solid, phase_liquid] = ...
      accumulate_phase_budget(ledger, xT, xf_ice, xf_liq, T, f_ice, f_liq, dz)
   %ACCUMULATE_PHASE_BUDGET Add one substep's phase-change storage increments.
   %
   %  [ledger, solid_p, liquid_p, phase_solid, phase_liquid] = ...
   %     icemodel.column.accumulate_phase_budget( ...
   %     ledger, xT, xf_ice, xf_liq, T, f_ice, f_liq, dz)
   %
   % Call this once per accepted substep, after the column enthalpy solve and
   % before any surface vapor exchange rewrites the top layer, so the ledger
   % separates thermodynamic phase change from vapor exchange.
   %
   % Inputs
   %   ledger              - Forcing-step ledger (see
   %                         icemodel.column.initialize_budget_state).
   %   xT, xf_ice, xf_liq  - Checkpointed column state entering the substep.
   %   T, f_ice, f_liq     - Column state after the accepted phase-change solve.
   %   dz                  - Control-volume thickness [m].
   %
   % Outputs
   %   ledger              - Ledger with the accumulated phase budgets.
   %   solid_p, liquid_p   - Post-phase storage [m w.e.], the reference state
   %                         for icemodel.column.accumulate_vapor_budget.
   %   phase_solid,        - Signed substep phase increments [m w.e.], needed to
   %   phase_liquid          close gross storage after remeshing.
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
   phase_solid = solid_p - solid_0;
   phase_liquid = liquid_p - liquid_0;

   ledger.mass_budget_phase_solid_mwe = ...
      ledger.mass_budget_phase_solid_mwe + phase_solid;
   ledger.mass_budget_phase_liquid_mwe = ...
      ledger.mass_budget_phase_liquid_mwe + phase_liquid;

   % Retain absolute substep activity so melt and refreezing within one forcing
   % step cannot cancel to a signed zero and hide the exchange that occurred.
   ledger.mass_budget_phase_solid_gross_mwe = ...
      ledger.mass_budget_phase_solid_gross_mwe + abs(phase_solid);
   ledger.mass_budget_phase_liquid_gross_mwe = ...
      ledger.mass_budget_phase_liquid_gross_mwe + abs(phase_liquid);
end
