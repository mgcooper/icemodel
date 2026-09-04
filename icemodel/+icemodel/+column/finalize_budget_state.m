function budget = finalize_budget_state(budget, T, f_ice, f_liq, dz)
   %FINALIZE_BUDGET_STATE Record the storage end endpoints for one forcing step.
   %
   %  budget = icemodel.column.finalize_budget_state(budget, T, f_ice, ...
   %     f_liq, dz)
   %
   % The driver calls this once per forcing step, after the substep loop.
   % Together with the start endpoints that initialize_budget_state
   % records, these verify the per-phase storage closures: endpoint delta
   % equals phase plus vapor plus remesh plus transport increments.
   %
   % See also: icemodel.column.initialize_budget_state,
   %  icemodel.column.integrate_column_budget
   %
   %#codegen

   % Integrate solid and liquid phase storage [mwe] at the end of the step.
   [solid_end, liquid_end] = ...
      icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz);
   budget.mass_budget_solid_end_mwe = solid_end;
   budget.mass_budget_liquid_end_mwe = liquid_end;
end
