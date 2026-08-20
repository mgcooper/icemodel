function budget = accumulate_redistribution_budget( ...
      budget, T, f_ice, f_liq, dz)
   %ACCUMULATE_REDISTRIBUTION_BUDGET Record interior vapor transport storage.
   %
   %  budget = icemodel.column.accumulate_redistribution_budget( ...
   %     budget, T, f_ice, f_liq, dz)
   %
   % icemodel.column.couple_vapor_step calls this once per accepted
   % substep, right after icemodel.column.apply_vapor_transfer realizes
   % the node-wise transport in the cells. The pre-transport baseline is
   % budget.substep (solid_v, liquid_v), written by
   % icemodel.column.accumulate_vapor_budget; the state is untouched
   % between that budget close and the transport, so no extra integration
   % is needed for the baseline.
   %
   % Interior transport moves vapor between cells and creates none. Both
   % sides of each face carry the donor phase, so when every increment
   % applies in full the per-phase totals are conserved and both channels
   % are zero. A per-cell storage limit in apply_vapor_transfer clamps
   % donor and receiver increments independently, so the channels record
   % the clamp-induced storage changes and their sum can be nonzero;
   % test_a_bound_interior_clamp_keeps_its_own_accounting exercises this
   % case.
   %
   % The transport runs after icemodel.column.accumulate_vapor_budget takes
   % its storage baseline, so no interior term reaches the surface vapor
   % channels, and the surface closure identity holds without a
   % redistribution correction:
   %
   %   vapor_potential = ro_liq * (Ls * vapor_solid + Lv * vapor_liquid
   %                     + Lv * condensation_overflow) + unapplied_vapor
   %
   % The per-phase storage closures use the increments recorded here. The
   % solid storage delta closes against the phase, surface-vapor, remesh, and
   % redistribution-solid terms. The liquid storage delta closes against the
   % corresponding liquid terms. Transport the per-cell limits reject is
   % visible through apply_vapor_transfer's unapplied outputs; no standing
   % channel records it.
   %
   % Inputs
   %   budget          - Forcing-step budget with the post-exchange baselines
   %                     in budget.substep.
   %   T, f_ice, f_liq - Column state after the redistribution step.
   %   dz              - Control-volume thickness [m].
   %
   % Output
   %   budget          - Budget with the per-phase redistribution increments
   %                     accumulated.
   %
   % See also: icemodel.column.apply_vapor_transfer,
   %  icemodel.column.accumulate_vapor_budget
   %
   %#codegen

   [solid_after, liquid_after] = ...
      icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz);

   % Add rather than assign: a forcing step can take several accepted
   % substeps and each one redistributes.
   budget.mass_budget_vapor_redistribution_solid_mwe = ...
      budget.mass_budget_vapor_redistribution_solid_mwe ...
      + (solid_after - budget.substep.solid_v);
   budget.mass_budget_vapor_redistribution_liquid_mwe = ...
      budget.mass_budget_vapor_redistribution_liquid_mwe ...
      + (liquid_after - budget.substep.liquid_v);
end
