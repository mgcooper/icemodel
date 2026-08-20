function fields = budgetoutputs(kind)
   %BUDGETOUTPUTS Return the canonical mass-budget output channels.
   %
   %  fields = icemodel.namelists.budgetoutputs()
   %  fields = icemodel.namelists.budgetoutputs(kind)
   %
   % KIND is 'all', 'first', 'last', or 'sum'. Diagnostic output configuration
   % and 15-minute to hourly retiming both use these aggregation classes.
   % Start storage uses the first sample, end storage uses the
   % last sample, and signed step increments, energies, heights, and counts sum.
   % Storage, phase, vapor, and remesh changes are positive into the column.
   % MWE channels use physical intrinsic phase densities and physical
   % liquid-water density as their fixed solver and reference basis.
   % Condensation overflow and merge export are positive out of the column;
   % unapplied vapor retains the sign of the rejected potential input.
   % mass_budget_unapplied_vapor_j_m2 is signed: positive is rejected
   % deposition energy and negative is unsatisfied sublimation demand.
   % mass_budget_vapor_redistribution_solid_mwe and _liquid_mwe are the
   % per-phase storage changes that production interior vapor transport
   % causes. When every increment applies in full, both sides of each
   % face carry the donor phase and both channels are zero. A per-cell
   % storage limit clamps donor and receiver increments independently, so
   % the channels record the clamp-induced storage changes and their sum
   % can be nonzero. The transport runs after
   % the surface vapor budget takes its storage baseline, so the surface
   % closure identity carries no redistribution term; the per-phase
   % storage closures consume these increments alongside the phase,
   % surface vapor, and remesh terms. Transport the per-cell limits
   % reject is visible through apply_vapor_transfer's unapplied outputs;
   % no standing channel records it.
   % Cross-phase transfer populates the redistribution channels.
   % Unconstrained same-phase redistribution can have nonzero face flux
   % while both channels remain zero.
   % mass_budget_interior_merge_count counts every non-top removal, including
   % a deepest-cell removal; only top removals contribute grid translation.
   % mass_budget_top_export_solid_mwe and mass_budget_top_export_liquid_mwe are
   % the mass a top removal exports from the column, split by phase. The
   % quantized top-deletion height is a grid-geometry counter and is never
   % mass. Neither export channel is a surface mass flux: a merge keeps the
   % mean of the pair, so the export over-counts what the removed cell held.
   %
   % Their SUM is always nonnegative, because the merged cell keeps exactly
   % half the pair's total water. The per-phase split has no such property.
   % merge_layers re-derives f_liq_C from the merged temperature, so the
   % liquid share is a solve result, not a mean. It is positive across the
   % realistic states checked, so there is no separate absolute-gross channel.
   % The liquid channel summed on its own has no fixed sign.

   if nargin == 0
      kind = 'all';
   end

   % Endpoint storage and step-ledger channels require different aggregation.
   first_fields = { ...
      'mass_budget_solid_start_mwe', ...
      'mass_budget_liquid_start_mwe'};
   last_fields = { ...
      'mass_budget_solid_end_mwe', ...
      'mass_budget_liquid_end_mwe'};
   sum_fields = { ...
      'mass_budget_phase_solid_mwe', ...
      'mass_budget_phase_liquid_mwe', ...
      'mass_budget_vapor_solid_mwe', ...
      'mass_budget_vapor_liquid_mwe', ...
      'mass_budget_remesh_solid_mwe', ...
      'mass_budget_remesh_liquid_mwe', ...
      'mass_budget_cloned_bottom_solid_mwe', ...
      'mass_budget_merge_export_solid_mwe', ...
      'mass_budget_vapor_potential_j_m2', ...
      'mass_budget_condensation_overflow_mwe', ...
      'mass_budget_unapplied_vapor_j_m2', ...
      'mass_budget_vapor_redistribution_solid_mwe', ...
      'mass_budget_vapor_redistribution_liquid_mwe', ...
      'mass_budget_top_deletion_count', ...
      'mass_budget_top_deletion_height_m', ...
      'mass_budget_top_export_solid_mwe', ...
      'mass_budget_top_export_liquid_mwe', ...
      'mass_budget_interior_merge_count'};

   switch lower(char(kind))
      case 'all'
         fields = [first_fields, last_fields, sum_fields];
      case 'first'
         fields = first_fields;
      case 'last'
         fields = last_fields;
      case 'sum'
         fields = sum_fields;
      otherwise
         error('icemodel:namelists:budgetoutputs:kind', ...
            'unsupported budget-output kind: %s', kind)
   end
end
