function fields = budgetoutputs(kind)
   %BUDGETOUTPUTS Return the mass-budget output channels.
   %
   %  fields = icemodel.namelists.budgetoutputs()
   %  fields = icemodel.namelists.budgetoutputs(kind)
   %
   % KIND is 'all', 'first', 'last', or 'sum'.
   %
   % Start storage uses the first sample, end storage uses the last sample, and
   % step increments, energy, height, and count channels sum.
   %
   % Storage, phase, vapor, and remesh changes are positive into the column.
   % Condensation overflow and merge export are positive out of the column.
   %
   % Surface vapor exchange can leave vapor demand unapplied when the column
   % reaches a phase-storage limit. The applied vapor energy is:
   %   ro_liq * (Ls * vapor_solid + Lv * vapor_liquid
   %             + Lv * condensation_overflow).
   % The residual is vapor_potential minus the applied vapor energy.
   %
   % mass_budget_vapor_transport_solid_mwe and
   % mass_budget_vapor_transport_liquid_mwe store the phase changes from
   % interior vapor transport. A complete face transfer has zero column total.
   % Per-cell storage limits can make either total nonzero.
   %
   % df_vap_liq and df_vap_ice store the applied vapor-driven change from
   % surface exchange and interior transport.
   %
   % mass_budget_interior_merge_count counts every non-top removal. Only top
   % removals lower the fixed grid.
   %
   % mass_budget_top_export_solid_mwe and mass_budget_top_export_liquid_mwe
   % record the phase mass removed with a top cell. They are remesh diagnostics,
   % not surface mass fluxes.
   %
   % See also: icemodel.column.initialize_budget_state,
   %  icemodel.postprocess

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
      'mass_budget_vapor_transport_solid_mwe', ...
      'mass_budget_vapor_transport_liquid_mwe', ...
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
