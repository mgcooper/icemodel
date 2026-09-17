function budget = accumulate_remesh_budget(budget, remesh)
   %ACCUMULATE_REMESH_BUDGET Add one substep's remesh events to the budget.
   %
   %  budget = icemodel.column.accumulate_remesh_budget(budget, remesh)
   %
   % Records the column storage changes caused by merging layers and cloning the
   % bottom layer. icemodel.column.merge_thin_layers calls this once per
   % accepted substep. These numerical mass changes close the column storage
   % budget.
   %
   % Inputs
   %   budget - Forcing-step budget carrying the phase and vapor budgets.
   %   remesh - Event ledger from icemodel.column.initialize_remesh_ledger.
   %
   % Outputs
   %   budget - Budget with the remesh, layer removal/insertion, and
   %            grid-translation changes accumulated.
   %
   % See also: icemodel.column.merge_thin_layers,
   %  icemodel.column.initialize_remesh_ledger
   %
   %#codegen

   % Net remeshing exchange: what the column gained or lost across the merge.
   budget.mass_budget_remesh_solid_mwe = ...
      budget.mass_budget_remesh_solid_mwe + remesh.solid_mwe;
   budget.mass_budget_remesh_liquid_mwe = ...
      budget.mass_budget_remesh_liquid_mwe + remesh.liquid_mwe;

   % Net solid: the mass from the inserted bottom layer minus the discarded
   % merged layer. remesh_solid = cloned_bottom_solid - merge_export_solid.
   budget.mass_budget_cloned_bottom_solid_mwe = ...
      budget.mass_budget_cloned_bottom_solid_mwe ...
      + remesh.cloned_bottom_solid_mwe;
   budget.mass_budget_merge_export_solid_mwe = ...
      budget.mass_budget_merge_export_solid_mwe ...
      + remesh.merge_export_solid_mwe;

   % A top-cell removal lowers the fixed grid by one layer thickness. The top
   % export records half of the merged pair's mass, which the merge removes from
   % column storage. The top deletion height records the surface displacement.
   budget.mass_budget_top_deletion_count = ...
      budget.mass_budget_top_deletion_count + remesh.top_deletion_count;
   budget.mass_budget_top_deletion_height_m = ...
      budget.mass_budget_top_deletion_height_m + remesh.top_deletion_height_m;
   budget.mass_budget_top_export_solid_mwe = ...
      budget.mass_budget_top_export_solid_mwe + remesh.top_export_solid_mwe;
   budget.mass_budget_top_export_liquid_mwe = ...
      budget.mass_budget_top_export_liquid_mwe + remesh.top_export_liquid_mwe;
   budget.mass_budget_interior_merge_count = ...
      budget.mass_budget_interior_merge_count + remesh.interior_merge_count;
end
