function budget = accumulate_remesh_budget(budget, remesh)
   %ACCUMULATE_REMESH_BUDGET Fold one substep's remesh events into the budget.
   %
   %  budget = icemodel.column.accumulate_remesh_budget(budget, remesh)
   %
   % icemodel.column.merge_thin_layers calls this once per accepted substep
   % after its merge loop. Remeshing is a numerical grid operation, not a
   % physical flux. Its storage exchange stays in dedicated channels, so the
   % scientific comparator can exclude it from physical solid loss.
   %
   % Inputs
   %   budget - Forcing-step budget carrying the phase and vapor channels.
   %   remesh - Event ledger from icemodel.column.initialize_remesh_ledger.
   %
   % Outputs
   %   budget - Budget with the remesh, domain-exchange, and grid-translation
   %            channels accumulated.
   %
   % See also: icemodel.column.merge_thin_layers,
   %  icemodel.column.accumulate_vapor_budget
   %
   %#codegen

   % Net remeshing exchange: what the column gained or lost across the merge.
   budget.mass_budget_remesh_solid_mwe = ...
      budget.mass_budget_remesh_solid_mwe + remesh.solid_mwe;
   budget.mass_budget_remesh_liquid_mwe = ...
      budget.mass_budget_remesh_liquid_mwe + remesh.liquid_mwe;

   % Domain exchange splits the solid net into two parts: the cloned bottom
   % reservoir the fixed-depth grid imports, and the merge export the grid
   % discards. Every event satisfies
   % remesh_solid = cloned_bottom_solid - merge_export_solid.
   budget.mass_budget_cloned_bottom_solid_mwe = ...
      budget.mass_budget_cloned_bottom_solid_mwe ...
      + remesh.cloned_bottom_solid_mwe;
   budget.mass_budget_merge_export_solid_mwe = ...
      budget.mass_budget_merge_export_solid_mwe ...
      + remesh.merge_export_solid_mwe;

   % Grid translation: only top removals lower the surface of the fixed grid.
   % The height is quantized grid geometry. The export is the mass that the
   % removal took out of the column. The export over-counts what the removed
   % cell held, so it is not a surface mass flux.
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
