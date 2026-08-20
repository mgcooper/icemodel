function [T, f_ice, f_liq, Sc, Sp, d_lyr, budget] = ...
      merge_thin_layers(T, f_ice, f_liq, Sc, Sp, dz_therm, d_pevp, ...
      d_lyr, f_ice_min, budget)
   %MERGE_THIN_LAYERS Merge layers that fall below the retained ice floor.
   %
   % merge_thin_layers combines control volumes whose ice fraction is already
   % below the allowable minimum. It also combines the top layer when one more
   % substep of the current surface vapor exchange would take that layer below
   % the same minimum. The caller applies that exchange before this call, so
   % the second test is a one-substep look-ahead, not a second application.
   %
   % Inputs
   %   T, f_ice, f_liq - Column thermodynamic state.
   %   Sc, Sp          - Column shortwave source-term linearization vectors.
   %   dz_therm        - Thermal control-volume thickness [m].
   %   d_pevp          - Potential surface vapor-driven liquid-fraction
   %                     change for the top layer [-].
   %   d_lyr           - Accumulated merge-export diagnostic, water-equivalent
   %                     fraction per cell; scale by dz for metres w.e.
   %   f_ice_min       - Minimum allowed surface ice fraction [-].
   %   budget          - Forcing-step budget from
   %                     icemodel.column.initialize_budget_state.
   %
   % Outputs
   %   T, f_ice, f_liq - Updated column state after any merges.
   %   Sc, Sp          - Updated source-term vectors after remeshing.
   %   d_lyr           - Updated merge-export diagnostic. Each merge adds the
   %                     mass it removed from the column, because combining two
   %                     cells into one retains their mean.
   %   budget          - Budget with the remesh channels accumulated. Solid
   %                     and liquid terms are metres water equivalent.
   %                     top_deletion_count counts removals of the surface
   %                     cell. Only these removals translate the fixed grid
   %                     downward. interior_merge_count counts every other
   %                     removal. top_export_* is the mass that a surface
   %                     removal exported. A merge keeps the mean of the pair,
   %                     so top_export_* over-counts what the removed cell
   %                     held, and it is not a surface mass flux. The
   %                     quantized top_deletion_height_m is grid geometry, not
   %                     mass. Domain exchange closes as
   %                     remesh_solid = cloned_bottom_solid - merge_export_solid.
   %
   % See also: icemodel.column.merge_layer_indices,
   %  icemodel.column.merge_layers,
   %  icemodel.column.accumulate_remesh_budget,
   %  icemodel.column.budget_surface_mass_balance,
   %  icemodel.column.potential_sublimation,
   %  icemodel.surface.apply_surface_vapor_exchange
   %
   %#codegen

   % Flag layers that already violate the allowed minimum ice fraction.
   merge_mask = f_ice <= f_ice_min;

   % Look one substep ahead: merge a layer the same vapor exchange would take
   % below f_ice_min next substep. The caller applies d_pevp first, so f_ice
   % holds the current exchange and its tendency estimates the next one. This
   % merges a spent layer at the end of this substep rather than the start of
   % the next.
   %
   % Index 1 only, because budget_surface_mass_balance applies d_pevp at
   % index 1. Broadcasting the scalar down the column would merge interior
   % cells for a mass change they never receive. This restriction changes
   % default-mode remeshing by design (bead icemodel-bhk.1, DesignSpec
   % decision 9), and test_mass_budget_bookkeeping holds a parity oracle for
   % the broadcast form. Do not broadcast it to make a default-mode diff go
   % away.
   %
   % Two biases. Charging the whole tendency to ice merges a wet layer one
   % substep early. Coupled interior transport gets no look-ahead, so a layer
   % its next substep would spend merges one substep late; bead icemodel-eea
   % measures whether that delay matters.
   merge_mask(1) = merge_mask(1) || (f_ice(1) ...
      + icemodel.column.potential_sublimation(d_pevp)) <= f_ice_min;

   % A substep with no merge leaves the budget unchanged: every remesh
   % channel is additive and the event ledger is zero.
   if ~any(merge_mask)
      return
   end

   % Start from a zeroed event ledger so the budget accumulates only
   % the events this substep realizes; callers never infer events from d_lyr.
   remesh = icemodel.column.initialize_remesh_ledger();

   % merge_mask records which layers this substep combines. The loop updates
   % do_merge so that it tracks the combined layers.
   do_merge = merge_mask;

   % ji tracks index drift relative to the original column while layers are
   % removed and replacement layers are appended at the bottom.
   ii = 0;
   for j = 1:numel(f_ice)
      ji = j + ii;

      if ~do_merge(ji)
         continue
      end

      % Take a separate snapshot for each event so the export closure uses
      % that event's own boundary state.
      [solid_before, liquid_before] = ...
         icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz_therm);

      % Flag merge layer indices.
      [j1, j2] = icemodel.column.merge_layer_indices(ji, f_ice);

      % Merge the flagged layer with its nearest eligible neighbor.
      [T(j2), f_ice(j2), f_liq(j2), Sc(j2), Sp(j2), d_lyr] = ...
         icemodel.column.merge_layers(T, f_ice, f_liq, Sc, Sp, j1, j2, ...
         d_lyr, dz_therm);

      % Count the boundary actually removed, not its pre-loop eligibility.
      if j1 == 1
         remesh.top_deletion_count = remesh.top_deletion_count + 1;
      else
         remesh.interior_merge_count = remesh.interior_merge_count + 1;
      end

      % Remove the merged layer, then preserve column length by repeating the
      % deepest remaining state.
      T = dropAndCloneBottom(T, j1);
      Sc = dropAndCloneBottom(Sc, j1);
      Sp = dropAndCloneBottom(Sp, j1);
      f_ice = dropAndCloneBottom(f_ice, j1);
      f_liq = dropAndCloneBottom(f_liq, j1);
      do_merge = dropAndCloneBottom(do_merge, j1);

      % Snapshot the bottom reservoir the clone imported.
      [bottom_solid, bottom_liquid] = ...
         icemodel.column.integrate_column_budget( ...
         T(end), f_ice(end), f_liq(end), dz_therm);

      % Close remesh_solid = cloned_bottom_solid - merge_export_solid for
      % this event before accumulating, so the ledger needs no
      % variable-length per-event arrays.
      [solid_after, liquid_after] = ...
         icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz_therm);
      remesh = accumulateMergeEvent(remesh, j1 == 1, ...
         solid_after - solid_before, liquid_after - liquid_before, ...
         bottom_solid, bottom_liquid);

      ii = ii - 1;
   end

   % Convert actual top-event count to current uniform-grid translation.
   remesh.top_deletion_height_m = remesh.top_deletion_count * dz_therm;

   % Fold this substep's events into the forcing-step budget channels.
   budget = icemodel.column.accumulate_remesh_budget(budget, remesh);
end


function remesh = accumulateMergeEvent(remesh, is_top_removal, ...
      event_solid, event_liquid, bottom_solid, bottom_liquid)
   %ACCUMULATEMERGEEVENT Add one completed merge to the remesh event ledger.
   %
   % IS_TOP_REMOVAL marks a surface removal, the only kind that translates the
   % fixed grid downward. EVENT_* are the signed column-storage changes across
   % the merge and BOTTOM_* are the reservoir the clone imported, so the
   % export the merge discarded closes as BOTTOM - EVENT for each phase.
   %
   %#codegen

   export_solid = bottom_solid - event_solid;
   export_liquid = bottom_liquid - event_liquid;

   % Record the mass a top removal exported. A merge keeps the mean of the
   % pair, so this over-counts what the removed cell held. It is kept separate
   % from interior merges, which move mass without lowering the grid.
   if is_top_removal
      remesh.top_export_solid_mwe = ...
         remesh.top_export_solid_mwe + export_solid;
      remesh.top_export_liquid_mwe = ...
         remesh.top_export_liquid_mwe + export_liquid;
   end

   remesh.solid_mwe = remesh.solid_mwe + event_solid;
   remesh.liquid_mwe = remesh.liquid_mwe + event_liquid;
   remesh.cloned_bottom_solid_mwe = ...
      remesh.cloned_bottom_solid_mwe + bottom_solid;
   remesh.merge_export_solid_mwe = ...
      remesh.merge_export_solid_mwe + export_solid;
end

function values = dropAndCloneBottom(values, j1)
   %DROPANDCLONEBOTTOM Remove cell j1 and repeat the deepest remaining cell.
   %
   % Assigning through values(:) keeps the column length fixed, so the arrays
   % never grow and no AGROW suppression is needed.
   %
   %#codegen

   kept = [values(1:j1 - 1); values(j1 + 1:end)];
   values(:) = [kept; kept(end)];
end
