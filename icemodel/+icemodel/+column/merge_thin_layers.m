function [T_ice, f_ice, f_liq, Sc, Sp, d_lyr, budget] = ...
      merge_thin_layers(T_ice, f_ice, f_liq, Sc, Sp, dz_therm, d_pevp, ...
      d_lyr, f_ice_min, budget)
   %MERGE_THIN_LAYERS Merge layers that fall below the minimum ice fraction.
   %
   % merge_thin_layers combines control volumes with ice fractions below the
   % allowable minimum. It also combines the top layer when the current
   % substep's sublimation would trigger a merge on the next substep.
   %
   % Inputs
   %   T_ice      - Column temperature vector.
   %   f_ice      - Column solid ice fraction vector.
   %   f_liq      - Column liquid water fraction vector.
   %   Sc, Sp     - Column shortwave source-term linearization vectors.
   %   dz_therm   - Thermal control-volume thickness [m].
   %   d_pevp     - Potential surface vapor-driven liquid-fraction
   %                change for the top layer [-].
   %   d_lyr      - Accumulated merge-export diagnostic, water-equivalent
   %                fraction per cell; scale by dz for metres w.e.
   %   f_ice_min  - Minimum allowed surface ice fraction [-].
   %   budget     - Forcing-step budget from
   %                icemodel.column.initialize_budget_state.
   %
   % Outputs
   %   T_ice      - Updated column temperature vector after remeshing.
   %   f_ice      - Updated column solid ice fraction vector after remeshing.
   %   f_liq      - Updated column liquid water fraction vector after remeshing.
   %   Sc, Sp     - Updated source-term vectors after remeshing.
   %   d_lyr      - Updated merge-export diagnostic. Each merge adds the
   %                mass it removed from the column, because combining two
   %                cells into one retains their mean.
   %   budget     - Budget with the remesh channels accumulated. Solid and
   %                liquid terms are metres water equivalent. top_deletion_count
   %                counts removals of the surface cell. Only these removals
   %                translate the fixed grid downward. interior_merge_count
   %                counts every other removal. top_export_* is the mass that a
   %                surface removal exported. A merge keeps the mean of the
   %                pair, so top_export_* over-counts what the removed cell
   %                held, and it is not a surface mass flux. The quantized
   %                top_deletion_height_m is grid geometry, not mass. Domain
   %                exchange closes as remesh_solid = cloned_bottom_solid -
   %                merge_export_solid.
   %
   % See also: icemodel.column.merge_layer_indices,
   %  icemodel.column.merge_layers, icemodel.column.accumulate_remesh_budget,
   %  icemodel.column.budget_surface_mass_balance,
   %  icemodel.column.potential_sublimation,
   %  icemodel.surface.apply_surface_vapor_exchange
   %
   %#codegen

   % Flag layers that are already below the minimum ice fraction.
   merge_mask = f_ice <= f_ice_min;

   % Use this substep's vapor exchange to predict a top-cell removal on the next
   % substep. Interior cell merges are not predicted here.
   merge_mask(1) = merge_mask(1) || (f_ice(1) ...
      + icemodel.column.potential_sublimation(d_pevp)) <= f_ice_min;

   % Return early if no cell requires merging.
   if ~any(merge_mask)
      return
   end

   % Initialize the remesh increments to zero.
   remesh = icemodel.column.initialize_remesh_ledger();

   % merge_mask records which layers this substep combined on the original
   % column index. The loop updates do_merge so it tracks the remeshed index.
   do_merge = merge_mask;

   % ji tracks index drift relative to the original column while layers are
   % removed and replacement layers are appended at the bottom.
   ii = 0;
   for j = 1:numel(f_ice)
      ji = j + ii;

      if ~do_merge(ji)
         continue
      end

      % Snapshot the pre-merge column-integrated state.
      [solid_before, liquid_before] = ...
         icemodel.column.integrate_column_budget(T_ice, f_ice, f_liq, dz_therm);

      % Flag merge layer indices.
      [j1, j2] = icemodel.column.merge_layer_indices(ji, f_ice);

      % Merge the flagged layer with its nearest eligible neighbor.
      [T_ice(j2), f_ice(j2), f_liq(j2), Sc(j2), Sp(j2), d_lyr] = ...
         icemodel.column.merge_layers(T_ice, f_ice, f_liq, Sc, Sp, j1, j2, ...
         d_lyr, dz_therm);

      % Count top layer removals and interior merges.
      if j1 == 1
         remesh.top_deletion_count = remesh.top_deletion_count + 1;
      else
         remesh.interior_merge_count = remesh.interior_merge_count + 1;
      end

      % Remove the merged layer and insert a new bottom layer.
      T_ice = dropAndCloneBottom(T_ice, j1);
      Sc = dropAndCloneBottom(Sc, j1);
      Sp = dropAndCloneBottom(Sp, j1);
      f_ice = dropAndCloneBottom(f_ice, j1);
      f_liq = dropAndCloneBottom(f_liq, j1);
      do_merge = dropAndCloneBottom(do_merge, j1);

      % Snapshot the new bottom layer.
      [bottom_solid, bottom_liquid] = ...
         icemodel.column.integrate_column_budget( ...
         T_ice(end), f_ice(end), f_liq(end), dz_therm);

      % Close remesh_solid = cloned_bottom_solid - merge_export_solid for this
      % merge event before accumulating.
      [solid_after, liquid_after] = ...
         icemodel.column.integrate_column_budget(T_ice, f_ice, f_liq, dz_therm);
      remesh = accumulateMergeEvent(remesh, j1 == 1, ...
         solid_after - solid_before, liquid_after - liquid_before, ...
         bottom_solid, bottom_liquid);

      ii = ii - 1;
   end

   % Convert top-layer removal to grid translation ("quantized ablation").
   remesh.top_deletion_height_m = remesh.top_deletion_count * dz_therm;

   % Add this substep's events to the forcing-step budget channels.
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
   %#codegen

   kept = [values(1:j1 - 1); values(j1 + 1:end)];
   values(:) = [kept; kept(end)];
end
