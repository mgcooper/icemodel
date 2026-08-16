function [T, f_ice, f_liq, Sc, Sp, d_lyr, merge_mask, remesh] = ...
      merge_thin_layers(T, f_ice, f_liq, Sc, Sp, dz_therm, d_pevp, ...
      d_lyr, f_ice_min)
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
   %
   % Outputs
   %   T, f_ice, f_liq - Updated column state after any merges.
   %   Sc, Sp          - Updated source-term vectors after remeshing.
   %   d_lyr           - Updated merge-export diagnostic. Each merge adds the
   %                     mass it removed from the column, because combining two
   %                     cells into one retains their mean.
   %   merge_mask      - Logical flag marking merge-eligible layers.
   %   remesh          - Optional ledger of events, signed storage exchange,
   %                     and absolute event gross. Solid and liquid terms are
   %                     metres water equivalent. Enthalpy terms are J m-2.
   %                     top_deletion_count counts removals of the surface
   %                     cell. Only these removals translate the fixed grid
   %                     downward. interior_merge_count counts every other
   %                     removal. top_export_* is the mass that a surface
   %                     removal exported. A merge keeps the mean of the pair,
   %                     so top_export_* over-counts what the removed cell
   %                     held, and it is not a surface mass flux. The
   %                     quantized top_deletion_height_m is grid geometry, not
   %                     mass. Domain exchange closes as
   %                     remesh = cloned_bottom - merge_export.
   %
   % See also: icemodel.column.merge_layer_indices,
   %  icemodel.column.merge_layers,
   %  icemodel.column.budget_surface_mass_balance,
   %  icemodel.column.potential_sublimation,
   %  icemodel.surface.apply_surface_vapor_exchange
   %
   %#codegen

   % The eighth output turns on the diagnostic ledger. A caller that asks for
   % seven outputs gets the state transition, and this function does not
   % allocate structs or integrate column storage on every substep.
   use_remesh_ledger = nargout > 7;

   if use_remesh_ledger
      % Return a complete zero ledger when no merge occurs so callers never
      % infer events from the eligibility mask or from d_lyr.
      remesh = icemodel.column.initialize_remesh_ledger();
   end

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

   if ~any(merge_mask)
      return
   end

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

      % Take a separate snapshot for each event. This keeps exchanges with
      % opposite signs from cancelling before the ledger records their absolute
      % gross.
      if use_remesh_ledger
         [solid_before, liquid_before, enthalpy_before] = ...
            icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz_therm);
      end

      % Flag merge layer indices.
      [j1, j2] = icemodel.column.merge_layer_indices(ji, f_ice);

      % Merge the flagged layer with its nearest eligible neighbor.
      [T(j2), f_ice(j2), f_liq(j2), Sc(j2), Sp(j2), d_lyr] = ...
         icemodel.column.merge_layers(T, f_ice, f_liq, Sc, Sp, j1, j2, ...
         d_lyr, dz_therm);

      % Count the boundary actually removed, not its pre-loop eligibility.
      if use_remesh_ledger
         if j1 == 1
            remesh.top_deletion_count = remesh.top_deletion_count + 1;
         else
            remesh.interior_merge_count = remesh.interior_merge_count + 1;
         end
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
      if use_remesh_ledger
         [bottom_solid, bottom_liquid, bottom_enthalpy] = ...
            icemodel.column.integrate_column_budget( ...
            T(end), f_ice(end), f_liq(end), dz_therm);
      end

      % Close remesh = cloned_bottom - merge_export for this event before
      % accumulating signed and absolute channels, so the ledger needs no
      % variable-length per-event arrays.
      if use_remesh_ledger
         [solid_after, liquid_after, enthalpy_after] = ...
            icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz_therm);
         remesh = accumulateMergeEvent(remesh, j1 == 1, ...
            solid_after - solid_before, liquid_after - liquid_before, ...
            enthalpy_after - enthalpy_before, ...
            bottom_solid, bottom_liquid, bottom_enthalpy);
      end

      ii = ii - 1;
   end

   % Convert actual top-event count to current uniform-grid translation.
   if use_remesh_ledger
      remesh.top_deletion_height_m = remesh.top_deletion_count * dz_therm;
   end
end


function remesh = accumulateMergeEvent(remesh, is_top_removal, ...
      event_solid, event_liquid, event_enthalpy, ...
      bottom_solid, bottom_liquid, bottom_enthalpy)
   %ACCUMULATEMERGEEVENT Add one completed merge to the remesh event ledger.
   %
   % IS_TOP_REMOVAL marks a surface removal, the only kind that translates the
   % fixed grid downward. EVENT_* are the signed column-storage changes across
   % the merge and BOTTOM_* are the reservoir the clone imported, so the export
   % that the merge discarded closes as BOTTOM - EVENT for each quantity.

   export_solid = bottom_solid - event_solid;
   export_liquid = bottom_liquid - event_liquid;
   export_enthalpy = bottom_enthalpy - event_enthalpy;

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
   remesh.enthalpy_j_m2 = remesh.enthalpy_j_m2 + event_enthalpy;
   remesh.cloned_bottom_solid_mwe = ...
      remesh.cloned_bottom_solid_mwe + bottom_solid;
   remesh.cloned_bottom_liquid_mwe = ...
      remesh.cloned_bottom_liquid_mwe + bottom_liquid;
   remesh.cloned_bottom_enthalpy_j_m2 = ...
      remesh.cloned_bottom_enthalpy_j_m2 + bottom_enthalpy;
   remesh.merge_export_solid_mwe = ...
      remesh.merge_export_solid_mwe + export_solid;
   remesh.merge_export_liquid_mwe = ...
      remesh.merge_export_liquid_mwe + export_liquid;
   remesh.merge_export_enthalpy_j_m2 = ...
      remesh.merge_export_enthalpy_j_m2 + export_enthalpy;

   % The absolute totals keep merges with opposite signs in one forcing step
   % from cancelling before the model aggregates them in time.
   remesh.solid_gross_mwe = ...
      remesh.solid_gross_mwe + abs(event_solid);
   remesh.liquid_gross_mwe = ...
      remesh.liquid_gross_mwe + abs(event_liquid);
   remesh.enthalpy_gross_j_m2 = ...
      remesh.enthalpy_gross_j_m2 + abs(event_enthalpy);
   remesh.cloned_bottom_solid_gross_mwe = ...
      remesh.cloned_bottom_solid_gross_mwe + abs(bottom_solid);
   remesh.cloned_bottom_liquid_gross_mwe = ...
      remesh.cloned_bottom_liquid_gross_mwe + abs(bottom_liquid);
   remesh.cloned_bottom_enthalpy_gross_j_m2 = ...
      remesh.cloned_bottom_enthalpy_gross_j_m2 + abs(bottom_enthalpy);
   remesh.merge_export_solid_gross_mwe = ...
      remesh.merge_export_solid_gross_mwe + abs(export_solid);
   remesh.merge_export_liquid_gross_mwe = ...
      remesh.merge_export_liquid_gross_mwe + abs(export_liquid);
   remesh.merge_export_enthalpy_gross_j_m2 = ...
      remesh.merge_export_enthalpy_gross_j_m2 + abs(export_enthalpy);
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
