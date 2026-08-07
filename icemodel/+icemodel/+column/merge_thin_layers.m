function [T, f_ice, f_liq, Sc, Sp, d_lyr, merge_mask, remesh] = ...
      merge_thin_layers( ...
      T, f_ice, f_liq, Sc, Sp, dz_therm, d_pevp, d_lyr, f_ice_min)
   %MERGE_THIN_LAYERS Merge layers that fall below the retained ice floor.
   %
   % merge_thin_layers combines control volumes whose ice fraction is already
   % below the allowable minimum or is predicted to fall below it after the
   % current surface vapor-mass exchange.
   %
   % Inputs
   %   T, f_ice, f_liq - Column thermodynamic state.
   %   Sc, Sp          - Column shortwave source-term linearization vectors.
   %   dz_therm        - Thermal control-volume thickness [m].
   %   d_pevp          - Potential surface vapor-driven liquid-fraction change.
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
   %   remesh          - Optional event, signed storage-exchange, and absolute
   %                     event-throughput ledger. Solid/liquid terms are metres
   %                     water equivalent and enthalpy terms are J m-2.
   %                     top_deletion_count counts removals of the surface cell,
   %                     the only ones that translate the fixed grid downward;
   %                     interior_merge_count counts every other removal.
   %                     top_export_* is the mass a surface removal actually
   %                     exports, which is the physical surface-loss quantity;
   %                     the quantized top_deletion_height_m is grid geometry
   %                     and is never mass. Domain exchange closes as
   %                     remesh = cloned_bottom - collapse_export.
   %
   % See also: icemodel.column.merge_layer_indices,
   %  icemodel.column.merge_layers,
   %  icemodel.column.budget_surface_mass_balance
   %
   %#codegen

   persistent Ls Lv ro_ice ro_liq
   if isempty(Ls)
      [Ls, Lv, ro_ice, ro_liq] = icemodel.physicalConstant( ...
         'Ls', 'Lv', 'ro_ice', 'ro_liq');
   end
   % The eighth output opts into the diagnostic ledger. Existing seven-output
   % solver callers keep the state transition without allocating structs or
   % integrating column storage on every substep.
   use_remesh_ledger = nargout > 7;
   if use_remesh_ledger
      % Return a complete zero ledger when no merge occurs so callers never
      % infer events from the eligibility mask or from d_lyr.
      remesh = zeroRemeshLedger();
   end

   % Flag layers that already violate the allowed minimum ice fraction or would
   % do so after the current surface vapor-driven mass change is applied.
   merge_mask = f_ice <= f_ice_min | ...
      (f_ice + d_pevp * (Lv * ro_liq) / (Ls * ro_ice)) <= f_ice_min;

   if ~any(merge_mask)
      return
   end

   % merge_mask records which layer(s) were combined this substep. do_merge is
   % updated in the loop to stay on track with the updated combined layers.
   do_merge = merge_mask;

   % ji tracks index drift relative to the original column while layers are
   % removed and replacement layers are appended at the bottom.
   ii = 0;
   for j = 1:numel(f_ice)
      ji = j + ii;

      if ~do_merge(ji)
         continue
      end

      % Snapshot each actual event independently so opposite-signed exchanges
      % cannot cancel before their absolute throughput is retained.
      if use_remesh_ledger
         [solid_before, liquid_before, enthalpy_before] = ...
            icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz_therm);
      end

      % Merge the flagged layer with its nearest eligible neighbor.
      [j1, j2] = icemodel.column.merge_layer_indices(ji, f_ice);
      [T(j2), f_ice(j2), f_liq(j2), Sc(j2), Sp(j2), d_lyr] = ...
         icemodel.column.merge_layers(T, f_ice, f_liq, Sc, Sp, j1, j2, ...
         d_lyr, dz_therm);

      % Count the boundary actually removed, not its pre-loop eligibility, and
      % snapshot the post-merge bottom reservoir the clone below imports.
      if use_remesh_ledger
         if j1 == 1
            remesh.top_deletion_count = remesh.top_deletion_count + 1;
         else
            remesh.interior_merge_count = remesh.interior_merge_count + 1;
         end
         [bottom_solid, bottom_liquid, bottom_enthalpy] = ...
            icemodel.column.integrate_column_budget( ...
            T(end), f_ice(end), f_liq(end), dz_therm);
      end

      % Remove the merged layer and preserve column length by repeating the
      % deepest remaining state at the bottom.
      T = vertcat(T, T(end)); T(j1) = []; %#ok<*AGROW>
      Sc = vertcat(Sc, Sc(end)); Sc(j1) = [];
      Sp = vertcat(Sp, Sp(end)); Sp(j1) = [];
      f_ice = vertcat(f_ice, f_ice(end)); f_ice(j1) = [];
      f_liq = vertcat(f_liq, f_liq(end)); f_liq(j1) = [];
      do_merge = vertcat(do_merge, do_merge(end)); do_merge(j1) = [];

      % Close remesh = cloned_bottom - collapse_export for this event before
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

function remesh = zeroRemeshLedger()
   %ZEROREMESHLEDGER Return the fixed all-zero remesh event ledger layout.

   % Declare every field literally so MATLAB Coder never constructs or extends
   % the event ledger from runtime field names.
   remesh = struct( ...
      'top_deletion_count', 0, ...
      'top_deletion_height_m', 0.0, ...
      'top_export_solid_mwe', 0.0, ...
      'top_export_liquid_mwe', 0.0, ...
      'interior_merge_count', 0, ...
      'solid_mwe', 0.0, ...
      'liquid_mwe', 0.0, ...
      'enthalpy_j_m2', 0.0, ...
      'cloned_bottom_solid_mwe', 0.0, ...
      'cloned_bottom_liquid_mwe', 0.0, ...
      'cloned_bottom_enthalpy_j_m2', 0.0, ...
      'collapse_export_solid_mwe', 0.0, ...
      'collapse_export_liquid_mwe', 0.0, ...
      'collapse_export_enthalpy_j_m2', 0.0, ...
      'solid_throughput_mwe', 0.0, ...
      'liquid_throughput_mwe', 0.0, ...
      'enthalpy_throughput_j_m2', 0.0, ...
      'cloned_bottom_solid_throughput_mwe', 0.0, ...
      'cloned_bottom_liquid_throughput_mwe', 0.0, ...
      'cloned_bottom_enthalpy_throughput_j_m2', 0.0, ...
      'collapse_export_solid_throughput_mwe', 0.0, ...
      'collapse_export_liquid_throughput_mwe', 0.0, ...
      'collapse_export_enthalpy_throughput_j_m2', 0.0);
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

   % The surface-loss comparator needs the mass a top removal actually exports,
   % separated from interior merges, which move mass without lowering the grid.
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
   remesh.collapse_export_solid_mwe = ...
      remesh.collapse_export_solid_mwe + export_solid;
   remesh.collapse_export_liquid_mwe = ...
      remesh.collapse_export_liquid_mwe + export_liquid;
   remesh.collapse_export_enthalpy_j_m2 = ...
      remesh.collapse_export_enthalpy_j_m2 + export_enthalpy;

   % Absolute event activity keeps opposite-signed merges within one forcing
   % step from cancelling before temporal aggregation.
   remesh.solid_throughput_mwe = ...
      remesh.solid_throughput_mwe + abs(event_solid);
   remesh.liquid_throughput_mwe = ...
      remesh.liquid_throughput_mwe + abs(event_liquid);
   remesh.enthalpy_throughput_j_m2 = ...
      remesh.enthalpy_throughput_j_m2 + abs(event_enthalpy);
   remesh.cloned_bottom_solid_throughput_mwe = ...
      remesh.cloned_bottom_solid_throughput_mwe + abs(bottom_solid);
   remesh.cloned_bottom_liquid_throughput_mwe = ...
      remesh.cloned_bottom_liquid_throughput_mwe + abs(bottom_liquid);
   remesh.cloned_bottom_enthalpy_throughput_j_m2 = ...
      remesh.cloned_bottom_enthalpy_throughput_j_m2 + abs(bottom_enthalpy);
   remesh.collapse_export_solid_throughput_mwe = ...
      remesh.collapse_export_solid_throughput_mwe + abs(export_solid);
   remesh.collapse_export_liquid_throughput_mwe = ...
      remesh.collapse_export_liquid_throughput_mwe + abs(export_liquid);
   remesh.collapse_export_enthalpy_throughput_j_m2 = ...
      remesh.collapse_export_enthalpy_throughput_j_m2 + abs(export_enthalpy);
end
