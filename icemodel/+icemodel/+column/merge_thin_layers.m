function [T, f_ice, f_liq, Sc, Sp, d_lyr, merge_mask] = merge_thin_layers( ...
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
   %   d_lyr           - Accumulated layer-change diagnostic [m].
   %   f_ice_min       - Minimum allowed surface ice fraction [-].
   %
   % Outputs
   %   T, f_ice, f_liq - Updated column state after any merges.
   %   Sc, Sp          - Updated source-term vectors after remeshing.
   %   d_lyr           - Updated cumulative layer-change diagnostic [m].
   %   merge_mask      - Logical flag marking merge-eligible layers.
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

      % Merge the flagged layer with its nearest eligible neighbor.
      [j1, j2] = icemodel.column.merge_layer_indices(ji, f_ice);
      [T(j2), f_ice(j2), f_liq(j2), Sc(j2), Sp(j2), d_lyr] = ...
         icemodel.column.merge_layers(T, f_ice, f_liq, Sc, Sp, j1, j2, ...
         d_lyr, dz_therm);

      % Remove the merged layer and preserve column length by repeating the
      % deepest remaining state at the bottom.
      T = vertcat(T, T(end)); T(j1) = []; %#ok<*AGROW>
      Sc = vertcat(Sc, Sc(end)); Sc(j1) = [];
      Sp = vertcat(Sp, Sp(end)); Sp(j1) = [];
      f_ice = vertcat(f_ice, f_ice(end)); f_ice(j1) = [];
      f_liq = vertcat(f_liq, f_liq(end)); f_liq(j1) = [];
      do_merge = vertcat(do_merge, do_merge(end)); do_merge(j1) = [];

      ii = ii - 1;
   end
end
