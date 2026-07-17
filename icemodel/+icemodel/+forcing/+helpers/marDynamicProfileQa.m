function qa = marDynamicProfileQa(dzsn1, rosn1, shsn3, ro1, outlay, kwargs)
   %MARDYNAMICPROFILEQA Diagnose MAR dynamic snow/firn layer consistency.
   %
   %  qa = icemodel.forcing.helpers.marDynamicProfileQa( ...
   %     dzsn1, rosn1, shsn3, ro1, outlay)
   %
   % DZSN1 and ROSN1 are the native dynamic-layer thickness and density
   % vectors at one MAR grid cell and daily snapshot. SHSN3 is permanent-ice
   % sector 1 total snow/firn thickness. RO1 and OUTLAY are the authoritative
   % public fixed-depth density product used only to quantify reconstruction
   % mismatch.
   %
   % MAR stores the dynamic layers numerically bottom-to-surface. This helper
   % reverses active source indices to reconstruct surface-down layer
   % midpoints, but it never changes RO1 and never declares the dynamic
   % profile publishable. The vertical-order interpretation remains an
   % internal source diagnostic until independently documented.
   %
   % Name-value
   %   thickness_tolerance_m : absolute DZSN1-versus-SHSN3 tolerance (2e-5 m)

   arguments
      dzsn1 double
      rosn1 double
      shsn3 double
      ro1 double
      outlay double
      kwargs.thickness_tolerance_m (1, 1) double {mustBeNonnegative} = 2e-5
   end

   % Initialize every field before early returns so callers can concatenate
   % diagnostics without normalizing status-specific schemas.
   qa = emptyDiagnostic(kwargs.thickness_tolerance_m);
   if any([isempty(dzsn1), isempty(rosn1), isempty(shsn3)])
      qa.detail = "DZSN1, ROSN1, or SHSN3 is unavailable";
      return
   end

   % Column orientation makes source-index reversal and cumulative depths
   % independent of the orientation returned by ncread.
   dzsn1 = double(dzsn1(:));
   rosn1 = double(rosn1(:));
   if numel(dzsn1) ~= numel(rosn1) || ~isscalar(shsn3)
      qa.status = "shape_mismatch";
      qa.detail = "DZSN1/ROSN1 sizes differ or SHSN3 is not scalar";
      return
   end

   % A layer is active only when both native thickness and density are
   % positive. Separate masks retain direct evidence of partial activation.
   qa.available = true;
   thickness_active = isfinite(dzsn1) & dzsn1 > 0;
   density_active = isfinite(rosn1) & rosn1 > 0;
   qa.activity_mask_match = isequal(thickness_active, density_active);
   qa.activity_pair_mismatch_count = nnz(xor( ...
      thickness_active, density_active));

   % Negative, nonfinite, or invalid total-thickness values are source errors,
   % even when the two positive-activity masks happen to agree.
   invalid_values = any(~isfinite(dzsn1) | ~isfinite(rosn1) ...
      | dzsn1 < 0 | rosn1 < 0) || ~isfinite(shsn3) || shsn3 < 0;
   active = thickness_active & density_active;
   qa.active_layer_count = nnz(active);
   qa.active_thickness_sum_m = sum(dzsn1(active));
   qa.shsn3_permanent_ice_m = double(shsn3);
   qa.thickness_residual_m = ...
      qa.active_thickness_sum_m - qa.shsn3_permanent_ice_m;
   qa.thickness_match = abs(qa.thickness_residual_m) ...
      <= qa.thickness_tolerance_m;

   % Reverse active source indices because the native numeric storage order is
   % bottom-to-surface; cumulative thickness then produces positive-down depth.
   source_indices = find(active);
   surface_indices = flipud(source_indices);
   qa.surface_down_source_indices = surface_indices(:)';
   if ~isempty(surface_indices)
      surface_thickness = dzsn1(surface_indices);
      surface_density = rosn1(surface_indices);
      qa.dynamic_depth_upper_m = cumsum(surface_thickness);
      qa.dynamic_depth_lower_m = [0; qa.dynamic_depth_upper_m(1:end-1)];
      qa.dynamic_depth_midpoint_m = ...
         (qa.dynamic_depth_lower_m + qa.dynamic_depth_upper_m) / 2;
      qa.dynamic_density_kg_m3 = surface_density;
      qa.dynamic_density_min_kg_m3 = min(surface_density);
      qa.dynamic_density_max_kg_m3 = max(surface_density);
   end

   % Interpolation is intentionally overlap-only and diagnostic. OUTLAY values
   % outside the dynamic midpoint support remain missing; no extrapolation or
   % RO1 correction is performed.
   ro1 = double(ro1(:));
   outlay = double(outlay(:));
   if numel(qa.dynamic_depth_midpoint_m) >= 2 ...
         && numel(ro1) == numel(outlay)
      reconstructed = interp1(qa.dynamic_depth_midpoint_m, ...
         qa.dynamic_density_kg_m3, outlay, 'linear', NaN);
      compare = isfinite(reconstructed) & isfinite(ro1);
      qa.reconstructed_ro1_kg_m3 = reconstructed;
      qa.reconstruction_overlap_count = nnz(compare);
      if any(compare)
         difference = reconstructed(compare) - ro1(compare);
         qa.reconstruction_bias_kg_m3 = mean(difference);
         qa.reconstruction_rmse_kg_m3 = sqrt(mean(difference .^ 2));
         qa.reconstruction_max_abs_kg_m3 = max(abs(difference));
      end
   end

   % Status precedence retains the most fundamental source inconsistency while
   % leaving all computed metrics available for diagnosis.
   if invalid_values
      qa.status = "invalid_values";
      qa.detail = "dynamic layers or SHSN3 contain invalid values";
   elseif ~qa.activity_mask_match
      qa.status = "activity_mask_mismatch";
      qa.detail = "positive DZSN1 and ROSN1 activity masks differ";
   elseif ~qa.thickness_match
      qa.status = "thickness_mismatch";
      qa.detail = "active DZSN1 thickness does not match SHSN3 sector 1";
   else
      qa.status = "ok";
      qa.detail = "dynamic layers pass activity and thickness diagnostics";
   end
end

function qa = emptyDiagnostic(tolerance)
   %EMPTYDIAGNOSTIC Return the stable diagnostic schema for all outcomes.
   qa = struct( ...
      'status', "not_available", ...
      'detail', "", ...
      'available', false, ...
      'publishable', false, ...
      'correction_applied', false, ...
      'storage_order', "bottom_to_surface; reverse active source indices", ...
      'activity_mask_match', false, ...
      'activity_pair_mismatch_count', 0, ...
      'active_layer_count', 0, ...
      'active_thickness_sum_m', NaN, ...
      'shsn3_permanent_ice_m', NaN, ...
      'thickness_residual_m', NaN, ...
      'thickness_tolerance_m', tolerance, ...
      'thickness_match', false, ...
      'surface_down_source_indices', zeros(1, 0), ...
      'dynamic_depth_lower_m', zeros(0, 1), ...
      'dynamic_depth_midpoint_m', zeros(0, 1), ...
      'dynamic_depth_upper_m', zeros(0, 1), ...
      'dynamic_density_kg_m3', zeros(0, 1), ...
      'dynamic_density_min_kg_m3', NaN, ...
      'dynamic_density_max_kg_m3', NaN, ...
      'reconstructed_ro1_kg_m3', zeros(0, 1), ...
      'reconstruction_overlap_count', 0, ...
      'reconstruction_bias_kg_m3', NaN, ...
      'reconstruction_rmse_kg_m3', NaN, ...
      'reconstruction_max_abs_kg_m3', NaN);
end
