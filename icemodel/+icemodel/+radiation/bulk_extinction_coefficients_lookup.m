function bulkcoefs = bulk_extinction_coefficients_lookup(ro_sno, lookup)
   %bulk_extinction_coefficients_lookup Bulk extinction from a quantized lookup.
   %
   % bulkcoefs = icemodel.radiation.bulk_extinction_coefficients_lookup( ...
   %    ro_sno, lookup)
   %
   %#codegen

   % Clamp densities to the lookup support before quantizing them.
   ro_sno = min(max(ro_sno(:), lookup.ro_min), lookup.ro_max);

   % Use direct integer-style indexing when the lookup grid is uniform, and
   % fall back to nearest-neighbor interpolation for irregular grids.
   if lookup.is_uniform
      idx = round((ro_sno - lookup.ro_min) / lookup.ro_step) + 1;
   else
      idx = interp1(double(lookup.ro_grid), ...
         transpose(1:numel(lookup.ro_grid)), double(ro_sno), ...
         'nearest', 'extrap');
      idx = round(idx);
   end
   idx = max(1, min(numel(lookup.ro_grid), idx));

   % Reuse the last spectral-layer density for the two padded boundaries.
   idx = vertcat(idx, idx(end), idx(end));
   row_idx = transpose(1:size(lookup.bulkcoefs, 1));
   bulkcoefs = lookup.bulkcoefs(sub2ind(size(lookup.bulkcoefs), row_idx, idx));
end
