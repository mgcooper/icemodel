function bulkcoefs = BULKEXTCOEFSLOOKUP(ro_sno, lookup)
   %BULKEXTCOEFSLOOKUP Approximate bulk extinction from a quantized lookup.
   %
   %  bulkcoefs = BULKEXTCOEFSLOOKUP(ro_sno, lookup)
   %
   %#codegen

   % Clamp densities to the lookup support before quantizing them.
   rho = max(ro_sno(:), lookup.rho_min);
   rho = min(rho, lookup.rho_max);

   % Use direct integer-style indexing when the lookup grid is uniform, and
   % fall back to nearest-neighbor interpolation for irregular grids.
   if lookup.is_uniform
      idx = round((rho - lookup.rho_min) ./ lookup.rho_step) + 1;
   else
      idx = interp1(double(lookup.rho_grid), transpose(1:numel(lookup.rho_grid)), ...
         double(rho), 'nearest', 'extrap');
      idx = round(idx);
   end
   idx = max(1, min(numel(lookup.rho_grid), idx));

   % Reuse the last spectral-layer density for the two padded boundaries.
   idx = vertcat(idx, idx(end), idx(end));
   rows = transpose(1:size(lookup.bulkcoefs, 1));
   bulkcoefs = lookup.bulkcoefs(sub2ind(size(lookup.bulkcoefs), rows, idx));
end
