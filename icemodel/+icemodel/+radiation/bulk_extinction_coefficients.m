function k_bulk = bulk_extinction_coefficients(dz_spect, ro_sno, tau_N, ...
      tau_S, solar_dwavel, lookup)
   %bulk_extinction_coefficients Compute bulk extinction coefficients.
   %
   %  k_bulk = icemodel.radiation.bulk_extinction_coefficients( ...
   %     dz_spect, ro_sno, tau_N, tau_S, solar_dwavel)
   %  k_bulk = icemodel.radiation.bulk_extinction_coefficients( ...
   %     [], ro_sno, [], [], [], lookup)
   %
   % Computes bulk (spectrally integrated) extinction coefficients.
   %
   % When LOOKUP is provided and non-empty, the fast density-table path is used
   % (see local subfunction bulk_extinction_coefficients_lookup). When LOOKUP is
   % absent or empty, the exact spectral transform is computed.
   %
   % Inputs correspond to one spectral-grid density profile and the fixed
   % spectral integration coefficients tau_N/S returned by
   % icemodel.radiation.initialize_spectral_model.
   %
   %#codegen

   if nargin > 5 && ~isempty(lookup)
      k_bulk = bulk_extinction_coefficients_lookup(ro_sno, lookup);
      return
   end

   % The transform assumes the spectral grid spacing is uniform, so it uses
   % DZ_SPECT(1) as the representative layer thickness. If DZ_SPECT is
   % nonuniform, this transform needs to be updated.
   ro_sno = max(ro_sno(:), 300.0);

   % Integrate the spectrally weighted extinction exactly on the spectral
   % grid, then pad the lower boundaries for the two-stream solve.
   k_bulk = -log((sum(solar_dwavel .* exp(tau_S .* ro_sno), 2)) ...
      ./ (sum(solar_dwavel .* exp(tau_N .* ro_sno), 2))) / dz_spect(1);

   k_bulk = vertcat(k_bulk, k_bulk(end), k_bulk(end));
end

function bulkcoefs = bulk_extinction_coefficients_lookup(ro_sno, lookup)
   %bulk_extinction_coefficients_lookup Bulk extinction from a quantized lookup.

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
