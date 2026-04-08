function [tau_N, tau_S, k_bulk_lookup, k_ext] = update_extinction_coefficients( ...
      qext, g, coalbedo, kabs, kice, wavel, radii, iradius, z_edges_spect, ...
      dz_spect, solar_dwavel, ro_ice, use_lookup)
   %update_extinction_coefficients Update spectral extinction coefficients for one grain radius.
   %
   % [tau_N, tau_S, k_bulk_lookup, k_ext] = ...
   %    icemodel.radiation.update_extinction_coefficients(qext, g, coalbedo, ...
   %    kabs, kice, wavel, radii, iradius, z_edges_spect, dz_spect, ...
   %    solar_dwavel, ro_ice, use_lookup)
   %
   % This helper computes k_ext for the requested optical grain-radius index,
   % or interpolated index, applies the optional impurity scaling, then
   % precomputes tau_N and tau_S for the exact bulk-extinction coefficient
   % transform. It is the natural entry point if a future grain-growth model
   % needs to refresh k_ext during the timestep loop without reloading the
   % optical tables.
   %
   % When USE_LOOKUP is true, the bulk-extinction lookup table is also built and
   % returned through K_BULK_LOOKUP. Otherwise K_BULK_LOOKUP is an empty struct,
   % which causes icemodel.column.shortwave_source_term to use the exact
   % transform.
   %
   %#codegen

   % Compute the spectral extinction coefficients.
   k_ext = icemodel.radiation.spectral_extinction_coefficients( ...
      qext, g, coalbedo, radii, iradius);
   if xor(isempty(kice), isempty(kabs))
      error('kice and kabs must both be empty or both be provided')
   end
   if ~isempty(kice)
      k_ext = icemodel.radiation.rescale_spectral_extinction_coefficients( ...
         k_ext, kabs, kice, wavel);
   end

   % Compute specific optical depth (precomputed arguments to the exponential
   % function in icemodel.radiation.bulk_extinction_coefficients)
   tau_edges = -(z_edges_spect * k_ext) / ro_ice;
   tau_N = tau_edges(1:end - 1, :);
   tau_S = tau_edges(2:end, :);

   % Build the bulk-extinction lookup table when requested.
   if use_lookup
      k_bulk_lookup = icemodel.makeBulkExtCoefsLookup( ...
         dz_spect, tau_N, tau_S, solar_dwavel);
   else
      k_bulk_lookup = struct([]);
   end
end
