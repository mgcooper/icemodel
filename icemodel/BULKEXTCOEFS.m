function k_bulk = BULKEXTCOEFS(dz_spect, ro_sno, tau_N, tau_S, solar_dwavel)
   %BULKEXTCOEFS Compute bulk (spectrally integrated) extinction coefficients.
   %
   %  k_bulk = BULKEXTCOEFS(dz_spect, ro_sno, tau_N, tau_S, solar_dwavel)
   %
   % Inputs correspond to one spectral-grid density profile and the fixed
   % spectral integration coefficients tau_N/S returned by EXTCOEFSINIT.
   %
   %#codegen

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
