function [Sc, chi] = UPDATEEXTCOEFSDECOMPOSEDCACHED(Qsi, albedo, Q0, ...
      dz_spect, spect_N, spect_S, solardwavl, Sc, dz_therm, ro_sno, ...
      grid_thermal, grid_spectral, z_walls)
   %UPDATEEXTCOEFSDECOMPOSEDCACHED Update spectral terms with cached geometry.
   %
   %  [Sc, chi] = UPDATEEXTCOEFSDECOMPOSEDCACHED(Qsi, albedo, Q0, ...
   %     dz_spect, spect_N, spect_S, solardwavl, Sc, dz_therm, ro_sno, ...
   %     grid_thermal, grid_spectral, z_walls)
   %
   % This helper mirrors UPDATEEXTCOEFS while assuming the caller already
   % cached the thermal/spectral grids and spectral control-volume walls.
   %
   %#codegen

   % Preserve the legacy dark/no-sun early exit exactly.
   if Qsi < 1e-3
      Sc = 0.0 * Sc;
      chi = 1.0;
      return
   end

   % Transform the thermal-grid density profile onto the spectral grid.
   ro_sno = max(GRIDFORWARD(ro_sno, grid_thermal, grid_spectral), 300);

   % Build the exact two-stream coefficients through the shared helpers.
   bulkcoefs = BULKEXTCOEFS(dz_spect, ro_sno, spect_N, spect_S, solardwavl);
   [a, r] = GETAANDR(bulkcoefs, albedo);
   [~, up, dn] = SOLVETWOSTREAM(a, r, bulkcoefs, Q0, albedo, z_walls);

   % Reconstruct the same net spectral flux profile used by UPDATEEXTCOEFS.
   M = numel(dz_spect) + 2;
   Q = vertcat(min(-Q0 * (1 - albedo), up(1) - dn(1)), ...
      (up(2:M-1) + up(3:M)) ./ 2.0 - (dn(2:M-1) + dn(3:M)) ./ 2.0);

   % Convert the net flux into the absorbed source term and chi partition.
   [Sc, chi] = SPECTRALSOURCETERM(Qsi, albedo, Q0, Q, dz_spect, dz_therm, ...
      Sc);
end
