function [Sc, chi] = UPDATEEXTCOEFSLOOKUP(Qsi, albedo, Q0, dz_spect, Sc, ...
      dz_therm, ro_sno, grid_thermal, grid_spectral, z_walls, lookup)
   %UPDATEEXTCOEFSLOOKUP Update spectral source terms via lookup bulkcoefs.
   %
   %  [Sc, chi] = UPDATEEXTCOEFSLOOKUP(Qsi, albedo, Q0, dz_spect, Sc, ...
   %     dz_therm, ro_sno, grid_thermal, grid_spectral, z_walls, lookup)
   %
   % This exploratory variant keeps the existing two-stream solve but swaps
   % the exact bulk-extinction transform for a precomputed density lookup.
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

   % Pull approximate bulk-extinction coefficients from the lookup.
   bulkcoefs = BULKEXTCOEFSLOOKUP(ro_sno, lookup);
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
