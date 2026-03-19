function [Sc, chi] = UPDATEEXTCOEFSDECOMPOSED(Qsi, albedo, Q0, dz_spect, ...
      spect_N, spect_S, solardwavl, Sc, dz_therm, ro_sno)
   %UPDATEEXTCOEFSDECOMPOSED Update spectral source terms via helper calls.
   %
   %  [Sc, chi] = UPDATEEXTCOEFSDECOMPOSED(Qsi, albedo, Q0, dz_spect, ...
   %     spect_N, spect_S, solardwavl, Sc, dz_therm, ro_sno)
   %
   % This mirrors UPDATEEXTCOEFS while routing the main spectral steps
   % through GRIDFORWARD, BULKEXTCOEFS, GETAANDR, and SOLVETWOSTREAM.
   %
   %#codegen

   % Build the geometry terms that the inlined UPDATEEXTCOEFS body computed
   % implicitly from dz_therm and dz_spect on every call.
   grid_thermal = cumsum(dz_therm) - dz_therm / 2;
   grid_spectral = cumsum(dz_spect) - dz_spect / 2;
   z_walls = [0; cumsum(dz_spect)];

   % Dispatch to the cached-geometry helper so the spectral algebra lives in
   % one canonical implementation.
   [Sc, chi] = UPDATEEXTCOEFSDECOMPOSEDCACHED(Qsi, albedo, Q0, dz_spect, ...
      spect_N, spect_S, solardwavl, Sc, dz_therm, ro_sno, ...
      grid_thermal, grid_spectral, z_walls);
end
