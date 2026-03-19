function [Sc, chi] = UPDATEEXTCOEFS(Qsi, albedo, Q0, dz_spect, spect_N, ...
      spect_S, solardwavl, Sc, dz_therm, ro_sno)
   %UPDATEEXTCOEFS Update the spectral source term and SEB chi partition.
   %
   %  [Sc, chi] = UPDATEEXTCOEFS(Qsi, albedo, Q0, dz_spect, spect_N, ...
   %     spect_S, solardwavl, Sc, dz_therm, ro_sno)
   %
   % This public kernel now routes through the decomposed helper path so the
   % spectral model can share BULKEXTCOEFS, GETAANDR, and SOLVETWOSTREAM.
   %
   %#codegen

   [Sc, chi] = UPDATEEXTCOEFSDECOMPOSED(Qsi, albedo, Q0, dz_spect, ...
      spect_N, spect_S, solardwavl, Sc, dz_therm, ro_sno);
end
