function bulkcoefs = BULKEXTCOEFS(dz_spect, ro_sno, spect_N, spect_S, solardwavl)
   %BULKEXTCOEFS Compute spectral bulk extinction coefficients.
   %
   %  bulkcoefs = BULKEXTCOEFS(dz_spect, ro_sno, spect_N, spect_S, solardwavl)
   %
   % Inputs correspond to one spectral-grid density profile and the fixed
   % spectral integration coefficients returned by EXTCOEFSINIT.
   %
   %#codegen

   % Match UPDATEEXTCOEFS by capping unrealistically light densities and
   % using the leading spectral grid spacing in the integrated transform.
   ro_sno = max(ro_sno, 300);
   dz_bulk = dz_spect(1);

   % Integrate the spectrally weighted extinction exactly on the spectral
   % grid, then pad the lower boundaries for the two-stream solve.
   bulkcoefs = -log((sum(solardwavl .* exp(spect_S .* ro_sno), 2)) ...
      ./ (sum(solardwavl .* exp(spect_N .* ro_sno), 2))) ./ dz_bulk;
   bulkcoefs = vertcat(bulkcoefs, bulkcoefs(end), bulkcoefs(end));
end
