function Q = SPECTRALNETFLUX(I0, albedo, up, dn, dz_spect)
   %SPECTRALNETFLUX Reconstruct the net spectral flux profile.
   %
   %  Q = SPECTRALNETFLUX(I0, albedo, up, dn, dz_spect)
   %
   % Q is the interface net spectral flux profile. In Schlatter's notation
   % this quantity is XYnet, formed from the staggered up/down fluxes and then
   % differenced to recover absorbed energy in the spectral control volumes.
   % The returned values correspond to the top surface of the top spectral
   % control volume and the lower interfaces beneath it.
   %
   %#codegen

   % Preserve the established top-boundary absorbed-flux correction.
   M = numel(dz_spect) + 2;
   Q = [min(-I0 * (1 - albedo), up(1) - dn(1)); ...
      (up(2:M - 1) + up(3:M)) / 2.0 - (dn(2:M - 1) + dn(3:M)) / 2.0];
end
