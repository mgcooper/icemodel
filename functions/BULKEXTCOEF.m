function bulkcoefs = BULKEXTCOEF(dz_spect,ro_sno,spect_lower,spect_upper,solardwavl)
   %BULKEXTCOEF Compute the spectrally integrated extinction coefficient.
   %
   %  bulkcoefs = BULKEXTCOEF(dz_spect,ro_sno,spect_lower,spect_upper,solardwavl)
   %
   % See also:

   % Compute the downward bulk extinction coefficient, scaled by the total ice
   %   equivalent thickness of each layer
   bulkcoefs = -log((sum(exp(spect_lower.*ro_sno).*solardwavl,2))./    ...
      (sum(exp(spect_upper.*fix(max(ro_sno,300))).*solardwavl,2)))./dz_spect;

   % Cast in a form that can be used in the two-stream computation (add the
   %   boundaries).  Here I have assumed that it is okay to call the
   %   boundaries equal to the value at the center of that grid cell (the
   %   value prevails thoughout the cell).
   bulkcoefs = vertcat(bulkcoefs,bulkcoefs(end),bulkcoefs(end));
end
