function valid = promiceAlbedoSourceValid(albedo)
   %PROMICEALBEDOSOURCEVALID Identify finite native PROMICE albedo samples.
   %
   % Native ingestion accepts the closed physical source interval [0, 1].
   % Reconstruction uses narrower post-fill bounds, so source provenance
   % must use this predicate rather than the reconstruction validator.

   arguments
      albedo (:, 1) double
   end

   valid = isfinite(albedo) & albedo >= 0 & albedo <= 1;
end
