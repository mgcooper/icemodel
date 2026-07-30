function bucket = gapDurationBucket(duration_hours, edges_hours)
   %GAPDURATIONBUCKET Assign positive durations to right-closed policy bins.
   %
   %  bucket = icemodel.forcing.reconstruct.gapDurationBucket(hours)
   %
   % Role
   %  Shared duration-bucket convention for census, held-out sampling, and
   %  composition. Policy boundaries belong to the shorter bucket, so exactly
   %  6 h is in the <=6 h stratum and exactly 24 h is in the 6--24 h stratum.

   arguments
      duration_hours double
      edges_hours (1, :) double = ...
         icemodel.forcing.reconstruct.bucketEdges()
   end

   % MATLAB defaults to left-closed bins; the policy is explicitly
   % right-closed at each finite upper duration boundary.
   bucket = discretize(duration_hours, edges_hours, ...
      'IncludedEdge', 'right');
end
