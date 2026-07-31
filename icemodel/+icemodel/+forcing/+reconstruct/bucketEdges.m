function edges = bucketEdges()
   %BUCKETEDGES Return the canonical gap-duration bucket edges in hours.
   %
   %  edges = icemodel.forcing.reconstruct.bucketEdges()
   %
   % Role
   %  Single source of the policy's gap-duration strata ([0 6 24 72 168
   %  Inf]) shared by the census, the synthetic-missingness sampler, the
   %  orchestrator, and the method planner. gapDurationBucket applies the
   %  policy's right-closed convention, so exact boundaries stay in the
   %  shorter-duration stratum.
   %
   % See also: icemodel.forcing.reconstruct.gapCensus,
   %  icemodel.forcing.reconstruct.syntheticMissingness

   % Bucket 1 is the tier-1 short-gap regime; 2..5 are the donor/proxy
   % strata the planner admits per bucket.
   edges = [0 6 24 72 168 Inf];
end
