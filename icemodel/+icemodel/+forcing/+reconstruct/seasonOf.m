function season = seasonOf(times)
   %SEASONOF Meteorological season label of each timestamp.
   %
   %  season = icemodel.forcing.reconstruct.seasonOf(times)
   %
   % Role
   %  Returns the DJF/MAM/JJA/SON season labels that the harness strata use.
   %  The census, the synthetic-missingness sampler, and the metrics all call
   %  this function, so they share the same stratum axes.
   %
   % Returns
   %  season : string column, one label per timestamp.
   %
   % See also: icemodel.forcing.reconstruct.gapCensus,
   %  icemodel.forcing.reconstruct.syntheticMissingness

   arguments
      times datetime
   end

   labels = ["DJF", "DJF", "MAM", "MAM", "MAM", "JJA", ...
      "JJA", "JJA", "SON", "SON", "SON", "DJF"];
   season = reshape(labels(month(times)), [], 1);
end
