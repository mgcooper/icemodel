function season = seasonOf(times)
   %SEASONOF Meteorological season label of each timestamp.
   %
   %  season = icemodel.forcing.reconstruct.seasonOf(times)
   %
   % Role
   %  Single source of the DJF/MAM/JJA/SON season labels the harness strata
   %  key on, shared by the census, the synthetic-missingness sampler, and
   %  the metrics so the stratum axes can never drift apart.
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
