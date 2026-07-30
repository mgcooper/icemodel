function [rainf, snowf] = partitionPrecipitation(ppt, tair, kwargs)
   %PARTITIONPRECIPITATION Split total precipitation by air temperature.
   %
   %  [rainf, snowf] = ...
   %     icemodel.forcing.reconstruct.partitionPrecipitation(ppt, tair)
   %
   % Total precipitation at or above the configured phase transition is
   % rain; colder precipitation is snow. Invalid or negative totals and
   % samples without air temperature remain unresolved rather than guessed.
   %
   % This is the RUNTIME phase kernel (POLICY A10/D-18): reconstruction
   % never partitions, and icemodel.resolvePrecipPhase calls this only
   % when the runtime phase-source option is 'threshold'.

   arguments
      ppt (:, 1) double
      tair (:, 1) double
      kwargs.transition_temperature_k (1, 1) double = ...
         icemodel.forcing.reconstruct.setopts(). ...
         rain_snow_transition_temperature_k
   end

   if numel(ppt) ~= numel(tair)
      error('icemodel:reconstruct:partitionPrecipitation:sizeMismatch', ...
         'ppt and tair must share one sample axis');
   end

   % Preserve missingness at invalid inputs, then assign every valid total
   % to exactly one phase so rainf + snowf equals ppt by construction.
   rainf = nan(size(ppt));
   snowf = nan(size(ppt));
   valid = isfinite(ppt) & ppt >= 0 & isfinite(tair);
   rain = valid & tair >= kwargs.transition_temperature_k;
   snow = valid & ~rain;
   rainf(valid) = 0;
   snowf(valid) = 0;
   rainf(rain) = ppt(rain);
   snowf(snow) = ppt(snow);
end
