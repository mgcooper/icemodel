function filled = persistenceEstimate(values, times, draw, draws)
   %PERSISTENCEESTIMATE Hold pre-gap values without held-out-data leakage.
   %
   %  filled = icemodel.forcing.reconstruct.persistenceEstimate( ...
   %     values, times, draw, draws)
   %
   % Every synthetic draw is evaluation data during the validation pass.
   % Mask their union before locating the active draw's pre-gap anchors.

   % Build one held-out mask for the validation pass, including draws from
   % other duration/season strata.
   heldout = false(size(values));
   for q = 1:numel(draws)
      if ~isempty(draws{q})
         heldout = heldout | draws{q}.mask;
      end
   end
   values(heldout) = NaN;

   % Hold the last remaining pre-gap observation across each active gap.
   estimate = nan(size(values));
   for g = 1:height(draw.gaps)
      first = find(times == draw.gaps.start_time(g), 1);
      last = find(times == draw.gaps.end_time(g), 1);
      if first > 1 && isfinite(values(first - 1))
         estimate(first:last) = values(first - 1);
      end
   end
   filled = estimate(draw.mask);
end
