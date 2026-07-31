function scale = stepScale(times, x)
   %STEPSCALE Per-season median absolute step of one observed channel.
   %
   %  scale = icemodel.forcing.reconstruct.stepScale(times, x)
   %
   % Role
   %  Single source of the boundary-jump scale (POLICY B6):
   %  the median absolute consecutive change of the observed channel at
   %  this station and season, with an all-season fallback where a season
   %  has too few observed steps for a stable median. Shared by the
   %  validation metrics (scoring jumps) and the engine tiers (rejecting
   %  fills that would create them) so the two can never disagree.
   %
   % Returns
   %  scale : struct with fields DJF, MAM, JJA, SON (finite nonnegative
   %     doubles; zero for a constant channel or when no finite adjacent
   %     step exists, so boundary checks fail closed rather than comparing
   %     against NaN).
   %
   % See also: icemodel.forcing.reconstruct.validationMetrics,
   %  icemodel.forcing.reconstruct.fillShortGaps

   arguments
      times datetime
      x (:, 1) double
   end

   % Steps between ADJACENT finite samples only: a difference spanning a
   % gap is not an hourly change and would inflate the scale. Zero observed
   % steps remain in the policy median; a flat/quantized channel must tighten
   % rather than inflate its seam threshold.
   adjacent = isfinite(x(1:end - 1)) & isfinite(x(2:end));
   step_from = find(adjacent);
   steps = abs(x(step_from + 1) - x(step_from));
   step_season = icemodel.forcing.reconstruct.seasonOf(times(step_from));
   overall = 0;
   if ~isempty(steps)
      overall = median(steps, 'omitnan');
   end

   scale = struct();
   for name = ["DJF", "MAM", "JJA", "SON"]
      in_season = steps(step_season == name);
      scale.(char(name)) = overall;
      if numel(in_season) >= 30
         scale.(char(name)) = median(in_season, 'omitnan');
      end
   end
end
