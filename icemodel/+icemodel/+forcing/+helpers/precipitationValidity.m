function valid = precipitationValidity(ppt, rainf, snowf)
   %PRECIPITATIONVALIDITY Validate complete and partial precipitation splits.
   %
   %  valid = icemodel.forcing.helpers.precipitationValidity( ...
   %     ppt, rainf, snowf)
   %
   % A missing component is allowed. Every finite value must be nonnegative.
   % A finite phase paired with a finite total must allow a nonnegative
   % complement. A complete split must conserve mass.

   arguments
      ppt (:, 1) double
      rainf (:, 1) double
      snowf (:, 1) double
   end

   % The three inputs must hold one value per sample on the same axis.
   if numel(rainf) ~= numel(ppt) || numel(snowf) ~= numel(ppt)
      error('icemodel:forcing:precipitationValidity:sizeMismatch', ...
         'ppt, rainf, and snowf must share one sample axis');
   end

   % Every finite total or component is a nonnegative accumulation rate.
   valid = ~(isfinite(ppt) & ppt < 0) ...
      & ~(isfinite(rainf) & rainf < 0) ...
      & ~(isfinite(snowf) & snowf < 0);

   % A finite phase cannot exceed its finite total. Call the central
   % mass-balance helper, which owns the tolerance.
   rain_pair = isfinite(ppt) & isfinite(rainf);
   if any(rain_pair)
      complement = max(ppt(rain_pair) - rainf(rain_pair), 0);
      valid(rain_pair) = valid(rain_pair) ...
         & icemodel.forcing.helpers.precipitationConsistency( ...
         ppt(rain_pair), rainf(rain_pair), complement);
   end
   snow_pair = isfinite(ppt) & isfinite(snowf);
   if any(snow_pair)
      complement = max(ppt(snow_pair) - snowf(snow_pair), 0);
      valid(snow_pair) = valid(snow_pair) ...
         & icemodel.forcing.helpers.precipitationConsistency( ...
         ppt(snow_pair), complement, snowf(snow_pair));
   end

   % When all three values exist, the two given phases must sum to the total.
   % A nonnegative complement is not enough here.
   complete = isfinite(ppt) & isfinite(rainf) & isfinite(snowf);
   if any(complete)
      valid(complete) = valid(complete) ...
         & icemodel.forcing.helpers.precipitationConsistency( ...
         ppt(complete), rainf(complete), snowf(complete));
   end
end
