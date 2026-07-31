function valid = precipitationValidity(ppt, rainf, snowf)
   %PRECIPITATIONVALIDITY Validate complete and partial precipitation splits.
   %
   %  valid = icemodel.forcing.helpers.precipitationValidity( ...
   %     ppt, rainf, snowf)
   %
   % Missing components are allowed. Every finite value must be
   % nonnegative, each finite phase paired with a finite total must permit
   % a nonnegative complement, and every complete split must conserve mass.

   arguments
      ppt (:, 1) double
      rainf (:, 1) double
      snowf (:, 1) double
   end

   % Refuse axes that cannot represent one samplewise product contract.
   if numel(rainf) ~= numel(ppt) || numel(snowf) ~= numel(ppt)
      error('icemodel:forcing:precipitationValidity:sizeMismatch', ...
         'ppt, rainf, and snowf must share one sample axis');
   end

   % Every finite total or component is a nonnegative accumulation rate.
   valid = ~(isfinite(ppt) & ppt < 0) ...
      & ~(isfinite(rainf) & rainf < 0) ...
      & ~(isfinite(snowf) & snowf < 0);

   % A finite phase cannot exceed its finite total. Reuse the central
   % mass-balance helper so its tolerance remains the single source.
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

   % When all three values exist, require the shipped phases themselves to
   % preserve the total rather than merely admitting hypothetical complements.
   complete = isfinite(ppt) & isfinite(rainf) & isfinite(snowf);
   if any(complete)
      valid(complete) = valid(complete) ...
         & icemodel.forcing.helpers.precipitationConsistency( ...
         ppt(complete), rainf(complete), snowf(complete));
   end
end
