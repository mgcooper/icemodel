function suffix = metTimestepSuffix(dt)
   %METTIMESTEPSUFFIX Return the canonical file tag for a model-met cadence.
   %
   %  suffix = icemodel.forcing.helpers.metTimestepSuffix(dt)
   %
   % Accepts the proven model-met cadences in seconds (900, 1800, 3600) or
   % their canonical tags ("15m", "30m", "1hr"). Centralizing this closed
   % registry keeps writer naming and runtime discovery in lockstep.

   % Preserve canonical text inputs while rejecting aliases that would create a
   % file name the runtime does not recognize.
   if ischar(dt) || isstring(dt)
      suffix = string(dt);
      if ~isscalar(suffix) || ~ismember(suffix, ["15m", "30m", "1hr"])
         unsupportedTimestep()
      end
      return
   end

   % Map only the source cadences proven by staged artifacts.
   if ~isnumeric(dt) || ~isscalar(dt) || ~isfinite(dt)
      unsupportedTimestep()
   end
   switch double(dt)
      case 900
         suffix = "15m";
      case 1800
         suffix = "30m";
      case 3600
         suffix = "1hr";
      otherwise
         unsupportedTimestep()
   end
end

function unsupportedTimestep()
   %UNSUPPORTEDTIMESTEP Raise the shared closed-registry validation error.
   error('icemodel:forcing:metTimestepSuffix:unsupportedTimestep', ...
      'model-met file cadence must be 900, 1800, or 3600 seconds')
end
