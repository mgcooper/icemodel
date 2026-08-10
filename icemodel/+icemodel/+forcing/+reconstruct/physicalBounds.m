function bounds = physicalBounds(channel)
   %PHYSICALBOUNDS Return the approved physical bounds for one channel.
   %
   %  bounds = icemodel.forcing.reconstruct.physicalBounds("tair")
   %
   % Role
   %  Holds the approved post-fill physical bounds (POLICY A15). The
   %  harness counts violations as hard method failures and the engine
   %  enforces the same limits. Bounds are inclusive [lower, upper] in the
   %  canonical met units. Channels with a data-dependent upper limit (swd
   %  vs top-of-atmosphere, swu vs swd) return Inf here; their relational
   %  checks live with the metrics, which have both channels in hand.
   %
   % Returns
   %  bounds : 1x2 double [lower, upper], inclusive.
   %
   % See also: icemodel.forcing.reconstruct.validationMetrics,
   %  icemodel.forcing.reconstruct.admissionGate

   arguments
      channel (1, 1) string
   end

   % Each approved channel is listed explicitly; an unlisted channel errors.
   switch channel
      case "tair"
         bounds = [193, 300];         % K
      case "rh"
         bounds = [5, 100];           % percent
      case "wspd"
         % metchecks reads the wind floor from here to clamp at runtime.
         % One shared floor keeps mean-preserving disaggregation from
         % synthesizing singular calm-air rows from valid hourly postings.
         bounds = [0.1, 60];          % m/s
      case "psfc"
         bounds = [60000, 108000];    % Pa
      case "swd"
         bounds = [0, Inf];           % W/m2; TOA-relative check is relational
      case "swu"
         bounds = [0, Inf];           % W/m2; swu <= swd check is relational
      case "albedo"
         bounds = [0.05, 0.98];       % fraction
      case "lwd"
         % The floor sits below real physics while still rejecting
         % garbage (POLICY A15/D-25): extreme-cold clear skies emit
         % 69-90 W/m2 in MAR over the interior, and a blackbody at the
         % tair floor (193 K) emits ~79 W/m2. The ceiling sits above real
         % physics for the same reason (POLICY D-26): warm-fjord stations
         % observe genuine 406-451 W/m2 under mild overcast, and
         % observations are never clamped or censored.
         bounds = [40, 470];          % W/m2
      case "ppt"
         bounds = [0, Inf];           % accumulation rate is nonnegative
      case "boom_height"
         bounds = [0, Inf];           % m; runtime additionally requires z > z0
      otherwise
         error('icemodel:reconstruct:physicalBounds:unknownChannel', ...
            'no approved physical bounds for channel: %s', channel)
   end
end
