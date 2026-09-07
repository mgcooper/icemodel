function bounds = physicalBounds(channel)
   %PHYSICALBOUNDS Return the approved physical bounds for one channel.
   %
   %  bounds = icemodel.forcing.reconstruct.physicalBounds("tair")
   %
   % Role
   %  Holds the approved post-fill physical bounds (POLICY A15). The
   %  harness counts violations as hard method failures and the engine
   %  enforces the same limits. Bounds are inclusive [lower, upper] in the
   %  met units. A channel with a data-dependent upper limit (swd
   %  against top-of-atmosphere, swu against swd) returns Inf here. The
   %  metrics hold both channels, so they run the relational check.
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
         % The floor is below every real value and still rejects invalid
         % data (POLICY A15/D-25). Extreme-cold clear skies emit
         % 69-90 W/m2 in MAR over the interior, and a blackbody at the
         % tair floor (193 K) emits about 79 W/m2. The ceiling is above
         % every real value for the same reason (POLICY D-26). Warm-fjord
         % stations record true values of 406-451 W/m2 under mild
         % overcast, and the pipeline never clamps or censors an
         % observation.
         bounds = [40, 470];          % W/m2
      case "ppt"
         bounds = [0, Inf];           % accumulation rate is nonnegative
      case "boom_height"
         bounds = [0, Inf];           % m; runtime also requires z > z0
      otherwise
         error('icemodel:reconstruct:physicalBounds:unknownChannel', ...
            'no approved physical bounds for channel: %s', channel)
   end
end
