function bounds = physicalBounds(channel)
   %PHYSICALBOUNDS Return the approved physical bounds for one channel.
   %
   %  bounds = icemodel.forcing.reconstruct.physicalBounds("tair")
   %
   % Role
   %  Single source of the approved post-fill physical bounds (POLICY
   %  A15). The harness counts violations as
   %  hard method failures and the engine enforces the same limits, so the
   %  registry lives once here. Bounds are inclusive [lower, upper] in the
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

   % Keep this explicit so a new channel is a deliberate policy change.
   switch channel
      case "tair"
         bounds = [193, 300];         % K
      case "rh"
         bounds = [5, 100];           % percent
      case "wspd"
         bounds = [0, 60];            % m/s
      case "psfc"
         bounds = [60000, 108000];    % Pa
      case "swd"
         bounds = [0, Inf];           % W/m2; TOA-relative check is relational
      case "swu"
         bounds = [0, Inf];           % W/m2; swu <= swd check is relational
      case "albedo"
         bounds = [0.05, 0.98];       % fraction
      case "lwd"
         % Extreme-cold clear skies genuinely emit below the old 100
         % floor: MAR produces 69-90 W/m2 over the interior, and a
         % blackbody at the tair floor (193 K) emits ~79 W/m2, so the
         % floor sits safely below real physics while still rejecting
         % garbage (POLICY A15/D-25). The ceiling moved 400 -> 470
         % because warm-fjord stations observe genuine 406-451 W/m2
         % under mild overcast; observations are never clamped or
         % censored, so the ceiling exists to reject garbage and must
         % sit above real physics (POLICY D-26).
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
