function x = elevationAdjust(channel, x, dz, kwargs)
   %ELEVATIONADJUST Adjust a donor channel across an elevation difference.
   %
   %  x = icemodel.forcing.reconstruct.elevationAdjust("tair", x, dz)
   %
   % Role
   %  The policy's donor elevation adjustments (POLICY B4): temperature by a
   %  lapse rate and pressure barometrically, applied to the DONOR series
   %  before transfer fitting whenever |dz| exceeds the adjustment
   %  threshold. Other channels pass through unchanged unless held-out
   %  validation later justifies channel-specific treatment.
   %
   % Inputs
   %  channel : canonical channel name.
   %  x : donor series values.
   %  dz : target elevation minus donor elevation, meters.
   %
   % Name-value
   %  lapse_rate : temperature lapse in K/m (default -0.0060; the policy's
   %     recorded fallback — a fitted overlap lapse may override it).
   %  threshold_m : |dz| below which no adjustment applies (default 100,
   %     a Section-C parameter).
   %  tair_for_pressure : coincident air temperature (K) for the
   %     barometric scale height; defaults to 255 K when omitted, the
   %     ice-sheet annual-mean order the policy uses for the fallback.
   %
   % See also: icemodel.forcing.reconstruct.fitDonorTransfer

   arguments
      channel (1, 1) string
      x (:, 1) double
      dz (1, 1) double
      kwargs.lapse_rate (1, 1) double = ...
         icemodel.forcing.reconstruct.setopts().lapse_rate
      kwargs.threshold_m (1, 1) double {mustBeNonnegative} = ...
         icemodel.forcing.reconstruct.setopts().elevation_threshold_m
      kwargs.tair_for_pressure (:, 1) double = ...
         icemodel.forcing.reconstruct.setopts().tair_for_pressure
   end

   if abs(dz) <= kwargs.threshold_m
      return
   end
   switch channel
      case "tair"
         % Linear lapse toward the target elevation.
         x = x + kwargs.lapse_rate * dz;
      case "psfc"
         % Barometric adjustment with scale height from the coincident (or
         % fallback) temperature: p_target = p_donor * exp(-dz g / (R T)).
         [Rd, g] = icemodel.physicalConstant('Rd', 'gravity');
         T = kwargs.tair_for_pressure;
         x = x .* exp(-dz * g ./ (Rd * T));
      otherwise
         % No approved elevation dependence for the remaining channels.
   end
end
