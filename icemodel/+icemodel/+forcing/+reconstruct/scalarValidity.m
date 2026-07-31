function valid = scalarValidity(channel, values)
   %SCALARVALIDITY Check finite samples against the A15 scalar registry.
   %
   %  valid = icemodel.forcing.reconstruct.scalarValidity(channel, values)
   %
   % This deliberately excludes relational checks such as swd versus TOA
   % and swu versus swd. Readiness grades completeness plus scalar bounds,
   % while physicalValidity owns the stricter fill-candidate rules.

   arguments
      channel (1, 1) string
      values (:, 1) double
   end

   % Bound values live only in physicalBounds so runtime and producer
   % readiness cannot drift when a policy parameter changes.
   bounds = icemodel.forcing.reconstruct.physicalBounds(channel);
   valid = isfinite(values) & values >= bounds(1) & values <= bounds(2);
end
