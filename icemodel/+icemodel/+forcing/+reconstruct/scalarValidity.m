function valid = scalarValidity(channel, values)
   %SCALARVALIDITY Check finite samples against the A15 scalar registry.
   %
   %  valid = icemodel.forcing.reconstruct.scalarValidity(channel, values)
   %
   % This check excludes relational checks such as swd versus TOA and swu
   % versus swd. Readiness grades completeness and scalar bounds.
   % physicalValidity applies the stricter fill-candidate rules.

   arguments
      channel (1, 1) string
      values (:, 1) double
   end

   % Bound values live only in physicalBounds, so runtime and producer
   % readiness always use the same policy parameters.
   bounds = icemodel.forcing.reconstruct.physicalBounds(channel);
   valid = isfinite(values) & values >= bounds(1) & values <= bounds(2);
end
