function [passed, lower_wall_s, upper_wall_s, reason] = ...
      performanceGate(current_wall_s, reference_wall_s, tolerance)
   %PERFORMANCEGATE Compare one runtime to a two-sided accepted band.
   %
   %  [passed, lower_wall_s, upper_wall_s, reason] = ...
   %     icemodel.test.helpers.performanceGate(12, 10, 0.2)
   %
   % A large apparent speedup can mean that work was skipped or that an
   % inflated reference was accepted. Treat both sides as review gates.

   arguments
      current_wall_s (1, 1) double {mustBeFinite, mustBePositive}
      reference_wall_s (1, 1) double {mustBeFinite, mustBePositive}
      tolerance (1, 1) double {mustBeFinite, mustBeNonnegative, ...
         mustBeLessThan(tolerance, 1)}
   end

   % Build one symmetric review band around the accepted median runtime.
   lower_wall_s = reference_wall_s * (1 - tolerance);
   upper_wall_s = reference_wall_s * (1 + tolerance);
   passed = current_wall_s >= lower_wall_s ...
      && current_wall_s <= upper_wall_s;
   reason = "";

   % Explain the failing side so an agent can distinguish slowdown from a
   % suspiciously fast measurement without retaining runner chatter.
   if current_wall_s < lower_wall_s
      reason = sprintf( ...
         'runtime %.6g s is below the accepted %.6g-%.6g s band', ...
         current_wall_s, lower_wall_s, upper_wall_s);
   elseif current_wall_s > upper_wall_s
      reason = sprintf( ...
         'runtime %.6g s is above the accepted %.6g-%.6g s band', ...
         current_wall_s, lower_wall_s, upper_wall_s);
   end
end
