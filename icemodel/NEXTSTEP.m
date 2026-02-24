function [metiter, subiter, dt_new] = NEXTSTEP(metiter, subiter, ...
      dt_new, dt_max, maxsubiter, OK, subfail, N)
   %NEXTSTEP Advance forcing index and adapt substep size for next full step.
   %
   % Syntax:
   % [metiter, subiter, dt_new] = NEXTSTEP(metiter, subiter, dt_new, ...
   %    dt_max, maxsubiter, OK, subfail, N)
   %
   % Inputs:
   % metiter    - current forcing iterator
   % subiter    - current substep divisor (dt_new = dt_max / subiter)
   % dt_new     - current substep length
   % dt_max     - full forcing-step length
   % maxsubiter - max allowed substep divisor (sets dt_min)
   % OK         - full-step success flag
   % subfail    - number of failed/retried substeps in this full step
   % N          - iterations used by the final subsurface solve
   %
   % Outputs:
   % metiter    - incremented forcing iterator
   % subiter    - updated substep divisor
   % dt_new     - next-step substep length, bounded by dt_min
   %
   %#codegen

   persistent cooldown
   if isempty(cooldown)
      cooldown = 0;
   end

   if nargin < 7
      subfail = 0;
   end
   if nargin < 8
      N = 0;
   end

   % Iteration thresholds for dt control.
   % N >= N_hi: solver work is high, shrink dt (increase subiter).
   % N <= N_lo: solver work is low, grow dt faster (decrease subiter by 2).
   N_hi = 12;
   N_lo = 2;

   % Start a short cooldown window after any failed substep/full-step fail
   % so dt does not re-grow immediately after a difficult step.
   if not(OK) || subfail > 0
      cooldown = 2;
   end

   % Adapt dt based on success/failure and solver work
   if OK

      if cooldown > 0
         cooldown = cooldown - 1;
         if subfail > 0 || N >= N_hi
            subiter = min(subiter + 1, maxsubiter);
         end

      else
         if subfail > 0 || N >= N_hi
            subiter = min(subiter + 1, maxsubiter);
         elseif N <= N_lo
            subiter = max(1, subiter - 2);
         else
            subiter = max(1, subiter - 1);
         end
      end

      dt_new = dt_max / subiter;

   else
      subiter = min(subiter + 1, maxsubiter);
      dt_new = dt_max / subiter;
   end

   % Enforce dt_min = dt_max / maxsubiter (if UPDATESUBSTEP shortens dt_new to
   % <dt_min to ensure the final substep exactly completes a full step).
   dt_new = max(dt_new, dt_max / maxsubiter);
   metiter = metiter + 1;
end
