function [metstep, substep, dt_new] = NEXTSTEP(metstep, substep, ...
      dt_new, dt_max, maxsubstep, ok, n_subfail, n_iters)
   %NEXTSTEP Advance forcing index and adapt substep size for next full step.
   %
   % Syntax:
   % [metstep, substep, dt_new] = NEXTSTEP(metstep, substep, dt_new, ...
   %    dt_max, maxsubstep, ok, n_subfail, n_iters)
   %
   % Inputs:
   % metstep    - current forcing step index
   % substep    - current substep divisor (dt_new = dt_max / substep)
   % dt_new     - current substep length
   % dt_max     - full forcing-step length
   % maxsubstep - max allowed substep divisor (sets dt_min)
   % ok         - full-step success flag
   % n_subfail  - number of failed/retried substeps in this full step
   % n_iters    - iterations used by the final subsurface solve
   %
   % Outputs:
   % metstep    - incremented forcing step index
   % substep    - updated substep divisor
   % dt_new     - next-step substep length, bounded by dt_min
   %
   %#codegen

   persistent cooldown
   if isempty(cooldown)
      cooldown = 0;
   end

   if nargin < 7
      n_subfail = 0;
   end
   if nargin < 8
      n_iters = 0;
   end

   % Iteration thresholds for dt control.
   % n_iters >= N_hi: solver work is high, shrink dt (increase substep).
   % n_iters <= N_lo: solver work is low, grow dt faster (decrease substep by 2).
   N_hi = 12;
   N_lo = 2;

   % Start a short cooldown window after any failed substep/full-step fail
   % so dt does not re-grow immediately after a difficult step.
   if not(ok) || n_subfail > 0
      cooldown = 2;
   end

   % Adapt dt based on success/failure and solver work
   if ok

      if cooldown > 0
         cooldown = cooldown - 1;
         if n_subfail > 0 || n_iters >= N_hi
            substep = min(substep + 1, maxsubstep);
         end

      else
         if n_subfail > 0 || n_iters >= N_hi
            substep = min(substep + 1, maxsubstep);
         elseif n_iters <= N_lo
            substep = max(1, substep - 2);
         else
            substep = max(1, substep - 1);
         end
      end

      dt_new = dt_max / substep;

   else
      substep = min(substep + 1, maxsubstep);
      dt_new = dt_max / substep;
   end

   % Enforce dt_min = dt_max / maxsubstep (if UPDATESUBSTEP shortens dt_new to
   % <dt_min to ensure the final substep exactly completes a full step).
   dt_new = max(dt_new, dt_max / maxsubstep);
   metstep = metstep + 1;
end
