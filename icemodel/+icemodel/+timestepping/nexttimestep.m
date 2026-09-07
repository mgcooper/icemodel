function [metstep, substep, dt_new] = nexttimestep(metstep, substep, ...
      ok, settings, diag)
   %NEXTTIMESTEP Advance the forcing index and adapt the next substep size.
   %
   % [metstep, substep, dt_new] = icemodel.timestepping.nexttimestep(...
   %    metstep, substep, ok, settings, diag)
   %
   % Update the forcing (met) index and timestep divisor, and ensure the new
   % substep size (dt_new) is valid to carry over to the next full forcing step.
   %
   % Inputs:
   %  metstep    - current forcing step index
   %  substep    - current timestep divisor (used to compute dt_new:
   %               dt_new = settings.dt_full_step / substep)
   %  ok         - full-step success flag
   %  settings   - solver settings; this function reads dt_full_step and
   %               maxsubstep (the max substep divisor, which sets dt_min)
   %  diag       - the forcing step's solver-diagnostics record; this
   %               function reads n_failed_substeps and the final solve
   %               attempt's inner iteration count diag.substep.n_iters
   %
   % Outputs:
   %  metstep    - incremented forcing step index
   %  substep    - updated substep divisor
   %  dt_new     - next-step substep length, bounded by dt_min
   %
   % See also: icemodel, skinmodel, icemodel.timestepping.newtimestep
   %
   %#codegen

   persistent cooldown
   if isempty(cooldown)
      cooldown = 0;
   end

   dt_max = settings.dt_full_step;
   maxsubstep = settings.maxsubstep;
   n_subfail = diag.n_failed_substeps;
   n_iters = diag.substep.n_iters;

   % Iteration thresholds for dt control.
   % n_iters >= N_hi: solver work is high, shrink dt (increase substep).
   % n_iters <= N_lo: solver work is low, grow dt faster (decrease substep by 2).
   N_hi = 12;
   N_lo = 2;

   % Start a short cooldown window after any failed substep/full-step fail
   % so dt does not re-grow immediately after a difficult step.
   if ~ok || n_subfail > 0
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

   % Enforce dt_min = dt_max / maxsubstep. acceptsubstep can shorten dt_new
   % below dt_min so the final substep exactly completes a full step.
   dt_new = max(dt_new, dt_max / maxsubstep);
   metstep = metstep + 1;
end
