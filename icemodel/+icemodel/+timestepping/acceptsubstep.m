function [T_sfc, T_ice, f_ice, f_liq, k_eff, dt_sum, dt_new, settings, diag] = ...
      acceptsubstep(T_sfc, T_ice, f_ice, f_liq, k_eff, dt_sum, dt_new, ...
      TINY, settings, settings0, diag)
   %ACCEPTSUBSTEP Accept the substep: checkpoint state and time bookkeeping.
   %
   %  [T_sfc, T_ice, f_ice, f_liq, k_eff, dt_sum, dt_new, settings, diag] = ...
   %     icemodel.timestepping.acceptsubstep( ...
   %     T_sfc, T_ice, f_ice, f_liq, k_eff, dt_sum, dt_new, TINY, ...
   %     settings, settings0, diag)
   %
   % Accepting a substep is one operation with four parts:
   %   1. Checkpoint the accepted state (T_sfc, T_ice, f_ice, f_liq, k_eff) for
   %      the next substep. The state passes through unchanged as the x-prefixed
   %      checkpoint names. k_eff is derived from the state and checkpointed in
   %      the main flow.
   %   2. Credit the substep duration to dt_sum and adjust dt_new to
   %      exactly complete the full step (settings.dt_full_step) without
   %      overshooting.
   %   3. Promote diag.substep into the step record when it holds an
   %      accepted solve (see icemodel.couplers.update_solver_diag). On a
   %      forced advance diag.substep holds the rejected solve, the
   %      promotion is skipped, and the step keeps its previous record.
   %   4. Restore the primary SETTINGS0 so the next substep never inherits
   %      the recovery settings.
   %
   % See also:
   %   icemodel.surface.update_surface_state,
   %   icemodel.couplers.update_solver_diag,
   %   icemodel.timestepping.checksubstep,
   %   icemodel.timestepping.nexttimestep
   %
   %#codegen

   % Allocate this substep to the timestep.
   dt_sum = dt_sum + dt_new;

   % Adjust dt to exactly complete the full step without going over.
   % The first condition is true if the full step is incomplete, the
   % second is true if the next substep will exceed the full step.
   if (settings.dt_full_step - dt_sum) > TINY ...
         && (dt_sum + dt_new - settings.dt_full_step) > TINY
      dt_new = settings.dt_full_step - dt_sum;
   end

   % Record the accepted solve, then return to the default solver settings.
   diag = icemodel.couplers.update_solver_diag( ...
      settings.cpl_recovery_active, diag);
   settings = settings0;
end
