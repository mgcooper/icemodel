function [Ts, T, f_ice, f_liq, k_eff, dt_sum, dt_new] = updatesubstep( ...
      Ts, T, f_ice, f_liq, k_eff, dt_FULL_STEP, dt_sum, dt_new, TINY)
   %UPDATESUBSTEP Checkpoint the accepted state and advance time within the
   % full step.
   %
   %  [Ts, T, f_ice, f_liq, k_eff, dt_sum, dt_new] = ...
   %     icemodel.timestepping.updatesubstep( ...
   %     Ts, T, f_ice, f_liq, k_eff, dt_FULL_STEP, dt_sum, dt_new, TINY)
   %
   % This function performs the timestepping bookkeeping after a
   % successful substep:
   %   1. Checkpoints the accepted column state and conductivity
   %      (Ts, T, f_ice, f_liq, k_eff) for the next substep.
   %   2. Accumulates the substep duration into dt_sum.
   %   3. Adjusts dt_new to exactly complete the full step without
   %      overshooting.
   %
   % This function does not compute the surface running state
   % (liqflag, ro_sfc, hv_atm, H_e, f_res_por).
   % icemodel.surface.update_surface_state derives those quantities at
   % substep entry, so they always match the current forcing step and
   % column state.
   %
   % See also:
   %   icemodel.surface.update_surface_state,
   %   icemodel.timestepping.checksubstep,
   %   icemodel.timestepping.nexttimestep
   %
   %#codegen

   % Allocate this substep to the timestep
   dt_sum = dt_sum + dt_new;

   % Adjust dt to exactly complete the full step without going over.
   % The first condition is true if the full step is incomplete, the
   % second is true if the next substep will exceed the full step.
   if (dt_FULL_STEP - dt_sum) > TINY ...
         && (dt_sum + dt_new - dt_FULL_STEP) > TINY
      dt_new = dt_FULL_STEP - dt_sum;
   end
end
