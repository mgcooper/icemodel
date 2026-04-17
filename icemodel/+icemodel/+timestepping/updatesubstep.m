function [Ts, T, f_ice, f_liq, dt_sum, dt_new] = updatesubstep( ...
      Ts, T, f_ice, f_liq, dt_FULL_STEP, dt_sum, dt_new, TINY)
   % Checkpoint the accepted state and advance time within the full step.
   %
   %  [Ts, T, f_ice, f_liq, dt_sum, dt_new] = ...
   %     icemodel.timestepping.updatesubstep( ...
   %     Ts, T, f_ice, f_liq, dt_FULL_STEP, dt_sum, dt_new, TINY)
   %
   % This function performs the timestepping bookkeeping after a
   % successful substep:
   %   1. Checkpoints the accepted column state (Ts, T, f_ice, f_liq)
   %      as the initial condition for the next substep.
   %   2. Accumulates the substep duration into dt_sum.
   %   3. Adjusts dt_new to exactly complete the full step without
   %      overshooting.
   %
   % Surface running state (liqflag, ro_sfc, hv_atm, H_e, f_res_por)
   % is no longer computed here. Those quantities are derived at
   % substep entry by icemodel.surface.update_surface_state so that
   % they are always consistent with the current forcing step and
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
