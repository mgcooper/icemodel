function [Ts, T, f_ice, f_liq, k_eff, dt_sum, dt_new] = acceptsubstep( ...
      Ts, T, f_ice, f_liq, k_eff, dt_FULL_STEP, dt_sum, dt_new, TINY)
   %ACCEPTSUBSTEP Accept the substep: checkpoint state and credit its time.
   %
   %  [Ts, T, f_ice, f_liq, k_eff, dt_sum, dt_new] = ...
   %     icemodel.timestepping.acceptsubstep( ...
   %     Ts, T, f_ice, f_liq, k_eff, dt_FULL_STEP, dt_sum, dt_new, TINY)
   %
   % Accepting a substep is one operation with two parts:
   %   1. Checkpoint the accepted column state and conductivity
   %      (Ts, T, f_ice, f_liq, k_eff) for the next substep. The state
   %      passes through unchanged; the caller binds these outputs to the
   %      x-prefixed checkpoint names.
   %   2. Credit the substep duration to dt_sum and adjust dt_new to
   %      exactly complete the full step without overshooting.
   %
   % A forced advance calls this with the restored checkpoint state: it
   % accepts only elapsed time, and the state pass-through returns the
   % unchanged checkpoint.
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
