function [Ts, T, f_ice, f_liq, n_subfail, substep, dt_new] = RESETSUBSTEP( ...
      Ts, T, f_ice, f_liq, dt_max, substep, maxsubstep, n_subfail, dt_sum)
   %RESETSUBSTEP Decrease timestep and reset variables on phase change overshoot
   %
   %#codegen
   if nargout > 4
      n_subfail = n_subfail + 1;
      substep = min(substep + 1, maxsubstep);

      % Activate this for aggressive timestep shortening.
      % substep = min(substep + min(n_subfail, 4), maxsubstep);

      % Note: this ensures dt_sum + dt_new does not exceed dt_max, but also
      % means dt gets off-track from dt_max / substep. NEXTSTEP performs a
      % similar check to ensure dt_new = dt_min if this shortens it <dt_min.
      dt_new = min(dt_max / substep, dt_max - dt_sum);
   end
end
