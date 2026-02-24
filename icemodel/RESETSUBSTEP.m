function [Ts, T, f_ice, f_liq, subfail, subiter, dt_new] = RESETSUBSTEP( ...
      Ts, T, f_ice, f_liq, dt_max, subiter, maxsubiter, subfail, dt_sum)
   %RESETSUBSTEP Decrease timestep and reset variables on phase change overshoot
   %
   %#codegen
   if nargout > 4
      subfail = subfail + 1;
      subiter = min(subiter + 1, maxsubiter);

      % Activate this for aggressive timestep shortening.
      % subiter = min(subiter + min(subfail, 4), maxsubiter);

      % Note: this ensures dt_sum + dt_new does not exceed dt_max, but also
      % means dt gets off-track from dt_max / subiter. NEXTSTEP performs a
      % similar check to ensure dt_new = dt_min if this shortens it <dt_min.
      dt_new = min(dt_max / subiter, dt_max - dt_sum);
   end
end
