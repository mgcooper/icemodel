function [T, Ts, f_ice, f_liq, subfail, subiter, dt_new] = RESETSUBSTEP( ...
      T, Ts, f_ice, f_liq, dt_max, subiter, maxsubiter, subfail, dt_sum)
   %RESETSUBSTEP Decrease timestep and reset variables on phase change overshoot
   if nargout > 4
      subfail = subfail + 1;
      subiter = min(subiter + 1, maxsubiter);

      % Note: this ensures dt_sum + dt_new does not exceed dt_max, but also
      % means dt gets off-track from dt_max / subiter. If that is problematic,
      % update subiter after this: subiter = dt_max / dt_new, but that likely
      % won't be an integer, so would need to fix it.
      dt_new = min(dt_max / subiter, dt_max - dt_sum);
   end
end
