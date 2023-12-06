function [T, Ts, f_ice, f_liq, subfail, subiter, dt_new] = RESETSUBSTEP( ...
      T, Ts, f_ice, f_liq, dt_max, subiter, maxsubiter, subfail)
   %RESETSUBSTEP Decrease timestep and reset variables on phase change overshoot
   if nargout > 4
      subfail = subfail + 1;
      subiter = min(subiter + 1, maxsubiter);
      dt_new = dt_max / subiter;
   end
end
