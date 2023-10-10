function [subfail, subiter, dt_new, T, Tsfc, f_ice, f_liq] = RESETSUBSTEP( ...
      xT, xTsfc, xf_ice, xf_liq, dt_max, subiter, maxsubiter, subfail)
   %RESETSUBSTEP Decrease timestep and reset variables on phase change overshoot

   subfail = subfail+1;
   subiter = min(subiter+1, maxsubiter);
   dt_new = dt_max/subiter;
   T = xT;
   Tsfc = xTsfc;
   f_ice = xf_ice;
   f_liq = xf_liq;
end
