function [T_sfc, T_ice, f_ice, f_liq, k_eff, n_subfail, substep, dt_new] = ...
      resetsubstep(T_sfc, T_ice, f_ice, f_liq, k_eff, dt_max, substep, ...
      maxsubstep, n_subfail, dt_sum)
   %RESETSUBSTEP Restore the accepted state and shorten the retry timestep.
   %
   % Call this after a failed substep. K_EFF travels with the temperature and
   % phase checkpoint so diagnostics never retain conductivity from a rejected
   % solve.
   %
   %#codegen
   if nargout > 5
      n_subfail = min(n_subfail + 1, maxsubstep);
      substep = min(substep + 1, maxsubstep);

      % Activate this for aggressive timestep shortening.
      % substep = min(substep + min(n_subfail, 4), maxsubstep);

      % This keeps dt_sum + dt_new at or below dt_max. It also means dt no
      % longer equals dt_max / substep. nexttimestep performs a similar check
      % to set dt_new = dt_min if this shortens dt_new below dt_min.
      dt_new = min(dt_max / substep, dt_max - dt_sum);
   end
end
