function [iter, metiter, subiter, dt_new] = NEXTSTEP(iter, metiter, subiter, ...
      dt_flag, dt_max, OK, dt_new)

   % i think this check can be deleted b/c next check gets it
   if dt_flag == true
      dt_new = dt_max/subiter;
   end

   % IF CONVERGENCE IS SUCCESSFUL, INCREASE THE TIMESTEP
   if OK == true
      subiter = max(1, subiter-1); % reverse the subiter decrement
      dt_new = dt_max/subiter;
   end

   iter = iter + 1;
   metiter = metiter + 1;
end
