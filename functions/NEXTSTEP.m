function [iter, metiter, subiter, dt_new] = NEXTSTEP(iter, metiter, subiter, ...
      dt_new, dt_max, dt_min, OK)

   % If convergence is successful, increase the timestep
   if OK == true
      subiter = max(1, subiter - 1); % reverse the subiter decrement
      dt_new = dt_max / subiter;
   else
      % If UPDATESUBSTEP shortens dt_new to less than dt_min so the final
      % substep exactly completes a full step, this ensures its set to dt_new.
      % However, as-is this only triggers if OK is false. If OK is true, the
      % step above will take care of this. 
      dt_new = max(dt_new, dt_min);
   end

   iter = iter + 1;
   metiter = metiter + 1;
end
