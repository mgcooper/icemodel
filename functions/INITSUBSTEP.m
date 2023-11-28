function [dt_sum, subfail, OK, d_liq, d_drn, d_evp] = INITSUBSTEP(f_liq)
   %INITSUBSTEP initialize new substep

   OK = false;             % assume sub-step failed
   d_liq = 0.0 * f_liq;    % reset the change in liq water content
   d_drn = 0.0 * f_liq;    % reset the drained water
   d_evp = 0.0 * f_liq;    % reset the drained water
   dt_sum = 0.0;
   subfail = 0;            % keep track of failed substeps
end
