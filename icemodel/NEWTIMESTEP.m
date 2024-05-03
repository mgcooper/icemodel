function [dt_sum, subfail, OK, d_liq, d_evp, d_lyr] = NEWTIMESTEP(f_liq)
   %NEWTIMESTEP initialize new timestep

   OK = false;             % assume sub-step failed
   d_liq = 0.0 * f_liq;    % reset the change in liq water content
   d_evp = 0.0 * f_liq;    % reset the evaporation change in water content
   d_lyr = 0.0 * f_liq;    % reset the layer change
   dt_sum = 0.0;
   subfail = 0;            % keep track of failed substeps
end
