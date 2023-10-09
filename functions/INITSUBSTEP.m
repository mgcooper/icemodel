function [dt_sum, subfail, dt_flag, OK, d_liq, d_drn, d_evp] = INITSUBSTEP(f_liq)
%INITSUBSTEP initialize new substep
dt_sum      =  0.0;
subfail     =  0;                % keep track of failed substeps
dt_flag     =  false;            % sub-step error flag
OK          =  false;            % assume sub-step failed
d_liq       =  0.0.*f_liq;       % reset the change in liq water content
d_drn       =  0.0.*f_liq;       % reset the drained water
d_evp       =  0.0.*f_liq;       % reset the drained water