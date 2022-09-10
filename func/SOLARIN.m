%--------------------------------------------------------------------------
%   Compute the incoming solar radiation.  
%--------------------------------------------------------------------------

function Qsi = SOLARIN(J_day_start,iter,dt,xlat,cloud_frac, ...
     slope_az,terrain_slope,ihrs_day,transmiss,ihour)
 
%--------------------------------------------------------------------------

% Here I am going to assume that we are using daily time steps.
% If you have some other time step, see MicroMet for ideas about what to do.

    J_day = iter + J_day_start - 1;
    Qsi_sum = 0.0;
    if (dt==86400.0)
        for ihour=1:ihrs_day
            xhour = real(ihour);
            Qsi = SOLAR_RAD(J_day,xlat,cloud_frac,xhour,slope_az, ...
                                        terrain_slope,transmiss);
            Qsi_sum = Qsi_sum + Qsi_tmp;
        end
        Qsi = Qsi_sum / real(ihrs_day);
    else
        J_day = iter/24 + J_day_start - 1;
        xhour = ihour;
        Qsi = SOLAR_RAD(J_day,xlat,cloud_frac,xhour,slope_az, ...
                                    terrain_slope,transmiss);
    end
