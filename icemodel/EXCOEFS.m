function De_h = EXCOEFS(wspd, wcoef)
   %EXCOEFS Compute the stability exchange coefficients
   De_h = wcoef * wspd;

   % wcoef is precomputed to improve speed:
   % wcoef = kappa ^ 2 / log(z_obs / z_0) ^ 2
   
   % function De_h = EXCOEFS(wspd, kappa, z_0, z_obs)
   %  De_h = kappa ^ 2 * wspd / log(z_obs / z_0) ^ 2; % [m s-1]
end
