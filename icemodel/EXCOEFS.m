function De_h = EXCOEFS(wspd, wcoef)
   %EXCOEFS Compute the stability exchange coefficients
   %
   %#codegen
   De_h = wcoef * wspd;

   % Note: wcoef is precomputed in WINDCOEFS to avoid recomputing here:
   % wcoef = kappa ^ 2 / log(z_wind / z_0) ^ 2
   %
   % Previous version of this function used the traditional formula:
   % function De_h = EXCOEFS(wspd, kappa, z_0, z_wind)
   %  De_h = kappa ^ 2 * wspd / log(z_wind / z_0) ^ 2; % [m s-1]
end
