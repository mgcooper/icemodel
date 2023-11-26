function De = WINDCOEF(wspd, z_0, z_obs)
   kappa = icemodel.physicalConstant('kappa');
   De = wspd .* (kappa / log(z_obs / z_0)) ^ 2;
end
