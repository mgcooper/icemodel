function [De, wcoef, scoef] = WINDCOEF(wspd, z_0, z_obs)
   %WINDCOEF Compute wind transfer coefficients
   
   [kappa, gravity] = icemodel.physicalConstant('kappa', 'gravity');

   wcoef = (kappa / log(z_obs / z_0)) ^ 2;
   De = wspd .* wcoef;

   if nargout > 2
      scoef = nan(1, 3);
      scoef(1) = 5.3 * 9.4 * wcoef * sqrt(z_obs / z_0); % gamma Eq. A15
      scoef(2) = 9.4 * gravity * z_obs;
      scoef(3) = scoef(1) * sqrt(gravity * z_obs);
   end
end
