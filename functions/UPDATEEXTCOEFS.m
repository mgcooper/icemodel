function [Sc, chi] = UPDATEEXTCOEFS(Qsi, albedo, z_spect, dz_spect, ...
      z_therm, dz_therm, JJ_therm, dz, ro_sno, total_solar, spect_lower, ...
      spect_upper, solardwavl)
   %UPDATEEXTCOEFS Update the bulk extinction coefficients
   
   persistent z_walls
   if isempty(z_walls)
      z_walls = round([0; z_spect + dz_spect/2; z_spect(end) + 3*dz_spect/2], 3);
   end

   % Transform the mass density to the spectral grid resolution
   ro_sno = interp1(z_therm, ro_sno, z_spect, 'nearest', 'extrap');
   ro_sno = max(ro_sno, 300);

   % Compute the downward bulk extinction coefficient.
   bulkcoefs = -log((sum(solardwavl .* exp(spect_lower .* ro_sno), 2))./ ...
      (sum(solardwavl .* exp(spect_upper .* ro_sno), 2))) ./ dz_spect;
   %  NOTE: solardwavl comes multiplied by dwavl, and sum/trapz are identical

   % Cast in a form that can be used in the two-stream computation (add the
   % boundaries).  Here I have assumed that it is okay to call the boundaries
   % equal to the value at the center of that grid cell (the value prevails
   % thoughout the cell).
   bulkcoefs = vertcat(bulkcoefs, bulkcoefs(end), bulkcoefs(end));

   % Compute the a and r coefficients from the surface albedo and the 
   % extinction coefficients.
   a = (1.0 - albedo) ./ (1.0 + albedo) .* bulkcoefs;
   r = 2.0 .* albedo .* bulkcoefs ./ (1.0 - albedo^2);

   % Solve the system of equations for the two-stream model
   M = numel(z_spect); % notation here roughly follows Schlatter
   M1 = M+1;
   M2 = M+2;

   % initialize the matrix
   e = zeros(M1, 1);
   f = zeros(M1, 1);
   g = zeros(M1, 1);
   b = zeros(M1, 1);

   % Account for the upper boundary condition
   alfa = 1.0 / (a(1) + r(1));
   e(1) = 0.0;
   f(1) = 1.0;
   g(1) = -alfa / (dz_spect + alfa);
   b(1) = r(1) * total_solar * dz_spect * alfa / (dz_spect + alfa);

   % Fill the vectors between the boundaries
   deltaz = z_walls(3:M2) - z_walls(2:M1);
   e(2:M1) = 1.0 + (r(3:M2) - r(1:M)) ./ (4.0 .* r(2:M1));
   f(2:M1) = deltaz ./ (2.0 .* r(2:M1)) .* (a(2:M1) .* (r(3:M2) - r(1:M)) ...
      - r(2:M1) .* (a(3:M2) - a(1:M))) - (2.0 + deltaz.^2 .* bulkcoefs(2:M1).^2);
   g(2:M1) = 1.0 - (r(3:M2) - r(1:M)) ./ (4.0 * r(2:M1));
   b(2:M1) = 0.0;

   % Account for the lower boundary condition
   g(M1) = 0.0;
   b(M1) = 0.0;

   % Solve the equation
   x = zeros(M1, 1);
   for k = 2:numel(f) % forward elimination
      f(k) = f(k) - e(k) / f(k-1) * g(k-1);
      b(k) = b(k) - e(k) / f(k-1) * b(k-1);
   end
   x(M1) = b(M1) / f(M1); % back substitution
   for k = M:-1:1
      x(k) = (b(k) - g(k) * x(k+1)) / f(k);
   end

   % Add the boundary conditions to up and reconstruct down. X = up, Y = down.
   up = vertcat(x(1:M1), 0); % Add the boundary conditions to rad.

   % Reconstruct y.
   dn = (a(2:M1) + r(2:M1)) ./ r(2:M1) .* up(2:M1) - (up(3:M2) - up(1:M2-2)) ...
      ./ (2.0 .* deltaz .*r (2:M1));
   dn = vertcat(total_solar, dn);
   dn = vertcat(dn, (a(M2) + r(M2)) / r(M2) * up(M2) - (up(M2) - up(M1)) / ...
      (deltaz(1) * r(M2)));

   % Smooth any small bumps in the up and down curve.  This will assist
   % in any delta up, delta down computations. Do the interior.
   dntmp = (dn(2:M1) + 0.5 .* (dn(1:M2-2) + dn(3:M2))) ./ 2.0;
   uptmp = (up(2:M1) + 0.5 .* (up(1:M2-2) + up(3:M2))) ./ 2.0;

   % Do the ends.
   dntmp = vertcat((dn(2) + 0.5 * dn(1)) / 1.5, dntmp);
   uptmp = vertcat((up(2) + 0.5 * up(1)) / 1.5, uptmp);
   dntmp = vertcat(dntmp, (dn(M1) + 0.5 * dn(M2)) / 1.5);
   uptmp = vertcat(uptmp, (up(M1) + 0.5 * up(M2)) / 1.5);

   % Rewrite the arrays. Adjust to get back the Qsi at the surface.
   dn = vertcat(dn(1), dntmp(2:M1), dn(end));
   up = vertcat(up(1),   uptmp(2:M1),   up(end));

   % Repeat the smoothing (deal with edge effects)
   dntmp = (dn(2:M1) + 0.5 * (dn(1:M) + dn(3:M2))) / 2.0;
   uptmp = (up(2:M1) + 0.5 * (up(1:M) + up(3:M2))) / 2.0;

   % Do the ends.
   dntmp = vertcat(dn(2) + 0.5 * dn(1) / 1.5, dntmp);
   uptmp = vertcat(up(2) + 0.5 * up(1) / 1.5, uptmp);
   dntmp = vertcat(dntmp, dn(M1) + 0.5 * dn(M2) / 1.5);
   uptmp = vertcat(uptmp, up(M1) + 0.5 * up(M2) / 1.5);

   % Rewrite the arrays.
   dn = vertcat(dn(1), dntmp(2:M1), dn(end));
   up = vertcat(up(1), uptmp(2:M1), up(end));

   % Build an array of source-term values on the c.v. boundaries.
   xynet = (up(2:M) + up(3:M1)) ./ 2.0 - (dn(2:M) + dn(3:M1)) ./ 2.0;
   xynet = vertcat(up(1) - dn(1), xynet);

   % Ensure xynet(1) is equal to the total absorbed radiation. This corrects
   % for any (very) small errors in the two-stream model, typically due to a
   % spectral grid that is too shallow to absorb all the radiation
   xynet(1) = min(-total_solar * (1 - albedo), xynet(1));

   % xynet is the net solar radiation at each level. In Schlatter it is
   % just Q. The source term (absorbed flux per unit volume) is dQ/dz.
   % Below, I call the absorbed flux dQ, so Sc = dQ/dz = d(xynet)/dz. In
   % Glen's original model, the quantity dQ was not defined, xynet was
   % converted directly to dQ/dz in ICE_HEAT by scaling d(xynet) by Qsip and
   % then to Sc.

   % Transform the pentrating radiation back to the thermal grid as dQ
   % extrapolate xynet to the bottom of the thermal grid
   
   JJnew = JJ_therm * dz_therm / dz_spect;
   JJext = round(JJnew-M, 0);
   delxy = xynet(M) - xynet(M-1);
   xyext = [xynet(M) + delxy; xynet(M) + 2*delxy; zeros(JJext-2, 1) + delxy];
   xynew = [xynet; xyext];
   
   % Upscale to the thermal grid
   dQ = xynew(1:JJnew-1) - xynew(2:JJnew);
   dQ(JJnew) = dQ(JJnew - 1);
   
   % Get the number of grid cells in the first 3 m segment on each grid
   rz = JJnew / JJ_therm;
   
   % Reshape the radiation into equal chunks along the thermal grid
   dQ = reshape(dQ, rz, JJ_therm);
   
   % Sum up those equal chunks to get the absorbed radiation in each thermal c.v.
   dQ = transpose(sum(dQ, 1));

   % Compute the amount of solar radiation absorbed in the top grid cell, to
   % be allocated to the SEB. Optionally enhance it by qsfactor to account
   % for impurities on the ice surface or other factors.
   if albedo > 0.65 % && snowd > 0.05
      chi = 0.9;
   else
      chi = dQ(1) / sum(dQ);
   end
   Sc = -(1.0 - chi) * Qsi / total_solar * dQ ./ dz; % [W/m3]
end
