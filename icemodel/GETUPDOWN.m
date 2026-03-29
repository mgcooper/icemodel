function [up, dn] = GETUPDOWN(a, r, x, I0, z_edges, M)
   %GETUPDOWN Reconstruct the up/down flux profiles from the two-stream solve.
   %
   % Compute up/down flux at each depth from the two-stream solution vector x
   %
   % See also:
   %
   %#codegen

   % X = up, Y = down
   N = M+2;

   % Add the boundary conditions to rad.
   up = vertcat(x(1:N-1), 0);

   % Reconstruct y.
   dz = z_edges(3:N) - z_edges(2:N-1);
   dn = (a(2:N-1) + r(2:N-1)) ./ r(2:N-1) .* up(2:N-1) - ...
      (up(3:N) - up(1:N-2)) ./ (2.0 * dz .* r(2:N-1));
   dn = vertcat(I0, dn);
   dn = vertcat(dn, (a(N) + r(N)) / r(N) * up(N) ...
      - (up(N) - up(N-1)) / (dz(1) * r(N)));

   % Smooth any small bumps in the up/down curves. This legacy double-pass
   % smoothing comes from Glen's original Fortran implementation.
   dntmp = (dn(2:N-1) + 0.5 * (dn(1:N-2) + dn(3:N))) / 2.0;
   uptmp = (up(2:N-1) + 0.5 * (up(1:N-2) + up(3:N))) / 2.0;

   % Do the ends.
   dntmp = vertcat((dn(2) + 0.5 * dn(1)) / 1.5, dntmp);
   uptmp = vertcat((up(2) + 0.5 * up(1)) / 1.5, uptmp);
   dntmp = vertcat(dntmp, (dn(N-1) + 0.5 * dn(N)) / 1.5);
   uptmp = vertcat(uptmp, (up(N-1) + 0.5 * up(N)) / 1.5);

   % Rewrite the arrays. Adjust to get back the Qsi at the surface.
   dn = vertcat(dn(1), dntmp(2:N-1), dn(end));
   up = vertcat(up(1), uptmp(2:N-1), up(end));

   % Repeat the smoothing once to reduce edge effects from the first pass.
   dntmp = (dn(2:N-1) + 0.5 * (dn(1:N-2) + dn(3:N))) / 2.0;
   uptmp = (up(2:N-1) + 0.5 * (up(1:N-2) + up(3:N))) / 2.0;

   % Do the ends.
   dntmp = vertcat((dn(2) + 0.5 * dn(1)) / 1.5, dntmp);
   uptmp = vertcat((up(2) + 0.5 * up(1)) / 1.5, uptmp);
   dntmp = vertcat(dntmp, (dn(N-1) + 0.5 * dn(N)) / 1.5);
   uptmp = vertcat(uptmp, (up(N-1) + 0.5 * up(N)) / 1.5);

   % Rewrite the arrays.
   dn = vertcat(dn(1), dntmp(2:N-1), dn(end));
   up = vertcat(up(1), uptmp(2:N-1), up(end));
end
