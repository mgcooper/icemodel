function [upwd, dnwd] = GETUPDOWN(a, r, x, I0, z_edges, M)
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
   upwd = vertcat(x(1:N-1), 0);

   % Reconstruct y.
   dz = z_edges(3:N) - z_edges(2:N-1);
   dnwd = (a(2:N-1) + r(2:N-1)) ./ r(2:N-1) .* upwd(2:N-1) - ...
      (upwd(3:N) - upwd(1:N-2)) ./ (2.0 * dz .* r(2:N-1));
   dnwd = vertcat(I0, dnwd);
   dnwd = vertcat(dnwd, (a(N) + r(N)) / r(N) * upwd(N) ...
      - (upwd(N) - upwd(N-1)) / (dz(1) * r(N)));

   % Smooth any small bumps in the up/down curves. This legacy double-pass
   % smoothing comes from Glen's original Fortran implementation. It is kept
   % because removing it tends to introduce a visible near-surface kink in
   % log-scale flux plots.
   dntmp = (dnwd(2:N-1) + 0.5 * (dnwd(1:N-2) + dnwd(3:N))) / 2.0;
   uptmp = (upwd(2:N-1) + 0.5 * (upwd(1:N-2) + upwd(3:N))) / 2.0;

   % Do the ends.
   dntmp = vertcat((dnwd(2) + 0.5 * dnwd(1)) / 1.5, dntmp);
   uptmp = vertcat((upwd(2) + 0.5 * upwd(1)) / 1.5, uptmp);
   dntmp = vertcat(dntmp, (dnwd(N-1) + 0.5 * dnwd(N)) / 1.5);
   uptmp = vertcat(uptmp, (upwd(N-1) + 0.5 * upwd(N)) / 1.5);

   % Rewrite the arrays. Adjust to get back the Qsi at the surface.
   dnwd = vertcat(dnwd(1), dntmp(2:N-1), dnwd(end));
   upwd = vertcat(upwd(1), uptmp(2:N-1), upwd(end));

   % Repeat the smoothing once to reduce edge effects from the first pass.
   dntmp = (dnwd(2:N-1) + 0.5 * (dnwd(1:N-2) + dnwd(3:N))) / 2.0;
   uptmp = (upwd(2:N-1) + 0.5 * (upwd(1:N-2) + upwd(3:N))) / 2.0;

   % Do the ends.
   dntmp = vertcat((dnwd(2) + 0.5 * dnwd(1)) / 1.5, dntmp);
   uptmp = vertcat((upwd(2) + 0.5 * upwd(1)) / 1.5, uptmp);
   dntmp = vertcat(dntmp, (dnwd(N-1) + 0.5 * dnwd(N)) / 1.5);
   uptmp = vertcat(uptmp, (upwd(N-1) + 0.5 * upwd(N)) / 1.5);

   % Rewrite the arrays.
   dnwd = vertcat(dnwd(1), dntmp(2:N-1), dnwd(end));
   upwd = vertcat(upwd(1), uptmp(2:N-1), upwd(end));

   % % Compute the net solar flux (don't call this from the main, just to test)
   % dz = z_edges(3:M-1) - z_edges(2:M-2);
   % xynet = (up(2:M-2)+up(3:M-1))./dz-(down(2:M-2)+down(3:M-1))./dz;
   % xynet = vertcat(up(1)-down(1), xynet, up(M) - down(M));
end
