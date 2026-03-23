function [up, dn, xynet] = SOLVETWOSTREAM(I0, albedo, k_bulk, z_edges)
   %SOLVETWOSTREAM Solve Schlatter's two-stream system.
   %
   %  [up, dn] = SOLVETWOSTREAM(I0, albedo, k_bulk, z_edges)
   %  [up, dn, xynet] = SOLVETWOSTREAM(I0, albedo, k_bulk, z_edges)
   %
   % The notation here stays close to Schlatter's formulation. The solved
   % quantity is the radiative transfer state on the staggered spectral grid,
   % while XYnet (returned as optional XYNET) is the derived net-flux profile.
   %
   %#codegen

   % Notation here roughly follows Schlatter.
   M = numel(z_edges) - 1;

   % Compute the absorption and reflection coefficients at each interface.
   a = ((1.0 - albedo) / (1.0 + albedo)) * k_bulk;
   r = (2.0 * albedo / (1.0 - albedo ^ 2)) * k_bulk;

   % e = A_sub = sub diagonal
   % f = A_main = main diagonal
   % g = A_super = super diagonal
   % b = b_vector = right hand side
   % x = rad = solution

   % Extend the spectral edge coordinates downward by one control volume so the
   % lower boundary matches the padded bulk-extinction coefficients.
   dz_bottom = z_edges(M + 1) - z_edges(M);
   z_edges(M + 2) = z_edges(M + 1) + dz_bottom;

   % BULKEXTCOEFS is parameterized on the spectral cell thickness, so use the
   % same top-edge spacing in the upper boundary condition.
   deltaz = z_edges(2) - z_edges(1);

   % Initialize the tridiagonal system.
   e = zeros(M + 1, 1);
   f = zeros(M + 1, 1);
   g = zeros(M + 1, 1);
   b = zeros(M + 1, 1);

   % Account for the upper boundary condition.
   alfa = 1.0 / (a(1) + r(1));
   e(1) = 0.0;
   f(1) = 1.0;
   g(1) = -alfa / (deltaz + alfa);
   b(1) = r(1) * I0 * deltaz * alfa / (deltaz + alfa);

   % Fill the system between the boundaries.
   deltaz = z_edges(3:M + 2) - z_edges(2:M + 1);
   e(2:M + 1) = 1.0 + (r(3:M + 2) - r(1:M)) ./ (4.0 * r(2:M + 1));
   g(2:M + 1) = 1.0 - (r(3:M + 2) - r(1:M)) ./ (4.0 * r(2:M + 1));
   b(2:M + 1) = 0.0;
   f(2:M + 1) = deltaz ./ (2.0 * r(2:M + 1)) .* ...
      (a(2:M + 1) .* (r(3:M + 2) - r(1:M)) ...
      - r(2:M + 1) .* (a(3:M + 2) - a(1:M))) ...
      - (2.0 + deltaz .^ 2 .* k_bulk(2:M + 1) .^ 2);

   % Account for the lower boundary condition.
   g(M + 1) = 0.0;
   b(M + 1) = 0.0;

   % Solve the tridiagonal system.
   x = TRISOLVE(e, f, g, b);

   % Reconstruct the up/down fluxes.
   [up, dn] = GETUPDOWN(a, r, x, I0, z_edges, M);

   % Optionally return Schlatter's XYnet profile using the same interface
   % indexing as the direct two-stream formulation.
   if nargout > 2
      xynet = (up(2:M) + up(3:M + 1)) / 2.0 - (dn(2:M) + dn(3:M + 1)) / 2.0;
      xynet = [up(1) - dn(1); xynet];

      % Ensure the surface value matches the total absorbed shortwave. This
      % corrects any very small residual from a spectral grid that is too
      % shallow to absorb all of the incoming radiation.
      xynet(1) = min(-I0 * (1 - albedo), xynet(1));
   end
end
