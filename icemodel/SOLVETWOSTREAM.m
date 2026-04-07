function [Qnet, Qup, Qdn] = SOLVETWOSTREAM(I0, albedo, k_bulk, z_edges)
   %SOLVETWOSTREAM Solve Schlatter's two-stream radiative transfer system.
   %
   %  Qnet = SOLVETWOSTREAM(I0, albedo, k_bulk, z_edges)
   %  [Qnet, Qup, Qdn] = SOLVETWOSTREAM(I0, albedo, k_bulk, z_edges)
   %
   %  Solves the two-stream equations on the staggered spectral grid for the
   %  upward (Qup) and downward (Qdn) diffuse flux profiles, then reconstructs
   %  the interface net flux Qnet (Schlatter's XYnet). The primary output Qnet
   %  gives the net spectral flux at each spectral control-volume interface.
   %
   %  Inputs:
   %     I0      - Incident spectral irradiance at the top surface [W m-2]
   %     albedo  - Surface albedo [1]
   %     k_bulk  - Bulk spectral extinction coefficients on the spectral grid
   %               [m-1], padded by BULKEXTCOEFS or BULKEXTCOEFSLOOKUP
   %     z_edges - Spectral control-volume edge depths [m] (M+1 values for M CVs)
   %
   %  Outputs:
   %     Qnet - Interface net spectral flux profile [W m-2] (M+1 values: one
   %            per spectral CV edge, from the top surface through the bottom).
   %            The absorbed flux in each spectral CV is Qnet(k) - Qnet(k+1).
   %            A min(-I0*(1-albedo),...) guard at the top boundary ensures the
   %            surface absorbed flux never exceeds the physical limit.
   %     Qup  - Upward diffuse flux on the staggered grid [W m-2] (M+2 values)
   %     Qdn  - Downward diffuse flux on the staggered grid [W m-2] (M+2 values)
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
   dz_bottom = z_edges(M+1) - z_edges(M);
   z_edges(M+2) = z_edges(M+1) + dz_bottom;

   % BULKEXTCOEFS is parameterized on the spectral cell thickness, so use the
   % same top-edge spacing in the upper boundary condition.
   deltaz = z_edges(2) - z_edges(1);

   % Initialize the tridiagonal system.
   e = zeros(M+1, 1);
   f = zeros(M+1, 1);
   g = zeros(M+1, 1);
   b = zeros(M+1, 1);

   % Account for the upper boundary condition.
   alfa = 1.0 / (a(1) + r(1));
   e(1) = 0.0;
   f(1) = 1.0;
   g(1) = -alfa / (deltaz + alfa);
   b(1) = r(1) * I0 * deltaz * alfa / (deltaz + alfa);

   % Fill the system between the boundaries.
   deltaz = z_edges(3:M+2) - z_edges(2:M+1);
   e(2:M+1) = 1.0 + (r(3:M+2) - r(1:M)) ./ (4.0 * r(2:M+1));
   g(2:M+1) = 1.0 - (r(3:M+2) - r(1:M)) ./ (4.0 * r(2:M+1));
   b(2:M+1) = 0.0;
   f(2:M+1) = deltaz ./ (2.0 * r(2:M+1)) .* (a(2:M+1) .* (r(3:M+2) - r(1:M)) ...
      - r(2:M+1) .* (a(3:M+2) - a(1:M))) ...
      - (2.0 + deltaz .^ 2 .* k_bulk(2:M+1) .^ 2);

   % Account for the lower boundary condition.
   g(M+1) = 0.0;
   b(M+1) = 0.0;

   % Solve the tridiagonal system.
   x = icemodel.numerics.trisolve(e, f, g, b);

   % Reconstruct the up/down fluxes.
   [Qup, Qdn] = GETUPDOWN(a, r, x, I0, z_edges, M);

   % Compute the net flux at each interface (Schlatter's XYnet). The staggered
   % up/dn have M+2 elements; averaging adjacent pairs gives M+1 interface
   % values (one per spectral CV edge, top through bottom). The min(-I0...
   % correction ensures the surface flux equals the total absorbed shortwave.
   Qnet = [min(-I0 * (1 - albedo), Qup(1) - Qdn(1)); ...
      (Qup(2:M+1) + Qup(3:M+2)) / 2.0 - (Qdn(2:M+1) + Qdn(3:M+2)) / 2.0];
end
