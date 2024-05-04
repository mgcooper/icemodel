function [Sc, chi] = UPDATEEXTCOEFS(Qsi, albedo, Q0, dz_spect, spect_N, ...
      spect_S, solardwavl, Sc, dz_therm, ro_sno)
   %UPDATEEXTCOEFS Update the bulk extinction coefficients
   %
   %#codegen

   if Qsi < 1e-3
      Sc = 0.0 * Sc;
      chi = 1.0;
      return
   end

   % Transform the mass density to the spectral grid
   ro_sno = max(interp1(cumsum(dz_therm)-dz_therm/2, ro_sno, ...
      cumsum(dz_spect) - dz_spect/2, 'nearest', 'extrap'), 300);

   % Compute the downward bulk extinction coefficient.
   bulkcoefs = -log((sum(solardwavl .* exp(spect_S .* ro_sno), 2)) ...
      ./ (sum(solardwavl .* exp(spect_N .* ro_sno), 2))) ./ dz_spect(1);
   % NOTE: solardwavl comes multiplied by dwavl

   % Add the boundaries for the two-stream computation.
   bulkcoefs = vertcat(bulkcoefs, bulkcoefs(end), bulkcoefs(end));

   % Compute a and r from the surface albedo and the extinction coefficients.
   a = (1.0 - albedo) / (1.0 + albedo) * bulkcoefs;
   r = 2.0 * albedo * bulkcoefs ./ (1.0 - albedo ^ 2);

   % Solve the two-stream system of equations for the up flux
   M = numel(dz_spect) + 2;

   % Initialize the matrix
   e = zeros(M-1, 1);
   f = zeros(M-1, 1);
   g = zeros(M-1, 1);
   b = zeros(M-1, 1);

   % Account for the upper boundary condition
   alfa = 1.0 / (a(1) + r(1));
   e(1) = 0.0;
   f(1) = 1.0;
   g(1) = -alfa / (dz_spect(1) + alfa);
   b(1) = r(1) * Q0 * dz_spect(1) * alfa / (dz_spect(1) + alfa);

   % Fill the vectors between the boundaries
   deltaz = [dz_spect(2:end); dz_spect(end)];
   e(2:M-1) = 1.0 + (r(3:M) - r(1:M-2)) ./ (4.0 * r(2:M-1));
   f(2:M-1) = deltaz ./ (2.0 * r(2:M-1)) ...
      .* (a(2:M-1) .* (r(3:M) - r(1:M-2)) - r(2:M-1) .* (a(3:M) - a(1:M-2))) ...
      - (2.0 + deltaz .^ 2 .* bulkcoefs(2:M-1) .^ 2);
   g(2:M-1) = 1.0 - (r(3:M) - r(1:M-2)) ./ (4.0 * r(2:M-1));
   b(2:M-1) = 0.0;

   % Account for the lower boundary condition
   g(M-1) = 0.0;
   b(M-1) = 0.0;

   % Solve the equation
   up = zeros(M-1, 1);
   for k = 2:numel(f)
      f(k) = f(k) - e(k) / f(k-1) * g(k-1);
      b(k) = b(k) - e(k) / f(k-1) * b(k-1);
   end
   up(M-1) = b(M-1) / f(M-1);
   for k = numel(f)-1:-1:1
      up(k) = (b(k) - g(k) * up(k+1)) / f(k);
   end

   % Add the boundary conditions to up and reconstruct down. X = up, Y = down.
   up = vertcat(up(1:M-1), 0);

   % Reconstruct down.
   dn = vertcat(Q0, ...
      (a(2:M-1) + r(2:M-1)) ./ r(2:M-1) .* up(2:M-1) - (up(3:M) - up(1:M-2)) ...
      ./ (2.0 * deltaz .* r(2:M-1)), ...
      (a(M) + r(M)) / r(M) * up(M) - (up(M) - up(M-1)) ...
      / (deltaz(M-2) * r(M))); % note: dz(M-1) for the bottom

   % Do the interior.
   dntmp = (dn(2:M-1) + 0.5 * (dn(1:M-2) + dn(3:M))) / 2.0;
   uptmp = (up(2:M-1) + 0.5 * (up(1:M-2) + up(3:M))) / 2.0;

   % Do the ends.
   dntmp = vertcat((dn(2) + 0.5 * dn(1)) / 1.5, dntmp);
   uptmp = vertcat((up(2) + 0.5 * up(1)) / 1.5, uptmp);
   dntmp = vertcat(dntmp, (dn(M-1) + 0.5 * dn(M)) / 1.5);
   uptmp = vertcat(uptmp, (up(M-1) + 0.5 * up(M)) / 1.5);

   % Rewrite the arrays, ensuring Qsi at the surface is conserved.
   dn = vertcat(dn(1), dntmp(2:M-1), dn(end));
   up = vertcat(up(1), uptmp(2:M-1), up(end));

   % Repeat the interior.
   dntmp = (dn(2:M-1) + 0.5 * (dn(1:M-2) + dn(3:M))) / 2.0;
   uptmp = (up(2:M-1) + 0.5 * (up(1:M-2) + up(3:M))) / 2.0;

   % Do the ends.
   dntmp = vertcat(dn(2) + 0.5 * dn(1) / 1.5, dntmp);
   uptmp = vertcat(up(2) + 0.5 * up(1) / 1.5, uptmp);
   dntmp = vertcat(dntmp, dn(M-1) + 0.5 * dn(M) / 1.5);
   uptmp = vertcat(uptmp, up(M-1) + 0.5 * up(M) / 1.5);

   % Rewrite the arrays.
   dn = vertcat(dn(1), dntmp(2:M-1), dn(end));
   up = vertcat(up(1), uptmp(2:M-1), up(end));

   % Compute the net flux at each level, ensuring Q(1) recovers the total flux.
   Q = vertcat(min(-Q0 * (1 - albedo), up(1) - dn(1)), ...
      (up(2:M-1) + up(3:M)) ./ 2.0 - (dn(2:M-1) + dn(3:M)) ./ 2.0);

   % Compute the absorbed flux at each level
   dQ = Q(1:end-1) - Q(2:end);
   dQ = transpose(sum(reshape(dQ, dz_therm(1) / dz_spect(1), []), 1));
   dQ = vertcat(dQ, zeros((sum(dz_therm) - sum(dz_spect)) / dz_therm(1), 1));

   % Compute the flux absorbed in the top layer, to be allocated to the SEB.
   if albedo > 0.65 % && snowd > 0.05
      chi = 0.9;
   else
      chi = dQ(1) / sum(dQ);
   end

   % Compute the source term (absorbed flux per unit volume), Sc = dQ/dz.
   Sc = -(1.0 - chi) * Qsi / Q0 * dQ ./ dz_therm; % [W m-3]
end
