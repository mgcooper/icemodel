function [Sc, chi] = UPDATEEXTCOEFSINLINELEGACY(Qsi, albedo, Q0, dz_spect, ...
      spect_N, spect_S, solardwavl, Sc, dz_therm, ro_sno)
   %UPDATEEXTCOEFSINLINELEGACY Historical inlined spectral source-term path.
   %
   %  [Sc, chi] = UPDATEEXTCOEFSINLINELEGACY(Qsi, albedo, Q0, dz_spect, ...
   %     spect_N, spect_S, solardwavl, Sc, dz_therm, ro_sno)
   %
   % This preserves the pre-refactor UPDATEEXTCOEFS implementation so the
   % benchmark suite can compare it directly against the helper-call paths.
   %
   %#codegen

   if Qsi < 1e-3
      Sc = 0.0 * Sc;
      chi = 1.0;
      return
   end

   % Transform the mass density to the spectral grid.
   ro_sno = max(interp1(cumsum(dz_therm) - dz_therm / 2, ro_sno, ...
      cumsum(dz_spect) - dz_spect / 2, 'nearest', 'extrap'), 300);

   % Compute the downward bulk extinction coefficient and pad the lower
   % boundaries expected by the two-stream solve.
   bulkcoefs = -log((sum(solardwavl .* exp(spect_S .* ro_sno), 2)) ...
      ./ (sum(solardwavl .* exp(spect_N .* ro_sno), 2))) ./ dz_spect(1);
   bulkcoefs = vertcat(bulkcoefs, bulkcoefs(end), bulkcoefs(end));

   % Compute the two-stream coefficients directly from the extinction profile.
   a = (1.0 - albedo) / (1.0 + albedo) * bulkcoefs;
   r = 2.0 * albedo * bulkcoefs ./ (1.0 - albedo ^ 2);

   % Solve the inlined two-stream system.
   M = numel(dz_spect) + 2;
   e = zeros(M - 1, 1);
   f = zeros(M - 1, 1);
   g = zeros(M - 1, 1);
   b = zeros(M - 1, 1);

   alfa = 1.0 / (a(1) + r(1));
   e(1) = 0.0;
   f(1) = 1.0;
   g(1) = -alfa / (dz_spect(1) + alfa);
   b(1) = r(1) * Q0 * dz_spect(1) * alfa / (dz_spect(1) + alfa);

   deltaz = [dz_spect(2:end); dz_spect(end)];
   e(2:M-1) = 1.0 + (r(3:M) - r(1:M-2)) ./ (4.0 * r(2:M-1));
   f(2:M-1) = deltaz ./ (2.0 * r(2:M-1)) .* ...
      (a(2:M-1) .* (r(3:M) - r(1:M-2)) ...
      - r(2:M-1) .* (a(3:M) - a(1:M-2))) ...
      - (2.0 + deltaz .^ 2 .* bulkcoefs(2:M-1) .^ 2);
   g(2:M-1) = 1.0 - (r(3:M) - r(1:M-2)) ./ (4.0 * r(2:M-1));
   b(2:M-1) = 0.0;
   g(M-1) = 0.0;
   b(M-1) = 0.0;

   up = zeros(M - 1, 1);
   for k = 2:numel(f)
      f(k) = f(k) - e(k) / f(k-1) * g(k-1);
      b(k) = b(k) - e(k) / f(k-1) * b(k-1);
   end
   up(M-1) = b(M-1) / f(M-1);
   for k = numel(f)-1:-1:1
      up(k) = (b(k) - g(k) * up(k+1)) / f(k);
   end
   up = vertcat(up(1:M-1), 0);

   % Reconstruct and smooth the downward and upward flux profiles.
   dn = vertcat(Q0, ...
      (a(2:M-1) + r(2:M-1)) ./ r(2:M-1) .* up(2:M-1) ...
      - (up(3:M) - up(1:M-2)) ./ (2.0 * deltaz .* r(2:M-1)), ...
      (a(M) + r(M)) / r(M) * up(M) ...
      - (up(M) - up(M-1)) / (deltaz(M-2) * r(M)));
   [up, dn] = smoothFluxProfiles(up, dn, M);

   % Convert the net spectral flux to the absorbed source term.
   Q = vertcat(min(-Q0 * (1 - albedo), up(1) - dn(1)), ...
      (up(2:M-1) + up(3:M)) ./ 2.0 - (dn(2:M-1) + dn(3:M)) ./ 2.0);
   [Sc, chi] = SPECTRALSOURCETERM(Qsi, albedo, Q0, Q, dz_spect, dz_therm, ...
      Sc);
end

function [up, dn] = smoothFluxProfiles(up, dn, M)
   %SMOOTHFLUXPROFILES Reproduce the legacy UPDATEEXTCOEFS smoothing pass.

   dntmp = (dn(2:M-1) + 0.5 * (dn(1:M-2) + dn(3:M))) / 2.0;
   uptmp = (up(2:M-1) + 0.5 * (up(1:M-2) + up(3:M))) / 2.0;

   dntmp = vertcat((dn(2) + 0.5 * dn(1)) / 1.5, dntmp);
   uptmp = vertcat((up(2) + 0.5 * up(1)) / 1.5, uptmp);
   dntmp = vertcat(dntmp, (dn(M-1) + 0.5 * dn(M)) / 1.5);
   uptmp = vertcat(uptmp, (up(M-1) + 0.5 * up(M)) / 1.5);

   dn = vertcat(dn(1), dntmp(2:M-1), dn(end));
   up = vertcat(up(1), uptmp(2:M-1), up(end));

   dntmp = (dn(2:M-1) + 0.5 * (dn(1:M-2) + dn(3:M))) / 2.0;
   uptmp = (up(2:M-1) + 0.5 * (up(1:M-2) + up(3:M))) / 2.0;

   dntmp = vertcat(dn(2) + 0.5 * dn(1) / 1.5, dntmp);
   uptmp = vertcat(up(2) + 0.5 * up(1) / 1.5, uptmp);
   dntmp = vertcat(dntmp, dn(M-1) + 0.5 * dn(M) / 1.5);
   uptmp = vertcat(uptmp, up(M-1) + 0.5 * up(M) / 1.5);

   dn = vertcat(dn(1), dntmp(2:M-1), dn(end));
   up = vertcat(up(1), uptmp(2:M-1), up(end));
end
