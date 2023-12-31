function [Sc, chi] = UPDATEEXTCOEFS(Qsi, albedo, I0, dz, z_spect, dz_spect, ...
      z_therm, dz_therm, spect_N, spect_S, solardwavl, ro_sno)
   %UPDATEEXTCOEFS Update the bulk extinction coefficients

   if Qsi < 1e-3
      Sc = 0.0 * dz;
      chi = 1.0;
      return
   end
   
   debug = true;
   debugplot = false;

   z_walls = [0; cumsum(dz_spect); sum(dz_spect) + dz_spect(end)];

   % Transform the mass density to the spectral grid
   ro_sno = interp1(z_therm, ro_sno, z_spect, 'nearest', 'extrap');
   ro_sno = max(ro_sno, 300);

   % Compute the downward bulk extinction coefficient.
   bulkcoefs = -log((sum(solardwavl .* exp(spect_S .* ro_sno), 2)) ...
      ./ (sum(solardwavl .* exp(spect_N .* ro_sno), 2))) ./ dz_spect(1);
   % NOTE: solardwavl comes multiplied by dwavl

   % Cast in a form that can be used in the two-stream computation (add the
   % boundaries).  Here I have assumed that it is okay to call the boundaries
   % equal to the value at the center of that grid cell (the value prevails
   % thoughout the cell).
   bulkcoefs = vertcat(bulkcoefs, bulkcoefs(end), bulkcoefs(end));

   % Compute the a and r coefficients from the surface albedo and the
   % extinction coefficients.
   a = (1.0 - albedo) / (1.0 + albedo) * bulkcoefs;
   r = 2.0 * albedo * bulkcoefs ./ (1.0 - albedo ^ 2);

   % Solve the system of equations for the two-stream model
   M = numel(z_spect) + 2; % notation here roughly follows Schlatter

   % initialize the matrix
   e = zeros(M-1, 1);
   f = zeros(M-1, 1);
   g = zeros(M-1, 1);
   b = zeros(M-1, 1);

   % Account for the upper boundary condition
   alfa = 1.0 / (a(1) + r(1));
   e(1) = 0.0;
   f(1) = 1.0;
   g(1) = -alfa / (dz_spect(1) + alfa);
   b(1) = r(1) * I0 * dz_spect(1) * alfa / (dz_spect(1) + alfa);

   % Fill the vectors between the boundaries
   deltaz = z_walls(3:M) - z_walls(2:M-1);
   e(2:M-1) = 1.0 + (r(3:M) - r(1:M-2)) ./ (4.0 * r(2:M-1));
   f(2:M-1) = deltaz ./ (2.0 * r(2:M-1)) ...
      .* (a(2:M-1) .* (r(3:M) - r(1:M-2)) - r(2:M-1) .* (a(3:M) - a(1:M-2))) ...
      - (2.0 + deltaz .^ 2 .* bulkcoefs(2:M-1) .^ 2);
   g(2:M-1) = 1.0 - (r(3:M) - r(1:M-2)) ./ (4.0 * r(2:M-1));
   b(2:M-1) = 0.0;

   % Account for the lower boundary condition
   g(M-1) = 0.0;
   b(M-1) = 0.0;
   
   %% TEST BEGIN
   
   % % e(M-1) = 0.0;
   % % f(M-1) = 1.0;
   % 
   % % keep the og ones
   % e1 = e;
   % f1 = f;
   % g1 = g;
   % b1 = b;
   % deltaz1 = deltaz;
   % 
   % % It should be possible to just compute the interior equations then append the
   % % upper and lower
   % deltaz = [dz_spect; dz_spect(end)];
   % e = 1.0 + (r(3:M) - r(1:M-2)) ./ (4.0 * r(2:M-1));
   % f = deltaz(2:M-1) ./ (2.0 * r(2:M-1)) ...
   %    .* (a(2:M-1) .* (r(3:M) - r(1:M-2)) - r(2:M-1) .* (a(3:M) - a(1:M-2))) ...
   %    - (2.0 + deltaz(2:M-1) .^ 2 .* bulkcoefs(2:M-1) .^ 2);
   % g = 1.0 - (r(3:M) - r(1:M-2)) ./ (4.0 * r(2:M-1));
   % b = 0 * g;
   % 
   % % Append the upper boundary condition
   % e = vertcat(0, e);
   % f = vertcat(1, f);
   % g = vertcat(-alfa / (dz_spect(1) + alfa), g);
   % b = vertcat(r(1) * I0 * dz_spect(1) * alfa / (dz_spect(1) + alfa), b);
   % 
   % % Account for the lower boundary condition
   % e(end) = 0.0;
   % f(end) = 1.0;
   % g(end) = 0.0;
   % 
   % %
   % isequal(e, e1) % no
   % isequal(f, f1) % no
   % isequal(g, g1) % yes
   % isequal(b, b1) % yes
   % 
   % isequal(e(1:end-1), e1(1:end-1)) % yes
   % isequal(f(1:end-1), f1(1:end-1)) % yes
   % 
   % e(end)
   % e1(end)
   % 
   % f(end)
   % f1(end)

   %% TEST END

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
   dn = vertcat(I0, ...
      (a(2:M-1) + r(2:M-1)) ./ r(2:M-1) .* up(2:M-1) - (up(3:M) - up(1:M-2)) ...
      ./ (2.0 * deltaz .* r(2:M-1)), ...
      (a(M) + r(M)) / r(M) * up(M) - (up(M) - up(M-1)) ...
      / (deltaz(1) * r(M))); % note: dz(M-1) for the bottom

   if debug == true; UP(:, 1) = up; DN(:, 1) = dn; end
   
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

   if debug == true; UP(:, 2) = up; DN(:, 2) = dn; end
   
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

   if debug == true; UP(:, 3) = up; DN(:, 3) = dn; end
   
   % Compute the net flux at each level, ensuring Q(1) recovers the total flux.
   Q = vertcat(min(-I0 * (1 - albedo), up(1) - dn(1)), ...
      (up(2:M-1) + up(3:M)) ./ 2.0 - (dn(2:M-1) + dn(3:M)) ./ 2.0);
   
   % Compute the absorbed flux at each level
   dQ = Q(1:end-1) - Q(2:end);
   
   if debug == true; dQspect = dQ; end
   
   % Upscale to the thermal grid.
   dQ = transpose(sum(reshape(dQ, dz_therm(1) / dz_spect(1), []), 1));
   
   % Compute the flux absorbed in the top layer, to be allocated to the SEB.
   % Optionally enhance it by a parameter to account for impurities on the ice
   % surface or other factors that enhance the surface absorption. 
   if albedo > 0.65 % && snowd > 0.05
      chi = 0.9;
   else
      chi = dQ(1) / sum(dQ);
   end
   
   % NEED TO ADD A 
   dQ = vertcat(dQ, zeros((sum(dz) - sum(dz_spect)) / dz_therm, 1));
   
   % Compute the source term (absorbed flux per unit volume), Sc = dQ/dz.
   Sc = -(1.0 - chi) * Qsi / I0 * dQ ./ dz; % [W m-3]
   
   if debug == true
      
      % sum(dQspect(1:40)) / sum(dQspect)
      % dQ(1) / sum(dQ)
      
      % This compares the different up/dn fluxes
      if debugplot == true
         icemodel.plot.twostream(UP, DN, Q, dQspect, z_walls, z_spect, false);
      
         % This compares the absorbed flux on the spectral grid to the thermal grid
         X1 = -dQspect ./ dz_spect;
         X2 = -dQ ./ dz;
         Y1 = z_spect;
         Y2 = z_therm;
         
         figure
         semilogx(X1, Y1, '-o'); hold on;
         semilogx(X2, Y2, '-o');
         formatPlotMarkers("markersize", 6)
         xlabel('dQ/dz')
         ylabel('depth')
         legend('spectral grid', 'thermal grid')
         set(gca, 'YDir', 'reverse', 'XDir', 'reverse')
         
         ratio = dz_therm(1)/dz_spect(1);
         itherm = 20;
         ispect = ratio * itherm;
         xlims = [min([X1(1:ispect); X2(1:itherm)]) max([X1(1:ispect); X2(1:itherm)])];
         ylims = [min(Y1(1), Y2(1)) max(Y2(itherm), Y1(ispect))];
         axis([0.98*xlims(1) 1.02*xlims(2) ylims(1) ylims(2)])
         
         
         figure
         semilogx(Sc, Y2, '-o'); hold on;
         semilogx(X2, Y2, '-o');
         formatPlotMarkers("markersize", 6)
         xlabel('dQ/dz')
         ylabel('depth')
         legend('S_c', 'dQ/dz')
         set(gca, 'YDir', 'reverse', 'XDir', 'reverse')
         
         close all
      end
   end
   
end
 