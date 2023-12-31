function twostream(up, dn, Q, dQ, z_walls, z_spect, zoom)

   if nargin < 7
      zoom = true;
   end
   zidx = 4; % this controls which index is used as the maximum zoom in plot depth

   %% Compute Q and dQ directly from up/dn for comparison with Q/dQ passed in

   % If up/dn have multiple columns, it's assumed the first one has not had any
   % centered averaging, the second was averaged once, the third twice.
   % Regardless of that, it's assumed that the centered averaging was applied to
   % up/dn when Q was computed and thus dQ. Therefore Q has effectively had
   % three centered averages applied to it. Here, Q and dQ are re-computed from
   % the columns of up/dn, using the value of Q(1) passed in to ensure the top
   % value recovers the incoming net shortwave. Thus by plotting Q and dQ passed
   % in along with the last column of Q0 and dQ0, we compare the effect of
   % applying the centered averaging when computing Q. If instead the first
   % column of Q0 and dQ0 is plotted, we compare the effect of not doing any
   % centered averaging.

   % The up(1) - dn(1) correction is essential to not get a negative dQ/dz in
   % the top layer, whereas the centered difference smoothing is a different
   % matter, so I brought this in so the comparison is between smoothing vs not
   Q0 = vertcat(repmat(Q(1), 1, size(up, 2)), up(2:end-1, :) - dn(2:end-1, :));
   dQ0 = -diff(Q0);

   % This would be the simplest method if not correcting for Q(1)
   % Q0 = up(1:end-1, :) - dn(1:end-1, :);

   %% Plot

   % If up/dn have multiple series, need to decide which one to plot. This uses
   % the last one. To see the improvement over the first one, change to idx = 1.
   idx = size(up, 2);

   % up / down
   figure
   subplot(1, 3, 1)
   semilogx(up(:, idx), z_walls, '-o'); hold on;
   semilogx(dn(:, idx), z_walls, '-o');
   formatPlotMarkers("markersize", 6)
   xlabel('up/dn flux')
   ylabel('depth')
   legend('up', 'down')
   set(gca, 'YDir', 'reverse', 'XDir', 'reverse')
   axis tight

   if zoom == true
      axis([0.98*up(zidx, idx) 1.02*dn(1, idx) z_walls(1) z_walls(zidx)])
   else
      axis tight
   end

   % Q
   subplot(1, 3, 2)
   semilogx(-Q0(:, idx), z_walls(1:end-1), '-o'); hold on;
   semilogx(-Q, z_walls(1:end-1), '-o');
   formatPlotMarkers("markersize", 6)
   xlabel('Q = up - dn')
   ylabel('depth')
   legend('uncorrected', 'corrected')
   set(gca, 'YDir', 'reverse', 'XDir', 'reverse')

   if zoom == true
      axis([0.98*min(-Q(1:zidx)) 1.02*max(-Q(1:zidx)) z_walls(1) z_walls(4)])
   else
      axis tight
   end

   % dQ
   subplot(1, 3, 3)
   semilogx(-dQ0(:, idx), z_spect, '-o'); hold on;
   semilogx(-dQ, z_spect, '-o');
   formatPlotMarkers("markersize", 6)
   xlabel('dQ')
   ylabel('depth')
   legend('uncorrected', 'corrected')
   set(gca, 'YDir', 'reverse', 'XDir', 'reverse')

   if zoom == true
      axis([0.98*min(-dQ(1:zidx)) 1.02*max(-dQ(1:zidx)) z_spect(1) z_spect(zidx)])
   else
      axis tight
   end
end
