function temperature(T, dz, z_nodes)
   %ICEMODEL.PLOT.TEMPERATURE
   %
   % Inputs
   %  T - vector of temperatures
   %  dz - vector of cv thicknesses
   %  z_nodes - optional precomputed thermal node coordinates

   % Reuse the caller-provided thermal node coordinates when available.
   if nargin < 3 || isempty(z_nodes)
      z_nodes = cumsum(dz) - dz / 2;
   end

   figure;
   plot(T, z_nodes); set(gca, 'YDir', 'Reverse')
   xlabel('Temperature [K]')
   ylabel('Depth [m]')

   % This was in a random file. It is the plot of subsurface T profile with
   % open circles for the skin and air temperature.
   %
   % figure
   % plot(T, cumsum(dz)); hold on
   % plot(Ta, 0, 'o')
   % plot(Ts, 0, 'o')
   % set(gca,'YDir','reverse');
   %
   % figure; plot(1:numel(numiter), numiter)
   % figure; plot(1:numel(Tflag), Tflag)
   % figure; plot(ice1.chf)
   % figure; plot(ice1.balance)
end
