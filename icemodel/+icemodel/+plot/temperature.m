function temperature(T, dz)
   % need this in a function, and an anonymous one to create z_therm:

   znodes = @(dz) cumsum(dz) - dz / 2;

   figure;
   plot(T, znodes(dz)); set(gca, 'YDir', 'Reverse')


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
