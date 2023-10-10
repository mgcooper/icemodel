function plotdiags

   % just a starter i copied out of the bottom of icemodel_region ... probably
   % better already in runoff//functions

   % compare the saved ice1 data with the simulated data here
   figure; scatterfit(ice1a.Tsfc,ice1.Tsfc)
   figure; scatterfit(ice1a.runoff,ice1.runoff)
   figure; plot(ice1a.Time,ice1a.runoff); hold on;
   plot(ice1.Time,ice1.runoff); legend('saved','new')
   figure; plot(ice1a.Time,ice1a.melt); hold on;
   plot(ice1.Time,ice1.melt); legend('saved','new')

   % compare the saved ice2 data with the simulated data here
   figure; plot(ice2a.Tice(1,:),ice2.Tice(1,:),'o');
   [max(ice2a.Tice(:)) min(ice2a.Tice(:))] % all zero
   figure; pcolor(ice2.Tice); shading flat; colorbar
   figure; pcolor(ice2a.Tice); shading flat; colorbar

   % cant compare saved ice1 with met b/c they don't share any common variables
end
