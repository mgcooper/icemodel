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
   [f_res, f_por] = snowphysics.residualWater(f_ice);
   
   [check.f_ice, check.f_liq, check.f_air, check.xice, check.xliq] = ...
      VOLBAL(f_ice, f_liq, 1.0 - f_ice - f_liq, f_res, ones(size(f_ice)));
   
   for n = 1:numel(f_ice)
      if check.xice(n) > 0.0
         fprintf('Layer %d extra ice: %.3f \n', n, check.xice(n))
      end
      if check.xliq(n) > 0.0
         fprintf('Layer %d extra liq: %5.3f \n', n, check.xliq(n))
      end
   end
   
   [min(check.f_ice) max(check.f_ice)]
   [min(check.f_liq) max(check.f_liq)]
   [min(check.f_air) max(check.f_air)]

   % compute the melt equivalent of a given temperature error
   max(cp_sno * 1e-13 / Lf * dz(1)) > eps

   
   
end
