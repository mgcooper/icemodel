function plotdiags(ice1a, ice1b, fields)

   % Compare ice1 data from two simulations.
   %
   % just a starter i copied out of the bottom of icemodel_region ... probably
   % better already in runoff//functions

   if isfile(ice1a)
      try
         load(ice1a, 'ice1')
         ice1a = ice1; clear ice1;
      catch e
         rethrow(e)
      end
   end

   if isfile(ice1b)
      try
         load(ice1b, 'ice1')
         ice1b = ice1; clear ice1;
      catch e
         rethrow(e)
      end
   end


   if nargin < 3
      fields = unique( ...
         [ice1a.Properties.VariableNames, ice1b.Properties.VariableNames] ...
         );
   end

   name1 = inputname(1);
   name2 = inputname(2);
   
   ice1a = rmleapinds(ice1a);
   ice1b = rmleapinds(ice1b);

   % Skip the first 1000 values when plotting timeseries 
   idx = 1000:height(ice1a);

   figure('Position', [100   100   900   450])
   for n = 0:numel(fields)
      
      % 1:1 plots
      subtight(1, 3, 1, "gap", [0.08 0.08], "marw", [0.06 0.02], "marh", [0.08 0.04])
      try
         % Custom melt-freeze
         if n == 0
            plot(ice1a.melt-ice1a.freeze, ice1b.melt-ice1b.freeze, 'o');
            title('melt - freeze')
         else
            plot(ice1a.(fields{n}), ice1b.(fields{n}), 'o');
            title((fields{n}))
         end
          addOnetoOne
      catch
      end
      xlabel(name1)
      ylabel(name2)
   
      % Timeseries plots
      subplot(1, 3, [2 3])
      hold on
      try
         % Custom melt-freeze
         if n == 0
            plot(ice1a.Time(idx), ice1a.melt(idx) - ice1a.freeze(idx), '--');
            plot(ice1b.Time(idx), ice1b.melt(idx) - ice1b.freeze(idx), ':');
            ylabel('melt - freeze')
         else
            plot(ice1a.Time(idx), ice1a.(fields{n})(idx), '--');
            plot(ice1b.Time(idx), ice1b.(fields{n})(idx), ':');
            ylabel((fields{n}))
         end
         legend(name1, name2)
      catch
      end
      hold off
      pause
      clf
   end


   % % compare the saved ice2 data with the simulated data here
   % figure
   % plot(ice2a.Tice(1,:), ice2.Tice(1,:),'o');
   % [max(ice2a.Tice(:)) min(ice2a.Tice(:))] % all zero
   % 
   % figure
   % pcolor(ice2.Tice); 
   % shading flat
   % colorbar
   % 
   % figure
   % pcolor(ice2a.Tice); 
   % shading flat
   % colorbar

   % Below here various debugging stuff 
   
   % % cant compare saved ice1 with met b/c they don't share any common variables
   % [f_res, f_por] = snowphysics.residualWater(f_ice);
   %
   % [check.f_ice, check.f_liq, check.f_air, check.xice, check.xliq] = ...
   %    VOLBAL(f_ice, f_liq, 1.0 - f_ice - f_liq, f_res, ones(size(f_ice)));
   %
   % for n = 1:numel(f_ice)
   %    if check.xice(n) > 0.0
   %       fprintf('Layer %d extra ice: %.3f \n', n, check.xice(n))
   %    end
   %    if check.xliq(n) > 0.0
   %       fprintf('Layer %d extra liq: %5.3f \n', n, check.xliq(n))
   %    end
   % end
   %
   % [min(check.f_ice) max(check.f_ice)]
   % [min(check.f_liq) max(check.f_liq)]
   % [min(check.f_air) max(check.f_air)]
   %
   % % compute the melt equivalent of a given temperature error
   % max(cp_sno * 1e-13 / Lf * dz(1)) > eps

end
