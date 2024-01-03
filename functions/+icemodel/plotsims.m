function plotsims(ice1a, ice1b, fields)

   % Compare ice1 data from two simulations.
   %
   % just a starter i copied out of the bottom of icemodel_region ... probably
   % better already in runoff//functions
   isskinmodel = false(2, 1);
   if isfile(ice1a)
      if contains(ice1a, 'skinmodel')
         isskinmodel(1) = true;
      end
      try
         load(ice1a, 'ice1')
         ice1a = ice1; clear ice1;
      catch e
         rethrow(e)
      end
   end

   if isfile(ice1b)
      if contains(ice1b, 'skinmodel')
         isskinmodel(2) = true;
      end
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

   %figure('Position', [100   100   900   450])
   figontop
   for n = 0:numel(fields)
      
      % 1:1 plots
      subtight(1, 3, 1, "gap", [0.08 0.08], "marw", [0.06 0.02], "marh", [0.08 0.04])
      try
         % Custom melt-freeze
         if n == 0 
            if ~all(isskinmodel)
               plot(ice1a.melt-ice1a.freeze, ice1b.melt-ice1b.freeze, 'o');
               title('melt - freeze')
            else
               continue
            end
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
end
