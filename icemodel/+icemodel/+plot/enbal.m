function varargout = enbal(ice1, met)
   %ENBAL plot the surface energy balance.
   %
   %  H = ICEMODEL.PLOT.ENBAL(ICE1, MET)
   %
   % See also

   parser = inputParser;
   parser.FunctionName = mfilename;
   parser.addRequired('ice1', @istimetable);
   parser.addRequired('met', @istimetable);
   parser.parse(ice1, met);

   % Retime to daily
   met = retime(met, 'daily', 'mean');
   ice1 = retime(ice1, 'daily', 'mean');

   % Moving mean window size
   win = 7;

   % Replicate figure 12 in van As
   vars = {'swd', 'swu', 'lwd', 'lwu', 'shf', 'lhf', 'netr'};
   colors = [
      0    0.0118    0.3569   % dark blue
      0.0118    0.2627    0.8745   % blue
      0.8980         0         0   % red
      0.9765    0.4510    0.0235   % orange
      0.0118    0.2078         0   % dark green
      0.0824    0.6902    0.1020   % green
      0.4941    0.1176    0.6118]; % purple

   % Create the figure
   figure; hold on
   for n = numel(vars):-1:1
      try
         h(n) = plot(met.Time, movmean(met.(vars{n}), win), ...
            'Color', colors(n, :));
         plot(ice1.Time, movmean(ice1.(vars{n}), win), ':', ...
            'Color', colors(n, :));
      catch
      end
   end
   set(gca,'XAxisLocation', 'origin'); axis tight
   ltext = {'SW\downarrow', 'SW\uparrow', 'LW\downarrow', 'LW\uparrow', ...
      'SHF', 'LHF', 'Rnet'};
   legend(h, ltext, 'Orientation', 'horizontal', 'Location', 'southoutside')
   title('icemodel (dotted) versus forcings (solid)');
   ylabel('heat flux (W m^{-2})')

   if nargout == 1
      varargout{1} = h;
   end
end
