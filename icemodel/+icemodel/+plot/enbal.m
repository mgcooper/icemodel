function varargout = enbal(ice1, met, varargin)
   %ENBAL Plot the surface energy balance.
   %
   %  H = icemodel.plot.enbal(ice1, met)
   %  H = icemodel.plot.enbal(ice1, met, ice1b)
   %  H = icemodel.plot.enbal(ice1, met, ice1b, 'labels', {'A', 'B'})
   %
   % Plots a 7-variable surface energy balance timeseries (7-day moving
   % mean of daily averages). Each variable is assigned a distinct color;
   % line style encodes the data source:
   %
   %   Observations (met):   solid  (–)
   %   Model A (ice1):       dashed (- -)   [dotted when ice1b is absent]
   %   Model B (ice1b):      dotted (···)
   %
   % Optional arguments:
   %   ice1b   - second model timetable to overlay (default: none)
   %   labels  - 1×2 cell of scheme names printed in the figure title
   %             when ice1b is provided (default: {'model A', 'model B'})

   import icemodel.helpers.rmttleapinds

   parser = inputParser;
   parser.FunctionName = mfilename;
   parser.addRequired('ice1', @istimetable);
   parser.addRequired('met', @istimetable);
   parser.addOptional('ice1b', [], @(x) isempty(x) || istimetable(x));
   parser.addParameter('labels', {}, @iscell);
   parser.addParameter('grouptitle', '', @ischar);
   parser.parse(ice1, met, varargin{:});

   ice1b = parser.Results.ice1b;
   labels = parser.Results.labels;
   grouptitle = parser.Results.grouptitle;
   two_model = ~isempty(ice1b);

   % Retime to daily and remove Feb 29.
   ice1 = rmttleapinds(retime(ice1, 'daily', 'mean'));
   met = rmttleapinds(retime(met, 'daily', 'mean'));
   if two_model
      ice1b = rmttleapinds(retime(ice1b, 'daily', 'mean'));
   end

   % 7-day moving-mean window.
   win = 7;

   % Variables and colors (replicates figure 12 in van As).
   vars = {'swd', 'swu', 'lwd', 'lwu', 'shf', 'lhf', 'netr'};
   colors = [
      0.0000  0.0118  0.3569   % dark blue
      0.0118  0.2627  0.8745   % blue
      0.8980  0.0000  0.0000   % red
      0.9765  0.4510  0.0235   % orange
      0.0118  0.2078  0.0000   % dark green
      0.0824  0.6902  0.1020   % green
      0.4941  0.1176  0.6118]; % purple

   nv = numel(vars);
   h = gobjects(nv, 1);

   figure('Position', [100 100 1400 900]); hold on
   for n = nv:-1:1
      try
         % Observations: solid.
         h(n) = plot(met.Time, movmean(met.(vars{n}), win), '-', ...
            'Color', colors(n, :));
         if two_model
            % Scheme A: dashed; scheme B: dotted.
            plot(ice1.Time, movmean(ice1.(vars{n}), win), '--', ...
               'Color', colors(n, :));
            plot(ice1b.Time, movmean(ice1b.(vars{n}), win), ':', ...
               'Color', colors(n, :), 'LineWidth', 1.2);
         else
            % Single model: dotted (original behavior).
            plot(ice1.Time, movmean(ice1.(vars{n}), win), ':', ...
               'Color', colors(n, :));
         end
      catch
      end
   end

   set(gca, 'XAxisLocation', 'origin'); axis tight
   ltext = {'SW\downarrow', 'SW\uparrow', 'LW\downarrow', 'LW\uparrow', ...
      'SHF', 'LHF', 'Rnet'};
   legend(h, ltext, 'Orientation', 'horizontal', 'Location', 'southoutside');
   ylabel('heat flux (W m^{-2})');

   if two_model
      if numel(labels) == 2
         la = labels{1};
         lb = labels{2};
      else
         la = 'model A';
         lb = 'model B';
      end
      scheme_str = sprintf('%s (dashed)  vs  %s (dotted)  vs  obs (solid)', la, lb);
   else
      scheme_str = 'icemodel (dotted) versus forcings (solid)';
   end

   % Two-line title when a group name is provided; single line otherwise.
   if ~isempty(grouptitle)
      title({grouptitle, scheme_str}, 'Interpreter', 'none');
   else
      title(scheme_str, 'Interpreter', 'none');
   end

   if nargout > 0
      varargout{1} = h;
   end
end
