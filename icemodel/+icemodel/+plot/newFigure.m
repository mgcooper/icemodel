function fig = newFigure(kwargs)
   %NEWFIGURE Create a white, export-ready figure with a stable pixel size.
   %
   %  fig = icemodel.plot.newFigure()
   %  fig = icemodel.plot.newFigure(width=1100, height=360, name="met")
   %
   % Role
   %  Builds the export-figure frame that every report and verification
   %  figure uses: a white background, hidden by default, and a stable
   %  pixel size so exported rasters are reproducible. Without the white
   %  background, a headless session exports on the default dark canvas.
   %  The caller must close the figure; pair it with onCleanup.
   %
   % See also: icemodel.plot.compareTimeseries, exportgraphics

   arguments
      kwargs.width (1, 1) double {mustBePositive} = 900
      kwargs.height (1, 1) double {mustBePositive} = 360
      kwargs.name (1, 1) string = ""
      kwargs.visible (1, 1) logical = false
   end

   % 'off'/'on' strings keep compatibility with the figure property API.
   state = 'off';
   if kwargs.visible
      state = 'on';
   end
   fig = figure('Name', char(kwargs.name), 'Visible', state, ...
      'Color', 'w', 'Position', [0 0 kwargs.width kwargs.height]);
end
