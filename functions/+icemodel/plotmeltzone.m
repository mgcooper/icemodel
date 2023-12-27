function plotmeltzone(ice2, iter, z_therm, TL, TH, zdepth, numsteps)

   if nargin < 6
      zdepth = 2;
   end
   if nargin < 7
      numsteps = 8;
   end
   
   zi = find(zdepth - z_therm < 0, 1, 'first');

   % Prepare data vertices
   if isstruct(ice2) || istimetable(ice2)
      Tice = ice2.Tice;
   elseif ismatrix(ice2)
      Tice = ice2;
   end
   
   X = Tice(1:zi, max(1, iter-1));
   Y = z_therm(1:zi);

   % Prepare melt zone vertices
   meltzoneT = (TL + TH)/2;
   meltzoneW = (TH - TL)/2;
   meltzoneX = [meltzoneT, meltzoneT];
   meltzoneY = [0 zdepth];
   meltzoneE = [meltzoneW, meltzoneW];
   
   figure; hold on

   % Plot the melt zone
   fillplot(meltzoneX, meltzoneY, meltzoneE, rgb('grey'), 'x', 'FaceAlpha', 0.7)

   % Plot the data
   plot(X, Y, '-o');
   set(gca, 'YDir', 'reverse');
   % legend(var, 'Location', 'best');
   xlabel('T [K]');
   ylabel('Depth [m]');
   
   for n = 1:numsteps
      plot(Tice(1:zi, max(1, iter-n+1)), Y, '-o');
   end
   formatPlotMarkers('markersize', 6)
   ylim([0 zdepth])
end
