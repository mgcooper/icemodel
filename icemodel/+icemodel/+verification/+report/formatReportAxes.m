function formatReportAxes(ax)
   %FORMATREPORTAXES Isolate exported graphics from interactive theme defaults.

   ax.Color = 'w';
   ax.XColor = 'k';
   ax.YColor = 'k';
   ax.GridColor = [0.65 0.65 0.65];
   ax.GridAlpha = 0.25;
   ax.FontSize = 11;
   ax.Box = 'off';
   ax.Title.Color = 'k';
   ax.XLabel.Color = 'k';
   ax.YLabel.Color = 'k';
end
