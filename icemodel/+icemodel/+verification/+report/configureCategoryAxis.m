function configureCategoryAxis(ax, labels)
   %CONFIGURECATEGORYAXIS Label horizontal evidence rows without clipping.

   yticks(ax, 1:numel(labels))
   yticklabels(ax, icemodel.verification.report.safeLabel(labels))
   ax.YDir = 'reverse';
   ax.TickLabelInterpreter = 'none';
   icemodel.verification.report.formatReportAxes(ax)
end
