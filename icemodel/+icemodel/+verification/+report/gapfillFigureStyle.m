function style = gapfillFigureStyle()
   %GAPFILLFIGURESTYLE Single source of the gap-fill report figure colors.
   %
   %  style = icemodel.verification.report.gapfillFigureStyle()
   %
   % Role
   %  One named registry for the gap-fill report's shared rendering style
   %  (STYLE.local SSOT rule), so overview figures, method-detail panels,
   %  and tests cannot drift apart. POLICY D-31: a method detail panel
   %  renders ONLY its own method in the accent color; fills from other
   %  methods inside the context window render in the muted context color
   %  so the panel can never visually claim a foreign fill.
   %
   % Returns
   %  style : struct —
   %     observed : RGB for observed (native-truth) samples.
   %     accent   : RGB for the panel's own method (and the overview's
   %                filled overlay, which aggregates all methods).
   %     context  : muted RGB for other methods' fills inside a method
   %                detail panel's context window.
   %     max_overview_points : per-channel point budget for the eight-axis
   %                full-period overview. Detail figures remain unthinned.
   %
   % See also: icemodel.verification.report.methodFillLayers

   % The observed dark and accent orange preserve the report's original
   % identity; the muted grey sits far from both so the three layers stay
   % separable in print and on screen.
   style = struct( ...
      'observed', [0.25 0.25 0.30], ...
      'accent', [0.85 0.45 0.10], ...
      'context', [0.62 0.62 0.66], ...
      'max_overview_points', 20000);
end
