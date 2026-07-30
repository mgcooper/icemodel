function layers = methodFillLayers(values, provenance, own_method_mask)
   %METHODFILLLAYERS Split one channel into observed/own-fill/other-fill layers.
   %
   %  layers = icemodel.verification.report.methodFillLayers( ...
   %     values, provenance, own_method_mask)
   %
   % Role
   %  The color logic behind the POLICY D-31 method detail panels: the
   %  filled product's per-sample provenance channel identifies which
   %  samples are reconstructions, and the plan audit's segment spans for
   %  the panel's method (OWN_METHOD_MASK) attribute each reconstructed
   %  sample to the panel's own method or to some other method. The
   %  caller plots the three returned layers so a panel accents only its
   %  own method and renders every foreign fill in the muted context
   %  color — a bounded_interp panel can never display a MAR fill in the
   %  same accent.
   %
   % Inputs
   %  values : channel sample values from the filled product.
   %  provenance : matching per-sample codes from the paired
   %     <channel>_provenance column (reconstruct.provenanceCodes()).
   %  own_method_mask : logical mask of samples inside the panel
   %     method's audited fill segments.
   %
   % Returns
   %  layers : struct —
   %     observed, own_fill, other_fill : VALUES with samples outside
   %        the layer replaced by NaN, ready to plot as separate lines.
   %     observed_color, own_color, other_color : the layer RGB rows
   %        from icemodel.verification.report.gapfillFigureStyle.
   %
   % See also: icemodel.verification.report.gapfillFigureStyle,
   %  icemodel.forcing.reconstruct.provenanceCodes

   arguments
      values (:, 1) double
      provenance (:, 1) {mustBeNumeric}
      own_method_mask (:, 1) logical
   end
   if numel(provenance) ~= numel(values) ...
         || numel(own_method_mask) ~= numel(values)
      error('icemodel:report:methodFillLayers:sizeMismatch', ...
         'values, provenance, and own_method_mask must share one axis');
   end

   % Raw-fallback and clamped shortwave (codes 13/14) are builder-selected
   % raw MEASUREMENTS (POLICY A7), so they belong to the observed layer
   % exactly as observedNative and the overview treat them.
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   is_observed = ismember(provenance, double([codes.observed, ...
      codes.raw_shortwave, codes.clamped_shortwave]));
   % A fill is a finite, stamped reconstruction: never an observed sample
   % and never the explicit missing sentinel.
   is_fill = isfinite(values) & ~is_observed ...
      & provenance ~= double(codes.missing);

   % Each layer keeps the shared time axis and NaN-masks the other layers
   % so the caller can overlay all three without reindexing.
   observed = values;
   observed(~is_observed) = NaN;
   own_fill = values;
   own_fill(~(is_fill & own_method_mask)) = NaN;
   other_fill = values;
   other_fill(~(is_fill & ~own_method_mask)) = NaN;

   % Colors ride along from the single style registry so the plotting
   % call sites and the layer split can never disagree about the accent.
   style = icemodel.verification.report.gapfillFigureStyle();
   layers = struct('observed', observed, 'own_fill', own_fill, ...
      'other_fill', other_fill, 'observed_color', style.observed, ...
      'own_color', style.accent, 'other_color', style.context);
end
