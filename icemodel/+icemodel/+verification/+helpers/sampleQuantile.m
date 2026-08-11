function value = sampleQuantile(values, probability)
   %SAMPLEQUANTILE Linear-interpolated sample quantile, no extra toolboxes.
   %
   %  value = icemodel.verification.helpers.sampleQuantile(values, probability)
   %
   % Sorts internally and returns NaN for empty input.
   %
   % Inputs
   %  values      - sample values, any shape
   %  probability - quantile in [0, 1]
   %
   % Outputs
   %  value - the interpolated quantile, or NaN when values is empty

   values = sort(double(values(:)));
   if isempty(values)
      value = NaN;
      return
   end
   position = 1 + (numel(values) - 1) * probability;
   lower_index = floor(position);
   upper_index = ceil(position);
   fraction = position - lower_index;
   value = values(lower_index) * (1 - fraction) ...
      + values(upper_index) * fraction;
end
