function metadata = marRefreezeMetadata(T, metadata)
   %MARREFREEZEMETADATA Describe signed native MAR RZ values in one artifact.
   %
   %  metadata = icemodel.forcing.helpers.marRefreezeMetadata(T)
   %  metadata = ... marRefreezeMetadata(T, metadata)
   %
   % MAR RZ is a native daily combined refreezing/deposition term. Rare
   % negative values are source-real and must not be clipped or relabelled as
   % pure refreezing. This helper stamps a canonical signed-policy token plus
   % cadence-independent strict-negative and material-negative statistics. The
   % 1e-8 mWE/h material threshold is reporting-only: it never clips, rounds,
   % accepts, rejects, or otherwise authorizes mutation of RZ.
   %
   % See also: icemodel.forcing.helpers.marDiagnosticMetadata

   arguments
      T timetable
      metadata (1, 1) struct = struct()
   end

   % Remove the superseded roundoff-only tolerance so one artifact cannot claim
   % both the old nonnegative policy and the current signed native contract.
   legacy = 'mar_diagnostic_refreeze_negative_tolerance_mwe_h';
   if isfield(metadata, legacy)
      metadata = rmfield(metadata, legacy);
   end

   % Count strict and material negatives on distinct UTC days rather than
   % postings so hourly Data and derived 15-minute met retain the same
   % source-native statistics.
   threshold = 1e-8;
   negative = false(height(T), 1);
   material_negative = false(height(T), 1);
   negative_minimum = NaN;
   material_minimum = NaN;
   if ismember("refreeze_deposition", string(T.Properties.VariableNames))
      values = T.refreeze_deposition;
      physical = isfinite(real(values)) & imag(values) == 0;
      negative = physical & real(values) < 0;
      material_negative = physical & real(values) < -threshold;
      if any(negative)
         negative_minimum = min(real(values(negative)));
      end
      if any(material_negative)
         material_minimum = min(real(values(material_negative)));
      end
   end
   % Use the row-time property instead of assuming the dimension is named Time;
   % valid timetable artifacts may preserve a source-specific row-dimension name.
   % Normalize the instants before grouping so the recorded statistic is UTC-day
   % based even when a valid in-memory timetable carries another display zone.
   times = T.Properties.RowTimes;
   times.TimeZone = 'UTC';
   negative_days = unique(dateshift(times(negative), 'start', 'day'));
   material_negative_days = unique(dateshift( ...
      times(material_negative), 'start', 'day'));

   metadata.mar_diagnostic_refreeze_deposition_sign = ...
      "native_signed_combined_term_preserved_not_pure_refreeze";
   metadata.mar_diagnostic_refreeze_negative_day_count = ...
      numel(negative_days);
   metadata.mar_diagnostic_refreeze_negative_minimum_mwe_h = negative_minimum;
   metadata.mar_diagnostic_refreeze_negative_statistics_basis = ...
      "distinct_utc_days_and_minimum_finite_negative_artifact_value";
   metadata.mar_diagnostic_refreeze_material_negative_threshold_mwe_h = ...
      threshold;
   metadata.mar_diagnostic_refreeze_material_negative_day_count = ...
      numel(material_negative_days);
   metadata.mar_diagnostic_refreeze_material_negative_minimum_mwe_h = ...
      material_minimum;
   metadata.mar_diagnostic_refreeze_material_negative_statistics_basis = ...
      "distinct_utc_days_and_minimum_below_declared_reporting_threshold";
end
