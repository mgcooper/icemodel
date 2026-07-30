function [Data, metadata, diagnostics] = applyMerraTimeSupport(Data, metadata)
   %APPLYMERRATIMESUPPORT Enforce MERRA interval-start/support semantics.
   %
   %  [Data, metadata, diagnostics] = ...
   %     icemodel.forcing.helpers.applyMerraTimeSupport(Data, metadata)
   %
   % MERRA tavg1 channels are already hourly means after the builder's native
   % center-to-start relabel. The glacier collection is tavg3, so its RUNOFF,
   % SNICEALB, SNOWDP_GL, and SNOMAS_GL derivatives must be held over each
   % three-hour interval rather than linearly interpolated. Exact source rows at
   % 00/03/.../21 UTC plus durable raw-grid proof recover legacy staged artifacts
   % without reopening NetCDF in this helper. The helper also stamps the
   % reader/application provenance required by QA.
   %
   % Inputs
   %  Data     - hourly MERRA Data timetable on clock-hour interval starts
   %  metadata - existing artifact metadata to preserve and extend
   %
   % Outputs
   %  Data        - support-corrected timetable
   %  metadata    - existing fields plus canonical MERRA time provenance
   %  diagnostics - replaced element count and metadata-change flag

   arguments
      Data timetable
      metadata (1, 1) struct = struct()
   end

   % Cached hourly artifacts must be regular application products. A native
   % center-stamped or clipped non-hourly input cannot be repaired reliably from
   % its derived table alone.
   if height(Data) > 1 && any(diff(Data.Time) ~= hours(1))
      error('icemodel:forcing:applyMerraTimeSupport:irregularTime', ...
         'MERRA application Data must have a uniform hourly axis')
   end
   if ~isempty(Data.Time) ...
         && any(minute(Data.Time) ~= 0 | second(Data.Time) ~= 0)
      error('icemodel:forcing:applyMerraTimeSupport:notIntervalStart', ...
         'MERRA application Data must use clock-hour interval starts')
   end

   % A current clipped artifact may begin at 01:00/02:00 inside a tavg3 support
   % block whose source row is outside the saved window. Its durable time and
   % native-grid contracts are sufficient proof; preserve every value instead of
   % trying to reconstruct it.
   diagnostics = struct('replaced_count', 0, 'metadata_changed', false);
   if icemodel.forcing.helpers.hasCanonicalMerraTimeSupport(metadata) ...
         && icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid( ...
         Data, metadata) ...
         && icemodel.forcing.helpers.hasConstantMerraTavg3Support(Data)
      Data.Properties.UserData = metadata;
      return
   end

   % Recover the exact tavg3 source rows and hold each over its declared support.
   % Missing channels are valid for reduced fixtures; their provenance remains
   % useful for the tavg1 channels that are present.
   channels = intersect(["runoff", "albedo", "snowd", "swe"], ...
      string(Data.Properties.VariableNames), 'stable');
   if ~isempty(channels) && ~isempty(Data.Time) ...
         && mod(hour(Data.Time(1)), 3) ~= 0
      error('icemodel:forcing:applyMerraTimeSupport:unrecoverableLeadingBlock', ...
         ['legacy MERRA window begins inside a tavg3 support block; the ' ...
          'leading source row is unavailable'])
   end
   if ~icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid(Data, metadata)
      error('icemodel:forcing:applyMerraTimeSupport:unprovenSourceGrid', ...
         ['legacy MERRA tavg3 repair requires an exact native glc timestamp ' ...
         'inventory; saved 3-hour rows alone cannot prove omitted source gaps'])
   end

   % Only inventory-proven source rows may seed a hold. Explicitly absent native
   % stamps leave their complete three-hour support missing.
   source_rows = mod(hour(Data.Time), 3) == 0;
   if ~isempty(channels)
      missing = metadata.merra_tavg3_missing_source_times;
      source_rows = source_rows & ~ismember(Data.Time, missing);
   end
   for channel = reshape(channels, 1, [])
      values = Data.(char(channel));
      corrected = nan(size(values));
      source_time = Data.Time(source_rows);
      source_values = values(source_rows, :);
      for offset = 0:2
         [present, target_rows] = ismember(source_time + hours(offset), Data.Time);
         corrected(target_rows(present), :) = source_values(present, :);
      end
      diagnostics.replaced_count = diagnostics.replaced_count ...
         + nnz(differentValues(values, corrected));
      Data.(char(channel)) = corrected;
   end

   % Preserve all unrelated provenance while making both source and application
   % timing decisions durable and machine-checkable.
   prior = metadata;
   metadata.merra_source_time_coordinate = 'native_at_reader';
   metadata.merra_time_relabel_policy = ...
      'time_averaged_center_to_interval_start';
   metadata.merra_time_upsample_policy = ...
      'zero_order_hold_over_declared_support';
   metadata.merra_collection_support_hours = ...
      struct('slv', 1, 'rad', 1, 'flx', 1, 'glc', 3);
   if diagnostics.replaced_count == 0 ...
         && isfield(prior, 'merra_time_replaced_value_count')
      metadata.merra_time_replaced_value_count = ...
         prior.merra_time_replaced_value_count;
   else
      metadata.merra_time_replaced_value_count = diagnostics.replaced_count;
   end
   diagnostics.metadata_changed = ~isequaln(prior, metadata);
   Data.Properties.UserData = metadata;
end

function different = differentValues(a, b)
   %DIFFERENTVALUES Elementwise inequality that treats paired NaNs as equal.
   different = ~(a == b | (isnan(a) & isnan(b)));
end
