function metadata = marDiagnosticMetadata(T, melt_daily_rate, metadata, kwargs)
   %MARDIAGNOSTICMETADATA Describe optional MAR mass-balance diagnostics.
   %
   %  metadata = ... marDiagnosticMetadata(T, melt_daily_rate)
   %  metadata = ... marDiagnosticMetadata(_, _, metadata, sector=1)
   %
   % T may carry the independently defined MAR diagnostics subl (SUH),
   % subl_evap (SU), and refreeze_deposition (RZ), all in mWE/h.
   % MELT_DAILY_RATE is the native daily ME value divided by 24 and held on
   % T.Time; it is used only to validate hourly MEH-derived melt and is not
   % promoted to a second public melt channel. An empty vector explicitly
   % records that the optional daily ME product was unavailable.
   %
   % SUH is hourly sublimation, whereas SU combines sublimation and
   % evaporation. RZ is a signed native combined meltwater-refreezing and
   % deposition term with rare source-real negatives. The metadata therefore
   % freezes distinct canonical names rather than claiming unsupported
   % equivalence with pure evaporation or nonnegative refreezing products, and
   % records both strict-negative statistics and a reporting-only material
   % subset without changing any RZ value.
   %
   % See also: icemodel.forcing.buildMarData,
   %  icemodel.forcing.helpers.dailyToHourly

   arguments
      T timetable
      melt_daily_rate (:, 1) double = zeros(0, 1)
      metadata (1, 1) struct = struct()
      kwargs.sector (1, 1) double {mustBeMember(kwargs.sector, [1 2])} = 1
      kwargs.melt_abs_tolerance_mwe_day (1, 1) double ...
         {mustBeNonnegative} = 1e-5
      kwargs.melt_rel_tolerance (1, 1) double {mustBeNonnegative} = 1e-4
      kwargs.abs_limit_mwe_h (1, 1) double {mustBePositive} = 0.05
   end

   % Infer availability from public channels plus the private daily-ME
   % reference. A source missing every optional field remains a valid reduced
   % source and is described as not_available rather than as forcing failure.
   names = string(T.Properties.VariableNames);
   public_names = ["subl", "subl_evap", "refreeze_deposition"];
   native_names = ["SUH", "SU", "RZ"];
   available = ismember(public_names, names);
   channels = public_names(available)';
   native = native_names(available)';
   has_daily_melt = ~isempty(melt_daily_rate);
   if has_daily_melt
      if numel(melt_daily_rate) ~= height(T)
         error('icemodel:forcing:marDiagnosticMetadata:heightMismatch', ...
            'MAR daily ME rate must have one value per timetable row')
      end
      native(end + 1, 1) = "ME";
   end
   if isempty(native)
      status = "not_available";
   elseif isequal(sort(native), sort(["SUH"; "SU"; "RZ"; "ME"]))
      status = "applied";
   else
      status = "partial";
   end

   % Daily ME is an independent validation of the public hourly MEH channel.
   % The compact ledger lets artifact QA repeat the comparison without reading
   % the large MAR archive and keeps incomplete boundary days explicit.
   [day_status, daily_reference, residual] = meltLedger( ...
      T, melt_daily_rate, kwargs.melt_abs_tolerance_mwe_day, ...
      kwargs.melt_rel_tolerance);
   if ~has_daily_melt
      melt_status = "not_available";
   elseif any(day_status == 2)
      melt_status = "mismatch";
   elseif any(day_status == 1)
      melt_status = "validated";
   else
      melt_status = "unverified";
   end
   finite_residual = residual(isfinite(residual));
   if isempty(finite_residual)
      maximum_residual = NaN;
   else
      maximum_residual = max(abs(finite_residual));
   end

   % Freeze source identities, temporal support, sign conventions, and the two
   % deliberately non-equivalent process relationships in durable metadata.
   metadata.mar_diagnostic_method = ...
      "native_hourly_and_daily_previous_hold";
   metadata.mar_diagnostic_status = status;
   metadata.mar_diagnostic_channels = channels;
   metadata.mar_diagnostic_native_variables = native;
   metadata.mar_diagnostic_units = "mWE/h";
   metadata.mar_diagnostic_subl_source = "SUH";
   metadata.mar_diagnostic_subl_evap_source = "SU";
   metadata.mar_diagnostic_refreeze_deposition_source = "RZ";
   metadata.mar_diagnostic_melt_hourly_source = "MEH";
   metadata.mar_diagnostic_melt_daily_source = "ME";
   metadata.mar_diagnostic_su_sector = kwargs.sector;
   metadata.mar_diagnostic_su_sector_name = sectorName(kwargs.sector);
   metadata.mar_diagnostic_rz_sector = 1;
   metadata.mar_diagnostic_me_sector = 1;
   metadata.mar_diagnostic_daily_distribution = ...
      "native_daily_divided_by_24_previous_hold";
   metadata.mar_diagnostic_subl_sign = ...
      "positive_loss_negative_deposition";
   metadata.mar_diagnostic_subl_evap_sign = ...
      "positive_loss_negative_deposition";
    metadata.mar_diagnostic_suh_su_relationship = ...
       "distinct_native_products_not_interchangeable";
   metadata.mar_diagnostic_rz_relationship = ...
      "combined_refreezing_and_deposition_not_pure_refreeze";
    metadata.mar_diagnostic_abs_limit_mwe_h = kwargs.abs_limit_mwe_h;
   metadata.mar_diagnostic_melt_validation_status = melt_status;
   metadata.mar_diagnostic_melt_abs_tolerance_mwe_day = ...
      kwargs.melt_abs_tolerance_mwe_day;
   metadata.mar_diagnostic_melt_rel_tolerance = kwargs.melt_rel_tolerance;
   metadata.mar_diagnostic_melt_day_status_codes = ...
      "1=validated;2=mismatch;3=unverified";
   metadata.mar_diagnostic_melt_day_status = day_status;
   metadata.mar_diagnostic_melt_daily_reference_mwe = daily_reference;
   metadata.mar_diagnostic_melt_residual_mwe_day = residual;
   metadata.mar_diagnostic_melt_validated_day_count = nnz(day_status == 1);
   metadata.mar_diagnostic_melt_mismatch_day_count = nnz(day_status == 2);
   metadata.mar_diagnostic_melt_unverified_day_count = nnz(day_status == 3);
   metadata.mar_diagnostic_melt_max_abs_error_mwe_day = maximum_residual;
    metadata.mar_diagnostic_basis = ...
       "MAR SUH is hourly sublimation; SU and RZ are native daily rates " ...
       + "divided by 24 and previous-held; daily ME validates but does not " ...
       + "replace hourly MEH";
    metadata = icemodel.forcing.helpers.marRefreezeMetadata(T, metadata);
end

function [status, reference, residual] = meltLedger( ...
      T, daily_rate, abs_tolerance, rel_tolerance)
   %MELTLEDGER Compare complete UTC-day MEH sums with native daily ME.
   if isempty(daily_rate)
      status = zeros(0, 1, 'uint8');
      reference = zeros(0, 1);
      residual = zeros(0, 1);
      return
   end

   % Only exact 24-row interval-start days have enough support for a daily
   % identity. Partial days retain an unverified code and NaN reference.
   times = T.Properties.RowTimes(:);
   days = dateshift(times, 'start', 'day');
   [groups, unique_days] = findgroups(days);
   ndays = numel(unique_days);
   status = repmat(uint8(3), ndays, 1);
   reference = nan(ndays, 1);
   residual = nan(ndays, 1);
   has_melt = ismember("melt", string(T.Properties.VariableNames));
   for day = 1:ndays
      rows = groups == day;
      actual_times = times(rows);
      expected_times = unique_days(day) + hours((0:23)');
      complete = isequal(actual_times, expected_times);
      source_daily = daily_rate(rows);
      if ~complete || ~all(isfinite(source_daily))
         continue
      end
      reference(day) = sum(source_daily);
      if ~has_melt
         continue
      end
      hourly = T.melt(rows);
      if ~all(isfinite(hourly))
         continue
      end
      residual(day) = sum(hourly) - reference(day);
      tolerance = abs_tolerance + rel_tolerance ...
         * max(abs([sum(hourly), reference(day)]));
      if abs(residual(day)) <= tolerance
         status(day) = uint8(1);
      else
         status(day) = uint8(2);
      end
   end
end

function name = sectorName(sector)
   %SECTORNAME Return MAR's documented two-surface-sector label.
   if sector == 1
      name = "permanent_ice";
   else
      name = "tundra";
   end
end
