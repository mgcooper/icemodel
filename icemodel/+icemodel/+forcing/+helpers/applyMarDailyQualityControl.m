function [T, metadata] = applyMarDailyQualityControl( ...
      T, replacements, metadata, kwargs)
   %APPLYMARDAILYQUALITYCONTROL Constrain MAR hourly mass data by daily totals.
   %
   %  [T, metadata] = ... applyMarDailyQualityControl(T, replacements)
   %  [T, metadata] = ... applyMarDailyQualityControl(_, _, metadata, ...
   %     sector=1)
   %
   % REPLACEMENTS is a scalar struct whose optional runoff and smb fields are
   % hourly rates derived from native MAR daily accumulations (daily / 24,
   % held over the UTC day). A complete hourly day is preserved when its sum
   % agrees with the native daily reference. Missing, partial, or inconsistent
   % days use the daily rate; days lacking a finite daily reference retain the
   % hourly source and are marked unverified.
   %
   % Metadata stores one compact status code and native daily reference per
   % UTC day and channel. This lets source-light artifact QA distinguish
   % preserved (1), replaced (2), and unverified (3) support and recheck every
   % constrained daily total without reopening the MAR archive.
   %
   % See also: icemodel.forcing.buildMarData,
   %  icemodel.forcing.helpers.dailyToHourly

   arguments
      T timetable
      replacements (1, 1) struct
      metadata (1, 1) struct = struct()
      kwargs.sector (1, 1) double {mustBeMember(kwargs.sector, [1 2])} = 1
      kwargs.abs_tolerance_mwe_day (1, 1) double ...
         {mustBeNonnegative} = 1e-5
      kwargs.rel_tolerance (1, 1) double {mustBeNonnegative} = 1e-4
   end

   % Accept only the two diagnostics covered by the source comparison. Stable
   % channel order keeps builder and repair metadata shape-identical.
   replacement_names = string(fieldnames(replacements));
   unexpected = setdiff(replacement_names, ["runoff", "smb"]);
   assert(isempty(unexpected), ...
      'MAR daily replacements may contain only runoff and smb')
   channels = intersect(["runoff", "smb"], replacement_names, 'stable');
   [groups, days, complete] = utcDayGroups(T.Time);
   ndays = numel(days);

   % Start both channel ledgers as unverified. Applied replacement fields
   % overwrite these defaults; a reduced source therefore has an explicit
   % per-day ledger without changing its native hourly values.
   ledgers = struct();
   for channel = ["runoff", "smb"]
      ledgers.(channel) = emptyLedger(ndays);
   end
   original = T;
   for channel = reshape(channels, 1, [])
      replacement = replacements.(channel);
      replacement = replacement(:);
      assert(numel(replacement) == height(T), ...
         'MAR %s replacement must have one value per timetable row', channel)
      if ismember(channel, string(T.Properties.VariableNames))
         source = T.(channel);
      else
         % A source-backed daily channel may complete an older reduced
         % artifact. Missing hourly input is treated as an incomplete day.
         source = nan(height(T), 1);
      end
      [T.(channel), ledgers.(channel)] = constrainChannel( ...
         source, replacement, groups, complete, ...
         kwargs.abs_tolerance_mwe_day, kwargs.rel_tolerance);
   end

   % Reapplying an unchanged current record is byte-idempotent. In particular,
   % a previously replaced constant day must not be relabelled as preserved on
   % the second pass merely because it now matches the daily reference.
   current = metadataCurrent(metadata, channels, ndays, kwargs);
   if current && isequaln(T, original)
      return
   end

   % Cumulative changed-sample counts retain repair history, while day-ledger
   % counts describe the current application and are reproducible from status.
   runoff_count = priorCount(metadata, ...
      'mar_qc_replaced_runoff_count', current) ...
      + nnz(unequalWithMissing(column(original, "runoff"), ...
      column(T, "runoff")));
   smb_count = priorCount(metadata, ...
      'mar_qc_replaced_smb_count', current) ...
      + nnz(unequalWithMissing(column(original, "smb"), ...
      column(T, "smb")));

   if isempty(channels)
      status = 'not_applicable';
      fallback = 'hourly_RUH_SMBH_retained_native_daily_unavailable';
      basis = 'MAR native daily RU/SMB unavailable; retained hourly RUH/SMBH';
   else
      status = 'applied';
      fallback = 'none';
      basis = [ ...
         'MAR hourly RUH/SMBH preserved where complete UTC-day sums agree ' ...
         'with native daily RU/SMB; missing, partial, or inconsistent days ' ...
         'use the native daily rate'];
   end
   if kwargs.sector == 1
      sector_name = 'permanent_ice';
   else
      sector_name = 'tundra';
   end

   % Stamp the complete daily-constrained provenance contract. The day axis is
   % derivable from T.Time, so only aligned status/reference vectors are saved.
   metadata.mar_qc_method = 'daily_constrained_hourly';
   metadata.mar_qc_status = status;
   metadata.mar_qc_fallback = fallback;
   metadata.mar_qc_channels = channels(:);
   metadata.mar_qc_sector = kwargs.sector;
   metadata.mar_qc_sector_name = sector_name;
   metadata.mar_qc_runoff_source = 'RU';
   metadata.mar_qc_smb_source = 'SMB';
   metadata.mar_qc_hourly_distribution = ...
      'preserve_matching_hourly_else_daily_divided_by_24';
   metadata.mar_qc_partial_day_policy = ...
      'native_daily_rate_applied_to_available_rows_marked_replaced';
   metadata.mar_qc_abs_tolerance_mwe_day = ...
      kwargs.abs_tolerance_mwe_day;
   metadata.mar_qc_rel_tolerance = kwargs.rel_tolerance;
   metadata.mar_qc_day_status_codes = ...
      '1=preserved_hourly;2=replaced_daily;3=unverified';
   metadata.mar_qc_daily_reference_units = 'mWE/day';
   metadata.mar_qc_complete_utc_day_count = nnz(complete);
   metadata.mar_qc_partial_utc_day_count = nnz(~complete);
   for channel = ["runoff", "smb"]
      ledger = ledgers.(channel);
      metadata.("mar_qc_" + channel + "_day_status") = ledger.status;
      metadata.("mar_qc_" + channel + "_daily_reference_mwe") = ...
         ledger.reference;
      metadata.("mar_qc_preserved_" + channel + "_day_count") = ...
         nnz(ledger.status == 1);
      metadata.("mar_qc_replaced_" + channel + "_day_count") = ...
         nnz(ledger.status == 2);
      metadata.("mar_qc_unverified_" + channel + "_day_count") = ...
         nnz(ledger.status == 3);
   end
   metadata.mar_qc_replaced_runoff_count = runoff_count;
   metadata.mar_qc_replaced_smb_count = smb_count;
   metadata.mar_qc_basis = basis;
end

function [values, ledger] = constrainChannel( ...
      source, replacement, groups, complete, abs_tolerance, rel_tolerance)
   %CONSTRAINCHANNEL Preserve matching hourly days and repair the remainder.
   values = source(:);
   ledger = emptyLedger(numel(complete));
   for day = 1:numel(complete)
      rows = groups == day;
      daily_rate = replacement(rows);
      valid_reference = all(isfinite(daily_rate));
      if valid_reference
         scale = max(1, max(abs(daily_rate)));
         valid_reference = max(daily_rate) - min(daily_rate) ...
            <= 16 * eps(scale);
      end
      if ~valid_reference
         % A missing or internally inconsistent native daily source cannot
         % constrain this day. Preserve the hourly source and expose status 3.
         continue
      end

      reference = 24 * daily_rate(1);
      ledger.reference(day) = reference;
      hourly = source(rows);
      matches = complete(day) && all(isfinite(hourly));
      if matches
         hourly_total = sum(hourly);
         tolerance = abs_tolerance + rel_tolerance ...
            * max(abs([reference, hourly_total]));
         numeric_slack = 16 * eps(max(1, ...
            max(abs([reference, hourly_total]))));
         matches = abs(hourly_total - reference) ...
            <= tolerance + numeric_slack;
      end
      if matches
         ledger.status(day) = 1;
      else
         values(rows) = daily_rate;
         ledger.status(day) = 2;
      end
   end
end

function ledger = emptyLedger(ndays)
   %EMPTYLEDGER Allocate an unverified per-day MAR provenance ledger.
   ledger = struct('status', uint8(3 * ones(ndays, 1)), ...
      'reference', nan(ndays, 1));
end

function [groups, days, complete] = utcDayGroups(times)
   %UTCDAYGROUPS Group rows and identify exact 00:00--23:00 UTC support.
   if isempty(times)
      groups = zeros(0, 1);
      days = times;
      complete = false(0, 1);
      return
   end
   [groups, days] = findgroups(dateshift(times(:), 'start', 'day'));
   complete = false(numel(days), 1);
   for day = 1:numel(days)
      complete(day) = isequal(times(groups == day), ...
         days(day) + hours((0:23)'));
   end
end

function tf = metadataCurrent(metadata, channels, ndays, kwargs)
   %METADATACURRENT True for the exact current hybrid-QC configuration.
   required = ["mar_qc_method", "mar_qc_channels", ...
      "mar_qc_abs_tolerance_mwe_day", "mar_qc_rel_tolerance", ...
      "mar_qc_sector", "mar_qc_sector_name", ...
      "mar_qc_runoff_day_status", "mar_qc_smb_day_status"];
   if kwargs.sector == 1
      sector_name = "permanent_ice";
   else
      sector_name = "tundra";
   end
   tf = all(isfield(metadata, required)) ...
      && string(metadata.mar_qc_method) == "daily_constrained_hourly" ...
      && isequal(string(metadata.mar_qc_channels), channels(:)) ...
      && double(metadata.mar_qc_sector) == kwargs.sector ...
      && string(metadata.mar_qc_sector_name) == sector_name ...
      && double(metadata.mar_qc_abs_tolerance_mwe_day) ...
      == kwargs.abs_tolerance_mwe_day ...
      && double(metadata.mar_qc_rel_tolerance) == kwargs.rel_tolerance ...
      && numel(metadata.mar_qc_runoff_day_status) == ndays ...
      && numel(metadata.mar_qc_smb_day_status) == ndays;
end

function values = column(T, channel)
   %COLUMN Return one table channel or an all-missing comparison vector.
   if ismember(channel, string(T.Properties.VariableNames))
      values = T.(channel);
   else
      values = nan(height(T), 1);
   end
end

function count = priorCount(metadata, field, current)
   %PRIORCOUNT Return a valid prior changed-sample count for restamping.
   count = 0;
   if current && isfield(metadata, field)
      candidate = double(metadata.(field));
      if isscalar(candidate) && isfinite(candidate) && candidate >= 0
         count = candidate;
      end
   end
end

function changed = unequalWithMissing(old, replacement)
   %UNEQUALWITHMISSING Compare numeric vectors with NaN treated as equal.
   old = old(:);
   replacement = replacement(:);
   same_missing = isnan(old) & isnan(replacement);
   changed = ~(old == replacement | same_missing);
end
