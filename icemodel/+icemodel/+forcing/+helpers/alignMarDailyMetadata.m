function metadata = alignMarDailyMetadata(metadata, source_days, retained_times)
   %ALIGNMARDAILYMETADATA Align MAR per-day provenance to a retained time axis.
   %
   %  metadata = icemodel.forcing.helpers.alignMarDailyMetadata( ...
   %     metadata, source_days, retained_times)
   %
   % SOURCE_DAYS is the exact UTC-day axis used to create the saved RU/SMB and
   % ME/MEH ledgers. RETAINED_TIMES is the hourly Data or derived-met axis kept
   % after a staging window is applied. Complete retained days preserve their
   % source ledger entries exactly. Partial boundary days become unverified
   % (status 3) with NaN references/residuals. Numeric payloads are not inputs
   % and therefore cannot be modified by this metadata-only operation.

   arguments
      metadata (1, 1) struct
      source_days (:, 1) datetime
      retained_times (:, 1) datetime
   end

   % Normalize both axes to UTC before comparing calendar days. A timezone-less
   % source is interpreted as UTC, matching every MAR reader/writer contract.
   source_days = utcColumn(source_days);
   retained_times = utcColumn(retained_times);
   if any(source_days ~= dateshift(source_days, 'start', 'day')) ...
         || any(diff(source_days) <= seconds(0))
      error('icemodel:forcing:alignMarDailyMetadata:sourceDayAxis', ...
         'MAR source days must be unique, increasing UTC midnights')
   end
   if any(diff(retained_times) <= seconds(0))
      error('icemodel:forcing:alignMarDailyMetadata:retainedTimeAxis', ...
         'MAR retained times must be strictly increasing')
   end

   % Validate each independently meaningful ledger as an atomic group. Empty
   % ME/MEH fields are the supported reduced-source representation; a partially
   % populated or wrong-length group is ambiguous and must fail closed.
   qc_runoff = ["mar_qc_runoff_day_status", ...
      "mar_qc_runoff_daily_reference_mwe"];
   qc_smb = ["mar_qc_smb_day_status", ...
      "mar_qc_smb_daily_reference_mwe"];
   melt = ["mar_diagnostic_melt_day_status", ...
      "mar_diagnostic_melt_daily_reference_mwe", ...
      "mar_diagnostic_melt_residual_mwe_day"];
   validateLedgerGroup(metadata, qc_runoff, numel(source_days));
   validateLedgerGroup(metadata, qc_smb, numel(source_days));
   validateLedgerGroup(metadata, melt, numel(source_days));

   % Select retained UTC days from the source ledger, rejecting times that have
   % no source-day identity rather than guessing a positional correspondence.
   retained_days = unique(dateshift(retained_times, 'start', 'day'), 'stable');
   [found, source_index] = ismember(retained_days, source_days);
   if ~all(found)
      error('icemodel:forcing:alignMarDailyMetadata:retainedDayOutsideSource', ...
         'MAR retained times include a UTC day absent from the source ledger')
   end
   complete = retainedDayCompleteness(retained_times, retained_days);
   partial = ~complete;

   % Subset every nonempty ledger. Status and numeric ledgers have distinct
   % partial-day sentinels but otherwise retain their exact source values.
   status_fields = [qc_runoff(1), qc_smb(1), melt(1)];
   value_fields = [qc_runoff(2), qc_smb(2), melt(2:3)];
   for field = status_fields
      if isfield(metadata, field) && ~isempty(metadata.(field))
         values = metadata.(field);
         values = values(source_index);
         values(partial) = cast(3, 'like', values);
         metadata.(field) = values(:);
      end
   end
   for field = value_fields
      if isfield(metadata, field) && ~isempty(metadata.(field))
         values = metadata.(field);
         values = values(source_index);
         values(partial) = NaN;
         metadata.(field) = values(:);
      end
   end

   % RU/SMB summaries describe only the retained artifact. Cumulative changed-
   % sample counts intentionally remain untouched as repair-history provenance.
   qc_status_fields = [qc_runoff(1), qc_smb(1)];
   if any(isfield(metadata, cellstr(qc_status_fields)))
      metadata.mar_qc_complete_utc_day_count = nnz(complete);
      metadata.mar_qc_partial_utc_day_count = nnz(partial);
      for channel = ["runoff", "smb"]
         field = "mar_qc_" + channel + "_day_status";
         if isfield(metadata, field)
            status = uint8(metadata.(field));
            metadata.("mar_qc_preserved_" + channel + "_day_count") = ...
               nnz(status == 1);
            metadata.("mar_qc_replaced_" + channel + "_day_count") = ...
               nnz(status == 2);
            metadata.("mar_qc_unverified_" + channel + "_day_count") = ...
               nnz(status == 3);
         end
      end
   end

   % Recompute the ME/MEH validation summary from the aligned ledger. Empty
   % reduced-source ledgers remain explicitly unavailable, not forcing failures.
   if isfield(metadata, melt(1))
      status = uint8(metadata.(melt(1)));
      if isempty(status)
         validation_status = "not_available";
      elseif any(status == 2)
         validation_status = "mismatch";
      elseif any(status == 1)
         validation_status = "validated";
      else
         validation_status = "unverified";
      end
      metadata.mar_diagnostic_melt_validation_status = validation_status;
      metadata.mar_diagnostic_melt_validated_day_count = nnz(status == 1);
      metadata.mar_diagnostic_melt_mismatch_day_count = nnz(status == 2);
      metadata.mar_diagnostic_melt_unverified_day_count = nnz(status == 3);
      residual = double(metadata.(melt(3)));
      residual = residual(isfinite(residual));
      if isempty(residual)
         maximum_residual = NaN;
      else
         maximum_residual = max(abs(residual));
      end
      metadata.mar_diagnostic_melt_max_abs_error_mwe_day = maximum_residual;
   end
end

function values = utcColumn(values)
   %UTCCOLUMN Normalize a datetime vector to one UTC column.
   values = values(:);
   values.TimeZone = 'UTC';
   if any(isnat(values))
      error('icemodel:forcing:alignMarDailyMetadata:missingTime', ...
         'MAR metadata alignment axes cannot contain NaT')
   end
end

function validateLedgerGroup(metadata, fields, source_count)
   %VALIDATELEDGERGROUP Require one complete numeric ledger or an empty group.
   present = isfield(metadata, cellstr(fields));
   populated = false(size(fields));
   for k = 1:numel(fields)
      if present(k)
         values = metadata.(fields(k));
         if ~isnumeric(values)
            error('icemodel:forcing:alignMarDailyMetadata:ledgerType', ...
               'MAR daily ledger %s must be numeric', fields(k))
         end
         populated(k) = ~isempty(values);
      end
   end
   if any(populated) && ~all(present & populated)
      error('icemodel:forcing:alignMarDailyMetadata:ledgerSchema', ...
         'MAR daily ledger group is only partially populated')
   end
   if any(populated)
      lengths = arrayfun(@(field) numel(metadata.(field)), fields);
      if any(lengths ~= source_count)
         error('icemodel:forcing:alignMarDailyMetadata:ledgerLength', ...
            'MAR daily ledger length does not match the source UTC-day axis')
      end
      status = double(metadata.(fields(1)));
      if any(~ismember(status, [1 2 3]))
         error('icemodel:forcing:alignMarDailyMetadata:statusCode', ...
            'MAR daily ledger status must use codes 1, 2, or 3')
      end
      for field = fields(2:end)
         if ~isfloat(metadata.(field))
            error('icemodel:forcing:alignMarDailyMetadata:ledgerType', ...
               'MAR daily references and residuals must support NaN')
         end
      end
   end
end

function complete = retainedDayCompleteness(times, days)
   %RETAINEDDAYCOMPLETENESS Identify full days at the retained regular cadence.
   complete = false(numel(days), 1);
   if isempty(times) || numel(times) < 2
      return
   end

   % Missing intervals may make a day partial, but every observed step must be
   % an integer multiple of one cadence that divides a UTC day exactly.
   steps = seconds(diff(times));
   cadence = min(steps);
   samples_per_day = 86400 / cadence;
   tolerance = 32 * eps(max(1, samples_per_day));
   if ~isfinite(cadence) || cadence <= 0 ...
         || abs(samples_per_day - round(samples_per_day)) > tolerance ...
         || any(abs(steps / cadence - round(steps / cadence)) > tolerance)
      error('icemodel:forcing:alignMarDailyMetadata:retainedCadence', ...
         'MAR retained times do not have one unambiguous regular cadence')
   end

   % Compare exact interval-start support within each retained UTC day. This is
   % cadence-neutral: hourly Data and 15-minute met are never conflated.
   samples_per_day = round(samples_per_day);
   offsets = seconds((0:samples_per_day - 1)' * cadence);
   groups = findgroups(dateshift(times, 'start', 'day'));
   for day = 1:numel(days)
      complete(day) = isequal(times(groups == day), days(day) + offsets);
   end
end
