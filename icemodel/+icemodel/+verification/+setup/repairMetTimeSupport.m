function report = repairMetTimeSupport(files, kwargs)
   %REPAIRMETTIMESUPPORT Repair legacy linear 15-minute met artifacts.
   %
   %  report = icemodel.verification.setup.repairMetTimeSupport(files)
   %  report = ... repairMetTimeSupport(files, dry_run=false)
   %
   % Reconstructs cached native rows from explicitly selected 15-minute met
   % files stamped by the legacy `linear_adjacent_finite_only` policy only when
   % their provenance proves a complete regular native grid, then delegates the
   % corrected interval-start hold to
   % icemodel.forcing.helpers.resampleMetTimestep. The default is read-only.
   % A write replaces only the named MAT file through a sibling temporary file;
   % no manifest discovery, raw-source read, or unrelated artifact mutation is
   % performed here.
   %
   % Inputs
   %  files   - explicit existing met MAT-file paths
   %
   % Name-value
   %  dry_run - report the exact repair without writing (default true)
   %
   % Output
   %  report  - deterministic per-file records plus changed/unchanged counts
   %
   % See also: icemodel.forcing.helpers.resampleMetTimestep,
   %  icemodel.verification.auditArtifacts

   arguments
      files string
      kwargs.dry_run (1, 1) logical = true
   end

   % Stable uniqueness makes dry-run inventories reproducible and prevents a
   % repeated input path from being replaced twice in one call.
   files = sort(unique(files(:)));
   records = repmat(emptyRecord(), numel(files), 1);
   for n = 1:numel(files)
      filename = files(n);
      [met, metadata] = loadMetArtifact(filename);
      record = emptyRecord();
      record.path = filename;
      record.old_policy = string(metadata.met_resample_policy);
      record.old_rows = height(met);
      record.old_end = timeString(met.Time(end));

      % Current hold/native artifacts are already idempotent. Only the exact
      % legacy policy is recoverable without reopening a raw source archive.
      if ismember(record.old_policy, ...
            ["interval_start_zero_order_hold", "native_15m_unchanged"])
         record.status = "unchanged";
         record.new_policy = record.old_policy;
         record.new_rows = record.old_rows;
         record.new_end = record.old_end;
         records(n) = record;
         continue
      end
      if record.old_policy ~= "linear_adjacent_finite_only"
         error('icemodel:verification:repairMetTimeSupport:unsupportedPolicy', ...
            'unsupported met resample policy in %s: %s', ...
            filename, record.old_policy)
      end

      % Linear interpolation preserved every native row at its original time.
      % Recover those rows on the recorded source cadence, then let the shared
      % helper rebuild all quarter-hour support and provenance.
      source = recoverSourceRows(met, metadata, filename);
      source = repairMerraSourceRows(source, metadata, filename);
      repaired = icemodel.forcing.helpers.resampleMetTimestep(source, "15m");
      verifyRecoveredRows(source, repaired, filename)
      record.new_policy = string( ...
         repaired.Properties.UserData.met_resample_policy);
      record.new_rows = height(repaired);
      record.new_end = timeString(repaired.Time(end));

      if kwargs.dry_run
         record.status = "would_repair";
      else
         saveAtomic(filename, repaired);
         record.status = "repaired";
      end
      records(n) = record;
   end

   % Keep the public result compact and directly usable as an acceptance log.
   statuses = string({records.status});
   report = struct();
   report.dry_run = kwargs.dry_run;
   report.files = records;
   report.summary = struct( ...
      'file_count', numel(records), ...
      'would_repair_count', nnz(statuses == "would_repair"), ...
      'repaired_count', nnz(statuses == "repaired"), ...
      'unchanged_count', nnz(statuses == "unchanged"));
end

function [met, metadata] = loadMetArtifact(filename)
   %LOADMETARTIFACT Validate the required staged-MAT payload and provenance.
   if ~isfile(filename)
      error('icemodel:verification:repairMetTimeSupport:fileNotFound', ...
         'met artifact does not exist: %s', filename)
   end
   names = string({whos('-file', filename).name});
   if ~ismember("met", names)
      error('icemodel:verification:repairMetTimeSupport:badArtifact', ...
         'met artifact lacks the met payload: %s', filename)
   end
   S = load(filename, 'met', '-mat');
   met = S.met;
   if ~istimetable(met)
      error('icemodel:verification:repairMetTimeSupport:badArtifact', ...
         'met artifact has inconsistent payload provenance: %s', filename)
   end
   userdata = met.Properties.UserData;
   if ~isstruct(userdata) || ~isfield(userdata, 'met_resample_policy')
      error('icemodel:verification:repairMetTimeSupport:badArtifact', ...
         'met artifact has inconsistent payload provenance: %s', filename)
   end
   metadata = userdata;
   if ismember("artifact_metadata", names)
      S = load(filename, 'artifact_metadata', '-mat');
      if ~isstruct(S.artifact_metadata) ...
            || ~isequaln(S.artifact_metadata, userdata)
         error('icemodel:verification:repairMetTimeSupport:badArtifact', ...
            'met artifact has inconsistent payload provenance: %s', filename)
      end
      metadata = S.artifact_metadata;
   end
end

function source = recoverSourceRows(met, metadata, filename)
   %RECOVERSOURCEROWS Select exact cached source stamps from a linear artifact.
   required = ["met_resample_source_cadence_seconds", ...
      "met_resample_source_row_count", ...
      "met_resample_source_time_gap_count"];
   if ~all(isfield(metadata, required)) || height(met) < 2 ...
         || any(diff(met.Time) ~= minutes(15))
      error('icemodel:verification:repairMetTimeSupport:badLegacyMetadata', ...
         'legacy met cadence/provenance is incomplete: %s', filename)
   end

   cadence_s = double(metadata.met_resample_source_cadence_seconds);
   source_count = double(metadata.met_resample_source_row_count);
   gap_count = double(metadata.met_resample_source_time_gap_count);
   if ~isscalar(cadence_s) || ~isfinite(cadence_s) || cadence_s < 900 ...
         || abs(mod(cadence_s, 900)) > 1e-6 ...
         || ~isscalar(source_count) || ~isfinite(source_count) ...
         || source_count < 1 || source_count ~= fix(source_count) ...
         || ~isscalar(gap_count) || ~isfinite(gap_count) ...
         || gap_count < 0 || gap_count ~= fix(gap_count)
      error('icemodel:verification:repairMetTimeSupport:badLegacyMetadata', ...
         'legacy met source cadence/row count is invalid: %s', filename)
   end

   % An omitted native timestamp cannot be distinguished safely from a finite
   % row invented by the old regular output grid. Fail closed unless both the
   % recorded gap count and exact output-row identity prove a complete grid.
   cadence_ratio = cadence_s / 900;
   expected_output_count = (source_count - 1) * cadence_ratio + 1;
   if gap_count ~= 0 || height(met) ~= expected_output_count
      error( ...
         'icemodel:verification:repairMetTimeSupport:ambiguousLegacyGrid', ...
         ['legacy met cannot prove an exact gap-free source reconstruction ' ...
         'in %s'], filename)
   end

   % The old uniform output ends at the final native row. With the proof above,
   % selecting timestamps on the recorded cadence recovers every source row and
   % no regular-grid row can stand in for an omitted source timestamp.
   elapsed = seconds(met.Time - met.Time(1));
   on_source_grid = abs(elapsed / cadence_s ...
      - round(elapsed / cadence_s)) <= 1e-9;
   source = met(on_source_grid, :);
   if height(source) ~= source_count
      error( ...
         'icemodel:verification:repairMetTimeSupport:ambiguousLegacyGrid', ...
         'legacy met source-row identity is ambiguous: %s', filename)
   end

   % The recovered timetable is native source state, not an already-derived
   % guarded output. Remove only legacy resample fields so the shared helper can
   % stamp the corrected policy while retaining every upstream source field.
   source_metadata = source.Properties.UserData;
   names = string(fieldnames(source_metadata));
   legacy = names(startsWith(names, "met_resample_"));
   if ~isempty(legacy)
      source_metadata = rmfield(source_metadata, cellstr(legacy));
   end
   source.Properties.UserData = source_metadata;
end

function verifyRecoveredRows(source, repaired, filename)
   %VERIFYRECOVEREDROWS Prove the repair changes only derived quarter-hour rows.
   [present, rows] = ismember(source.Time, repaired.Time);
   if ~all(present) || ~isequaln(source{:, :}, repaired{rows, :})
      error('icemodel:verification:repairMetTimeSupport:sourceChanged', ...
         'repaired met changed a cached source row: %s', filename)
   end
end

function source = repairMerraSourceRows(source, metadata, filename)
   %REPAIRMERRASOURCEROWS Re-hold proven hourly glc rows before 15m resampling.
   channels = intersect(["runoff", "albedo", "snowd", "swe"], ...
      string(source.Properties.VariableNames), 'stable');
   fields = string(fieldnames(metadata));
   [~, stem] = fileparts(filename);
   lower_filename = lower(string(stem));
   is_merra = contains(lower_filename, "_merra2_") ...
      || contains(lower_filename, "_merra_") ...
      || any(startsWith(fields, "merra_"));
   if ~is_merra || isempty(channels)
      return
   end

   % A legacy MERRA met artifact is unsafe unless its recovered hourly rows carry
   % the exact raw glc timestamp inventory introduced by the Data builder/repair.
   if height(source) > 1 && any(diff(source.Time) ~= hours(1))
      error('icemodel:verification:repairMetTimeSupport:badMerraCadence', ...
         'legacy MERRA met source rows are not uniformly hourly: %s', filename)
   end
   source_metadata = source.Properties.UserData;
   if ~icemodel.forcing.helpers.hasProvenMerraTavg3SourceGrid( ...
         source, source_metadata)
      error('icemodel:verification:repairMetTimeSupport:unprovenMerraSourceGrid', ...
         ['legacy MERRA met lacks an exact native glc timestamp inventory; ' ...
         'regenerate it from corrected hourly MERRA Data: %s'], filename)
   end
   [source, source_metadata] = ...
      icemodel.forcing.helpers.applyMerraTimeSupport( ...
      source, source_metadata);
   source.Properties.UserData = source_metadata;
end

function saveAtomic(filename, repaired)
   %SAVEATOMIC Replace one artifact only after its sibling save succeeds.
   temporary = string(tempname(fileparts(filename))) + ".mat";
   cleanup = onCleanup(@() deleteIfPresent(temporary));
   [copied, message] = copyfile(filename, temporary, 'f');
   if ~copied
      error('icemodel:verification:repairMetTimeSupport:copyFailed', ...
         'could not copy %s before repair: %s', filename, message)
   end
   met = repaired;
   artifact_metadata = repaired.Properties.UserData;
   save(temporary, 'met', 'artifact_metadata', '-append')
   [moved, message] = movefile(temporary, filename, 'f');
   if ~moved
      error('icemodel:verification:repairMetTimeSupport:replaceFailed', ...
         'could not replace %s: %s', filename, message)
   end
   clear cleanup
end

function deleteIfPresent(filename)
   %DELETEIFPRESENT Remove an abandoned sibling temporary artifact.
   if isfile(filename)
      delete(filename)
   end
end

function value = timeString(time)
   %TIMESTRING Serialize one UTC datetime without display-format dependence.
   value = string(time, 'yyyy-MM-dd HH:mm:ss');
end

function record = emptyRecord()
   %EMPTYRECORD Stable per-file dry-run/write result shape.
   record = struct('path', "", 'status', "", 'old_policy', "", ...
      'new_policy', "", 'old_rows', 0, 'new_rows', 0, ...
      'old_end', "", 'new_end', "");
end
