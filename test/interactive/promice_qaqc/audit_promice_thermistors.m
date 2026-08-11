%AUDIT_PROMICE_THERMISTORS Reproduce the all-site PROMICE tice10m QC audit.
%
% Run from the repository root after staging the pypromice L3 hourly bundle:
%   run('test/interactive/promice_qaqc/audit_promice_thermistors.m')
%
% The script writes one row per >1 K exactly-hourly source discontinuity to
% artifacts/promice_thermistor_audit.csv. It relates each event to the native
% depth-tagged thermistor string and verifies the production QC classification.

% Resolve the canonical source and output paths so the audit is independent of
% the caller's current directory after the project is on the MATLAB path.
source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
   'verification', 'promice'));
script_dir = string(fileparts(mfilename('fullpath')));
artifact_dir = fullfile(script_dir, 'artifacts');
if ~isfolder(artifact_dir)
   mkdir(artifact_dir)
end
output_file = fullfile(artifact_dir, 'promice_thermistor_audit.csv');

% Accumulate the complete cross-site distribution and a compact event ledger.
files = dir(fullfile(source_dir, 'hour', '*_hour.nc'));
all_jumps_K = zeros(0, 1);
record_prototype = struct( ...
   'site', "", 'time', NaT(1, 1, TimeZone="UTC"), ...
   'jump_K', NaN, 'source_before_K', NaN, 'source_after_K', NaN, ...
   'qc_code', NaN, 'native_relationship', "", ...
   'comparable_sensor_count', 0, 'near_10m_sensor_count', 0, ...
   'max_jump_sensor', "", 'max_sensor_jump_K', NaN, ...
   'max_jump_depth_before_m', NaN, 'max_jump_depth_after_m', NaN, ...
   'max_depth_jump_m', NaN, ...
   'sensor_support_changed', false, ...
   'mask_start', NaT(1, 1, TimeZone="UTC"), ...
   'mask_end', NaT(1, 1, TimeZone="UTC"), ...
   'masked_sample_count', 0);
records = repmat(record_prototype, 0, 1);

% One jump block and one event block per source file. Sizing the block buffers
% to the file list keeps the cross-site accumulators out of the per-sample and
% per-event growth path; both are concatenated once after the scan.
jump_blocks = repmat({zeros(0, 1)}, numel(files), 1);
record_blocks = repmat({records}, numel(files), 1);
sites_with_tice10m = 0;
failed_sample_count = 0;
unreviewed_endpoint_count = 0;
persistent_sample_count = 0;

for file_index = 1:numel(files)
   % Use the production reader so this audit proves both the native evidence and
   % the exact source/mask/flag contract that will be staged.
   site = erase(string(files(file_index).name), '_hour.nc');
   aws = icemodel.forcing.readPromiceAws(site, source_dir=source_dir);
   names = string(aws.Properties.VariableNames);
   if ~ismember("tice10m_source", names)
      continue
   end
   sites_with_tice10m = sites_with_tice10m + 1;
   failed_sample_count = failed_sample_count ...
      + nnz(aws.tice10m_qc_flag == 1);
   unreviewed_endpoint_count = unreviewed_endpoint_count ...
      + nnz(aws.tice10m_qc_flag == 2);
   persistent_sample_count = persistent_sample_count ...
      + nnz(aws.tice10m_qc_flag == 3);

   source = aws.tice10m_source;
   consecutive = abs(seconds(diff(aws.Time)) - 3600) <= 1;
   jump = abs(diff(source));
   finite_pair = isfinite(source(1:end-1)) & isfinite(source(2:end));
   jump_blocks{file_index} = jump(consecutive & finite_pair);
   events = find(consecutive & finite_pair & jump > 1);
   thermistors = names(~cellfun('isempty', ...
      regexp(cellstr(names), '^tice\d+$', 'once')));

   % Every event contributes exactly one ledger row, so this file's block is
   % sized to its event list and filled by index.
   n_thermistors = numel(thermistors);
   event_records = repmat(record_prototype, numel(events), 1);
   n_records = 0;

   for first = reshape(events, 1, [])
      % Quantify native temperature/depth continuity across the same endpoint
      % pair. The 8..12 m band matches pypromice's documented 2 m target support.
      pair = first:first + 1;
      % Each thermistor contributes at most one entry to every per-event array,
      % so they are sized to the thermistor string once and trimmed to their
      % realized counts after the loop.
      sensor_jumps = zeros(n_thermistors, 1);
      depth_jumps = zeros(n_thermistors, 1);
      comparable_names = strings(n_thermistors, 1);
      depths_before = zeros(n_thermistors, 1);
      depths_after = zeros(n_thermistors, 1);
      comparable = 0;
      near_10m = 0;
      before_support = strings(n_thermistors, 1);
      after_support = strings(n_thermistors, 1);
      n_before_support = 0;
      n_after_support = 0;
      for name = reshape(thermistors, 1, [])
         depth_name = "d" + name;
         if ~ismember(depth_name, names)
            continue
         end
         temperature = aws.(name);
         depth = aws.(depth_name);
         before_valid = isfinite(temperature(first)) ...
            && isfinite(depth(first)) && depth(first) > 0;
         after_valid = isfinite(temperature(first + 1)) ...
            && isfinite(depth(first + 1)) && depth(first + 1) > 0;
         if before_valid
            n_before_support = n_before_support + 1;
            before_support(n_before_support) = name;
         end
         if after_valid
            n_after_support = n_after_support + 1;
            after_support(n_after_support) = name;
         end
         if before_valid && after_valid
            comparable = comparable + 1;
            comparable_names(comparable) = name;
            sensor_jumps(comparable) = abs(diff(temperature(pair)));
            depth_jumps(comparable) = abs(diff(depth(pair)));
            depths_before(comparable) = depth(first);
            depths_after(comparable) = depth(first + 1);
            if any(depth(pair) >= 8 & depth(pair) <= 12)
               near_10m = near_10m + 1;
            end
         end
      end

      % Drop the unused tails so the support comparison and the max-jump
      % lookups below see only the thermistors that actually reported.
      sensor_jumps = sensor_jumps(1:comparable);
      depth_jumps = depth_jumps(1:comparable);
      comparable_names = comparable_names(1:comparable);
      depths_before = depths_before(1:comparable);
      depths_after = depths_after(1:comparable);
      before_support = before_support(1:n_before_support);
      after_support = after_support(1:n_after_support);

      % Classify the evidence without changing the production result. Every
      % event remains visible even if native neighbors are insufficient.
      relationship = "derived_target_only";
      if comparable < 2
         relationship = "neighbor_insufficient";
      elseif any(sensor_jumps > 1)
         relationship = "native_sensor_jump";
      elseif ~isequal(before_support, after_support) ...
            || any(depth_jumps > 0.5)
         relationship = "depth_or_sensor_transition";
      end
      if isempty(sensor_jumps)
         max_jump_sensor = "";
         max_sensor_jump = NaN;
         max_jump_depth_before = NaN;
         max_jump_depth_after = NaN;
         max_depth_jump = NaN;
      else
         [max_sensor_jump, max_index] = max(sensor_jumps);
         max_jump_sensor = comparable_names(max_index);
         max_jump_depth_before = depths_before(max_index);
         max_jump_depth_after = depths_after(max_index);
         max_depth_jump = max(depth_jumps);
      end

      % Resolve the complete mask segment containing the event. Code 3 may
      % extend far beyond the two transition endpoints until a depth reset.
      qc_code = max(aws.tice10m_qc_flag(pair));
      mask_first = first;
      mask_last = first + 1;
      if qc_code == 3
         before = find(aws.tice10m_qc_flag(1:first) ~= 3, 1, 'last');
         after = find(aws.tice10m_qc_flag(first + 1:end) ~= 3, 1);
         if isempty(before)
            mask_first = 1;
         else
            mask_first = before + 1;
         end
         if isempty(after)
            mask_last = height(aws);
         else
            mask_last = first + after - 1;
         end
      end

      n_records = n_records + 1;
      event_records(n_records) = struct( ...
         'site', site, 'time', aws.Time(first + 1), ...
         'jump_K', jump(first), ...
         'source_before_K', source(first), ...
         'source_after_K', source(first + 1), ...
         'qc_code', qc_code, ...
         'native_relationship', relationship, ...
         'comparable_sensor_count', comparable, ...
         'near_10m_sensor_count', near_10m, ...
         'max_jump_sensor', max_jump_sensor, ...
         'max_sensor_jump_K', max_sensor_jump, ...
         'max_jump_depth_before_m', max_jump_depth_before, ...
         'max_jump_depth_after_m', max_jump_depth_after, ...
         'max_depth_jump_m', max_depth_jump, ...
         'sensor_support_changed', ...
         ~isequal(before_support, after_support), ...
         'mask_start', aws.Time(mask_first), ...
         'mask_end', aws.Time(mask_last), ...
         'masked_sample_count', mask_last - mask_first + 1);
   end
   record_blocks{file_index} = event_records(1:n_records);
end

% Stack the per-file blocks once so the ledger and the cross-site distribution
% match what the per-event growth loop used to accumulate.
all_jumps_K = vertcat(all_jumps_K, jump_blocks{:});
records = vertcat(records, record_blocks{:});

% Persist the deterministic event ledger and print the cross-site threshold
% evidence needed to review future product versions.
event_table = struct2table(records);
event_table = sortrows(event_table, {'jump_K', 'site'}, {'descend', 'ascend'});
writetable(event_table, output_file)
fprintf('PROMICE tice10m sites: %d\n', sites_with_tice10m)
fprintf('Finite consecutive hourly pairs: %d\n', numel(all_jumps_K))
fprintf('P50/P95/P99/P99.9/max |delta T10| [K]: %.6g %.6g %.6g %.6g %.6g\n', ...
   prctile(all_jumps_K, 50), prctile(all_jumps_K, 95), ...
   prctile(all_jumps_K, 99), prctile(all_jumps_K, 99.9), max(all_jumps_K))
fprintf('Events >1 K: %d (failed=%d, unreviewed=%d, persistent=%d)\n', ...
   height(event_table), nnz(event_table.qc_code == 1), ...
   nnz(event_table.qc_code >= 2), nnz(event_table.qc_code == 3))
fprintf(['Unique masked samples: failed=%d, neighbor-insufficient=%d, ' ...
   'persistent=%d\n'], failed_sample_count, unreviewed_endpoint_count, ...
   persistent_sample_count)
fprintf('Audit ledger: %s\n', output_file)
