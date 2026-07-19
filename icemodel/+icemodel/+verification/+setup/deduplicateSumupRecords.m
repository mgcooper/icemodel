function [record, raw_rows, unique_rows, removed_rows] = ...
      deduplicateSumupRecords(record, variable)
   %DEDUPLICATESUMUPRECORDS Keep one row per SUMup scientific identity.
   %
   %  [record, raw_rows, unique_rows, removed_rows] = ...
   %     icemodel.verification.setup.deduplicateSumupRecords(record, variable)
   %
   % RECORD is one selected, numeric, pre-datetime SUMup density, temperature,
   % or SMB table. VARIABLE chooses the source-release-specific scientific
   % identity columns. Provenance aliases such as name_key/name/elevation and
   % measurement_id are deliberately excluded, so stable retention preserves
   % their first-row values without treating aliases as new measurements.
   %
   % MATLAB releases before R2026a treat repeated missing table values as
   % distinct in UNIQUE. This helper normalizes only a comparison copy: each
   % numeric missing value becomes zero beside a separate logical missingness
   % flag. Genuine zero therefore remains distinct from missing, while repeated
   % missing values compare equal without sentinel assumptions.

   arguments
      record table
      variable (1, 1) string
   end

   % Identity is variable-local; the source variable token is therefore not an
   % additional column in the comparison table.
   switch lower(variable)
      case "density"
         identity_names = ["latitude", "longitude", "start_depth", ...
            "stop_depth", "midpoint", "timestamp", "reference_key", ...
            "method_key", "density", "error"];
      case "temperature"
         identity_names = ["latitude", "longitude", "depth", "duration", ...
            "timestamp", "reference_key", "method_key", "temperature", ...
            "error"];
      case "smb"
         identity_names = ["latitude", "longitude", "start_date", ...
            "end_date", "start_year", "end_year", "reference_key", ...
            "method_key", "smb", "error"];
      otherwise
         error('icemodel:verification:deduplicateSumupRecords:badVariable', ...
            'variable must be density, temperature, or SMB')
   end

   % Fail before subsetting so a changed SUMup schema cannot silently weaken the
   % scientific identity definition.
   present = string(record.Properties.VariableNames);
   missing_names = identity_names(~ismember(identity_names, present));
   if ~isempty(missing_names)
      error('icemodel:verification:deduplicateSumupRecords:missingIdentity', ...
         'SUMup %s table lacks identity column(s): %s', ...
         variable, strjoin(missing_names, ', '))
   end

   raw_rows = height(record);
   unique_rows = raw_rows;
   removed_rows = 0;
   if raw_rows == 0
      return
   end

   % Work only on the identity projection. The returned rows are always selected
   % from RECORD itself, leaving first-row provenance bytes untouched.
   identity = record(:, cellstr(identity_names));
   for name = identity_names
      values = identity.(name);
      missing = ismissing(values);
      if ~any(missing)
         continue
      end

      % A paired flag makes missing unequal to genuine numeric zero but equal to
      % another missing value on R2024b/R2025b table UNIQUE.
      identity.(name + "_is_missing") = missing;
      values(missing) = 0;
      identity.(name) = values;
   end

   % UNIQUE supplies stable first-row indices and MATLAB's native numeric
   % equality, including the required +0/-0 equivalence.
   [~, first_rows] = unique(identity, "rows", "stable");
   record = record(first_rows, :);
   unique_rows = height(record);
   removed_rows = raw_rows - unique_rows;
end
