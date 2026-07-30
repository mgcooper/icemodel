function [data, metadata] = readRetmipProtocolTable(filename, kwargs)
   %READRETMIPPROTOCOLTABLE Read a RetMIP tab-delimited protocol time series.
   %
   %  [data, metadata] = ...
   %     icemodel.verification.setup.readRetmipProtocolTable(filename)
   %
   % Role
   %  Thin RetMIP protocol parser. The official surface protocol variables are
   %  userdata/comparison inputs; this helper validates the time cadence before
   %  any importer records them in a manifest.

   arguments
      filename (1, 1) string
      kwargs.expected_hours (1, 1) double {mustBePositive} = 3
   end

   % RetMIP source files use semicolons and quote each data row, while tests and
   % ad hoc fixtures may use tabs. Normalize those two simple delimited forms
   % before validating cadence.
   data = readSimpleDelimitedTable(filename);

   % Convert the first time-like column to UTC datetimes; RetMIP protocol files
   % are time series, so a missing time axis is a malformed source file.
   [time_name, time_values] = readTimeColumn(data);
   data.(time_name) = time_values;

   % Cadence validation is intentionally simple and strict: the protocol series
   % is expected at a uniform 3-hour posting.
   dt_hours = hours(diff(time_values));
   if ~isempty(dt_hours) && any(abs(dt_hours - kwargs.expected_hours) > 1e-9)
      error('icemodel:verification:readRetmipProtocolTable:badCadence', ...
         'RetMIP protocol table is not %.3g-hourly', kwargs.expected_hours);
   end

   % Normalize protocol userdata to the comparison names used by manifests while
   % retaining the source variable inventory in metadata.
   raw_variables = string(data.Properties.VariableNames);
   data = canonicalizeProtocolVariables(data, kwargs.expected_hours);

   metadata = struct( ...
      'filename', filename, ...
      'time_variable', time_name, ...
      'timestep_hours', kwargs.expected_hours, ...
      'mass_flux_policy', ...
         "source mm w.e. timestep amounts converted to canonical mWE/h", ...
      'mass_flux_conversion_factor', 1e-3 / kwargs.expected_hours, ...
      'variables', string(data.Properties.VariableNames), ...
      'raw_variables', raw_variables);
end

function data = canonicalizeProtocolVariables(data, timestep_hours)
   %CANONICALIZEPROTOCOLVARIABLES Map RetMIP protocol names to eval variables.
   names = string(data.Properties.VariableNames);
   data = copyAs(data, names, "Tsurf_K", "tsfc", 0, 1);
   data = copyAs(data, names, "Tsurf_degC", "tsfc", 273.15, 1);

   % RetMIP posts mass as millimetres water equivalent accumulated over each
   % three-hour interval. Canonical melt-style diagnostics are rates [mWE/h].
   mass_flux_scale = 1e-3 / timestep_hours;
   data = copyAs(data, names, "melt_mmweq", "melt", 0, mass_flux_scale);
   data = copyAs(data, names, ["acc_subl_mmweq", "net_acc_mmweq"], ...
      "snowf_subl", 0, mass_flux_scale);

   % The protocol does not split precipitation phase, so keep explicit empty
   % placeholders for consumers that expect the standard met-like channels.
   for name = ["ppt", "rainf", "snowf"]
      if ~ismember(name, string(data.Properties.VariableNames))
         data.(name) = NaN(height(data), 1);
      end
   end
end

function data = copyAs(data, names, sources, target, offset, scale)
   %COPYAS Copy the first present source variable to a canonical target.
   if ismember(target, string(data.Properties.VariableNames))
      return
   end

   hit = sources(find(ismember(sources, names), 1));
   if isempty(hit)
      return
   end

   data.(target) = (data.(hit) + offset) .* scale;
   if hit ~= target
      data.(hit) = [];
   end
end

function data = readSimpleDelimitedTable(filename)
   %READSIMPLEDELIMITEDTABLE Parse RetMIP's small quoted .tab tables.
   lines = strip(readlines(filename));
   lines = lines(strlength(lines) > 0);
   if isempty(lines)
      error('icemodel:verification:readRetmipProtocolTable:emptyFile', ...
         'RetMIP protocol table is empty');
   end

   % Choose the delimiter from the header, then remove RetMIP's row-level quotes
   % before splitting the compact numeric/time table.
   delimiter = sprintf('\t');
   if contains(lines(1), ";")
      delimiter = ";";
   end
   lines = erase(lines, '"');
   names = split(lines(1), delimiter).';
   raw = strings(numel(lines) - 1, numel(names));
   for n = 2:numel(lines)
      raw(n - 1, :) = split(lines(n), delimiter).';
   end

   % Keep the first column as strings for datetime parsing; numeric protocol
   % variables become double columns.
   data = table();
   for k = 1:numel(names)
      values = raw(:, k);
      nums = str2double(values);
      if k > 1 && all(~isnan(nums))
         data.(matlab.lang.makeValidName(names(k))) = nums;
      else
         data.(matlab.lang.makeValidName(names(k))) = values;
      end
   end
   data.Properties.VariableNames = cellstr(names);
end

function [name, values] = readTimeColumn(data)
   %READTIMECOLUMN Convert the first time-like table column to datetime.
   names = string(data.Properties.VariableNames);
   idx = find(ismember(lower(names), ["time", "datetime", "date"]), 1);
   if isempty(idx)
      error('icemodel:verification:readRetmipProtocolTable:noTime', ...
         'RetMIP protocol table needs a time/datetime/date column');
   end

   name = char(names(idx));
   raw = data.(name);
   if isdatetime(raw)
      values = raw;
   else
      try
         values = datetime(string(raw), 'TimeZone', 'UTC');
      catch
         values = datetime(string(raw), 'InputFormat', 'dd-MMM-uuuu HH:mm:ss', ...
            'Locale', 'en_US', 'TimeZone', 'UTC');
      end
   end
end
