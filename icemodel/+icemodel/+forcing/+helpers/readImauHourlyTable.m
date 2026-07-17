function [data, metadata] = readImauHourlyTable(filename)
   %READIMAUHOURLYTABLE Parse one IMAU hourly PANGAEA AWS table.
   %
   %  [data, metadata] = icemodel.forcing.helpers.readImauHourlyTable(filename)
   %
   % Role
   %  Source-specific parser for the Van Tiggelen et al. PANGAEA hourly IMAU
   %  S21/S22/S23 files. It skips the PANGAEA metadata block, maps corrected
   %  meteorological channels to icemodel-native names, and preserves source
   %  metadata needed by staging/import manifests.

   arguments
      filename (1, 1) string
   end

   lines = readlines(filename);
   header_idx = find(startsWith(lines, "Date/Time"), 1, 'first');
   if isempty(header_idx)
      error('icemodel:forcing:readImauHourlyTable:missingHeader', ...
         'No Date/Time table header found in %s.', filename);
   end

   header = split(lines(header_idx), sprintf('\t')).';
   raw = readRawTable(filename, header_idx + 1, numel(header));

   % S22/S23 encode a missing derived surface temperature as -1273.05 degC.
   % Normalize that source sentinel before converting the physical values to K;
   % otherwise it makes metchecks misclassify the whole Kelvin series as Celsius.
   tsfc = raw.col26;
   tsfc(tsfc == -1273.05) = NaN;

   % Map duplicate PANGAEA column labels by stable position. The source tables
   % use repeated Height/Surface elevation names, so sanitized headers alone are
   % not enough to identify variables safely.
   time = datetime(raw.col1, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss", ...
      'TimeZone', 'UTC');
   data = timetable(time, ...
      raw.col2, raw.col3, ...
      raw.col5 + 273.15, ...
      raw.col9, ...
      raw.col11, ...
      raw.col12, ...
      raw.col13 * 100, ...
      raw.col14, raw.col15, raw.col16, raw.col17, ...
      raw.col18, raw.col19, raw.col20, raw.col21, raw.col22, ...
      raw.col23, raw.col24, raw.col25, tsfc + 273.15, ...
      'VariableNames', {'height_1', 'height_2', 'tair', 'rh', ...
      'wspd', 'wdir', 'psfc', 'swd', 'swu', 'lwd', 'lwu', ...
      'surface_elevation_raw', 'surface_elevation', ...
      'surface_height_change_sonic', 'surface_height_change', ...
      'boom_height', 'lat_wgs84', 'lon_wgs84', 'albedo', 'tsfc'});
   data.Properties.DimensionNames{1} = 'time';

   % S22/S23 include an ablation draw-wire length column; S21 does not.
   if width(raw) >= 27
      data.wire_length = raw.col27;
   end

   % Preserve source metadata and row-derived coordinate summaries for later
   % staging without pretending row-varying coordinates are one exact point.
   metadata = sourceMetadata(filename, lines, header, data);
end

function raw = readRawTable(filename, first_data_line, n_columns)
   %READRAWTABLE Read PANGAEA rows with stable positional column names.
   opts = delimitedTextImportOptions('NumVariables', n_columns);
   opts.Delimiter = '\t';
   opts.DataLines = [first_data_line Inf];
   opts.VariableNames = "col" + string(1:n_columns);
   opts.VariableTypes = ["string", repmat("double", 1, n_columns - 1)];
   opts.ExtraColumnsRule = 'ignore';
   opts.EmptyLineRule = 'read';
   opts.ConsecutiveDelimitersRule = 'split';
   opts.LeadingDelimitersRule = 'ignore';
   raw = readtable(filename, opts);
end

function metadata = sourceMetadata(filename, lines, header, data)
   %SOURCEMETADATA Collect citation, DOI, event, and row-summary metadata.
   text = strjoin(lines(1:find(startsWith(lines, "*/"), 1, 'first')), newline);
   event = eventMetadata(text);
   metadata = struct( ...
      'filename', string(filename), ...
      'site_id', event.site_id, ...
      'citation', icemodel.forcing.helpers.regexpOnce( ...
      text, 'Citation:\s*(.*?)\n'), ...
      'doi', icemodel.forcing.helpers.regexpOnce(text, ...
      'Citation:.*?https://doi.org/([0-9.]+/PANGAEA\.[0-9]+)'), ...
      'bundle_doi', icemodel.forcing.helpers.regexpOnce(text, ...
      'bundled publication\].*?https://doi.org/([0-9.]+/PANGAEA\.[0-9]+)'), ...
      'raw_headers', header, ...
      'event', event, ...
      'row_summary', rowSummary(data), ...
      'variables', string(data.Properties.VariableNames));
end

function event = eventMetadata(text)
   %EVENTMETADATA Parse the PANGAEA Event(s) station line.
   line = icemodel.forcing.helpers.regexpOnce( ...
      text, 'Event\(s\):\s*(.*?)\n');
   event = struct('site_id', "", 'event_id', "", 'lat_wgs84', NaN, ...
      'lon_wgs84', NaN, 'elev_m', NaN, 'start', NaT, 'end', NaT);
   event.start.TimeZone = 'UTC';
   event.end.TimeZone = 'UTC';
   if line == ""
      return
   end

   event.event_id = icemodel.forcing.helpers.regexpOnce(line, '^([^\s]+)');
   event.site_id = icemodel.forcing.helpers.regexpOnce( ...
      line, '\((S[0-9]+)\)');
   event.lat_wgs84 = str2double(icemodel.forcing.helpers.regexpOnce( ...
      line, 'LATITUDE:\s*([-0-9.]+)'));
   event.lon_wgs84 = str2double(icemodel.forcing.helpers.regexpOnce( ...
      line, 'LONGITUDE:\s*([-0-9.]+)'));
   event.elev_m = str2double(icemodel.forcing.helpers.regexpOnce( ...
      line, 'ELEVATION:\s*([-0-9.]+)'));

   stamps = regexp(line, 'DATE/TIME (?:START|END):\s*([0-9T:-]+)', ...
      'tokens');
   if numel(stamps) >= 1
      event.start = parseTimestamp(stamps{1}{1});
   end
   if numel(stamps) >= 2
      event.end = parseTimestamp(stamps{2}{1});
   end
end

function summary = rowSummary(data)
   %ROWSUMMARY Capture row-derived coordinate and period ranges.
   summary = struct( ...
      'period', struct('start', data.time(1), 'end', data.time(end)), ...
      'lat_min', min(data.lat_wgs84, [], 'omitnan'), ...
      'lat_median', median(data.lat_wgs84, 'omitnan'), ...
      'lat_max', max(data.lat_wgs84, [], 'omitnan'), ...
      'lon_min', min(data.lon_wgs84, [], 'omitnan'), ...
      'lon_median', median(data.lon_wgs84, 'omitnan'), ...
      'lon_max', max(data.lon_wgs84, [], 'omitnan'));
end

function time = parseTimestamp(text)
   %PARSETIMESTAMP Convert PANGAEA ISO timestamps to UTC datetimes.
   time = datetime(text, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss", ...
      'TimeZone', 'UTC');
end
