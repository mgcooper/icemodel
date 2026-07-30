function [raw, header, metadata] = readPangaeaTab(filename, kwargs)
   %READPANGAEATAB Read one PANGAEA tab-delimited dataset export.
   %
   %  [raw, header, metadata] = icemodel.forcing.helpers.readPangaeaTab(filename)
   %  [raw, header, metadata] = icemodel.forcing.helpers.readPangaeaTab( ...
   %     filename, site_id_pattern='\((S[0-9]+)\)')
   %
   % Role
   %  Shared, source-agnostic ingest for PANGAEA .tab exports. Every export
   %  carries the same skeleton: a "/* ... */" metadata block (citation with the
   %  per-dataset DOI, an optional bundle/series DOI, and an Event(s) station
   %  line), a tab-delimited header starting with "Date/Time", and numeric data
   %  rows. This helper owns that skeleton once so per-source readers
   %  (readImauHourlyTable, readKtransectTable, ...) own only their column
   %  mapping and unit conversions.
   %
   % Name-value
   %  site_id_pattern : string. Regex with one capture group applied to the
   %     Event(s) line to recover the station identifier, e.g.
   %     '\((S[0-9]+)\)' for the IMAU hourly network or
   %     '^K-transect_(AWS[0-9]+)' for the K-transect series. Empty leaves
   %     metadata.site_id and event.site_id blank.
   %  bundle_pattern : string. Regex with one capture group recovering the
   %     bundle or publication-series DOI from the metadata block. The default
   %     matches the "bundled publication]" phrasing used by the IMAU hourly
   %     exports; series-style exports pass their own phrasing.
   %  missing_header_error_id : string. Error identifier raised when the file
   %     has no "Date/Time" table header, so per-source readers keep their
   %     established error ids for callers and tests.
   %
   % Returns
   %  raw : table. Tabular section with stable positional names col1..colN
   %     (col1 is the timestamp string; remaining columns are double). Blank
   %     cells become NaN, never zero.
   %  header : string row. The source column labels, split on tabs.
   %  metadata : struct with fields filename, site_id, citation, doi,
   %     bundle_doi, and event (event_id, site_id, lat_wgs84, lon_wgs84,
   %     elev_m, start, end as UTC datetimes).
   %
   % See also: icemodel.forcing.helpers.readImauHourlyTable,
   %  icemodel.forcing.helpers.readKtransectTable

   arguments
      filename (1, 1) string
      kwargs.site_id_pattern (1, 1) string = ""
      kwargs.bundle_pattern (1, 1) string = ...
         "bundled publication\].*?https://doi.org/([0-9.]+/PANGAEA\.[0-9]+)"
      kwargs.missing_header_error_id (1, 1) string = ...
         "icemodel:forcing:readPangaeaTab:missingHeader"
   end

   % Locate the tabular section; every PANGAEA export labels it "Date/Time".
   lines = readlines(filename);
   header_idx = find(startsWith(lines, "Date/Time"), 1, 'first');
   if isempty(header_idx)
      error(kwargs.missing_header_error_id, ...
         'No Date/Time table header found in %s.', filename);
   end

   % Split the source labels and read the rows with positional names so
   % duplicate or unit-decorated source labels cannot misidentify columns.
   header = split(lines(header_idx), sprintf('\t')).';
   raw = readRawTable(filename, header_idx + 1, numel(header));

   % Parse citation, DOIs, and the Event(s) station line from the metadata
   % block; the per-dataset DOI is the first DOI after "Citation:".
   text = strjoin(lines(1:find(startsWith(lines, "*/"), 1, 'first')), newline);
   event = eventMetadata(text, kwargs.site_id_pattern);
   metadata = struct( ...
      'filename', string(filename), ...
      'site_id', event.site_id, ...
      'citation', icemodel.forcing.helpers.regexpOnce( ...
      text, 'Citation:\s*(.*?)\n'), ...
      'doi', icemodel.forcing.helpers.regexpOnce(text, ...
      'Citation:.*?https://doi.org/([0-9.]+/PANGAEA\.[0-9]+)'), ...
      'bundle_doi', icemodel.forcing.helpers.regexpOnce(text, ...
      kwargs.bundle_pattern), ...
      'event', event);
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

function event = eventMetadata(text, site_id_pattern)
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

   % The event id is the leading token; the station id needs the caller's
   % source-specific pattern because networks encode it differently.
   event.event_id = icemodel.forcing.helpers.regexpOnce(line, '^([^\s]+)');
   if site_id_pattern ~= ""
      event.site_id = icemodel.forcing.helpers.regexpOnce( ...
         line, site_id_pattern);
   end
   event.lat_wgs84 = str2double(icemodel.forcing.helpers.regexpOnce( ...
      line, 'LATITUDE:\s*([-0-9.]+)'));
   event.lon_wgs84 = str2double(icemodel.forcing.helpers.regexpOnce( ...
      line, 'LONGITUDE:\s*([-0-9.]+)'));
   event.elev_m = str2double(icemodel.forcing.helpers.regexpOnce( ...
      line, 'ELEVATION:\s*([-0-9.]+)'));

   % Both stamps are optional; a start without an end stays NaT-terminated.
   % Series-style events (e.g. K-transect) post one bare "DATE/TIME:" stamp,
   % which maps to the start.
   stamps = regexp(line, 'DATE/TIME(?: START| END)?:\s*([0-9T:-]+)', ...
      'tokens');
   if numel(stamps) >= 1
      event.start = parseTimestamp(stamps{1}{1});
   end
   if numel(stamps) >= 2
      event.end = parseTimestamp(stamps{2}{1});
   end
end

function time = parseTimestamp(text)
   %PARSETIMESTAMP Convert PANGAEA ISO timestamps to UTC datetimes.
   time = datetime(text, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss", ...
      'TimeZone', 'UTC');
end
