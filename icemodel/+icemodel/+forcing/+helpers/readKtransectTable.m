function [data, metadata] = readKtransectTable(filename)
   %READKTRANSECTTABLE Parse one K-transect annual PANGAEA AWS table.
   %
   %  [data, metadata] = icemodel.forcing.helpers.readKtransectTable(filename)
   %
   % Role
   %  Source-specific parser for the Smeets et al. (2022) PANGAEA.947483
   %  K-transect annual files (AWS5/AWS6/AWS9/AWS10, 30-minute cadence). The
   %  shared PANGAEA ingest skips the metadata block and reads positional rows;
   %  this parser maps channels to icemodel-native names by source label
   %  because the column set varies across files: type-0 AWS10 files omit the
   %  "T tech" column and type-1 AWS5/AWS6 files add an "Ice melt" draw-wire
   %  column, so fixed positions cannot identify variables safely.
   %
   % Channel policy
   %  The AWS-generation flag (source "ID", 0=modular 2003-era type, 1=compact
   %  integrated type) ships as the aws_type channel because the acoustic
   %  height-ranger semantics differ by generation: type 0 records surface
   %  melt/snow height, type 1 records sensor-plus-snow height. The record is
   %  therefore kept source-faithful as height_rel rather than renamed to a
   %  canonical surface-height channel. Instrument diagnostics ("T body",
   %  "T tech") are not carried; the raw cache retains them. Battery voltage
   %  ("Vlog") is not carried either, but its once-yearly station-visit marker
   %  (value == 100) is preserved as metadata.visits because visits time the
   %  generation switch and height-pole re-mounts.
   %
   % See also: icemodel.forcing.helpers.readPangaeaTab,
   %  icemodel.forcing.buildKtransectData

   arguments
      filename (1, 1) string
   end

   % Shared PANGAEA ingest. The K-transect series line reads "[dataset
   % publication series]. PANGAEA, https://doi.org/..." rather than the IMAU
   % "bundled publication]" phrasing, and the station id rides the event id
   % (e.g. "K-transect_AWS9").
   [raw, header, source] = icemodel.forcing.helpers.readPangaeaTab( ...
      filename, site_id_pattern='^K-transect_(AWS[0-9]+)', ...
      bundle_pattern=['publication series\]\.\s*PANGAEA,\s*' ...
      'https://doi.org/([0-9.]+/PANGAEA\.[0-9]+)'], ...
      missing_header_error_id='icemodel:forcing:readKtransectTable:missingHeader');

   % Identify columns by their label prefix (text before the units bracket) so
   % the varying column sets and the non-ASCII degree sign in temperature
   % labels cannot break the mapping.
   labels = strtrim(extractBefore(header + " [", " ["));
   required = ["Date/Time", "dd", "ff", "SWD", "SWU", "LWD", "LWU", ...
      "TTT", "RH", "PPPP", "Height rel", "Vlog", "ID"];
   missing = required(~ismember(required, labels));
   if ~isempty(missing)
      error('icemodel:forcing:readKtransectTable:missingColumns', ...
         'missing required column(s) %s in %s', ...
         strjoin(missing, ", "), filename);
   end
   col = @(label) raw.("col" + find(labels == label, 1, 'first'));

   % Map source channels onto icemodel-native names and SI-consistent units
   % (degC -> K, hPa -> Pa). Radiative fluxes are already W m-2.
   time = parseRowTime(col("Date/Time"), filename);
   wspd = col("ff");
   wdir = col("dd");
   wdir(wspd < 0.1) = NaN;
   data = timetable(time, ...
      col("TTT") + 273.15, ...
      col("RH"), ...
      wspd, ...
      wdir, ...
      col("PPPP") * 100, ...
      col("SWD"), col("SWU"), col("LWD"), col("LWU"), ...
      col("Height rel"), ...
      col("ID"), ...
      'VariableNames', {'tair', 'rh', 'wspd', 'wdir', 'psfc', ...
      'swd', 'swu', 'lwd', 'lwu', 'height_rel', 'aws_type'});
   data.Properties.DimensionNames{1} = 'time';

   % Type-1 AWS5/AWS6 files add the draw-wire cumulative ice-melt column.
   if ismember("Ice melt", labels)
      data.ice_melt = col("Ice melt");
   end

   % Preserve the once-yearly maintenance-visit marker before dropping the
   % battery-voltage diagnostic channel.
   vlog = col("Vlog");
   visits = time(isfinite(vlog) & vlog == 100);

   % Preserve source metadata for staging manifests and DOI pinning; every
   % annual child carries its own citation and DOI in the header block.
   % Parse only the annual child's citation line; the generic metadata
   % extractor can otherwise fall through to the following series DOI.
   child_doi = icemodel.forcing.helpers.regexpOnce( ...
      string(source.citation), ...
      'https://doi.org/(10\.1594/PANGAEA\.[0-9]+)');
   valid_doi = ~isempty(regexp(char(child_doi), ...
      '^10\.1594/PANGAEA\.[0-9]+$', 'once'));
   if ~valid_doi
      error('icemodel:forcing:readKtransectTable:missingChildDoi', ...
         'missing or malformed annual child DOI in %s', filename);
   end
   metadata = struct( ...
      'filename', source.filename, ...
      'site_id', source.site_id, ...
      'citation', source.citation, ...
      'doi', child_doi, ...
      'bundle_doi', source.bundle_doi, ...
      'raw_headers', header, ...
      'event', source.event, ...
      'visits', visits, ...
      'variables', string(data.Properties.VariableNames));
end

function time = parseRowTime(text, filename)
   %PARSEROWTIME Parse K-transect timestamps with and without seconds.
   % Annual files post minute-resolution stamps ("2010-01-01T00:00"); accept a
   % seconds-bearing variant defensively since PANGAEA event lines carry one.
   time = NaT(size(text), 'TimeZone', 'UTC');
   nonempty = strlength(text) > 0;
   with_minutes = nonempty & strlength(text) == 16;
   with_seconds = nonempty & strlength(text) == 19;
   try
      time(with_minutes) = datetime(text(with_minutes), ...
         'InputFormat', "yyyy-MM-dd'T'HH:mm", 'TimeZone', 'UTC');
      time(with_seconds) = datetime(text(with_seconds), ...
         'InputFormat', "yyyy-MM-dd'T'HH:mm:ss", 'TimeZone', 'UTC');
   catch
      error('icemodel:forcing:readKtransectTable:badTimestamp', ...
         'unparseable Date/Time value(s) in %s', filename);
   end
   if any(ismissing(time) & nonempty)
      error('icemodel:forcing:readKtransectTable:badTimestamp', ...
         'unparseable Date/Time value(s) in %s', filename);
   end
end
