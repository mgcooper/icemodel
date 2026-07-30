function [Data, metadata] = buildKtransectData(station, kwargs)
   %BUILDKTRANSECTDATA Build native K-transect 30-minute AWS Data.
   %
   %  [Data, metadata] = icemodel.forcing.buildKtransectData("AWS9")
   %  [Data, metadata] = ... buildKtransectData("AWS9", source_dir=...)
   %
   % Reads every cached Smeets et al. (2022) PANGAEA.947483 annual tab file for
   % one K-transect station (AWS5/AWS6/AWS9/AWS10), merges the annual children
   % without duplicate timestamps or cadence conversion, and maps the channels
   % onto icemodel's canonical forcing/userdata names at the native 30-minute
   % cadence. Every annual child's DOI is preserved in metadata.children so
   % staging manifests can pin the exact upstream datasets.
   %
   % Albedo policy: the source publishes no albedo channel, so albedo is the
   % shared screened radiometer ratio (swd >= 10 W m-2, solar elevation > 20
   % degrees) with conservative recovered daily-collapse episodes removed; the
   % input radiation channels are never modified.
   %
   % Height policy: height_rel stays source-faithful with the aws_type
   % generation flag because its meaning differs by AWS type (surface
   % melt/snow height for type 0; sensor-plus-snow height for type 1), and
   % type-1 AWS5/AWS6 files add the ice_melt draw-wire channel.
   %
   % Source precipitation policy: the K-transect files have no precipitation
   % channel. `rainf` and `snowf` are explicit all-NaN placeholders so
   % `data2met` derives missing `ppt` rather than fabricating zero
   % precipitation.
   %
   % See also: icemodel.forcing.buildKtransectMet,
   %  icemodel.forcing.helpers.readKtransectTable

   arguments
      station (1, 1) string
      kwargs.source_dir (1, 1) string = ""
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.fillgaps (1, 1) logical = false
   end

   % Reject malformed public windows before resolving or reading any source file.
   % The shared mask validates again defensively when it applies the bounds.
   [window_start, window_end] = ...
      icemodel.internal.pairedWindow(kwargs.startdate, kwargs.enddate);

   station = canonicalStation(station);
   source_dir = icemodel.forcing.helpers.verificationSourceDir( ...
      kwargs.source_dir, "ktransect");
   filenames = icemodel.forcing.helpers.locateKtransectFiles( ...
      source_dir, station);

   % Read every annual child and merge on the native axis. Ordering by first
   % timestamp keeps the merge stable regardless of directory listing order.
   parts = cell(1, numel(filenames));
   children = repmat(childRecord(), 1, numel(filenames));
   for k = 1:numel(filenames)
      [parts{k}, part_metadata] = ...
         icemodel.forcing.helpers.readKtransectTable(filenames(k));
      if string(part_metadata.site_id) ~= station
         error('icemodel:forcing:buildKtransectData:stationMismatch', ...
            'child %s identifies %s, expected %s', filenames(k), ...
            string(part_metadata.site_id), station);
      end
      children(k) = childRecord(filenames(k), part_metadata, parts{k});
   end
   [~, order] = sort([children.start_time]);
   parts = parts(order);
   children = children(order);
   raw = mergeAnnualParts(parts, station);
   visits = mergeVisits(children);
   location = siteLocation(children);

   % Attach the published per-visit sensor-height records when the series
   % workbook is cached; donor evaluation and height-aware methods need the
   % wind/temperature sensor heights alongside the acoustic height_rel record.
   sensor_heights = stationHeights(source_dir, station);

   Time = raw.time;
   keep = icemodel.forcing.helpers.timeWindowMask( ...
      Time, window_start, window_end);
   if ~any(keep)
      error('icemodel:forcing:buildKtransectData:emptyWindow', ...
         'requested window does not overlap the %s record', station)
   end

   % Build a fresh timetable so the saved Data uses the standard `Time` row
   % dimension even though the source parser preserves PANGAEA's lower-case
   % `time` row label.
   Data = timetable(Time(keep));
   Data.tair = raw.tair(keep);
   Data.rh = raw.rh(keep);
   Data.wspd = raw.wspd(keep);
   Data.wdir = raw.wdir(keep);
   Data.psfc = raw.psfc(keep);
   Data.swd = raw.swd(keep);
   Data.swu = raw.swu(keep);
   Data.lwd = raw.lwd(keep);
   Data.lwu = raw.lwu(keep);
   Data.height_rel = raw.height_rel(keep);
   Data.aws_type = raw.aws_type(keep);
   if ismember("ice_melt", string(raw.Properties.VariableNames))
      Data.ice_melt = raw.ice_melt(keep);
   end
   Data.rainf = nan(height(Data), 1);
   Data.snowf = nan(height(Data), 1);

   % The source has no albedo channel; publish only the screened radiometer
   % ratio so low-sun and low-light noise never masquerades as albedo.
   [Data.albedo, albedo_qc_counts] = icemodel.forcing.helpers.sourceAlbedo( ...
      Data.swd, Data.swu, Time=Data.Time, ...
      latitude=location.lat_wgs84, longitude=location.lon_wgs84);

   % Radiative diagnostics are useful userdata/evaluation columns, while the
   % minimal met contract still comes from data2met at the met boundary.
   Data.swn = Data.swd - Data.swu;
   Data.lwn = Data.lwd - Data.lwu;
   Data.netr = Data.swn + Data.lwn;

   % The record can have real outage gaps. Keep the axis regular for icemodel
   % while preserving each outage as missing data rather than interpolating
   % across it.
   Data = regularizeHalfHourlyData(Data);

   % Preserve measured swd/swu but remove conservative recovered-collapse
   % episodes from the derived albedo and energy-balance channels. Running after
   % regularization gives the shared detector the canonical UTC timestamp grid.
   [transient_rows, transient_report] = ...
      icemodel.forcing.helpers.dailyAlbedoAnomalyFlags( ...
      Data.Time, Data.swd, Data.swu);
   albedo_transient_qc = rmfield(transient_report, 'diagnostics');
   if any(transient_rows)
      Data.albedo(transient_rows) = NaN;
      Data.swn(transient_rows) = NaN;
      Data.netr(transient_rows) = NaN;
   end

   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, ...
      fillgaps=kwargs.fillgaps);
   Data = orderDataColumns(Data);
   Data = icemodel.forcing.helpers.stampMetadata(Data);

   metadata = icemodel.forcing.helpers.columnizeMetadata( ...
      sourceMetadata(station, children, visits, location, sensor_heights, ...
      checks, albedo_qc_counts, albedo_transient_qc));
   Data = icemodel.forcing.helpers.attachLocationMetadata( ...
      Data, metadata.site_location);
   Data.Properties.UserData = metadata;
end

function Data = regularizeHalfHourlyData(Data)
   %REGULARIZEHALFHOURLYDATA Insert missing rows on the native 30-minute axis.
   if height(Data) < 2
      return
   end
   Data = sortrows(Data);
   Data = retime(Data, 'regular', 'fillwithmissing', ...
      'TimeStep', minutes(30));
end

function station = canonicalStation(station)
   %CANONICALSTATION Normalize K-transect station ids.
   token = upper(strrep(string(station), "_", ""));
   switch token
      case "AWS5"
         station = "AWS5";
      case "AWS6"
         station = "AWS6";
      case "AWS9"
         station = "AWS9";
      case "AWS10"
         station = "AWS10";
      otherwise
         error('icemodel:forcing:buildKtransectData:unknownStation', ...
            'unknown K-transect station: %s', station)
   end
end

function record = childRecord(filename, metadata, data)
   %CHILDRECORD Build one annual-child provenance record.
   if nargin == 0
      % Prototype for preallocation.
      record = struct('filename', "", 'year', NaN, 'doi', "", ...
         'citation', "", 'start_time', NaT, 'end_time', NaT, ...
         'lat_wgs84', NaN, 'lon_wgs84', NaN, 'elev_m', NaN, ...
         'visits', datetime.empty(0, 1), 'site_id', "");
      record.start_time.TimeZone = 'UTC';
      record.end_time.TimeZone = 'UTC';
      record.visits.TimeZone = 'UTC';
      return
   end
   record = childRecord();
   record.filename = string(filename);
   record.year = str2double(icemodel.forcing.helpers.regexpOnce( ...
      filename, '_([0-9]{4})\.tab$'));
   record.doi = metadata.doi;
   record.citation = metadata.citation;
   record.start_time = data.time(1);
   record.end_time = data.time(end);
   record.lat_wgs84 = metadata.event.lat_wgs84;
   record.lon_wgs84 = metadata.event.lon_wgs84;
   record.elev_m = metadata.event.elev_m;
   record.visits = metadata.visits;
   record.site_id = metadata.site_id;
end

function raw = mergeAnnualParts(parts, station)
   %MERGEANNUALPARTS Concatenate annual children without hidden edits.
   % Column sets differ across the AWS-generation boundary (ice_melt appears
   % only in type-1 AWS5/AWS6 files), so align on the union with explicit
   % missing values before concatenation.
   names = strings(1, 0);
   for k = 1:numel(parts)
      names = union(names, string(parts{k}.Properties.VariableNames), ...
         'stable');
   end
   for k = 1:numel(parts)
      absent = setdiff(names, string(parts{k}.Properties.VariableNames));
      for name = reshape(absent, 1, [])
         parts{k}.(name) = nan(height(parts{k}), 1);
      end
      parts{k} = parts{k}(:, cellstr(names));
   end
   raw = vertcat(parts{:});

   % Overlapping annual children would silently double-count samples through
   % retime; refuse rather than deduplicate.
   if numel(unique(raw.time)) ~= height(raw)
      error('icemodel:forcing:buildKtransectData:overlappingAnnualFiles', ...
         'duplicate timestamps across %s annual files', station)
   end
end

function visits = mergeVisits(children)
   %MERGEVISITS Pool the once-yearly maintenance-visit markers.
   visits = sort(vertcat(children.visits));
end

function location = siteLocation(children)
   %SITELOCATION Summarize per-year event coordinates for a moving station.
   % Ablation-zone stations ride the ice flow, so the annual event coordinates
   % drift; the median is the staged point and the ranges stay in metadata.
   lat = [children.lat_wgs84];
   lon = [children.lon_wgs84];
   elev = [children.elev_m];
   location = struct( ...
      'lat_wgs84', median(lat, 'omitnan'), ...
      'lon_wgs84', median(lon, 'omitnan'), ...
      'elev_m', median(elev, 'omitnan'));
   location = icemodel.forcing.helpers.projectLocation(location);
end

function heights = stationHeights(source_dir, station)
   %STATIONHEIGHTS Read this station's sensor-height records when cached.
   % The workbook is a series-level PANGAEA attachment; its absence is a
   % recorded condition, not an error, so tab-only caches still build.
   patterns = [ ...
      fullfile(source_dir, "metadata_GRL_AWS_*.xlsx")
      fullfile(source_dir, '**', "metadata_GRL_AWS_*.xlsx")];
   filename = "";
   for pattern = reshape(patterns, 1, [])
      hits = dir(pattern);
      hits = hits(~[hits.isdir]);
      if ~isempty(hits)
         filename = string(fullfile(hits(1).folder, hits(1).name));
         break
      end
   end
   if filename == ""
      heights = struct('present', false, 'file', "", ...
         'records', struct([]), 'notes', strings(0, 1), ...
         'reason', "sensor-height workbook not present in cache");
      return
   end
   [records, notes] = icemodel.forcing.helpers.readKtransectHeights( ...
      filename, station=station);
   heights = struct('present', true, 'file', filename, ...
      'records', records, 'notes', notes, 'reason', "");
end

function metadata = sourceMetadata(station, children, visits, location, ...
      sensor_heights, checks, albedo_qc_counts, albedo_transient_qc)
   %SOURCEMETADATA Assemble compact K-transect source provenance.
   % DOI and source basenames are durable identity; cache roots are local.
   for k = 1:numel(children)
      [~, name, ext] = fileparts(children(k).filename);
      children(k).filename = string(name) + string(ext);
   end
   if sensor_heights.present
      [~, name, ext] = fileparts(sensor_heights.file);
      sensor_heights.file = string(name) + string(ext);
   end
   metadata = struct( ...
      'source_file', strjoin([children.filename], ";"), ...
      'source_family', "ktransect", ...
      'source_package', "Smeets et al. (2022) PANGAEA K-transect AWS", ...
      'station', station, ...
      'doi', strjoin(unique([children.doi], 'stable'), ";"), ...
      'bundle_doi', "10.1594/PANGAEA.947483", ...
      'citation', children(1).citation, ...
      'license', "CC-BY-4.0", ...
      'children', children, ...
      'visits', visits, ...
      'sensor_heights', sensor_heights, ...
      'coordinate_summary', coordinateSummary(children), ...
      'site_location', location, ...
      'checks', checks, ...
      'albedo_qc_counts', albedo_qc_counts, ...
      'albedo_transient_qc', albedo_transient_qc, ...
      'channel_map', channelMap(), ...
      'albedo_policy', ...
         "no source albedo channel; albedo is the screened swu/swd radiometer ratio (swd >= 10 W m-2, solar elevation > 20 degrees) with conservative recovered daily collapse episodes removed", ...
      'height_policy', ...
         "height_rel kept source-faithful with the aws_type flag because acoustic ranger semantics differ by AWS generation; ice_melt draw-wire channel exists only for type-1 AWS5/AWS6", ...
      'precip_policy', ...
         "rainf = NaN placeholder; snowf = NaN placeholder; ppt = NaN placeholder via data2met because no precipitation channel exists");
end

function summary = coordinateSummary(children)
   %COORDINATESUMMARY Record per-year event coordinate ranges.
   lat = [children.lat_wgs84];
   lon = [children.lon_wgs84];
   elev = [children.elev_m];
   summary = struct( ...
      'lat_min', min(lat, [], 'omitnan'), ...
      'lat_max', max(lat, [], 'omitnan'), ...
      'lon_min', min(lon, [], 'omitnan'), ...
      'lon_max', max(lon, [], 'omitnan'), ...
      'elev_min', min(elev, [], 'omitnan'), ...
      'elev_max', max(elev, [], 'omitnan'));
end

function map = channelMap()
   %CHANNELMAP Record canonical-name to PANGAEA column mapping.
   map = struct( ...
      'tair', "TTT", ...
      'rh', "RH", ...
      'wspd', "ff", ...
      'wdir', "dd", ...
      'psfc', "PPPP", ...
      'swd', "SWD", ...
      'swu', "SWU", ...
      'lwd', "LWD", ...
      'lwu', "LWU", ...
      'height_rel', "Height rel", ...
      'ice_melt', "Ice melt", ...
      'aws_type', "ID");
end

function Data = orderDataColumns(Data)
   %ORDERDATACOLUMNS Keep K-transect Data files stable and metadata-mapped.
   preferred = ["tair", "rh", "wspd", "wdir", "psfc", "swd", "swu", ...
      "lwd", "lwu", "swn", "lwn", "netr", "albedo", "height_rel", ...
      "ice_melt", "aws_type", "rainf", "snowf"];
   names = string(Data.Properties.VariableNames);
   ordered = preferred(ismember(preferred, names));
   Data = Data(:, cellstr(ordered));
end
