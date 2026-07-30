function [Data, metadata] = buildImauHourlyData(station, kwargs)
   %BUILDIMAUHOURLYDATA Build native IMAU hourly AWS Data.
   %
   %  [Data, metadata] = icemodel.forcing.buildImauHourlyData("S21")
   %  [Data, metadata] = ... buildImauHourlyData("S21", source_dir=...)
   %
   % Reads the Van Tiggelen et al. PANGAEA hourly S21/S22/S23 tab files and
   % maps the corrected meteorological channels onto icemodel's canonical
   % forcing/userdata names. The first-pass IMAU staging family keeps these
   % hourly sites separate from the daily 19-station SEB QA product.
   %
   % Source precipitation policy: the hourly PANGAEA files have no rain/snow
   % split. `rainf` and `snowf` are explicit all-NaN placeholders so `data2met`
   % derives missing `ppt` rather than fabricating zero precipitation.
   %
   % See also: icemodel.forcing.buildImauHourlyMet,
   %  icemodel.forcing.helpers.readImauHourlyTable

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
      kwargs.source_dir, "imau");
   filename = icemodel.forcing.helpers.locateImauHourlyFile(source_dir, ...
      station);
   [raw, parsed] = icemodel.forcing.helpers.readImauHourlyTable(filename);
   location = siteLocation(parsed);

   Time = raw.time;
   keep = icemodel.forcing.helpers.timeWindowMask( ...
      Time, window_start, window_end);
   if ~any(keep)
      error('icemodel:forcing:buildImauHourlyData:emptyWindow', ...
         'requested window does not overlap %s', filename)
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
   source_albedo = raw.albedo(keep);
   Data.albedo = source_albedo;
   Data.tsfc = raw.tsfc(keep);
   Data.boom_height = raw.boom_height(keep);
   Data.surface_height = raw.surface_height_change(keep);
   Data.rainf = nan(height(Data), 1);
   Data.snowf = nan(height(Data), 1);

   % The source uses exactly 0.2 as a repeated lower-bound code, including
   % periods where reflected shortwave has collapsed. Screen that code and
   % low-light ratios without replacing the source's processed albedo; values
   % above 0.2 remain observational, including S23 summer dark-ice values.
   [ratio_validity, ratio_qc_counts] = ...
      icemodel.forcing.helpers.sourceAlbedo( ...
      Data.swd, Data.swu, Time=Data.Time, ...
      latitude=location.lat_wgs84, longitude=location.lon_wgs84);
   source_floor = isfinite(source_albedo) & source_albedo <= 0.2;
   Data.albedo(~isfinite(ratio_validity) | source_floor) = NaN;

   % Radiative diagnostics are useful userdata/evaluation columns, while the
   % minimal met contract still comes from data2met at the met boundary. Keep
   % raw swd/swu intact, but do not publish derived balances for bright samples
   % where both the invalid source floor and the raw flux ratio show a collapsed
   % reflected-shortwave channel. A separately plausible raw balance survives.
   Data.swn = Data.swd - Data.swu;
   raw_shortwave_ratio = Data.swu ./ Data.swd;
   invalid_shortwave_balance = source_floor & Data.swd >= 10 ...
      & (~isfinite(raw_shortwave_ratio) | raw_shortwave_ratio <= 0.2);
   Data.swn(invalid_shortwave_balance) = NaN;
   Data.lwn = Data.lwd - Data.lwu;
   Data.netr = Data.swn + Data.lwn;
   albedo_qc_counts = ratio_qc_counts;
   albedo_qc_counts.source_floor = nnz(source_floor);
   albedo_qc_counts.total = nnz(isfinite(source_albedo) ...
      & (~isfinite(ratio_validity) | source_floor));
   albedo_qc_counts.derived_balance = nnz(invalid_shortwave_balance);

   % The PANGAEA hourly record can have real outage gaps. Keep the met axis
   % regular for icemodel while preserving the outage as missing data rather
   % than interpolating across it.
   Data = regularizeHourlyData(Data);

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
      sourceMetadata(filename, station, parsed, checks, albedo_qc_counts, ...
      albedo_transient_qc));
   Data = icemodel.forcing.helpers.attachLocationMetadata( ...
      Data, metadata.site_location);
   Data.Properties.UserData = metadata;
end

function Data = regularizeHourlyData(Data)
   %REGULARIZEHOURLYDATA Insert missing rows on the native hourly axis.
   if height(Data) < 2
      return
   end
   Data = sortrows(Data);
   Data = retime(Data, 'regular', 'fillwithmissing', ...
      'TimeStep', hours(1));
end

function station = canonicalStation(station)
   %CANONICALSTATION Normalize first-pass IMAU station ids.
   token = upper(strrep(string(station), "_", ""));
   switch token
      case "S21"
         station = "S21";
      case "S22"
         station = "S22";
      case "S23"
         station = "S23";
      otherwise
         error('icemodel:forcing:buildImauHourlyData:unknownStation', ...
            'unknown IMAU hourly station: %s', station)
   end
end

function Data = orderDataColumns(Data)
   %ORDERDATACOLUMNS Keep IMAU Data files stable and metadata-mapped.
   preferred = ["tair", "rh", "wspd", "wdir", "psfc", "swd", "swu", ...
      "lwd", "lwu", "swn", "lwn", "netr", "albedo", "tsfc", ...
      "boom_height", "surface_height", "rainf", "snowf"];
   names = string(Data.Properties.VariableNames);
   ordered = preferred(ismember(preferred, names));
   Data = Data(:, cellstr(ordered));
end

function metadata = sourceMetadata(filename, station, parsed, checks, ...
      albedo_qc_counts, albedo_transient_qc)
   %SOURCEMETADATA Assemble compact IMAU source provenance.
   location = siteLocation(parsed);
   metadata = struct( ...
      'source_file', filename, ...
      'source_family', "imau", ...
      'source_package', "Van Tiggelen et al. PANGAEA hourly AWS", ...
      'station', station, ...
      'doi', string(parsed.doi), ...
      'bundle_doi', string(parsed.bundle_doi), ...
      'citation', string(parsed.citation), ...
      'event', parsed.event, ...
      'source_variables', string(parsed.raw_headers(:)), ...
      'canonical_variables', string(parsed.variables(:)), ...
      'channel_map', channelMap(), ...
      'site_location', location, ...
      'checks', checks, ...
      'albedo_qc_counts', albedo_qc_counts, ...
      'albedo_transient_qc', albedo_transient_qc, ...
      'albedo_policy', ...
         "source Alb frac retained only where solar elevation > 20 degrees, swd >= 10 W m-2, swu > 0, and albedo > 0.2; repeated 0.2 lower-bound code and conservative recovered daily collapse episodes remain NaN", ...
      'shortwave_balance_policy', ...
         "source swd/swu preserved; derived swn/netr remain NaN where raw swu/swd <= 0.2 confirms a source-floor collapse and on dailyAlbedoAnomalyFlags episodes", ...
      'precip_policy', ...
         "rainf = NaN placeholder; snowf = NaN placeholder; ppt = NaN placeholder via data2met because no precipitation channel exists");
end

function location = siteLocation(parsed)
   %SITELOCATION Prefer row-median station coordinates over event nominal coords.
   lat = parsed.row_summary.lat_median;
   lon = parsed.row_summary.lon_median;
   if ~isfinite(lat)
      lat = parsed.event.lat_wgs84;
   end
   if ~isfinite(lon)
      lon = parsed.event.lon_wgs84;
   end
   location = struct('lat_wgs84', lat, 'lon_wgs84', lon, ...
      'elev_m', parsed.event.elev_m);
   location = icemodel.forcing.helpers.projectLocation(location);
end

function map = channelMap()
   %CHANNELMAP Record canonical-name to PANGAEA hourly column mapping.
   map = struct( ...
      'tair', "TTT corrected at 2m height", ...
      'rh', "RH corrected at 2m height", ...
      'wspd', "FF10 corrected at 10m height", ...
      'wdir', "dd", ...
      'psfc', "PPPP", ...
      'swd', "SWD", ...
      'swu', "SWU", ...
      'lwd', "LWD", ...
      'lwu', "LWU", ...
      'albedo', "Alb frac", ...
      'tsfc', "Surf temp", ...
      'boom_height', "Height of sensor AWS boom above surface", ...
      'surface_height', "Surf elev change sonic and ablation");
end
