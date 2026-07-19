function [Data, metadata] = buildSamimiDye2Data(kwargs)
   %BUILDSAMIMIDYE2DATA Build RetMIP Dye-2 2016 native Data from Samimi AWS.
   %
   %  [Data, metadata] = icemodel.forcing.buildSamimiDye2Data()
   %  [Data, metadata] = ... buildSamimiDye2Data(source_dir=...)
   %
   % Reads the Samimi/Marshall Dye-2 summer 2016 AWS workbook and maps the
   % native columns onto icemodel's canonical forcing/userdata names. The source
   % record is half-hourly and remains native here. Shared artifact writers own
   % the public 15-minute met and hourly userdata output cadences.
   %
   % Source precipitation policy: this workbook has no precipitation channel.
   % `rainf`, `snowf`, and the `ppt` derived by `data2met` are explicit NaN
   % placeholders rather than fabricated zero precipitation.
   %
   % See also: icemodel.forcing.buildSamimiDye2Met,
   %  icemodel.forcing.data2met

   arguments
      kwargs.source_dir (1, 1) string = ""
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.fillgaps (1, 1) logical = false
   end

   % Reject malformed public windows before resolving or reading any workbook.
   % The shared mask validates again defensively when it applies the bounds.
   [window_start, window_end] = ...
      icemodel.internal.pairedWindow(kwargs.startdate, kwargs.enddate);

   source_dir = icemodel.forcing.helpers.verificationSourceDir( ...
      kwargs.source_dir, ["retmip", "samimi"]);
   filename = icemodel.forcing.helpers.locateSamimiDye2Workbook(source_dir);
   [raw, source_names] = readAwsSheet(filename);

   Time = normalizeTime(raw.(source_names(1)));
   keep = icemodel.forcing.helpers.timeWindowMask( ...
      Time, window_start, window_end);
   if ~any(keep)
      error('icemodel:forcing:buildSamimiDye2Data:emptyWindow', ...
         'requested window does not overlap %s', filename)
   end

   % Map source units to canonical icemodel units at native half-hourly cadence.
   % The workbook stores air temperature in degC, pressure in hPa, and its
   % derived continuous snow-depth channel `dsnow` in centimetres. The source
   % `surface` column is also centimetres, not surface temperature.
   Data = timetable(Time(keep));
   Data.tair = celsiusToKelvin(raw.airT(keep));
   Data.rh = raw.RH(keep);
   Data.psfc = raw.P(keep) .* 100;
   native_swd = raw.("SW in")(keep);
   native_swu = raw.("SW out")(keep);
   shortwave_negative_counts = struct( ...
      'swd', nnz(native_swd < 0), 'swu', nnz(native_swu < 0));
   Data.swd = max(native_swd, 0);
   Data.swu = max(native_swu, 0);
   Data.lwd = raw.("LW in")(keep);
   Data.lwu = raw.("LW out")(keep);
   Data.wspd = raw.wind(keep);
   Data.wdir = raw.winddir(keep);
   Data.rainf = nan(height(Data), 1);
   Data.snowf = nan(height(Data), 1);
   Data.snow_depth = raw.dsnow(keep) ./ 100;
   Data.swn = Data.swd - Data.swu;
   Data.lwn = Data.lwd - Data.lwu;
   Data.netr = Data.swn + Data.lwn;
   location = stationLocation();
   [Data.albedo, albedo_qc_counts] = ...
      icemodel.forcing.helpers.sourceAlbedo( ...
      Data.swd, Data.swu, minimum=0.3, Time=Data.Time, ...
      latitude=location.lat_wgs84, longitude=location.lon_wgs84);

   % Apply source-level checks before any artifact-writer resampling.
   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, ...
      fillgaps=kwargs.fillgaps);
   Data = icemodel.forcing.helpers.stampMetadata(Data);

   metadata = icemodel.forcing.helpers.columnizeMetadata( ...
      sourceMetadata(filename, source_names, checks, ...
      shortwave_negative_counts, albedo_qc_counts));
   Data = icemodel.forcing.helpers.attachLocationMetadata( ...
      Data, metadata.site_location);
   Data.Properties.UserData = metadata;
end

function [raw, source_names] = readAwsSheet(filename)
   %READAWSSHEET Read the Samimi AWS table preserving workbook names.
   opts = detectImportOptions(filename, 'VariableNamingRule', 'preserve');
   raw = readtable(filename, opts);
   source_names = string(raw.Properties.VariableNames);
   required = ["airT", "RH", "P", "SW in", "SW out", "LW in", ...
      "LW out", "wind", "winddir", "dsnow"];
   missing = required(~ismember(required, source_names));
   if ~isempty(missing)
      error('icemodel:forcing:buildSamimiDye2Data:missingVariables', ...
         'Samimi workbook missing required variable(s): %s', ...
         strjoin(missing, ', '))
   end
end

function Time = normalizeTime(value)
   %NORMALIZETIME Convert workbook timestamps to UTC datetimes.
   if isa(value, 'datetime')
      Time = value;
   elseif isnumeric(value)
      Time = datetime(value, 'ConvertFrom', 'excel');
   else
      Time = datetime(value);
   end
   Time = Time(:);
   if isempty(Time.TimeZone)
      Time.TimeZone = 'UTC';
   end
end

function value = celsiusToKelvin(value)
   %CELSIUSTOKELVIN Convert degC source temperatures to Kelvin.
   value = value + 273.15;
end

function metadata = sourceMetadata(filename, source_names, checks, ...
      shortwave_negative_counts, albedo_qc_counts)
   %SOURCEMETADATA Assemble compact Samimi source provenance.
   metadata = struct( ...
      'source_file', filename, ...
      'source_family', "samimi", ...
      'source_package', "Samimi Dye-2 summer 2016 AWS workbook", ...
      'station', "Dye-2", ...
      'source_variables', source_names(:), ...
      'channel_map', channelMap(), ...
       'site_location', icemodel.forcing.helpers.projectLocation(stationLocation()), ...
       'time_step_policy', ...
       "native half-hourly workbook; artifact writers select output cadence", ...
       'shortwave_policy', ...
       "negative SW in/SW out sensor offsets clipped to physical zero", ...
        'shortwave_negative_counts', shortwave_negative_counts, ...
        'albedo_qc_counts', albedo_qc_counts, ...
        'albedo_policy', ...
        "swu/swd only where solar elevation > 20 degrees, swd >= 10 W m-2, swu > 0, and albedo >= 0.3; rejected accumulation-site ratios remain NaN", ...
        'snow_depth_policy', "dsnow centimetres converted to metres", ...
       'checks', checks, ...
       'precip_policy', ...
          "rainf = NaN placeholder; snowf = NaN placeholder; " + ...
          "ppt = NaN placeholder via data2met because no precipitation channel exists");
end

function map = channelMap()
   %CHANNELMAP Record canonical-name to workbook-name mapping.
   map = struct( ...
      'tair', "airT", ...
      'rh', "RH", ...
      'psfc', "P", ...
      'swd', "SW in", ...
      'swu', "SW out", ...
      'lwd', "LW in", ...
      'lwu', "LW out", ...
      'wspd', "wind", ...
      'wdir', "winddir", ...
      'snow_depth', "dsnow");
end

function location = stationLocation()
   %STATIONLOCATION Return RetMIP Dye-2 coordinates.
   location = struct('lat_wgs84', 66.48001, 'lon_wgs84', -46.27889, ...
      'elev_m', 2165);
end
