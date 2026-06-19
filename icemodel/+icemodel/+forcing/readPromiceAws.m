function [aws, metadata] = readPromiceAws(site, kwargs)
   %READPROMICEAWS Read a PROMICE v3 AWS NetCDF into icemodel-named channels.
   %
   %  [aws, metadata] = icemodel.forcing.readPromiceAws(site)
   %  [aws, metadata] = ... readPromiceAws(site, source_dir=..., ...
   %     timescale="hourly", startdate=..., enddate=...)
   %
   % Reads one PROMICE automatic-weather-station NetCDF (the v3 bundle
   % files <STATION>_hour_v03.nc / <STATION>_day_v03.nc) and returns a
   % timetable with icemodel-standard channel names. The v3 NetCDF files
   % are CF-converted upstream (air/surface temperature in kelvin,
   % pressure in pascal, radiation in W m-2, humidity and cloud fraction
   % as fractions), so the only unit change applied here is relative
   % humidity fraction -> percent, the icemodel convention. This reader
   % supersedes the legacy 46-column .txt parser (runoff readPromiceAWS)
   % and the per-station read_save scripts it replaced.
   %
   % Channel mapping (NetCDF -> output):
   %    ta -> tair [K]            pa  -> psfc [Pa]
   %    rh -> rh [%]              wspd/wdir -> wspd [m s-1] / wdir [deg]
   %    fsds_cor (fallback fsds) -> swd [W m-2]
   %    fsus_cor (fallback fsus) -> swu [W m-2]
   %    flds -> lwd, flus -> lwu [W m-2]
   %    alb -> albedo [-]         clt -> cfrac [-]
   %    ts -> tsfc [K]            hfss/hfls -> shf/lhf [W m-2]
   %    height_sensor_boom -> boom_height [m]
   %    height_stakes -> stake_height [m]
   %    depth_pressure_transducer_cor -> transducer_depth [m]
   %    tice1..tice8 -> tice1..tice8 [K]
   %    orog -> elev [m] (GPS elevation time series)
   %
   % Inputs
   %  site - station name. Accepts the canonical id ("KAN_M") or the
   %         compact lowercase alias ("kanm"); matching ignores case and
   %         underscores against the station files in source_dir.
   %
   % Name-value
   %  source_dir : directory containing the station NetCDF files, or its
   %      parent holding hourly/ and daily/ subdirectories. Defaults to
   %      the gitignored cache data/forcing/promice under the icemodel
   %      data root. The reference external-drive layout is
   %      /Volumes/S03/DATA/greenland/geus/aws/v3.
   %  timescale : "hourly" (default) or "daily".
   %  startdate, enddate : optional datetime window; default full range.
   %
   % Outputs
   %  aws      - timetable, UTC time axis, channels above (a channel
   %             missing from the file is omitted from the output)
   %  metadata - provenance struct: source file, station name, lat, lon,
   %             median GPS elevation, row count, time bounds
   %
   % See also: icemodel.forcing.buildPromiceMet,
   %  icemodel.forcing.buildPromiceData, icemodel.forcing.fillPromiceAlbedo

   arguments
      site (1, 1) string
      kwargs.source_dir (1, 1) string = ""
      kwargs.timescale (1, 1) string ...
         {mustBeMember(kwargs.timescale, ["hourly", "daily"])} = "hourly"
      kwargs.startdate = ""
      kwargs.enddate = ""
   end

   filename = locateStationFile(site, kwargs.source_dir, kwargs.timescale);

   % Time axis: v3 files encode seconds since the unix epoch, stamped at
   % the bin center (00:30 for the 00:00-01:00 hourly average; 12:00 for
   % daily). Snap stamps to the interval start, the icemodel met-file
   % convention (and the legacy .txt HourOfDay labeling).
   t = ncread(filename, 'time');
   time_units = ncreadatt(filename, 'time', 'units');
   assert(startsWith(time_units, 'seconds since 1970'), ...
      'unexpected time units in %s: %s', filename, time_units)
   Time = datetime(1970, 1, 1, 'TimeZone', 'UTC') + seconds(double(t));
   if kwargs.timescale == "hourly"
      Time = dateshift(Time, 'start', 'hour');
   else
      Time = dateshift(Time, 'start', 'day');
   end

   % Channel map: output name -> primary nc variable (and optional
   % fallback for the tilt-corrected radiation channels).
   tice = arrayfun(@(n) sprintf('tice%d', n), 1:8, 'UniformOutput', false)';
   channels = [{
      'tair',             'ta',                            ''
      'psfc',             'pa',                            ''
      'rh',               'rh_wrtwater',                   ''
      'wspd',             'wspd',                          ''
      'wdir',             'wdir',                          ''
      'swd',              'fsds_cor',                      'fsds'
      'swu',              'fsus_cor',                      'fsus'
      'lwd',              'flds',                          ''
      'lwu',              'flus',                          ''
      'albedo',           'alb',                           ''
      'cfrac',            'clt',                           ''
      'tsfc',             'ts',                            ''
      'shf',              'hfss',                          ''
      'lhf',              'hfls',                          ''
      'boom_height',      'height_sensor_boom',            ''
      'stake_height',     'height_stakes',                 ''
      'transducer_depth', 'depth_pressure_transducer_cor', ''
      'elev',             'orog',                          ''
      }; [tice, tice, repmat({''}, 8, 1)]];

   info = ncinfo(filename);
   available = string({info.Variables.Name});
   aws = timetable(Time);
   for n = 1:size(channels, 1)
      outname = channels{n, 1};
      ncname = channels{n, 2};
      fallback = string(channels{n, 3});
      if ~ismember(ncname, available)
         continue
      end
      data = double(ncread(filename, ncname));
      if fallback ~= "" && ismember(fallback, available)
         data = preferPrimary(data, double(ncread(filename, fallback)));
      end
      aws.(outname) = data;
   end

   % Relative humidity: the v3 rh_wrtwater channel is already in percent
   % despite its units="1" attribute (verified against the legacy met
   % artifacts: raw rh_wrtwater matches to ~0.4% mean), and it is the
   % standard RH-with-respect-to-water series every legacy artifact used,
   % consistent with the icemodel.vapor (w.r.t. water above Tf) convention.
   % The separate nc "rh" channel is RH with respect to ICE, which
   % legitimately exceeds 100% in subfreezing conditions (full-record
   % median ~119%) - a valid but different quantity, deliberately not used
   % here. (Both channels carry units="1" but are stored as percent.)

   % Optional window subset.
   keep = true(height(aws), 1);
   if ~strcmp(string(kwargs.startdate), "")
      keep = keep & aws.Time >= ensureUtc(kwargs.startdate);
   end
   if ~strcmp(string(kwargs.enddate), "")
      keep = keep & aws.Time <= ensureUtc(kwargs.enddate);
   end
   if ~any(keep)
      error('icemodel:forcing:readPromiceAws:emptyWindow', ...
         'no samples in the requested window for station %s', site)
   end
   aws = aws(keep, :);

   metadata = struct( ...
      'site', site, ...
      'source_file', string(filename), ...
      'station_name', string(deblank(ncread(filename, 'station_name')')), ...
      'lat', double(ncread(filename, 'latitude')), ...
      'lon', double(ncread(filename, 'longitude')), ...
      'elev', median(aws.elev, 'omitnan'), ...
      'n_rows', height(aws), ...
      'window_start', aws.Time(1), ...
      'window_end', aws.Time(end));
end

%% Local functions
function filename = locateStationFile(site, source_dir, timescale)
   %LOCATESTATIONFILE Resolve the station NetCDF, matching site loosely.
   if source_dir == ""
      source_dir = string(fullfile(icemodel.getpath('data'), ...
         'forcing', 'promice'));
   end
   % Accept either the directory holding the files or its parent with
   % hourly/ and daily/ subdirectories.
   if isfolder(fullfile(source_dir, timescale))
      source_dir = fullfile(source_dir, timescale);
   end
   if ~isfolder(source_dir)
      error('icemodel:forcing:readPromiceAws:sourceNotFound', ...
         ['PROMICE source directory not found: %s. Pass source_dir or ' ...
         'stage the v3 NetCDF bundle (reference layout: ' ...
         '/Volumes/S03/DATA/greenland/geus/aws/v3).'], source_dir)
   end

   files = dir(fullfile(source_dir, '*_v03.nc'));
   wanted = lower(erase(site, "_"));
   names = string({files.name});
   stations = lower(erase(extractBefore(names, "_hour"), "_"));
   if timescale == "daily"
      stations = lower(erase(extractBefore(names, "_day"), "_"));
   end
   idx = find(stations == wanted, 1);
   if isempty(idx)
      error('icemodel:forcing:readPromiceAws:stationNotFound', ...
         'no %s station file matching "%s" in %s', timescale, site, ...
         source_dir)
   end
   filename = fullfile(files(idx).folder, files(idx).name);
end

function primary = preferPrimary(primary, secondary)
   %PREFERPRIMARY Fill missing primary samples from the secondary series.
   missing = ~isfinite(primary);
   primary(missing) = secondary(missing);
end

function t = ensureUtc(t)
   %ENSUREUTC Coerce a datetime (or text) to a UTC-zoned datetime.
   t = datetime(t);
   if isempty(t.TimeZone)
      t.TimeZone = 'UTC';
   else
      t = datetime(t, 'TimeZone', 'UTC');
   end
end
