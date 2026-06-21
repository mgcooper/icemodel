function [aws, metadata] = readPromiceAws(site, kwargs)
   %READPROMICEAWS Read a pypromice L3 AWS NetCDF into icemodel channels.
   %
   %  [aws, metadata] = icemodel.forcing.readPromiceAws(site)
   %  [aws, metadata] = ... readPromiceAws(site, source_dir=..., ...
   %     timescale="hourly", startdate=..., enddate=...)
   %
   % Reads one GEUS/PROMICE automatic-weather-station NetCDF from the
   % pypromice Level-3 product (the hour/<STATION>_hour.nc bundle) and
   % returns a timetable with icemodel-standard channel names. Channel
   % names, units, levels, and physical ranges follow the product variable
   % dictionary (data/verification/promice/AWS_variables.csv).
   %
   % Forcing channels (NetCDF -> output, with the unit change applied):
   %    t_u [degC] -> tair [K]         p_u [hPa] -> psfc [Pa]
   %    rh_u [%] -> rh [%]             wspd_u -> wspd [m s-1]
   %    wdir_u -> wdir [deg]           dsr -> swd [W m-2]
   %    usr -> swu [W m-2]             dlr -> lwd [W m-2]
   %    ulr -> lwu [W m-2]             albedo -> albedo [-]
   %    cc [%] -> cfrac [-]            t_surf [C] -> tsfc [K]
   %    dshf_u -> shf, dlhf_u -> lhf [W m-2]
   %    rainfall_cor_u [mm] -> rainf [mm per timestep]
   %
   % Evaluation channels (the QC'd L3 surface vars, read not derived):
   %    snow_height -> snow_height [m]   (snow surface rel. ice surface)
   %    z_ice_surf -> z_ice_surf [m]     (ice surface rel. installation)
   %    z_surf_combined -> z_surf_combined [m] (multi-sensor total surface)
   %    z_boom_cor_u (fallback z_boom_u) -> boom_height [m]
   %    z_pt_cor (fallback z_pt) -> transducer_depth [m]
   %    z_stake_cor (fallback z_stake) -> stake_height [m]
   %    t_i_1..t_i_N [degC] -> tice1..ticeN [K], clamped to the dictionary
   %       physical range (-80..1 C) per sensor
   %    alt -> elev [m] (smoothed postprocessed GPS altitude)
   %
   % The upper-boom channels (the dictionary "all" / "_u" set) are used so
   % one-boom and two-boom station configurations resolve identically; the
   % lower-boom (_l) channels are ignored.
   %
   % Inputs
   %  site - station name. Accepts the canonical id ("KAN_L") or the compact
   %         lowercase alias ("kanl"); matching ignores case and underscores
   %         against the station files in source_dir. The product carries the
   %         full GEUS AWS network (~150 stations).
   %
   % Name-value
   %  source_dir : directory holding the station NetCDF files, or its parent
   %      holding hour/ (and day/, month/) subdirectories. Defaults to the
   %      staged product under the repo-root data/verification/promice.
   %  timescale : "hourly" (default; hour/) or "daily" (day/).
   %  startdate, enddate : optional datetime window; default full range.
   %
   % Outputs
   %  aws      - timetable, UTC time axis, channels above (a channel missing
   %             from the file is omitted from the output)
   %  metadata - provenance struct: source file, station id, lat, lon,
   %             elevation, boom count, row count, time bounds
   %
   % See also: icemodel.forcing.buildPromiceMet,
   %  icemodel.forcing.buildPromiceData

   arguments
      site (1, 1) string
      kwargs.source_dir (1, 1) string = ""
      kwargs.timescale (1, 1) string ...
         {mustBeMember(kwargs.timescale, ["hourly", "daily"])} = "hourly"
      kwargs.startdate = ""
      kwargs.enddate = ""
   end

   filename = locateStationFile(site, kwargs.source_dir, kwargs.timescale);

   % Time axis: pypromice L3 files encode the bin-start stamp as hours (or
   % days) since a station-specific epoch carried in the time:units attribute.
   Time = readTimeAxis(filename, kwargs.timescale);

   info = ncinfo(filename);
   available = string({info.Variables.Name});

   % Forcing channel map: output name -> nc variable (and optional fallback).
   % Upper-boom (_u) channels are the dictionary "all" set, so one-boom and
   % two-boom stations resolve through the same names.
   channels = {
      'tair',      't_u',            ''
      'psfc',      'p_u',            ''
      'rh',        'rh_u',           ''
      'wspd',      'wspd_u',         ''
      'wdir',      'wdir_u',         ''
      'swd',       'dsr',            ''
      'swu',       'usr',            ''
      'lwd',       'dlr',            ''
      'lwu',       'ulr',            ''
      'albedo',    'albedo',         ''
      'cfrac',     'cc',             ''
      'tsfc',      't_surf',         ''
      'shf',       'dshf_u',         ''
      'lhf',       'dlhf_u',         ''
      'rainf',     'rainfall_cor_u', ''
      'snow_height',     'snow_height',     ''
      'z_ice_surf',      'z_ice_surf',      ''
      'z_surf_combined', 'z_surf_combined', ''
      'boom_height',     'z_boom_cor_u',    'z_boom_u'
      'transducer_depth','z_pt_cor',        'z_pt'
      'stake_height',    'z_stake_cor',     'z_stake'
      'elev',            'alt',             'gps_alt'
      };

   aws = timetable(Time);
   for n = 1:size(channels, 1)
      outname = channels{n, 1};
      ncname = channels{n, 2};
      fallback = string(channels{n, 3});
      if ~ismember(ncname, available)
         if fallback ~= "" && ismember(fallback, available)
            ncname = char(fallback);
         else
            continue
         end
      end
      data = double(ncread(filename, ncname));
      if fallback ~= "" && ismember(fallback, available) ...
            && ~strcmp(ncname, fallback)
         data = preferPrimary(data, double(ncread(filename, fallback)));
      end
      aws.(outname) = data;
   end

   % Unit conversions to icemodel standard.
   Tf = icemodel.physicalConstant('Tf');
   if ismember('tair', string(aws.Properties.VariableNames))
      aws.tair = aws.tair + Tf;                   % degC -> K
   end
   if ismember('psfc', string(aws.Properties.VariableNames))
      aws.psfc = aws.psfc * 100;                  % hPa -> Pa
   end
   if ismember('tsfc', string(aws.Properties.VariableNames))
      aws.tsfc = aws.tsfc + Tf;                   % C -> K
   end
   if ismember('cfrac', string(aws.Properties.VariableNames))
      aws.cfrac = aws.cfrac / 100;                % percent -> fraction
   end

   % Ice-temperature string: t_i_1..t_i_N [degC] -> tice1..ticeN [K], clamped
   % to the dictionary physical range (-80..1 C) before the K offset, so out-
   % of-range thermistor spikes never reach the evaluation axis.
   icerange = [-80, 1];   % [degC], from AWS_variables.csv (t_i_*)
   nice = 0;
   while ismember(sprintf('t_i_%d', nice + 1), available)
      nice = nice + 1;
      v = double(ncread(filename, sprintf('t_i_%d', nice)));
      v(v < icerange(1) | v > icerange(2)) = NaN;
      aws.(sprintf('tice%d', nice)) = v + Tf;
   end

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
      'station_name', readGlobalString(filename, 'site_id', site), ...
      'lat', readGlobalNumber(filename, 'latitude'), ...
      'lon', readGlobalNumber(filename, 'longitude'), ...
      'elev', readGlobalNumber(filename, 'altitude'), ...
      'n_booms', detectBooms(available), ...
      'n_tice', nice, ...
      'n_rows', height(aws), ...
      'window_start', aws.Time(1), ...
      'window_end', aws.Time(end));
end

%% Local functions
function filename = locateStationFile(site, source_dir, timescale)
   %LOCATESTATIONFILE Resolve the station NetCDF, matching site loosely.
   if source_dir == ""
      source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
         'verification', 'promice'));
   end
   subdir = "hour";
   if timescale == "daily"
      subdir = "day";
   end
   % Accept either the directory holding the files or its parent with the
   % hour/ (day/) subdirectories.
   if isfolder(fullfile(source_dir, subdir))
      source_dir = fullfile(source_dir, subdir);
   end
   if ~isfolder(source_dir)
      error('icemodel:forcing:readPromiceAws:sourceNotFound', ...
         ['PROMICE source directory not found: %s. Pass source_dir or ' ...
         'stage the pypromice L3 product (data/verification/promice).'], ...
         source_dir)
   end

   files = dir(fullfile(source_dir, '*.nc'));
   wanted = lower(erase(site, "_"));
   names = string({files.name});
   suffix = "_" + subdir + ".nc";
   stations = lower(erase(erase(names, suffix), "_"));
   idx = find(stations == wanted, 1);
   if isempty(idx)
      error('icemodel:forcing:readPromiceAws:stationNotFound', ...
         'no %s station file matching "%s" in %s', timescale, site, ...
         source_dir)
   end
   filename = fullfile(files(idx).folder, files(idx).name);
end

function Time = readTimeAxis(filename, timescale)
   %READTIMEAXIS Decode the pypromice L3 time axis onto a UTC bin-start axis.
   t = double(ncread(filename, 'time'));
   units = ncreadatt(filename, 'time', 'units');
   tok = regexp(units, ...
      '(\w+) since (\d{4}-\d{2}-\d{2})[ T](\d{2}:\d{2}:\d{2})', ...
      'tokens', 'once');
   if isempty(tok)
      error('icemodel:forcing:readPromiceAws:timeUnits', ...
         'unexpected time units in %s: %s', filename, units)
   end
   epoch = datetime([tok{2} ' ' tok{3}], ...
      'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', 'UTC');
   switch lower(tok{1})
      case {'hours', 'hour'}
         Time = epoch + hours(t);
      case {'days', 'day'}
         Time = epoch + days(t);
      case {'seconds', 'second'}
         Time = epoch + seconds(t);
      otherwise
         error('icemodel:forcing:readPromiceAws:timeUnits', ...
            'unsupported time unit "%s" in %s', tok{1}, filename)
   end
   if timescale == "hourly"
      Time = dateshift(Time, 'start', 'hour');
   else
      Time = dateshift(Time, 'start', 'day');
   end
end

function n = detectBooms(available)
   %DETECTBOOMS Station boom count from the presence of lower-boom channels.
   n = 1;
   if any(ismember(["t_l", "z_boom_l", "wspd_l"], available))
      n = 2;
   end
end

function primary = preferPrimary(primary, secondary)
   %PREFERPRIMARY Fill missing primary samples from the secondary series.
   missing = ~isfinite(primary);
   primary(missing) = secondary(missing);
end

function value = readGlobalNumber(filename, attname)
   %READGLOBALNUMBER Read a numeric global attribute stored as a string.
   raw = ncreadatt(filename, '/', attname);
   value = str2double(string(raw));
end

function value = readGlobalString(filename, attname, fallback)
   %READGLOBALSTRING Read a string global attribute, with a fallback.
   try
      value = string(ncreadatt(filename, '/', attname));
   catch
      value = string(fallback);
   end
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
