function [times, record] = stationTransitionTimes(composing_stations, kwargs)
   %STATIONTRANSITIONTIMES Within-record AWS handover times for a PROMICE site.
   %
   %  times = icemodel.forcing.helpers.stationTransitionTimes(stations)
   %  [times, record] = ... stationTransitionTimes(stations, ...
   %     window_start=..., window_end=..., source_dir=...)
   %
   % A PROMICE "site" can merge several AWS ("stations") over time. The
   % per-station install dates live in AWS_stations_metadata.csv (downloaded from
   % the GEUS thredds server); the L3 NetCDF and AWS_sites_metadata.csv carry only
   % the composing-station NAMES and a single SITE-level install date, so the
   % per-station handover TIMES are recovered here from the stations CSV.
   %
   % A composing station's install date is a within-record HANDOVER only when it
   % falls strictly AFTER the site record start (the FOUNDING station's install
   % coincides with the record start - it begins the record, it is not a handover
   % within it) and at/before the record end. Stations whose name is absent from
   % the CSV (legacy GC-Net names like GITS / DYE-2 / SwissCamp that predate the
   % modern v3 ids) contribute no date and are reported in the record.
   %
   % Inputs
   %  composing_stations  - string array of the site's composing AWS names
   %                        (from promiceSiteCatalog / AWS_sites_metadata.csv)
   %
   % Name-value
   %  window_start, window_end : the site record window (datetime). When set, a
   %                        station install is a handover only when it falls
   %                        strictly after window_start (+tol_days) and at/before
   %                        window_end. Default: NaT (no record clamp; every
   %                        non-founding install date is returned).
   %  source_dir          : directory holding AWS_stations_metadata.csv (default
   %                        the staged data/verification/promice tree)
   %  tol_days            : how far after the record start an install must fall to
   %                        count as a handover, not the founding install (1)
   %
   % Outputs
   %  times  - datetime column of within-record handover times (UTC, may be empty)
   %  record - struct array, one entry per composing station:
   %           .station .install_date .in_csv (logical) .is_handover (logical)
   %
   % See also: icemodel.forcing.helpers.surfaceFlags,
   %  icemodel.forcing.destepSurface, icemodel.forcing.buildPromiceData,
   %  icemodel.verification.setup.promiceSiteCatalog

   arguments
      composing_stations (1, :) string
      kwargs.window_start (1, 1) datetime = NaT
      kwargs.window_end (1, 1) datetime = NaT
      kwargs.source_dir (1, 1) string = ""
      kwargs.tol_days (1, 1) double {mustBeNonnegative} = 1
   end

   times = NaT(0, 1, 'TimeZone', 'UTC');
   record = emptyRecord();
   record = record([]);

   csv = locateStationsCsv(kwargs.source_dir);
   if csv == "" || ~isfile(csv)
      return
   end
   T = readtable(csv, 'TextType', 'string');
   if ~all(ismember({'station_id', 'date_installation'}, ...
         T.Properties.VariableNames))
      return
   end

   % Record window as a UTC bound (the install dates are date-only).
   w0 = kwargs.window_start;
   w1 = kwargs.window_end;
   if ~isnat(w0) && isempty(w0.TimeZone); w0.TimeZone = 'UTC'; end
   if ~isnat(w1) && isempty(w1.TimeZone); w1.TimeZone = 'UTC'; end

   record = repmat(emptyRecord(), 1, numel(composing_stations));
   keep = false(1, numel(composing_stations));
   for k = 1:numel(composing_stations)
      stn = composing_stations(k);
      idx = find(T.station_id == stn, 1);
      entry = emptyRecord();
      entry.station = stn;
      if isempty(idx)
         % Legacy GC-Net name not in the modern stations CSV: no install date.
         record(k) = entry;
         continue
      end
      di = parseInstall(T.date_installation(idx));
      entry.install_date = di;
      entry.in_csv = true;
      if isnat(di)
         record(k) = entry;
         continue
      end
      % Handover: strictly after the record start (founding install begins the
      % record, not a handover within it) and at/before the record end.
      after_start = isnat(w0) || di > w0 + days(kwargs.tol_days);
      before_end = isnat(w1) || di <= w1;
      entry.is_handover = after_start && before_end;
      keep(k) = entry.is_handover;
      record(k) = entry;
   end

   if any(keep)
      handover_dates = [record(keep).install_date];
      times = sort(handover_dates(:));
   end
end

%% Local functions
function di = parseInstall(s)
   %PARSEINSTALL Parse a yyyy-MM-dd install date string to a UTC datetime.
   s = strtrim(string(s));
   if s == "" || ismissing(s)
      di = NaT('TimeZone', 'UTC');
      return
   end
   di = datetime(s, 'InputFormat', 'yyyy-MM-dd', 'TimeZone', 'UTC');
end

function entry = emptyRecord()
   %EMPTYRECORD One station-transition record entry (canonical field order).
   entry = struct('station', "", 'install_date', NaT('TimeZone', 'UTC'), ...
      'in_csv', false, 'is_handover', false);
end

function csv = locateStationsCsv(source_dir)
   %LOCATESTATIONSCSV Resolve AWS_stations_metadata.csv under the staged product.
   if source_dir == ""
      source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
         'verification', 'promice'));
   end
   % Accept a selected product root, or the parent only when the caller
   % explicitly selected that product's hour/day subdirectory.
   source_dir = strip(source_dir, 'right', filesep);
   candidates = fullfile(source_dir, 'AWS_stations_metadata.csv');
   [parent, leaf] = fileparts(char(source_dir));
   if ismember(lower(string(leaf)), ["hour", "day"])
      candidates(end + 1) = fullfile(parent, 'AWS_stations_metadata.csv');
   end
   csv = "";
   for c = string(candidates(:)')
      if isfile(c)
         csv = c;
         return
      end
   end
end
