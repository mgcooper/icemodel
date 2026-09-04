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
   % The mapping covers more than the minimal met subset. Every L3 channel
   % useful for snow, firn, and ice model forcing or evaluation maps to one
   % icemodel name. The mapping omits the housekeeping and
   % diagnostic channels: battery voltage, fan current, raw per-timestep GPS
   % latitude and longitude, and radiation-sensor temperature.
   %
   % Forcing channels (NetCDF -> output, with the unit change applied):
   %    t_u [degC] -> tair [K]         p_u [hPa] -> psfc [Pa]
   %    rh_u [%] -> rh [%]             wspd_u -> wspd [m s-1]
   %    wdir_u -> wdir [deg]           dsr -> swd [W m-2]
   %    usr -> swu [W m-2]             dlr -> lwd [W m-2]
   %    ulr -> lwu [W m-2]             albedo -> albedo [-]
   %    cc [0..1 or %] -> cfrac [-]    t_surf [C] -> tsfc [K]
   %    dshf_u -> shf, dlhf_u -> lhf [W m-2]
   %    rainfall_cor_u [mm] -> rainf [mm per timestep]
   %    qh_u [g/kg] -> shum [kg/kg]    (specific humidity, /1000; the L3
   %       label "kg/kg" is wrong, the values are g/kg)
   %    rh_u_wrt_ice_or_water -> rh_ice [%] (RH wrt ice/water; correct
   %       subfreezing humidity for sublimation/condensation)
   %    wspd_x_u -> uwind, wspd_y_u -> vwind [m s-1] (wind components)
   %    dsr_cor -> swd_cor, usr_cor -> swu_cor [W m-2] (tilt/bias-corrected
   %       shortwave; preferred over raw dsr/usr for SEB evaluation)
   %    tilt_x -> tilt_x, tilt_y -> tilt_y, rot -> rot [deg] (platform tilt
   %       and azimuth; the radiation-correction geometry)
   %
   % Evaluation channels (the QC'd L3 surface vars, read not derived):
   %    snow_height -> snow_height [m]   (snow surface rel. ice surface)
   %    z_ice_surf -> z_ice_surf [m]     (ice surface rel. installation)
   %    z_surf_combined -> z_surf_combined [m] (multi-sensor total surface)
   %    z_boom_cor_u (fallback z_boom_u) -> boom_height [m]
   %    z_pt_cor (fallback z_pt) -> transducer_depth [m]
   %    z_stake_cor (fallback z_stake) -> stake_height [m]
   %    t_i_1..t_i_N [degC] -> tice1..ticeN [K], clamped to the dictionary
   %       physical range (-80..1 C) per sensor, with SURFACED thermistors
   %       (depth <= 0) discarded so only subsurface readings survive
   %    d_t_i_1..d_t_i_N [m] -> dtice1..dticeN (thermistor depths below the
   %       surface, positive down, so the tice string can be interpreted
   %       vertically and the surfaced-sensor discard applied)
   %    t_i_10m [degC] -> tice10m_source [K] (unmodified GEUS value) and
   %       tice10m [K] (PRIMARY subsurface-temperature evaluation channel,
   %       masked at physically discontinuous consecutive-hour endpoints).
   %       tice10m_qc_flag records 0=accepted, 1=failed by the source-range
   %       screen or by a discontinuity with native-neighbor support,
   %       2=unreviewed because fewer than two depth-tagged native thermistors
   %       span the endpoint pair, and 3=an unresolved isolated-sensor epoch.
   %       The raw depth-tagged tice string remains diagnostic.
   %    alt -> elev [m] (smoothed postprocessed GPS altitude)
   %
   % The upper-boom channels (the dictionary "all" / "_u" set) are used so
   % one-boom and two-boom station configurations resolve identically; the
   % lower-boom (_l) channels are ignored.
   %
   % Inputs
   %  site - station name. Accepts the official id ("KAN_L") or the compact
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
   %             elevation, boom count, row count, time bounds, and whole-file
   %             raw/corrected shortwave availability
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

   % Validate and normalize the optional interval before source discovery so
   % malformed public input cannot be obscured by an unrelated file error.
   [window_start, window_end, has_window] = ...
      icemodel.pairedWindow(kwargs.startdate, kwargs.enddate);

   filename = locateStationFile(site, kwargs.source_dir, kwargs.timescale);

   % TIME CONVENTION. This rule governs all model-versus-observation alignment
   % downstream. The pypromice L3 hourly timestamp is the START of the averaged
   % hour (product readme, "Temporal averaging": "the timestamp of the hourly
   % averages indicate the start of the averaged hour"). The file encodes the
   % time as integer hours (0,1,2,...) since a station epoch given in the
   % time:units attribute, so epoch+hours(t) reproduces the bin-START stamp
   % exactly. The dateshift('start','hour') call below snaps sub-hour epoch
   % jitter to the hour. It repeats without effect and does not re-bin the
   % data. The icemodel met and Data axis uses this same bin-START hourly grid
   % in UTC. icemodel.loadmet and the timestepping loop read a met row's Time
   % as the forcing valid AT that timestamp, then integrate forward over
   % [t, t+dt). A START-of-hour averaged forcing is the correct mean for the
   % [t, t+1h) step, so this code applies no half-hour shift. When you compare
   % a simulation with these observations, align on this START-of-hour stamp
   % and do not recentre to the middle of the hour.
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
      'rh_ice',    'rh_u_wrt_ice_or_water', ''
      'shum',      'qh_u',           ''
      'wspd',      'wspd_u',         ''
      'wdir',      'wdir_u',         ''
      'uwind',     'wspd_x_u',       ''
      'vwind',     'wspd_y_u',       ''
      'swd',       'dsr',            ''
      'swu',       'usr',            ''
      'swd_cor',   'dsr_cor',        ''
      'swu_cor',   'usr_cor',        ''
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
      'tilt_x',          'tilt_x',          ''
      'tilt_y',          'tilt_y',          ''
      'rot',             'rot',             ''
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
      % Current pypromice L3 files store cc as a 0..1 fraction despite the
      % percent label, while older products may contain 0..100 percentages.
      % Mask impossible values before inferring the whole-series scale so a
      % fill sentinel cannot make a fractional product divide a second time.
      invalid = ~isnan(aws.cfrac) & ...
         (~isfinite(aws.cfrac) | aws.cfrac < 0 | aws.cfrac > 100);
      aws.cfrac(invalid) = NaN;
      finite = isfinite(aws.cfrac);
      if any(aws.cfrac(finite) > 1)
         aws.cfrac(finite) = aws.cfrac(finite) / 100;
      end
   end
   if ismember('shum', string(aws.Properties.VariableNames))
      % qh_u is specific humidity reported in g/kg: the L3 magnitudes are
      % O(0.1..6) (e.g. KAN_M median 1.68 g/kg), implausible for kg/kg. The
      % product dictionary (AWS_variables.csv) and the NetCDF units attribute
      % both LABEL it "kg/kg" but the values are g/kg (the dictionary lo/hi
      % 0..100 is the g/kg range, confirming the label is wrong). Convert to
      % kg/kg so shum matches the MAR/MERRA convention consumed by the vapor
      % kernel (icemodel.vapor.relative_humidity_from_specific_humidity).
      aws.shum = aws.shum / 1000;                 % g/kg -> kg/kg
   end

   % Ice-temperature string: t_i_1..t_i_N [degC] -> tice1..ticeN [K]. Clamp to
   % the dictionary physical range (-80..1 C) before the K offset, so an
   % out-of-range thermistor spike never reaches the evaluation axis. The raw
   % string carries a depth tag d_t_i_N [m], the time-dependent depth BELOW
   % the surface, positive downward. A SURFACED thermistor has depth <= 0,
   % which means the sensor sits at or above the surface as the firn melts out
   % at an ablation site. The product readme says to discard it, so its tice
   % sample becomes NaN and only subsurface temperatures reach the evaluation
   % axis. GEUS already discards most surfaced readings upstream, so this
   % depth gate is a second check. It runs here so the icemodel channel is
   % depth-clean for every source version. When a depth tag is NaN, the code
   % keeps the temperature, because the sample is still a subsurface reading
   % with an unknown depth, and dticeN carries NaN for that sample.
   icerange = [-80, 1];   % [degC], from AWS_variables.csv (t_i_*)
   nice = 0;
   while ismember(sprintf('t_i_%d', nice + 1), available)
      nice = nice + 1;
      v = double(ncread(filename, sprintf('t_i_%d', nice)));
      v(v < icerange(1) | v > icerange(2)) = NaN;
      % Thermistor depth d_t_i_N [m] -> dticeN, so the tice string can be
      % interpreted vertically (sensors settle as the surface changes).
      if ismember(sprintf('d_t_i_%d', nice), available)
         depth = double(ncread(filename, sprintf('d_t_i_%d', nice)));
         v(depth <= 0) = NaN;   % discard surfaced thermistors
         aws.(sprintf('dtice%d', nice)) = depth;
      end
      aws.(sprintf('tice%d', nice)) = v + Tf;
   end

   % Keep GEUS's derived t_i_10m unchanged apart from degC -> K, and add one
   % conservative tice10m target. The source processor documents that
   % noisy-thermistor filtering is off. A 45-site audit of 2.30 million
   % consecutive hourly pairs found a 99.9th-percentile change of 0.350 C,
   % while known KAN_U sensor failures jump 5.23 and 6.90 C in one hour. A
   % change above 1 C between consecutive hours is therefore an impossible
   % 10 m thermal response, not seasonal evolution. This code masks both
   % endpoints so no plot or scorer draws a line across the discontinuity. The
   % source channel and a review-status flag stay in the output.
   if ismember('t_i_10m', available)
      aws.tice10m_source = double(ncread(filename, 't_i_10m')) + Tf;
      [aws.tice10m, aws.tice10m_qc_flag] = ...
         qualityControlTice10m(aws, kwargs.timescale);
   end

   % Record radiation support before applying the requested window. Builders
   % need the whole-file status, so a build of one outage window makes the
   % same geometry-derived darkness decision as a wider build that contains
   % finite source samples. This also separates an absent channel, or a
   % channel that is missing across the whole file, from a window of an
   % observed record that happens to be all missing.
   full_names = string(aws.Properties.VariableNames);
   swd_source_file_present = ismember("swd", full_names);
   swd_corrected_source_file_present = ismember("swd_cor", full_names);
   swd_source_file_observations_present = false;
   if swd_source_file_present
      swd_source_file_observations_present = any(isfinite(aws.swd), 'all');
   end
   if swd_corrected_source_file_present
      swd_source_file_observations_present = ...
         swd_source_file_observations_present ...
         || any(isfinite(aws.swd_cor), 'all');
   end
   swu_source_file_present = ismember("swu", full_names);
   swu_corrected_source_file_present = ismember("swu_cor", full_names);
   swu_source_file_observations_present = false;
   if swu_source_file_present
      swu_source_file_observations_present = any(isfinite(aws.swu), 'all');
   end
   if swu_corrected_source_file_present
      swu_source_file_observations_present = ...
         swu_source_file_observations_present ...
         || any(isfinite(aws.swu_cor), 'all');
   end

   % Optional window subset.
   keep = true(height(aws), 1);
   if has_window
      keep = aws.Time >= window_start & aws.Time <= window_end;
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
      'swd_source_file_present', swd_source_file_present, ...
      'swd_corrected_source_file_present', ...
      swd_corrected_source_file_present, ...
      'swd_source_file_observations_present', ...
      swd_source_file_observations_present, ...
      'swu_source_file_present', swu_source_file_present, ...
      'swu_corrected_source_file_present', ...
      swu_corrected_source_file_present, ...
      'swu_source_file_observations_present', ...
      swu_source_file_observations_present, ...
      'n_rows', height(aws), ...
      'window_start', aws.Time(1), ...
      'window_end', aws.Time(end));
   if ismember('tice10m_qc_flag', string(aws.Properties.VariableNames))
      metadata.tice10m_qc_status = "applied";
      metadata.tice10m_qc_method = ...
         "mask_gt_1K_hourly_endpoints_and_large_isolated_sensor_epochs";
      metadata.tice10m_qc_source_variable = "t_i_10m";
      metadata.tice10m_qc_source_channel = "tice10m_source";
      metadata.tice10m_qc_jump_threshold_K = 1;
      metadata.tice10m_qc_persistent_jump_threshold_K = 4;
      metadata.tice10m_qc_other_sensor_median_threshold_K = 0.25;
      metadata.tice10m_qc_target_depth_tolerance_m = 2;
      metadata.tice10m_qc_recovery_window_hours = 24;
      metadata.tice10m_qc_depth_reset_threshold_m = 0.5;
      metadata.tice10m_qc_flag_codes = struct( ...
         'accepted', 0, 'failed', 1, 'unreviewed', 2, ...
         'persistent_unreviewed', 3);
      metadata.tice10m_qc_flagged_sample_count = ...
         nnz(aws.tice10m_qc_flag > 0);
      metadata.tice10m_qc_failed_sample_count = ...
         nnz(aws.tice10m_qc_flag == 1);
      metadata.tice10m_qc_unreviewed_sample_count = ...
         nnz(aws.tice10m_qc_flag >= 2);
      metadata.tice10m_qc_persistent_sample_count = ...
         nnz(aws.tice10m_qc_flag == 3);
      metadata.tice10m_qc_basis = ...
         "GEUS t_i_10m preserved in tice10m_source; canonical tice10m " ...
         + "masks source values outside -80..1 degC and both endpoints of " ...
         + ">1 K changes across exactly one hour. " ...
         + "Depth-tagged native thermistors classify events as reviewed " ...
         + "or neighbor-insufficient. An isolated >4 K jump in a sensor " ...
         + "within 2 m of the 10 m target, with no 24-hour recovery, is " ...
         + "masked as unreviewed until that sensor's next >0.5 m depth " ...
         + "reset; gaps are not treated as jumps.";
   end
end

%% Local functions
function [clean, flag] = qualityControlTice10m(aws, timescale)
   %QUALITYCONTROLTICE10M Mask impossible hourly 10 m temperature jumps.
   %
   % This screen runs only across contiguous hourly source samples. The native
   % thermistor temperatures and their time-varying depths decide whether an
   % event has at least two comparable subsurface sensors. An event with fewer
   % sensors is still masked, and its flag is unreviewed (code 2).

   source = aws.tice10m_source;
   clean = source;
   flag = zeros(size(source));

   % A finite value outside the physical ice-temperature envelope is invalid
   % even when there is no neighboring sample. Apply this cadence-independent
   % screen before the hourly discontinuity rule and preserve only the source.
   out_of_range = isfinite(source) ...
      & (source < icemodel.physicalConstant('Tf') - 80 ...
      | source > icemodel.physicalConstant('Tf') + 1);
   flag(out_of_range) = 1;
   clean(out_of_range) = NaN;

   % The cross-site jump threshold was validated for the hourly product. Daily
   % files skip that rule; buildPromiceData always reads the hourly source.
   if timescale ~= "hourly" || numel(source) < 2
      return
   end

   % Only exact adjacent-hour pairs can prove a temporal discontinuity. A large
   % change across an omitted interval is missing evidence, not a sensor jump.
   contiguous = abs(seconds(diff(aws.Time)) - 3600) <= 1;
   jump = abs(diff(source));
   candidates = find(contiguous & isfinite(source(1:end-1)) ...
      & isfinite(source(2:end)) & jump > 1);
   names = string(aws.Properties.VariableNames);
   thermistors = names(~cellfun('isempty', ...
      regexp(cellstr(names), '^tice\d+$', 'once')));

   closed_through = 0;
   for first = reshape(candidates, 1, [])
      % Count thermistors whose temperature and positive depth are known at
      % both endpoints. Two provide the minimum native-profile context needed
      % to call the target discontinuity reviewed rather than merely suspect.
      comparable = 0;
      % Each ledger holds at most one entry per thermistor, so size them at
      % the thermistor count and trim to the filled portion after the loop.
      n_thermistors = numel(thermistors);
      comparable_names = strings(n_thermistors, 1);
      sensor_jumps = zeros(n_thermistors, 1);
      depth_jumps = zeros(n_thermistors, 1);
      depths_before = zeros(n_thermistors, 1);
      depths_after = zeros(n_thermistors, 1);
      before_support = strings(n_thermistors, 1);
      after_support = strings(n_thermistors, 1);
      n_before = 0;
      n_after = 0;
      for name = reshape(thermistors, 1, [])
         depth_name = "d" + name;
         if ~ismember(depth_name, names)
            continue
         end
         temperature = aws.(name);
         depth = aws.(depth_name);
         pair = first:first + 1;
         before_valid = isfinite(temperature(first)) ...
            && isfinite(depth(first)) && depth(first) > 0;
         after_valid = isfinite(temperature(first + 1)) ...
            && isfinite(depth(first + 1)) && depth(first + 1) > 0;
         if before_valid
            n_before = n_before + 1;
            before_support(n_before) = name;
         end
         if after_valid
            n_after = n_after + 1;
            after_support(n_after) = name;
         end
         if before_valid && after_valid
            comparable = comparable + 1;
            comparable_names(comparable) = name;
            sensor_jumps(comparable) = abs(diff(temperature(pair)));
            depth_jumps(comparable) = abs(diff(depth(pair)));
            depths_before(comparable) = depth(first);
            depths_after(comparable) = depth(first + 1);
         end
      end
      before_support = before_support(1:n_before);
      after_support = after_support(1:n_after);
      comparable_names = comparable_names(1:comparable);
      sensor_jumps = sensor_jumps(1:comparable);
      depth_jumps = depth_jumps(1:comparable);
      depths_before = depths_before(1:comparable);
      depths_after = depths_after(1:comparable);
      code = 1;
      if comparable < 2
         code = 2;
      end
      flag(first:first + 1) = max(flag(first:first + 1), code);

      % A large isolated sensor jump can shift the derived target level for a
      % long time, past its first missing endpoint. Extend the mask only when
      % every one of these conditions holds: at least three comparable native
      % sensors, exactly one jump above 1 K from a sensor within 2 m of the
      % 10 m target, no change in the other sensors, their support, or their
      % depths, and no return to the pre-jump level within 24 hours. The next
      % depth reset of that sensor above 0.5 m starts a new epoch.
      offender_index = NaN;
      offender_near_target = false;
      if comparable > 0
         [~, offender_index] = max(sensor_jumps);
         offender_depths = [depths_before(offender_index), ...
            depths_after(offender_index)];
         offender_near_target = any(abs(offender_depths - 10) <= 2);
      end
      isolated = jump(first) > 4 && comparable >= 3 ...
         && nnz(sensor_jumps > 1) == 1 ...
         && median(sensor_jumps) <= 0.25 ...
         && isequal(before_support, after_support) ...
         && ~any(depth_jumps > 0.5) && offender_near_target;
      event_time = aws.Time(first + 1);
      recovery_window = find(aws.Time > event_time ...
         & aws.Time <= event_time + hours(24));
      recovery_offset = find(isfinite(source(recovery_window)) ...
         & abs(source(recovery_window) - source(first)) <= 1, 1);
      recovered = ~isempty(recovery_offset);
      if recovered
         recovery = recovery_window(recovery_offset);
         closed_through = max(closed_through, recovery);
      end
      % A recovery jump belongs to the transient segment that it closes; do not
      % reinterpret the same opposite edge as the start of a persistent epoch.
      if isolated && ~recovered && first >= closed_through
         offender = comparable_names(offender_index);
         offender_depth = aws.("d" + offender);
         depth_step = [false; abs(diff(offender_depth)) > 0.5 ...
            & isfinite(offender_depth(1:end-1)) ...
            & isfinite(offender_depth(2:end)) ...
            & abs(seconds(diff(aws.Time)) - 3600) <= 1];
         reset = find((1:height(aws))' > first + 1 & depth_step, 1);
         if isempty(reset)
            reset = height(aws);
         end
         flag(first:reset) = 3;
      end
   end

   % Apply every discontinuity decision to the tice10m target.
   clean(flag > 0) = NaN;
end

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
