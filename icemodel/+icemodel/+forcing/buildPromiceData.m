function [Data, metadata] = buildPromiceData(site, kwargs)
   %BUILDPROMICEDATA Build PROMICE data for model evaluation.
   %
   %  [Data, metadata] = icemodel.forcing.buildPromiceData(site)
   %  [Data, metadata] = icemodel.forcing.buildPromiceData(site, ...
   %     startdate=..., enddate=..., frequency="daily")
   %
   % Reads one pypromice L3 hourly station file. The output timetable contains
   % surface energy fluxes, surface-height observations, ice temperatures,
   % station coordinates, and variable units.
   %
   % Site types
   % A z_ice_surf channel identifies an ablation site. The function records the
   % site type in metadata.site_surface_type.
   %
   %  Ablation sites (z_ice_surf present):
   %    ablation   = -(z_ice_surf - z_ice_surf(window start)), positive down
   %                 (cumulative ice-surface lowering relative to installation).
   %    snow_depth = max(snow_height, 0). metadata.snow_depth_negatives counts
   %                 negative source samples.
   %
   %  Accumulation, percolation, and bedrock sites:
   %    surface_height = z_surf_combined - z_surf_combined(window start), a
   %                     net surface-height change, positive upward.
   %
   % Surface flags
   % Per-sample flags describe the ablation or surface_height channel.
   %  - surface_height_flag (0/1): marks a slope-interpolated surface value
   %    when all available surface-ranging sensors are NaN. Exclude flagged
   %    samples from per-timestep rate diagnostics. Use the full series for
   %    cumulative and visual comparisons.
   %  - station_transition_flag (0/1): marks station handover windows from
   %    AWS_stations_metadata.csv. A detected step inside such a window adds a
   %    'station_transition' evidence line in
   %    icemodel.forcing.destepSurface. metadata.is_multistation records the
   %    merge fact. metadata.station_transition_times and
   %    metadata.station_transition_record contain the dates.
   %  - step_detected_flag / step_correctable_flag / step_magnitude: the staged
   %    de-stepping detection from icemodel.forcing.destepSurface. Analysis can
   %    apply the correction to unambiguous steps.
   %
   % Ice temperatures include tice1..ticeN from the L3 t_i_* string and the
   % standardized tice10m evaluation channel. dtice1..dticeN contain the
   % changing sensor depths. tice10m_source contains the GEUS-derived value.
   % tice10m_qc_flag marks discontinuities, and tice10m masks flagged endpoints.
   %
   % tice10m comparison
   % tice10m is the standardized GEUS temperature 10 m below the changing
   % surface. GEUS tracks the time-dependent depth of each thermistor below the
   % surface (d_t_i_*), discards surfaced thermistors, and depth-interpolates
   % the remaining string to 10 m below the current surface at each timestep.
   % Compare tice10m with the model temperature 10 m below the model's current
   % surface. Use the raw tice1..N channels with their d_t_i_N sensor depths for
   % secondary diagnostics.
   %
   % Radiation
   % SWD and SWU use the corrected dsr_cor and usr_cor values when available,
   % then use the raw measurements. Finite negative values are set to zero.
   % Metadata records the number of values selected from each source channel.
   %
   % Shortwave gaps that lie entirely below deep civil night are set to zero.
   % Other missing observations remain missing. metchecks applies the physical
   % limits, and the surface-height channel carries the gap flag above.
   %
   % Inputs
   %  site - station id ("KAN_M" or compact alias "kanm").
   %
   % Name-value
   %  source_dir - PROMICE NetCDF directory. See readPromiceAws.
   %  startdate  - first requested sample. The default uses the record start.
   %  enddate    - last requested sample. The default uses the record end.
   %  frequency  - "hourly" (default) or "daily" (daily means).
   %
   % Outputs
   %  Data     - timetable with CustomProperties (X, Y, Lat, Lon, Elev,
   %             Slope, ScalarUnits). Common channels: tair, tsfc [K]; swd, swu,
   %             lwd, lwu, swn, lwn, netr, shf, lhf, thf [W m-2]; albedo, cfrac
   %             [-]; rh [%]; wspd [m s-1]; wdir [deg]; psfc [Pa]; tice1..ticeN,
   %             tice10m, tice10m_source [K], dtice1..dticeN [m], and
   %             tice10m_qc_flag [1]. Surface-height channels are site-type
   %             dependent: ablation sites carry ablation + snow_depth [m];
   %             accumulation sites carry surface_height [m]. Both carry the
   %             surface flag channels: surface_height_flag (0/1 gap-bridged
   %             mask), station_transition_flag (0/1 handover window),
   %             step_detected_flag / step_correctable_flag (0/1) and
   %             step_magnitude [m] (the staged de-stepping detection).
   %  metadata - provenance: source file, station, site_surface_type, surface
   %             channel + source, gap-flag counts, station-transition counts,
   %             composing_stations / is_multistation, steps_detected /
   %             steps_correctable / step_record, snow-depth negative count,
   %             QA/QC gap counts.
   %
   % See also: icemodel.forcing.readPromiceAws,
   %  icemodel.forcing.buildPromiceMet, icemodel.forcing.helpers.writeuserdata

   arguments
      site (1, 1) string
      kwargs.source_dir (1, 1) string = ""
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.frequency (1, 1) string ...
         {mustBeMember(kwargs.frequency, ["hourly", "daily"])} = "hourly"
   end

   [aws, source_meta] = icemodel.forcing.readPromiceAws(site, ...
      source_dir=kwargs.source_dir, timescale="hourly", ...
      startdate=kwargs.startdate, enddate=kwargs.enddate);

   % Use one corrected-first shortwave selection for observational Data and the
   % met builder. Missing radiation becomes zero only when the complete source
   % hour lies below deep civil twilight. Daylight and twilight outages stay
   % missing. The builder adds only channels that the source product carries, so
   % an absent observational channel stays absent.
   [public_swd, public_swu, shortwave_meta] = ...
      icemodel.forcing.helpers.promiceShortwave(aws, fill_darkness=true, ...
      latitude=source_meta.lat, longitude=source_meta.lon, ...
      swd_source_file_observations_present= ...
      source_meta.swd_source_file_observations_present, ...
      swu_source_file_observations_present= ...
      source_meta.swu_source_file_observations_present);
   if shortwave_meta.swd_source_present ...
         || shortwave_meta.swd_corrected_source_present
      aws.swd = public_swd;
   end
   if shortwave_meta.swu_source_present ...
         || shortwave_meta.swu_corrected_source_present
      aws.swu = public_swu;
   end

   % Record whether the source supplied finite albedo observations. Unlike the
   % forcing builder, observational Data never fills an absent/all-missing
   % channel or a temporal gap.
   has_albedo_source = ismember("albedo", ...
      string(aws.Properties.VariableNames));
   has_albedo_observations = has_albedo_source ...
      && any(isfinite(aws.albedo), 'all');
   if has_albedo_observations
      albedo_policy = ...
         "albedo retains PROMICE L3 observations; missing samples remain " + ...
         "missing; physical-range clamp only";
   elseif has_albedo_source
      albedo_policy = ...
         "albedo = NaN placeholder (PROMICE L3 albedo is all missing); no observations invented";
   else
      albedo_policy = ...
         "albedo = NaN placeholder (PROMICE L3 albedo source channel " + ...
         "absent); no observations invented";
   end

   % Derived net fluxes. Sparse stations may ship no radiation or turbulent-
   % flux channels, so each derived term is computed only when its inputs are
   % present (a missing input drops the derived channel rather than erroring).
   % has() is re-bound after channels are added so terms that build on a derived
   % channel (netr on swn/lwn) see it.
   has = @(v) ismember(v, string(aws.Properties.VariableNames));
   swn_negative_invalid_count = 0;
   if has("swd") && has("swu")
      aws.swn = aws.swd - aws.swu;

      % Net shortwave cannot be negative at the surface. Keep the selected
      % component observations, but set their derived net flux to missing, so
      % evaluation totals do not use physically inconsistent radiometry.
      invalid_swn = isfinite(aws.swn) & aws.swn < 0;
      swn_negative_invalid_count = nnz(invalid_swn);
      aws.swn(invalid_swn) = NaN;
   end
   if has("lwd") && has("lwu")
      aws.lwn = aws.lwd - aws.lwu;
   end
   has = @(v) ismember(v, string(aws.Properties.VariableNames));
   if has("swn") && has("lwn")
      aws.netr = aws.swn + aws.lwn;
   end
   if has("shf") && has("lhf")
      aws.thf = aws.shf + aws.lhf;
   end

   % Surface-height channels, branched on site type (see header). The
   % discriminator is whether the L3 file ships z_ice_surf (ablation sites).
   awsnames0 = string(aws.Properties.VariableNames);
   is_ablation = ismember("z_ice_surf", awsnames0);

   % The underlying L3 surface-ranging sensors (transducer/boom/stake) give the
   % gap flag: a sample is gap-bridged (slope-interpolated, not measured) when
   % every one of these sensors is NaN and the surface series is finite. The
   % sensors present at this station go into one matrix for surfaceFlags. GEUS
   % recomputes z_ice_surf from z_pt_cor, the transducer, and falls back to the
   % stake sonic ranger. The boom adds the accumulation-site surface ranging.
   sensor_names = intersect(["transducer_depth", "boom_height", ...
      "stake_height"], awsnames0, 'stable');
   if isempty(sensor_names)
      sensors = [];
   else
      sensors = cell2mat(arrayfun(@(v) aws.(v), sensor_names, ...
         'UniformOutput', false));
   end

   % Known station-handover times for the station-transition flag. The
   % composing-station NAMES come from the curated catalog. Their per-station
   % install DATES come from AWS_stations_metadata.csv, the GEUS thredds product
   % staged beside the L3 NetCDFs. stationTransitionTimes maps each composing
   % station to its install date and keeps only WITHIN-RECORD handovers, that
   % is, an install strictly after the record start. The install of the founding
   % station begins the record, so it is not a handover within the record. With
   % the CSV present, the flag carries those dates. Without the CSV,
   % transition_times stays empty and the flag is all zero. The merge FACT still
   % reaches metadata, so destepSurface can recover a transition as a step at
   % the same time. The window clamp uses the L3 record bounds (source_meta), so
   % an install outside this station's record is excluded.
   info = icemodel.verification.setup.promiceSiteCatalog(site, ...
      source_dir=kwargs.source_dir);
   composing_stations = info.stations;
   [transition_times, transition_record] = ...
      icemodel.forcing.helpers.stationTransitionTimes(composing_stations, ...
      window_start=source_meta.window_start, ...
      window_end=source_meta.window_end, source_dir=kwargs.source_dir);

   surface_meta = struct();
   if is_ablation
      surface_meta.site_surface_type = "ablation";

      % Cumulative ice-surface ablation from the QC'd L3 ice-surface height.
      % z_ice_surf is the ice surface relative to installation (decreases as the
      % surface lowers); ablation = -(z - z(start)) is positive downward, zeroed
      % at the first finite sample of the window.
      aws.ablation = surfaceLowering(aws.z_ice_surf);
      surface_meta.surface_channel = "ablation";
      surface_meta.surface_source = "L3 z_ice_surf";

      % Snow depth: the L3 snow_height channel, clamped >= 0 per the readme.
      % Negatives are counted (provenance), then clamped, never deleted.
      sd = aws.snow_height;
      surface_meta.snow_depth_negatives = nnz(sd < -1e-6);
      sd(sd < 0) = 0;
      aws.snow_depth = sd;
      surface_meta.snow_depth_source = "L3 snow_height (clamped >= 0)";

      % Flags derived from the underlying sensors + the ablation series. The
      % step flags are detected on the ablation channel (a +down series).
      sflags = icemodel.forcing.helpers.surfaceFlags(aws.z_ice_surf, sensors, ...
         aws.Time, transition_times=transition_times);
      [~, step_record, step_flags] = icemodel.forcing.destepSurface( ...
         aws.Time, aws.ablation, mode="detect", gap_flag=sflags.gap, ...
         transition_times=transition_times, season="ablation");
      surf_channels = ["ablation", "snow_depth", "surface_height_flag", ...
         "station_transition_flag", "step_detected_flag", ...
         "step_correctable_flag", "step_magnitude"];
   else
      surface_meta.site_surface_type = "accumulation";

      % Accumulation, percolation, and bedrock sites ship no z_ice_surf. Emit
      % z_surf_combined as a NET surface-height channel (positive up), NOT an
      % ablation channel. Add no snow_depth channel, because snow_height here is
      % not a true snow-over-ice depth.
      if ismember("z_surf_combined", awsnames0)
         z = aws.z_surf_combined;
         first = find(isfinite(z), 1);
         if isempty(first)
            aws.surface_height = nan(size(z));
         else
            aws.surface_height = z - z(first);
         end
         sflags = icemodel.forcing.helpers.surfaceFlags(z, sensors, aws.Time, ...
            transition_times=transition_times);
         [~, step_record, step_flags] = icemodel.forcing.destepSurface( ...
            aws.Time, aws.surface_height, mode="detect", gap_flag=sflags.gap, ...
            transition_times=transition_times, season="accumulation");
         surface_meta.surface_channel = "surface_height";
         surface_meta.surface_source = "L3 z_surf_combined";
      else
         aws.surface_height = nan(height(aws), 1);
         sflags = struct('gap', ones(height(aws), 1), ...
            'station_transition', zeros(height(aws), 1));
         [~, step_record, step_flags] = icemodel.forcing.destepSurface( ...
            aws.Time, aws.surface_height, mode="detect", season="accumulation");
         surface_meta.surface_channel = "surface_height";
         surface_meta.surface_source = "none (no L3 surface-height channel)";
      end
      surface_meta.snow_depth_source = "n/a (accumulation site: no snow_height)";
      surf_channels = ["surface_height", "surface_height_flag", ...
         "station_transition_flag", "step_detected_flag", ...
         "step_correctable_flag", "step_magnitude"];
   end

   % Attach the per-sample flag channels. They are masks, and this code changes
   % no GEUS data. surface_height_flag is the sensor-derived gap-bridged mask.
   % station_transition_flag marks station-handover windows. The step_* channels
   % are the staged de-stepping detection (step_detected, step_correctable, and
   % the signed magnitude). The staged .mat therefore matches the source, and a
   % caller applies de-stepping at analysis time through
   % icemodel.forcing.destepSurface.
   aws.surface_height_flag = sflags.gap;
   aws.station_transition_flag = sflags.station_transition;
   aws.step_detected_flag = step_flags.step_detected;
   aws.step_correctable_flag = step_flags.step_correctable;
   aws.step_magnitude = step_flags.step_magnitude;

   surface_meta.gap_flagged_samples = nnz(aws.surface_height_flag == 1);
   surface_meta.station_transition_samples = ...
      nnz(aws.station_transition_flag == 1);
   surface_meta.composing_stations = composing_stations;
   surface_meta.is_multistation = numel(composing_stations) > 1;
   surface_meta.station_transition_times = transition_times;
   surface_meta.station_transition_record = transition_record;
   surface_meta.steps_detected = numel(step_record);
   if isempty(step_record)
      surface_meta.steps_correctable = 0;
   else
      surface_meta.steps_correctable = ...
         nnz([step_record.classification] == "unambiguous");
   end
   surface_meta.step_record = step_record;

   % Order the output channels (the ice-temperature string is variable length).
   % Thermistor channels with no finite samples (sensors absent on this
   % station's record) are dropped; tice10m is the primary subsurface evaluation
   % channel and is kept first in the string block. Its source and QC channels
   % remain adjacent diagnostics, not members of the sensor string.
   awsnames = string(aws.Properties.VariableNames);
   tice_diagnostics = ["tice10m", "tice10m_source", "tice10m_qc_flag"];
   tice = awsnames(startsWith(awsnames, "tice") ...
      & ~ismember(awsnames, tice_diagnostics));
   tice = tice(arrayfun(@(v) any(isfinite(aws.(v))), tice));
   dtice = awsnames(~cellfun('isempty', ...
      regexp(cellstr(awsnames), '^dtice\d+$', 'once')));
   dtice = dtice(arrayfun(@(v) any(isfinite(aws.(v))), dtice));
   % Keep tice10m even when a mask covers a whole window. Dropping it would
   % break the source/target/flag alignment when QC finds no tice10m samples
   % usable.
   tice10m = awsnames(awsnames == "tice10m");
   tice10m_diagnostics = intersect(["tice10m_source", ...
      "tice10m_qc_flag"], awsnames, 'stable');
   channels = ["tair", "tsfc", "swd", "swu", "lwd", "lwu", "swn", ...
      "lwn", "netr", "shf", "lhf", "thf", "albedo", "cfrac", "rh", ...
      "wspd", "wdir", "psfc", surf_channels, tice10m, ...
      tice10m_diagnostics, tice, dtice];
   channels = channels(ismember(channels, awsnames));
   Data = aws(:, cellstr(channels));

   % Apply physical-range clamps without interpolation. Observational gaps stay
   % missing except the explicit whole-hour deep-night shortwave zeros selected
   % above. The gap flag is a 0/1 mask: metchecks leaves it untouched (not a
   % clamp var, and fillgaps=false), so it stays a faithful per-sample mask.
   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, ...
      fillgaps=false);

   if kwargs.frequency == "daily"
      % Daily means for the physical channels. The 0/1 flag channels aggregate
      % by MAX, so a day with any flagged hour stays flagged. A mean would turn
      % the binary mask into a fraction. step_magnitude keeps the signed jump of
      % largest magnitude in the day (max over |.|), so the daily series still
      % reports the size of any step that day.
      flag_channels = intersect(["surface_height_flag", ...
         "station_transition_flag", "step_detected_flag", ...
         "step_correctable_flag", "tice10m_qc_flag"], ...
         string(Data.Properties.VariableNames));
      flags_daily = retime(Data(:, cellstr(flag_channels)), 'daily', 'max');
      if ismember("step_magnitude", string(Data.Properties.VariableNames))
         mag_daily = retime(Data(:, "step_magnitude"), 'daily', ...
            @(x) maxAbs(x));
      end
      Data = retime(Data, 'daily', 'mean');
      for fc = flag_channels
         Data.(fc) = flags_daily.(fc);
      end
      if ismember("step_magnitude", string(Data.Properties.VariableNames))
         Data.step_magnitude = mag_daily.step_magnitude;
      end
      % A daily mean must not revive tice10m on a day containing an invalid
      % hourly discontinuity; the source daily mean remains available.
      if all(ismember(["tice10m", "tice10m_qc_flag"], ...
            string(Data.Properties.VariableNames)))
         Data.tice10m(Data.tice10m_qc_flag > 0) = NaN;
      end
   end

   % Per-variable units come from the shared unit map.
   Data.Properties.VariableUnits = icemodel.forcing.helpers.variableUnits( ...
      string(Data.Properties.VariableNames));

   % Attach the same location metadata used by every forcing-family Data
   % builder.
   location_metadata = struct( ...
      'lat_wgs84', source_meta.lat, 'lon_wgs84', source_meta.lon, ...
      'elev_m', source_meta.elev);
   Data = icemodel.forcing.helpers.attachLocationMetadata( ...
      Data, location_metadata);

   metadata = source_meta;
   metadata.checks = checks;
   metadata.frequency = kwargs.frequency;
   if ismember("tice10m_qc_flag", string(Data.Properties.VariableNames))
      % Reader counts describe hourly source samples; artifact counts must match
      % the emitted cadence (daily flags aggregate by max).
      metadata.tice10m_qc_flagged_sample_count = ...
         nnz(Data.tice10m_qc_flag > 0);
      metadata.tice10m_qc_failed_sample_count = ...
         nnz(Data.tice10m_qc_flag == 1);
      metadata.tice10m_qc_unreviewed_sample_count = ...
         nnz(Data.tice10m_qc_flag >= 2);
      metadata.tice10m_qc_persistent_sample_count = ...
         nnz(Data.tice10m_qc_flag == 3);
   end

   % Carry source-faithful radiation provenance into both userdata and the
   % observations.mat timetable that reuses this exact Data product.
   fields = fieldnames(shortwave_meta);
   for k = 1:numel(fields)
      metadata.(fields{k}) = shortwave_meta.(fields{k});
   end
   metadata.albedo_source_present = has_albedo_source;
   metadata.albedo_observations_present = has_albedo_observations;
   metadata.albedo_policy = albedo_policy;
   metadata.swn_negative_invalid_count = swn_negative_invalid_count;
   metadata.swn_policy = "swd - swu; finite negative derived net shortwave " ...
      + "is physically inconsistent and remains NaN; component observations " ...
      + "remain source-faithful";
   metadata.site_surface_type = surface_meta.site_surface_type;
   metadata.surface_channel = surface_meta.surface_channel;
   metadata.surface_source = surface_meta.surface_source;
   % Back-compat alias for callers/tests that read metadata.ablation_source.
   metadata.ablation_source = surface_meta.surface_source;
   metadata.snow_depth_source = surface_meta.snow_depth_source;
   if isfield(surface_meta, 'snow_depth_negatives')
      metadata.snow_depth_negatives = surface_meta.snow_depth_negatives;
   end
   metadata.gap_flagged_samples = surface_meta.gap_flagged_samples;
   metadata.station_transition_samples = ...
      surface_meta.station_transition_samples;
   metadata.composing_stations = surface_meta.composing_stations;
   metadata.is_multistation = surface_meta.is_multistation;
   metadata.station_transition_times = surface_meta.station_transition_times;
   metadata.station_transition_record = surface_meta.station_transition_record;
   metadata.steps_detected = surface_meta.steps_detected;
   metadata.steps_correctable = surface_meta.steps_correctable;
   metadata.step_record = surface_meta.step_record;
   % Flagging rule: keep the GEUS series as delivered and ATTACH per-sample
   % flags. This code never edits the data. The gap-bridged surface height keeps
   % a valid CUMULATIVE and visual trend (readme: "the surface height trend over
   % the entire period should be unaffected by the gaps"). Only per-timestep
   % RATE diagnostics through a gap are unreliable, so RATE-based scoring flags
   % or excludes surface_height_flag==1 segments. Cumulative and visual
   % comparison use the FULL series. station_transition_flag marks known
   % AWS-handover windows. The step_* channels stage de-stepping DETECTION only.
   % The de-stepping CORRECTION is opt-in at analysis time through
   % icemodel.forcing.destepSurface (default: correct UNAMBIGUOUS steps only).
   metadata.gap_policy = ["no temporal interpolation (observational); " ...
      "missing shortwave is zero only for whole-hour deep civil night; " ...
      "other source gaps remain missing; clamps applied; " ...
      "surface_height_flag=1 marks gap-bridged samples (RATE-based " ...
      "diagnostics exclude them; cumulative/visual comparison uses the full " ...
      "series). station_transition_flag marks AWS-handover windows. step_* " ...
      "channels stage de-stepping DETECTION; correction is opt-in via " ...
      "icemodel.forcing.destepSurface (default: unambiguous steps only)."];
   metadata = icemodel.forcing.helpers.columnizeMetadata(metadata);
   Data.Properties.UserData = metadata;
end

%% Local functions
function m = maxAbs(x)
   %MAXABS The signed value of largest magnitude in x (0 if x is empty/all-NaN).
   %
   % Aggregates the signed step_magnitude channel to daily resolution. The step
   % of one day is its jump of largest magnitude, and the sign stays.
   x = x(isfinite(x));
   if isempty(x)
      m = 0;
   else
      [~, i] = max(abs(x));
      m = x(i);
   end
end

function lowering = surfaceLowering(z)
   %SURFACELOWERING Cumulative surface lowering from an L3 ice-surface height.
   %
   % z (z_ice_surf) is the ice surface relative to installation and decreases as
   % the surface lowers, so lowering = -(z - z(first finite)) is positive
   % downward, zeroed at the first finite sample of the window.
   first = find(isfinite(z), 1);
   if isempty(first)
      lowering = nan(size(z));
   else
      lowering = -(z - z(first));
   end
end
