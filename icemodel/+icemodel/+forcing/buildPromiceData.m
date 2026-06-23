function [Data, metadata] = buildPromiceData(site, kwargs)
   %BUILDPROMICEDATA Build a PROMICE evaluation/userdata Data timetable.
   %
   %  [Data, metadata] = icemodel.forcing.buildPromiceData(site)
   %  [Data, metadata] = ... buildPromiceData(site, source_dir=..., ...
   %     startdate=..., enddate=..., frequency="daily")
   %
   % Reads the station's pypromice L3 hourly NetCDF and assembles the
   % observational channels used to evaluate icemodel output and to feed the
   % met-swap (userdata) mechanism: the surface energy balance terms, derived
   % net fluxes, the QC'd L3 surface-height channels, and the ice-temperature
   % string. Location metadata attaches as table CustomProperties so the
   % result can be written with icemodel.forcing.helpers.writeuserdata.
   %
   % SITE-TYPE BRANCHING (the core correctness rule, per the product readme):
   % z_ice_surf and snow_height (as a true snow-over-ice depth) are provided
   % ONLY at ABLATION sites; accumulation-zone sites ship only z_surf_combined
   % (a cumulative net-surface height relative to installation that never
   % returns to zero). The builder branches on the PRESENCE of z_ice_surf in
   % the L3 file (the operational ablation-site signal), which is recorded in
   % metadata.site_surface_type:
   %
   %  ABLATION sites (z_ice_surf present):
   %    ablation   = -(z_ice_surf - z_ice_surf(window start)), positive down
   %                 (cumulative ice-surface lowering rel. installation).
   %    snow_depth = snow_height, clamped >= 0 per the readme ("strictly
   %                 positive"); any negative input sample is counted in
   %                 metadata.snow_depth_negatives and clamped, not deleted.
   %
   %  ACCUMULATION / percolation / bedrock sites (no z_ice_surf):
   %    surface_height = z_surf_combined - z_surf_combined(window start), a
   %                 net surface-height change (positive UP); this is NOT an
   %                 ablation channel and is NOT relabeled as one. No ablation
   %                 or snow_depth channel is fabricated for these sites.
   %
   % SURFACE FLAGS (we flag, never silently fix): companion per-sample flag
   % channels attach to whichever surface-height-derived channel the site carries
   % (ablation at ablation sites, surface_height at accumulation sites). Data are
   % FLAGGED, never deleted, and the staged values are the raw GEUS series.
   %  - surface_height_flag (0/1): gap-bridged mask derived from the underlying
   %    L3 sensors (a sample is gap-bridged when the surface value is finite but
   %    every surface-ranging sensor is NaN -> slope-interpolated, not measured;
   %    see icemodel.forcing.helpers.surfaceFlags). Per the readme the cumulative
   %    TREND is preserved through gaps; only per-timestep RATE diagnostics should
   %    exclude flag==1 samples (cumulative/visual comparison uses the full
   %    series).
   %  - station_transition_flag (0/1): marks station-handover windows (a step,
   %    not a NaN, so the gap flag never sees it). Populated from the per-station
   %    install dates in AWS_stations_metadata.csv (within-record handovers only;
   %    see icemodel.forcing.helpers.stationTransitionTimes); all-zero when the
   %    CSV is absent. A detected step coincident with such a window also gains
   %    the 'station_transition' evidence line in icemodel.forcing.destepSurface.
   %    metadata.is_multistation records the merge fact; metadata.
   %    station_transition_times / station_transition_record carry the dates.
   %  - step_detected_flag / step_correctable_flag / step_magnitude: the staged
   %    de-stepping DETECTION (icemodel.forcing.destepSurface in detect mode).
   %    The CORRECTION is opt-in at analysis time (default: unambiguous only),
   %    never baked into the staged data.
   %
   % Ice temperatures: tice1..ticeN from the L3 t_i_* string (surfaced
   % thermistors discarded, clamped to the dictionary physical range in
   % icemodel.forcing.readPromiceAws) plus tice10m, the PRIMARY standardized
   % 10 m subsurface-temperature evaluation channel.
   %
   % tice10m COMPARISON PROTOCOL (the primary subsurface channel). tice10m is
   % GEUS's standardized 10 m-BELOW-the-EVOLVING-SURFACE temperature: GEUS
   % tracks each thermistor's time-dependent depth below the CURRENT surface
   % (d_t_i_*), discards surfaced thermistors, and depth-interpolates the
   % surviving string to 10 m below the current surface at each time step. It is
   % a MOVING (Lagrangian) 10 m depth, not a fixed 10 m from installation. To
   % compare, the model must be sampled at 10 m BELOW ITS OWN CURRENT SURFACE -
   % a moving Lagrangian depth at ablation sites where the surface lowers - NOT
   % at a fixed 10 m from the run start. tice10m is the PRIMARY subsurface
   % comparison channel; the raw tice1..N string needs the per-sensor d_t_i_N
   % depths to place each reading and is SECONDARY / diagnostic.
   %
   % Gap policy: observational channels are NOT gap-filled (missing data
   % stays missing so evaluations are honest); the physical-range clamps of
   % metchecks are still applied. The surface-height channel additionally
   % carries the gap flag described above.
   %
   % Inputs
   %  site - station id ("KAN_M" or compact alias "kanm")
   %
   % Name-value
   %  source_dir : PROMICE NetCDF directory (see readPromiceAws)
   %  startdate, enddate : optional window; default full station record
   %  frequency : "hourly" (default) or "daily" (daily means)
   %
   % Outputs
   %  Data     - timetable with CustomProperties (X, Y, Lat, Lon, Elev,
   %             Slope, ScalarUnits). Common channels: tair, tsfc [K]; swd,
   %             swu, lwd, lwu, swn, lwn, netr, shf, lhf, thf [W m-2]; albedo,
   %             cfrac [-]; rh [%]; wspd [m s-1]; wdir [deg]; psfc [Pa];
   %             tice1..ticeN, tice10m [K]. Surface-height channels are
   %             site-type dependent: ABLATION sites carry ablation +
   %             snow_depth [m]; ACCUMULATION sites carry surface_height [m].
   %             Either carries the surface flag channels: surface_height_flag
   %             (0/1 gap-bridged mask), station_transition_flag (0/1 handover
   %             window), step_detected_flag / step_correctable_flag (0/1) and
   %             step_magnitude [m] (the staged de-stepping detection).
   %  metadata - provenance: source file, station, site_surface_type, surface
   %             channel + source, gap-flag counts, station-transition counts,
   %             composing_stations / is_multistation, steps_detected /
   %             steps_correctable / step_record, snow-depth negative count,
   %             QA/QC gap counts
   %
   % See also: icemodel.forcing.readPromiceAws,
   %  icemodel.forcing.buildPromiceMet,
   %  icemodel.forcing.helpers.writeuserdata

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

   % Derived net fluxes. Sparse stations may ship no radiation or turbulent-
   % flux channels, so each derived term is computed only when its inputs are
   % present (a missing input drops the derived channel rather than erroring).
   % has() is re-bound after channels are added so terms that build on a
   % derived channel (netr on swn/lwn) see it.
   has = @(v) ismember(v, string(aws.Properties.VariableNames));
   if has("swd") && has("swu")
      aws.swn = aws.swd - aws.swu;
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

   % The underlying L3 surface-ranging sensors (transducer/boom/stake) used to
   % derive the gap flag from FIRST PRINCIPLES: a sample is gap-bridged (slope-
   % interpolated, not measured) when every one of these is NaN yet the surface
   % series is finite. Whichever are present on this station are gathered into a
   % matrix for surfaceFlags (z_ice_surf is re-computed by GEUS from z_pt_cor,
   % the transducer, falling back to the stake sonic ranger; the boom adds the
   % accumulation-site surface ranging).
   sensor_names = intersect(["transducer_depth", "boom_height", ...
      "stake_height"], awsnames0, 'stable');
   if isempty(sensor_names)
      sensors = [];
   else
      sensors = cell2mat(arrayfun(@(v) aws.(v), sensor_names, ...
         'UniformOutput', false));
   end

   % Known station-handover times for the station-transition flag. The
   % composing-station NAMES come from the curated catalog; their per-station
   % install DATES come from AWS_stations_metadata.csv (the GEUS thredds product,
   % staged alongside the L3 NetCDFs). stationTransitionTimes maps each composing
   % station to its install date and keeps only WITHIN-RECORD handovers (an
   % install strictly after the record start; the founding station's install
   % begins the record, it is not a handover within it). With the CSV present the
   % flag is now populated; without it transition_times stays empty (all-zero
   % flag) and the merge FACT is still recorded in metadata so destepSurface can
   % recover a transition as a coincident step. The window clamp uses the L3
   % record bounds (source_meta), so an install outside this station's record is
   % excluded.
   info = icemodel.verification.helpers.promicesiteinfo(site);
   composing_stations = info.stations;
   [transition_times, transition_record] = ...
      icemodel.forcing.helpers.stationTransitionTimes(composing_stations, ...
      window_start=source_meta.window_start, ...
      window_end=source_meta.window_end, source_dir=kwargs.source_dir);

   surface_meta = struct();
   if is_ablation
      surface_meta.site_surface_type = "ablation";

      % Cumulative ice-surface ablation from the QC'd L3 ice-surface height.
      % z_ice_surf is the ice surface relative to installation (decreases as
      % the surface lowers); ablation = -(z - z(start)) is positive downward,
      % zeroed at the first finite sample of the window.
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

      % Accumulation / percolation / bedrock sites ship no z_ice_surf: emit
      % z_surf_combined as a NET surface-height channel (positive up), NOT an
      % ablation channel. No snow_depth is fabricated (snow_height here is not
      % a true snow-over-ice depth).
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

   % Attach the per-sample flag channels (faithful masks; we modify no GEUS
   % data). surface_height_flag is the gap-bridged mask (now sensor-derived);
   % station_transition_flag marks station-handover windows; the step_* channels
   % are the staged de-stepping detection (step_detected/correctable + signed
   % magnitude), so the staged .mat is faithful and de-stepping is applied opt-in
   % at analysis time via icemodel.forcing.destepSurface.
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

   % Order the output channels (the ice-temperature string is variable
   % length). Thermistor channels with no finite samples (sensors absent on
   % this station's record) are dropped; tice10m is the primary subsurface
   % evaluation channel and is kept first in the string block.
   awsnames = string(aws.Properties.VariableNames);
   tice = awsnames(startsWith(awsnames, "tice") & awsnames ~= "tice10m");
   tice = tice(arrayfun(@(v) any(isfinite(aws.(v))), tice));
   tice10m = awsnames(awsnames == "tice10m");
   tice10m = tice10m(arrayfun(@(v) any(isfinite(aws.(v))), tice10m));
   channels = ["tair", "tsfc", "swd", "swu", "lwd", "lwu", "swn", ...
      "lwn", "netr", "shf", "lhf", "thf", "albedo", "cfrac", "rh", ...
      "wspd", "wdir", "psfc", surf_channels, tice10m, tice];
   channels = channels(ismember(channels, awsnames));
   Data = aws(:, cellstr(channels));

   % Physical-range clamps only; observational gaps stay missing. The gap
   % flag is a 0/1 mask: metchecks leaves it untouched (not a clamp var, and
   % fillgaps=false), so it stays a faithful per-sample quality mask.
   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, ...
      fillgaps=false);

   if kwargs.frequency == "daily"
      % Daily means for the physical channels; the 0/1 flag channels aggregate
      % by MAX so a day touched by any flagged hour stays flagged (a mean would
      % blur the binary mask into a meaningless fraction). step_magnitude keeps
      % the day's largest-magnitude signed jump (max over |.|), so the daily
      % series still reports the size of any step that day.
      flag_channels = intersect(["surface_height_flag", ...
         "station_transition_flag", "step_detected_flag", ...
         "step_correctable_flag"], string(Data.Properties.VariableNames));
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
   end

   % Per-variable units from the shared canonical map.
   Data.Properties.VariableUnits = icemodel.forcing.helpers.variableUnits( ...
      string(Data.Properties.VariableNames));

   % Location metadata for the Data-file (userdata) contract.
   proj = icemodel.forcing.helpers.psnProjection();
   [x, y] = projfwd(proj, source_meta.lat, source_meta.lon);
   Data = addprop(Data, ...
      {'X', 'Y', 'Lat', 'Lon', 'Elev', 'Slope', 'ScalarUnits'}, ...
      repmat({'table'}, 1, 7));
   Data.Properties.CustomProperties.X = x;
   Data.Properties.CustomProperties.Y = y;
   Data.Properties.CustomProperties.Lat = source_meta.lat;
   Data.Properties.CustomProperties.Lon = source_meta.lon;
   Data.Properties.CustomProperties.Elev = source_meta.elev;
   Data.Properties.CustomProperties.Slope = NaN;
   Data.Properties.CustomProperties.ScalarUnits = ...
      ["m", "m", "degrees N", "degrees W", "m asl", "m/m"];

   metadata = source_meta;
   metadata.checks = checks;
   metadata.frequency = kwargs.frequency;
   metadata.site_surface_type = surface_meta.site_surface_type;
   metadata.surface_channel = surface_meta.surface_channel;
   metadata.surface_source = surface_meta.surface_source;
   % Back-compat aliases for callers/tests reading the legacy metadata names.
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
   % Flagging philosophy: we preserve the authoritative GEUS series and ATTACH
   % per-sample flags; we never silently edit data. The gap-bridged surface
   % height keeps a valid CUMULATIVE/visual trend (readme: "the surface height
   % trend over the entire period should be unaffected by the gaps"); only per-
   % timestep RATE diagnostics through a gap are unreliable, so RATE-based
   % scoring flags/excludes surface_height_flag==1 segments while cumulative and
   % visual comparison use the FULL series. station_transition_flag marks known
   % AWS-handover windows. The step_* channels stage de-stepping DETECTION only;
   % the de-stepping CORRECTION is opt-in at analysis time via
   % icemodel.forcing.destepSurface (default: correct UNAMBIGUOUS steps only).
   metadata.gap_policy = ["no gap fill (observational); clamps applied; " ...
      "surface_height_flag=1 marks gap-bridged samples (RATE-based " ...
      "diagnostics exclude them; cumulative/visual comparison uses the full " ...
      "series). station_transition_flag marks AWS-handover windows. step_* " ...
      "channels stage de-stepping DETECTION; correction is opt-in via " ...
      "icemodel.forcing.destepSurface (default: unambiguous steps only)."];
end

%% Local functions
function m = maxAbs(x)
   %MAXABS The signed value of largest magnitude in x (0 if x is empty/all-NaN).
   %
   % Used to aggregate the signed step_magnitude channel to daily resolution: a
   % day's step is summarised by its largest-magnitude jump, keeping the sign.
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
   % z (z_ice_surf) is the ice surface relative to installation and decreases
   % as the surface lowers, so lowering = -(z - z(first finite)) is positive
   % downward, zeroed at the first finite sample of the window.
   first = find(isfinite(z), 1);
   if isempty(first)
      lowering = nan(size(z));
   else
      lowering = -(z - z(first));
   end
end
