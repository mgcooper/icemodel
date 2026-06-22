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
   % GAP FLAG: surface-height change ACROSS DATA GAPS is not a direct
   % observation (the readme bridges gaps by manual slope). A companion
   % per-sample flag channel, surface_height_flag (0 = direct observation,
   % 1 = gap-bridged / interpolated), marks samples whose underlying L3
   % surface sensor is missing so model-vs-observation comparison can EXCLUDE
   % gap-bridged segments. Data are FLAGGED, never deleted. The flag attaches
   % to whichever surface-height-derived channel the site carries (ablation
   % at ablation sites, surface_height at accumulation sites).
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
   %             Either carries surface_height_flag (0/1 gap-bridged mask).
   %  metadata - provenance: source file, station, site_surface_type,
   %             surface channel + source, gap-flag counts, snow-depth
   %             negative count, QA/QC gap counts
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

   surface_meta = struct();
   if is_ablation
      surface_meta.site_surface_type = "ablation";

      % Cumulative ice-surface ablation from the QC'd L3 ice-surface height.
      % z_ice_surf is the ice surface relative to installation (decreases as
      % the surface lowers); ablation = -(z - z(start)) is positive downward,
      % zeroed at the first finite sample of the window. The gap flag marks
      % samples whose underlying z_ice_surf is missing (gap-bridged).
      [aws.ablation, ~] = surfaceLowering(aws.z_ice_surf);
      aws.surface_height_flag = gapFlag(aws.z_ice_surf);
      surface_meta.surface_channel = "ablation";
      surface_meta.surface_source = "L3 z_ice_surf";

      % Snow depth: the L3 snow_height channel, clamped >= 0 per the readme.
      % Negatives are counted (provenance), then clamped, never deleted.
      sd = aws.snow_height;
      surface_meta.snow_depth_negatives = nnz(sd < -1e-6);
      sd(sd < 0) = 0;
      aws.snow_depth = sd;
      surface_meta.snow_depth_source = "L3 snow_height (clamped >= 0)";

      surf_channels = ["ablation", "snow_depth", "surface_height_flag"];
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
         aws.surface_height_flag = gapFlag(z);
         surface_meta.surface_channel = "surface_height";
         surface_meta.surface_source = "L3 z_surf_combined";
      else
         aws.surface_height = nan(height(aws), 1);
         aws.surface_height_flag = ones(height(aws), 1);
         surface_meta.surface_channel = "surface_height";
         surface_meta.surface_source = "none (no L3 surface-height channel)";
      end
      surface_meta.snow_depth_source = "n/a (accumulation site: no snow_height)";
      surf_channels = ["surface_height", "surface_height_flag"];
   end
   surface_meta.gap_flagged_samples = nnz(aws.surface_height_flag == 1);

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
      % Daily means for the physical channels; the gap flag aggregates by MAX
      % so a day touched by any gap-bridged hour stays flagged (a mean would
      % blur the 0/1 mask into a meaningless fraction).
      flag = retime(Data(:, "surface_height_flag"), 'daily', 'max');
      Data = retime(Data, 'daily', 'mean');
      Data.surface_height_flag = flag.surface_height_flag;
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
   metadata.gap_policy = ["no gap fill (observational); clamps applied; " ...
      "surface_height_flag=1 marks gap-bridged surface-height samples " ...
      "(exclude before model comparison)"];
end

%% Local functions
function [lowering, source] = surfaceLowering(z)
   %SURFACELOWERING Cumulative surface lowering from an L3 ice-surface height.
   %
   % z (z_ice_surf) is the ice surface relative to installation and decreases
   % as the surface lowers, so lowering = -(z - z(first finite)) is positive
   % downward, zeroed at the first finite sample of the window.
   source = "L3 z_ice_surf";
   first = find(isfinite(z), 1);
   if isempty(first)
      lowering = nan(size(z));
   else
      lowering = -(z - z(first));
   end
end

function flag = gapFlag(z)
   %GAPFLAG Per-sample quality flag for a surface-height series.
   %
   % 0 = direct observation (the underlying L3 surface sensor reports a finite
   % value at this sample); 1 = gap-bridged (the sensor is missing here, so a
   % surface-height change spanning this sample is interpolated/slope-bridged,
   % not a direct observation per the product readme). Leading/trailing NaNs
   % before the first / after the last finite sample are also flagged 1. The
   % flag is returned for the full series so it aligns with the timetable.
   flag = double(~isfinite(z));
end
