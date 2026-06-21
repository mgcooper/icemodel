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
   % net fluxes, the QC'd L3 surface-height channels (snow depth and ice-
   % surface ablation), and the ice-temperature string. Location metadata
   % attaches as table CustomProperties so the result can be written with
   % icemodel.forcing.helpers.writeuserdata.
   %
   % Snow depth: the L3 snow_height channel (snow surface height relative to
   % the ice surface), read directly. It is already QC'd, so no homegrown
   % boom-reference derivation is applied.
   %
   % Ablation: cumulative surface lowering from the L3 ice-surface-height
   % channel z_ice_surf (ice surface relative to its installation height),
   % as -(z_ice_surf - z_ice_surf at the window start), positive downward.
   % Accumulation-zone stations that ship no z_ice_surf fall back to
   % z_surf_combined (the multi-sensor total-surface height); the channel
   % used is recorded in metadata.ablation_source. The L3 channels are QC'd,
   % so no despiking or service-window correction is applied.
   %
   % Ice temperatures: tice1..ticeN from the L3 t_i_* string, read and
   % clamped to the dictionary physical range in icemodel.forcing.
   % readPromiceAws.
   %
   % Gap policy: observational channels are NOT gap-filled (missing data
   % stays missing so evaluations are honest); the physical-range clamps of
   % metchecks are still applied.
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
   %             Slope, ScalarUnits). Channels: tair, tsfc [K]; swd, swu,
   %             lwd, lwu, swn, lwn, netr, shf, lhf, thf [W m-2]; albedo,
   %             cfrac [-]; rh [%]; wspd [m s-1]; wdir [deg]; psfc [Pa];
   %             ablation, snow_depth [m]; tice1..ticeN [K].
   %  metadata - provenance: source file, station, ablation source channel,
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

   % Derived net fluxes.
   aws.swn = aws.swd - aws.swu;
   aws.lwn = aws.lwd - aws.lwu;
   aws.netr = aws.swn + aws.lwn;
   aws.thf = aws.shf + aws.lhf;

   % Snow depth: the QC'd L3 snow_height channel, read directly.
   aws.snow_depth = aws.snow_height;

   % Cumulative ice-surface ablation from the QC'd L3 ice-surface height.
   % z_ice_surf is the ice surface relative to its installation height
   % (decreases as the surface lowers); accumulation-zone stations ship only
   % z_surf_combined. Ablation is the surface lowering, positive downward,
   % zeroed at the first finite sample of the requested window.
   [aws.ablation, ablation_source] = cumulativeAblation(aws);

   % Order the output channels (the ice-temperature string is variable
   % length). Thermistor channels with no finite samples (sensors absent on
   % this station's record) are dropped.
   awsnames = string(aws.Properties.VariableNames);
   tice = awsnames(startsWith(awsnames, "tice"));
   tice = tice(arrayfun(@(v) any(isfinite(aws.(v))), tice));
   channels = ["tair", "tsfc", "swd", "swu", "lwd", "lwu", "swn", ...
      "lwn", "netr", "shf", "lhf", "thf", "albedo", "cfrac", "rh", ...
      "wspd", "wdir", "psfc", "ablation", "snow_depth", tice];
   channels = channels(ismember(channels, awsnames));
   Data = aws(:, cellstr(channels));

   % Physical-range clamps only; observational gaps stay missing.
   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, ...
      fillgaps=false);

   if kwargs.frequency == "daily"
      Data = retime(Data, 'daily', 'mean');
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
   metadata.ablation_source = ablation_source;
   metadata.snow_depth_source = "L3 snow_height";
   metadata.gap_policy = "no gap fill (observational); clamps applied";
end

%% Local functions
function [ablation, source] = cumulativeAblation(aws)
   %CUMULATIVEABLATION Cumulative surface lowering from the L3 surface height.
   %
   % Prefers z_ice_surf (ice surface relative to installation); falls back to
   % z_surf_combined (multi-sensor total surface) when z_ice_surf is absent
   % (accumulation-zone stations). The L3 channel decreases as the surface
   % lowers, so ablation = -(z - z(first finite)) is positive downward,
   % zeroed at the first finite sample of the window.
   names = string(aws.Properties.VariableNames);
   if ismember("z_ice_surf", names)
      z = aws.z_ice_surf;
      source = "L3 z_ice_surf";
   elseif ismember("z_surf_combined", names)
      z = aws.z_surf_combined;
      source = "L3 z_surf_combined";
   else
      ablation = nan(height(aws), 1);
      source = "none (no L3 surface-height channel)";
      return
   end
   first = find(isfinite(z), 1);
   if isempty(first)
      ablation = nan(size(z));
   else
      ablation = -(z - z(first));
   end
end
