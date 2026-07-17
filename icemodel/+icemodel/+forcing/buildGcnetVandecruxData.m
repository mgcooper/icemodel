function [Data, metadata] = buildGcnetVandecruxData(station, kwargs)
   %BUILDGCNETVANDECRUXDATA Build GC-Net/Vandecrux surface forcing Data.
   %
   %  [Data, metadata] = icemodel.forcing.buildGcnetVandecruxData(station)
   %  [Data, metadata] = ... buildGcnetVandecruxData(station, source_dir=...)
   %
   % Reads the Vandecrux et al. GC-Net surface/SEB NetCDF product
   % (<station>_surface.nc) and maps the source channels onto icemodel's
   % canonical forcing/userdata names. The target RetMIP stations are Dye-2
   % long ("DYE_2", aliases "dye2", "dye2_long") and Summit ("Summit",
   % aliases "sum", "summit").
   %
   % Source precipitation policy: the files carry snowfall estimates but no
   % rain channel. The builder converts source snowfall amounts [m_weq per
   % source timestep] to the canonical snowfall rate `snowf` [m s-1], adds
   % `rainf` as all-NaN, and leaves total precipitation to data2met. This keeps
   % missing precipitation missing; no absent channel is zero-filled.
   %
   % Longwave policy: LRin is source-filled regional-climate-model longwave, not
   % an observed GC-Net longwave sensor. Metadata records both the settled
   % RetMIP/Vandecrux HIRHAM5 context and the local surface-package RACMO2.3p2
   % attribute text so downstream comparisons do not treat lwd as observed.
   %
   % Outputs
   %  Data     - timetable with canonical channel names, userdata location
   %             CustomProperties, and source provenance in Properties.UserData.
   %  metadata - provenance, channel mapping, unit/policy notes, and QA checks.
   %
   % See also: icemodel.forcing.buildGcnetVandecruxMet,
   %  icemodel.forcing.data2met, icemodel.forcing.helpers.writeuserdata

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

   [source_dir, station] = icemodel.forcing.helpers.gcnetVandecruxInputs( ...
      kwargs.source_dir, station);
   filename = locateSurfaceFile(source_dir, station);
   station_info = icemodel.forcing.helpers.gcnetVandecruxStationMetadata( ...
      station);

   info = ncinfo(filename);
   source_names = string({info.Variables.Name});
   if ~ismember("time", source_names)
      error('icemodel:forcing:buildGcnetVandecruxData:missingTime', ...
         'GC-Net/Vandecrux surface file has no time coordinate: %s', filename)
   end

   Time = regularizeHourlyAxis(icemodel.forcing.helpers.gcnetTime( ...
      double(ncread(filename, 'time')), ...
      icemodel.forcing.helpers.readNetcdfAttribute(filename, "time", ...
      "units")));
   keep = icemodel.forcing.helpers.timeWindowMask( ...
      Time, window_start, window_end);
   if ~any(keep)
      error('icemodel:forcing:buildGcnetVandecruxData:emptyWindow', ...
         'requested window does not overlap %s', filename)
   end

   % The source posts hourly timestep amounts for mass terms. Convert snowfall
   % to the met-contract m s-1 rate and mass diagnostics to mWE/h rates. The
   % actual files are hourly, but deriving dt from the coordinate makes the
   % policy explicit for fixtures and future source revisions.
   dt_hours = timeStepHours(Time);
   dt_seconds = dt_hours * 3600;

   Data = timetable(Time(keep));
   Data = mapIfPresent(Data, filename, source_names, "Ta_2m", "tair", keep, 1);
   Data = mapIfPresent(Data, filename, source_names, "Tsurf", "tsfc", keep, 1);
   Data = mapIfPresent(Data, filename, source_names, "SRin", "swd", keep, 1);
   Data = mapIfPresent(Data, filename, source_names, "SRout", "swu", keep, 1);
   Data = mapIfPresent(Data, filename, source_names, "LRin", "lwd", keep, 1);
   Data = mapIfPresent(Data, filename, source_names, "LRout", "lwu", keep, 1);
   Data = mapIfPresent(Data, filename, source_names, "SHF", "shf", keep, 1);
   Data = mapIfPresent(Data, filename, source_names, "LHF", "lhf", keep, 1);
   Data = mapIfPresent(Data, filename, source_names, "RH_2m", "rh", keep, 1);
   Data = mapIfPresent(Data, filename, source_names, "WS_10m", "wspd", keep, 1);
   Data = mapIfPresent(Data, filename, source_names, "Pres", "psfc", keep, 100);
   Data = mapIfPresent(Data, filename, source_names, "snowfall", "snowf", ...
      keep, 1 / dt_seconds);
   Data = mapIfPresent(Data, filename, source_names, "melt", "melt", ...
      keep, 1 / dt_hours);
   % Vandecrux declares sublimation as negative and deposition as positive.
   % Flip once at the source boundary to the canonical positive-loss sign.
   Data = mapIfPresent(Data, filename, source_names, "sublimation", "subl", ...
      keep, -1 / dt_hours);
   Data = mapIfPresent(Data, filename, source_names, "SMB", "smb", ...
      keep, 1 / dt_hours);

   % The source has no rain channel. Keep the split explicit so data2met derives
   % ppt as NaN, not as snowfall-only precipitation or fabricated zero rain.
   Data.rainf = nan(height(Data), 1);

   % Derived radiative/turbulent diagnostics stay alongside the mapped source
   % fluxes for userdata and evaluation, while downstream processing can still
   % recompute the minimal met set from swd/lwd/albedo/tsfc.
   albedo_qc_counts = struct( ...
      'low_light', 0, ...
      'low_solar_elevation', 0, ...
      'nonpositive_swu', 0, ...
      'below_minimum', 0, ...
      'total', 0);
   albedo_transient_qc = struct();
   albedo_minimum = 0.3;
   if station == "Summit"
      % Summit is a dry-snow accumulation site. Ratios below the conservative
      % dry-snow floor are not credible surface states there, but the measured
      % SWD/SWU pair remains staged for source diagnosis.
      albedo_minimum = 0.7;
   end
   names = string(Data.Properties.VariableNames);
   if all(ismember(["swd", "swu"], names))
      Data.swn = Data.swd - Data.swu;
      [Data.albedo, albedo_qc_counts] = ...
         icemodel.forcing.helpers.sourceAlbedo( ...
         Data.swd, Data.swu, minimum=albedo_minimum, Time=Data.Time, ...
         latitude=station_info.site_location.lat_wgs84, ...
         longitude=station_info.site_location.lon_wgs84);
   end
   names = string(Data.Properties.VariableNames);
   if all(ismember(["lwd", "lwu"], names))
      Data.lwn = Data.lwd - Data.lwu;
   end
   names = string(Data.Properties.VariableNames);
   if all(ismember(["swn", "lwn"], names))
      Data.netr = Data.swn + Data.lwn;
   end
   names = string(Data.Properties.VariableNames);
   if all(ismember(["shf", "lhf"], names))
      Data.thf = Data.shf + Data.lhf;
   end

   % Preserve the measured radiation but exclude conservative, recovered
   % reflected-shortwave collapse episodes from derived forcing/evaluation
   % channels. The shared helper owns the cross-day policy and provenance.
   names = string(Data.Properties.VariableNames);
   if all(ismember(["swd", "swu"], names))
      [transient_rows, transient_report] = ...
         icemodel.forcing.helpers.dailyAlbedoAnomalyFlags( ...
         Data.Time, Data.swd, Data.swu);
      albedo_transient_qc = rmfield(transient_report, 'diagnostics');
      if any(transient_rows)
         Data.albedo(transient_rows) = NaN;
         Data.swn(transient_rows) = NaN;
         if ismember("netr", names)
            Data.netr(transient_rows) = NaN;
         end
      end
   end
   if ismember("H_surf_obs", source_names)
      surface_height = double(ncread(filename, 'H_surf_obs'));
      Data.surface_height = relativeToFirstFinite(surface_height(keep));
      surface_height_source = "H_surf_obs";
   elseif ismember("H_surf_mod", source_names)
      surface_height = double(ncread(filename, 'H_surf_mod'));
      Data.surface_height = relativeToFirstFinite(surface_height(keep));
      surface_height_source = "H_surf_mod";
   else
      surface_height_source = "";
   end

   [Data, checks] = icemodel.forcing.helpers.metchecks(Data, ...
      fillgaps=kwargs.fillgaps);

   % Keep the output stable and labelled through the canonical metadata map.
   preferred = ["tair", "tsfc", "swd", "swu", "lwd", "lwu", "swn", ...
      "lwn", "netr", "shf", "lhf", "thf", "albedo", "rh", "wspd", ...
      "psfc", "rainf", "snowf", "melt", "subl", "smb", "surface_height"];
   names = string(Data.Properties.VariableNames);
   ordered = preferred(ismember(preferred, names));
   Data = Data(:, cellstr(ordered));
   Data = icemodel.forcing.helpers.stampMetadata(Data);

   metadata = icemodel.forcing.helpers.columnizeMetadata( ...
      sourceMetadata(filename, station, source_names, info, ...
      dt_hours, surface_height_source, checks, albedo_qc_counts, ...
      albedo_transient_qc, albedo_minimum));
   Data = icemodel.forcing.helpers.attachLocationMetadata( ...
      Data, metadata.site_location);
   Data.Properties.UserData = metadata;
end

%% Local functions
function filename = locateSurfaceFile(source_dir, station)
   %LOCATESURFACEFILE Find a flat or package-subfolder surface NetCDF.
   filename = icemodel.forcing.helpers.locateNormalizedFile(source_dir, ...
      station + "_surface.nc", ...
      error_id="icemodel:forcing:buildGcnetVandecruxData:fileNotFound", ...
      error_label="GC-Net/Vandecrux surface file");
end

function Data = mapIfPresent(Data, filename, source_names, src, dst, keep, scale)
   %MAPIFPRESENT Copy one source variable to a scaled canonical Data channel.
   if ismember(src, source_names)
      data = double(ncread(filename, src)) .* scale;
      Data.(dst) = data(keep);
   end
end

function data = relativeToFirstFinite(data)
   %RELATIVETOFIRSTFINITE Zero a surface-height series at its first sample.
   first = find(isfinite(data), 1);
   if isempty(first)
      data(:) = NaN;
   else
      data = data - data(first);
   end
end

function metadata = sourceMetadata(filename, station, source_names, info, ...
      dt_hours, surface_height_source, checks, albedo_qc_counts, ...
      albedo_transient_qc, albedo_minimum)
   %SOURCEMETADATA Assemble policy and provenance metadata for the Data file.
   station_info = icemodel.forcing.helpers.gcnetVandecruxStationMetadata( ...
      station);
   metadata = struct( ...
      'source_file', filename, ...
      'source_family', "gcnet_vandecrux", ...
      'source_package', "Vandecrux GC-Net gap-filled surface/SEB", ...
      'station', station, ...
      'aliases', station_info.aliases, ...
      'source_variables', source_names(:), ...
      'source_variable_attributes', ...
         icemodel.forcing.helpers.sourceVariableAttributes(info), ...
      'channel_map', channelMap(), ...
      'site_location', station_info.site_location, ...
      'time_step_hours', dt_hours, ...
      'time_axis_policy', ...
         "source fractional-day coordinate regularized to continuous hourly row-index axis", ...
      'surface_height_source', surface_height_source, ...
      'checks', checks, ...
      'albedo_qc_counts', albedo_qc_counts, ...
      'albedo_transient_qc', albedo_transient_qc, ...
      'mass_flux_policy', ...
         "source m_weq timestep amounts converted to canonical rates", ...
      'gcnet_subl_native_sign_convention', ...
         "negative_loss_positive_deposition", ...
      'gcnet_subl_sign_convention', ...
         "positive_loss_negative_deposition", ...
      'albedo_policy', ...
         sprintf("swu/swd only where solar elevation > 20 degrees, swd >= 10 W m-2, swu > 0, and albedo >= %.1f; conservative recovered daily reflected-shortwave collapse episodes remain NaN", albedo_minimum), ...
      'shortwave_balance_policy', ...
         "source swd/swu preserved; derived albedo/swn/netr are NaN on dailyAlbedoAnomalyFlags episodes", ...
      'precip_policy', ...
         "snowfall -> snowf [m s-1]; rainf = NaN placeholder; ppt = NaN placeholder via data2met because no rain channel exists", ...
      'lwd_policy', ...
         "LRin is source-filled regional-climate-model longwave, not observed; RetMIP/Vandecrux context records HIRHAM5 lineage and local surface XML identifies RACMO2.3p2 for LRin");
end

function map = channelMap()
   %CHANNELMAP Record the canonical-name to Vandecrux-name mapping.
   map = struct( ...
      'tair', "Ta_2m", ...
      'tsfc', "Tsurf", ...
      'swd', "SRin", ...
      'swu', "SRout", ...
      'lwd', "LRin", ...
      'lwu', "LRout", ...
      'shf', "SHF", ...
      'lhf', "LHF", ...
      'rh', "RH_2m", ...
      'wspd', "WS_10m", ...
      'psfc', "Pres", ...
      'snowf', "snowfall", ...
      'melt', "melt", ...
      'subl', "sublimation", ...
      'smb', "SMB", ...
      'surface_height', "H_surf_obs");
end

function time = regularizeHourlyAxis(time)
   %REGULARIZEHOURLYAXIS Use the documented hourly row-index time convention.
   %
   % The real Vandecrux surface files encode time as fractional days, but the
   % coordinate drifts inside leap years and jumps at several year boundaries
   % while the row count and endpoint still describe a continuous hourly series.
   % The product is hourly, so met-building uses the unambiguous row-index axis.
   time = time(:);
   if numel(time) < 2
      return
   end
   idx = 0:numel(time) - 1;
   time = time(1) + hours(idx(:));
end

function dt_hours = timeStepHours(Time)
   %TIMESTEPHOURS Infer the regular source timestep in hours.
   if numel(Time) < 2
      error('icemodel:forcing:buildGcnetVandecruxData:tooFewSamples', ...
         'GC-Net/Vandecrux surface files must contain at least two samples')
   end
   dt_hours = median(hours(diff(Time)), 'omitnan');
end
