function [met, metadata] = buildPromiceMet(site, kwargs)
   %BUILDPROMICEMET Build an icemodel met timetable from PROMICE AWS data.
   %
   %  [met, metadata] = icemodel.forcing.buildPromiceMet(site)
   %  [met, metadata] = ... buildPromiceMet(site, source_dir=..., ...
   %     startdate=..., enddate=..., fillwinter=true, fillwithmissing=true)
   %
   % Reads the station's PROMICE v3 hourly NetCDF (icemodel.forcing.
   % readPromiceAws), selects pypromice's corrected shortwave product with a
   % physical raw fallback, derives zero only for missing whole-hour deep-civil-
   % night radiation, applies the albedo winter-fill policy, runs the standard
   % source-faithful QA/QC pass (physical clamps without unbounded gap fill), and
   % returns a timetable
   % satisfying the icemodel met contract for any station and any window in the
   % source record. Save it with
   % icemodel.forcing.helpers.writemet.
   %
   % PROMICE rainfall_cor_u is corrected liquid precipitation from a tipping-
   % bucket gauge when that channel exists. The builder converts its hourly
   % timestep amount to canonical rainf [m s-1]. PROMICE supplies neither
   % reliable solid precipitation nor total precipitation, so snowf and ppt
   % remain explicit NaN placeholders for a later fill/swap step.
   %
   % Inputs
   %  site - station name ("KAN_M" or compact alias "kanm"); see
   %         icemodel.forcing.readPromiceAws
   %
   % Name-value
   %  source_dir : PROMICE NetCDF directory (see readPromiceAws)
   %  startdate, enddate : optional window; default full station record
   %  fillwinter : fill polar-night albedo with the dry-snow constant
   %               (default true); see icemodel.forcing.fillPromiceAlbedo
   %  fillwithmissing : add absent required channels as NaN placeholders
   %                    instead of rejecting the met build up front (default true)
   %
   % Outputs
   %  met      - timetable: tair [K], swd, swu, lwd [W m-2], albedo [-],
   %             wspd [m s-1], rh [%], psfc [Pa], rainf, snowf, ppt
   %             [m s-1], plus pass-through boom_height, wdir, tsfc, cfrac,
   %             shf, lhf
   %  metadata - provenance: source file, station, lat/lon/elev, QA/QC
   %             gap counts, policies applied
   %
   % See also: icemodel.forcing.readPromiceAws,
   %  icemodel.forcing.buildPromiceData,
   %  icemodel.forcing.helpers.writemet, icemodel.loadmet

   arguments
      site (1, 1) string
      kwargs.source_dir (1, 1) string = ""
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.fillwinter (1, 1) logical = true
      kwargs.fill_lwd (1, 1) logical = false
      kwargs.fillwithmissing (1, 1) logical = true
   end

   [aws, source_meta] = icemodel.forcing.readPromiceAws(site, ...
      source_dir=kwargs.source_dir, timescale="hourly", ...
      startdate=kwargs.startdate, enddate=kwargs.enddate);

   % Preserve known within-record station handovers with the met artifact
   % as provenance for surface-height QC flags and for the future
   % maintenance-visit registry refinement of the runtime interpolation
   % rung (icemodel-1ps.16); the A3 chain itself may interpolate across
   % handovers by design.
   site_info = icemodel.verification.setup.promiceSiteCatalog(site, ...
      source_dir=kwargs.source_dir);
   [transition_times, transition_record] = ...
      icemodel.forcing.helpers.stationTransitionTimes(site_info.stations, ...
      window_start=source_meta.window_start, ...
      window_end=source_meta.window_end, source_dir=kwargs.source_dir);

   % Public forcing uses the corrected pypromice shortwave product where
   % available and a nonnegative raw fallback otherwise. Missing values become
   % zero only when the complete source hour is below deep civil twilight;
   % readPromiceAws itself remains source-faithful and exposes both source pairs.
   [public_swd, public_swu, shortwave_meta] = ...
      icemodel.forcing.helpers.promiceShortwave(aws, fill_darkness=true, ...
      latitude=source_meta.lat, longitude=source_meta.lon, ...
      swd_source_file_observations_present= ...
      source_meta.swd_source_file_observations_present, ...
      swu_source_file_observations_present= ...
      source_meta.swu_source_file_observations_present);

   % lwd source. Normally the observed channel. OPT-IN fallback (fill_lwd, off
   % by default) for stations with no longwave sensor (older GC-Net/firn sites,
   % where dlr is absent or all-NaN in the window): estimate lwd from air
   % temperature + vapor pressure via the legacy empirical relation, FLAGGED as
   % estimated in metadata. The default met-builder contract leaves an absent
   % lwd as an explicit NaN placeholder; fillwithmissing=false is the opt-in
   % strict path for callers that want missing required channels to abort.
   has_lwd = ismember('lwd', aws.Properties.VariableNames) ...
      && any(isfinite(aws.lwd));
   lwd_estimated = false;
   lwd_placeholder = false;
   if has_lwd
      lwd = aws.lwd;
      lwd_policy = "lwd from the observed dlr channel";
   elseif kwargs.fill_lwd
      ea = icemodel.vapor.saturation_vapor_pressure(aws.tair, true) ...
         .* aws.rh ./ 100;
      lwd = icemodel.surface.empirical_incoming_longwave_radiation(aws.tair, ea);
      lwd_estimated = true;
      lwd_policy = ...
         "lwd ESTIMATED from tair + vapor pressure (no longwave sensor; " + ...
         "icemodel.surface.empirical_incoming_longwave_radiation, fill_lwd=true)";
   elseif kwargs.fillwithmissing && ismember('lwd', aws.Properties.VariableNames)
      lwd = aws.lwd;                  % present but all-NaN -> explicit placeholder
      lwd_placeholder = true;
      lwd_policy = "lwd = NaN placeholder (source dlr channel is all missing)";
   elseif kwargs.fillwithmissing
      lwd = nan(height(aws), 1);      % absent -> explicit placeholder
      lwd_placeholder = true;
      lwd_policy = "lwd = NaN placeholder (source dlr channel absent)";
   else
      error('icemodel:forcing:buildPromiceMet:missingForcing', ...
         'required forcing channel(s) absent at %s: lwd (no met built)', site);
   end

   % In strict mode, a station that lacks a required forcing channel aborts with
   % a clean, named error instead of a cryptic "Unrecognized table variable
   % name" from the assembly below. Normal met-builder mode keeps those absent
   % channels as explicit placeholders so a runtime source swap can fill them.
   % lwd is excluded here - it has its own observed/estimated/absent handling
   % above.
   forcing_channels = ["tair", "albedo", "wspd", "rh", "psfc"];
   absent = forcing_channels(~ismember(forcing_channels, ...
      string(aws.Properties.VariableNames)));
   if ~(shortwave_meta.swd_source_present ...
         || shortwave_meta.swd_corrected_source_present)
      absent(end + 1) = "swd";
   end
   if ~isempty(absent) && ~kwargs.fillwithmissing
      error('icemodel:forcing:buildPromiceMet:missingForcing', ...
         'required forcing channel(s) absent at %s: %s (no met built)', ...
         site, strjoin(absent, ', '));
   end
   if ~kwargs.fillwithmissing
      error('icemodel:forcing:buildPromiceMet:missingForcing', ...
         ['required forcing channel(s) absent at %s: ppt ' ...
         '(PROMICE does not supply total precipitation; no met built)'], site);
   end

   % Assemble the met variables from the required forcings. PROMICE does not
   % supply total precipitation, so ppt is an explicit NaN placeholder. A
   % missing OPTIONAL (non-forcing) channel must not block met building.
   Time = aws.Time;
   met = timetable(Time);
   met.tair = channelOrNan(aws, "tair");
   met.swd = public_swd;
   met.swu = public_swu;
   met.lwd = lwd;

   % Distinguish an observed albedo series from an intentional all-missing or
   % absent placeholder. The placeholder stays NaN for a later forcing swap;
   % no constant observation is invented merely to satisfy the met schema.
   has_albedo_source = ismember("albedo", ...
      string(aws.Properties.VariableNames));
   has_albedo_observations = has_albedo_source ...
      && any(isfinite(aws.albedo), 'all');
   if has_albedo_observations
      met.albedo = icemodel.forcing.fillPromiceAlbedo(aws.albedo, Time, ...
         fillwinter=kwargs.fillwinter);
      albedo_policy = sprintf( ...
         "albedo from PROMICE L3 albedo; fillPromiceAlbedo(fillwinter=%d)", ...
         kwargs.fillwinter);
   elseif has_albedo_source
      met.albedo = nan(height(aws), 1);
      albedo_policy = ...
         "albedo = NaN placeholder (PROMICE L3 albedo is all missing); no observations invented";
   else
      met.albedo = nan(height(aws), 1);
      albedo_policy = ...
         "albedo = NaN placeholder (PROMICE L3 albedo source channel " + ...
         "absent); no observations invented";
   end
   met.wspd = channelOrNan(aws, "wspd");
   met.rh = channelOrNan(aws, "rh");
   met.psfc = channelOrNan(aws, "psfc");
   met.ppt = nan(height(aws), 1);

   % Pass through the optional (non-forcing) diagnostics the station happens to
   % carry (time-varying boom height, wind direction, surface temperature,
   % cloud fraction, turbulent fluxes, ...). These are not required scalar
   % meteorological channels, so a station missing any of them (e.g. KAN_B has
   % no cfrac) must still yield a native met file. Boom height passes
   % through as optional geometry: the runtime A3 fallback chain
   % guarantees usable heights and never blocks a run (POLICY A3).
   % The canonical optional set is the single source of truth for pass-through;
   % rainf is handled explicitly below because it requires a unit conversion.
   [~, optional] = icemodel.forcing.helpers.metvariables();
   for v = setdiff(optional, "rainf", 'stable')
      if ismember(v, string(aws.Properties.VariableNames))
         met.(v) = aws.(v);
      end
   end

   % PROMICE rainfall_cor_u is a corrected LIQUID amount in millimetres per
   % hourly timestep. Convert it at the builder boundary to the canonical
   % water-equivalent rainf rate [m s-1]. The tipping-bucket gauge is not a
   % reliable solid-precipitation measurement, so it must not populate snowf or
   % total ppt. Stations without the source channel retain a NaN placeholder.
   has_rainf_source = ismember("rainf", string(aws.Properties.VariableNames));
   has_rainf_observations = has_rainf_source && any(isfinite(aws.rainf));
   if has_rainf_source
      met.rainf = aws.rainf * 1e-3 / 3600;
   end
   if has_rainf_observations
      rainf_policy = ...
         "rainf converted from rainfall_cor_u [mm per hourly timestep] to [m s-1]";
   elseif has_rainf_source
      rainf_policy = ...
         "rainf = NaN placeholder (rainfall_cor_u is all missing in the requested window)";
   else
      rainf_policy = ...
         "rainf = NaN placeholder (rainfall_cor_u source channel absent)";
   end

   % Physical QA/QC must not replace station outages. The final 15-minute writer
   % may interpolate only inside contiguous finite hourly support.
   [met, checks] = icemodel.forcing.helpers.metchecks(met, fillgaps=false);
   if kwargs.fillwithmissing
      met = icemodel.forcing.helpers.completeMetVariables(met, ...
         include_split_precip=true);
   end

   % Per-variable units from the shared canonical map. rainf carries observed
   % liquid precipitation where available; snowf and total ppt stay all-NaN
   % placeholders unless a runtime source swap fills them.
   met.Properties.VariableUnits = icemodel.forcing.helpers.variableUnits( ...
      string(met.Properties.VariableNames));

   metadata = source_meta;
   metadata.checks = checks;
   % Pin the raw bytes that define builder-derived support masks so later
   % reconstruction cannot replay a changed file with coincidentally equal
   % aggregate counts.
   source_info = dir(string(source_meta.source_file));
   metadata.source_size_bytes = source_info.bytes;
   metadata.source_sha256 = ...
      icemodel.verification.setup.fileSha256(source_meta.source_file);

   % Carry the corrected/raw shortwave selection and exact replacement counts
   % into the saved artifact contract, including the channel-specific swd
   % placeholder policy used by sparse stations.
   fields = fieldnames(shortwave_meta);
   for k = 1:numel(fields)
      metadata.(fields{k}) = shortwave_meta.(fields{k});
   end
   metadata.albedo_source_present = has_albedo_source;
   metadata.albedo_observations_present = has_albedo_observations;
   metadata.albedo_policy = albedo_policy;
   metadata.precip_policy = ...
      "rainfall_cor_u -> rainf [m s-1] when present (corrected liquid " + ...
      "tipping-bucket precipitation; not reliable solid precipitation); " + ...
      "snowf and ppt = NaN placeholders for runtime fill/swap";
   metadata.rainf_source_present = has_rainf_source;
   metadata.rainf_observations_present = has_rainf_observations;
   metadata.rainf_source = ...
      "PROMICE L3 rainfall_cor_u [mm per hourly timestep]";
   metadata.rainf_policy = rainf_policy;
   metadata.humidity_policy = ...
      "rh from the PROMICE v3 rh_wrtwater channel (percent, w.r.t. water)";
   metadata.lwd_estimated = lwd_estimated;
   metadata.lwd_placeholder = lwd_placeholder;
   metadata.lwd_policy = lwd_policy;
   metadata.composing_stations = site_info.stations;
   metadata.station_transition_times = transition_times;
   metadata.station_transition_record = transition_record;
   metadata.gap_policy = ...
      "shortwave missing values become zero only for whole-hour deep civil " + ...
      "night; other source gaps preserved; no metchecks gap interpolation";
   % Record the exact delivered met contract just like every Data-backed met
   % builder while preserving PROMICE's source-specific leg assembly above.
   metadata.met_variables = string(met.Properties.VariableNames);
   metadata.fillwithmissing = kwargs.fillwithmissing;
   metadata = icemodel.forcing.helpers.columnizeMetadata(metadata);
   met.Properties.UserData = metadata;
end

function data = channelOrNan(tt, name)
   %CHANNELORNAN Return a source channel or an all-NaN placeholder.
   if ismember(name, string(tt.Properties.VariableNames))
      data = tt.(name);
   else
      data = nan(height(tt), 1);
   end
end
