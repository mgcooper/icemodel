function [met, metadata] = buildPromiceMet(site, kwargs)
   %BUILDPROMICEMET Build an icemodel met timetable from PROMICE AWS data.
   %
   %  [met, metadata] = icemodel.forcing.buildPromiceMet(site)
   %  [met, metadata] = ... buildPromiceMet(site, source_dir=..., ...
   %     startdate=..., enddate=..., fillwinter=true, fillwithmissing=true)
   %
   % Reads the station's PROMICE v3 hourly NetCDF (icemodel.forcing.
   % readPromiceAws), applies the albedo winter-fill policy, runs the
   % standard QA/QC pass (gap fill + physical clamps), and returns a
   % timetable satisfying the icemodel met contract for any station and
   % any window in the source record. Save it with
   % icemodel.forcing.helpers.writemet.
   %
   % PROMICE stations carry no precipitation sensor, so ppt is zero with
   % a metadata note; in practice precipitation forcing comes from a
   % gridded source, either as the met file itself or through the
   % met-swap (userdata) mechanism.
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
   %                    instead of rejecting the met build up front
   %
   % Outputs
   %  met      - timetable: tair [K], swd, lwd [W m-2], albedo [-],
   %             wspd [m s-1], rh [%], psfc [Pa], ppt [m] (zero), plus
   %             pass-through wdir, tsfc, cfrac, shf, lhf
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
      kwargs.fillwithmissing (1, 1) logical = false
   end

   [aws, source_meta] = icemodel.forcing.readPromiceAws(site, ...
      source_dir=kwargs.source_dir, timescale="hourly", ...
      startdate=kwargs.startdate, enddate=kwargs.enddate);

   % lwd source. Normally the observed channel. OPT-IN fallback (fill_lwd, off
   % by default) for stations with no longwave sensor (older GC-Net/firn sites,
   % where dlr is absent or all-NaN in the window): estimate lwd from air
   % temperature + vapor pressure via the legacy empirical relation, FLAGGED as
   % estimated in metadata. The strict met contract still requires a real lwd
   % for genuine forcing sources, so by default an absent lwd fails validation.
   has_lwd = ismember('lwd', aws.Properties.VariableNames) ...
      && any(isfinite(aws.lwd));
   lwd_estimated = false;
   if has_lwd
      lwd = aws.lwd;
   elseif kwargs.fill_lwd
      ea = icemodel.vapor.saturation_vapor_pressure(aws.tair, true) ...
         .* aws.rh ./ 100;
      lwd = icemodel.surface.empirical_incoming_longwave_radiation(aws.tair, ea);
      lwd_estimated = true;
   elseif kwargs.fillwithmissing && ismember('lwd', aws.Properties.VariableNames)
      lwd = aws.lwd;                  % present but all-NaN -> explicit placeholder
   elseif kwargs.fillwithmissing
      lwd = nan(height(aws), 1);      % absent -> explicit placeholder
   else
      error('icemodel:forcing:buildPromiceMet:missingForcing', ...
         'required forcing channel(s) absent at %s: lwd (no met built)', site);
   end

   % A station that lacks a required FORCING channel (e.g. LYN_L has no swd
   % sensor) genuinely cannot drive the model. Detect that up front and throw a
   % clean, named error so the caller's met-leg guard records a readable skip
   % reason (the eval observations still stage) instead of a cryptic
   % "Unrecognized table variable name" from the assembly below. lwd is excluded
   % here - it has its own observed/estimated/absent handling above.
   forcing_channels = ["tair", "swd", "albedo", "wspd", "rh", "psfc"];
   absent = forcing_channels(~ismember(forcing_channels, ...
      string(aws.Properties.VariableNames)));
   if ~isempty(absent) && ~kwargs.fillwithmissing
      error('icemodel:forcing:buildPromiceMet:missingForcing', ...
         'required forcing channel(s) absent at %s: %s (no met built)', ...
         site, strjoin(absent, ', '));
   end

   % Assemble the met variables from the required forcings. PROMICE has no
   % precipitation channel; ppt is explicitly zero (see header note). A missing
   % OPTIONAL (non-forcing) channel must not block met building.
   Time = aws.Time;
   met = timetable(Time);
   met.tair = channelOrNan(aws, "tair");
   met.swd = channelOrNan(aws, "swd");
   met.lwd = lwd;
   if ismember("albedo", string(aws.Properties.VariableNames))
      met.albedo = icemodel.forcing.fillPromiceAlbedo(aws.albedo, Time, ...
         fillwinter=kwargs.fillwinter);
   else
      met.albedo = nan(height(aws), 1);
   end
   met.wspd = channelOrNan(aws, "wspd");
   met.rh = channelOrNan(aws, "rh");
   met.psfc = channelOrNan(aws, "psfc");
   met.ppt = zeros(height(aws), 1);

   % Pass through the optional (non-forcing) diagnostics the station happens to
   % carry (wind direction, surface temperature, cloud fraction, turbulent
   % fluxes, ...). These are NOT forcings, so a station missing any of them
   % (e.g. KAN_B has no cfrac) must still yield a met file. The canonical
   % optional set is the single source of truth for what is a pass-through.
   [~, optional] = icemodel.forcing.helpers.metvariables();
   for v = optional
      if ismember(v, string(aws.Properties.VariableNames))
         met.(v) = aws.(v);
      end
   end

   % Standard QA/QC: gap fill + physical clamps.
   [met, checks] = icemodel.forcing.helpers.metchecks(met);
   if kwargs.fillwithmissing
      met = icemodel.forcing.helpers.completeMetVariables(met);
   end

   % Per-variable units from the shared canonical map. ppt is the canonical
   % m s-1 water-equivalent rate (PROMICE has no precipitation sensor, so the
   % channel is zero, but the label stays consistent with the gridded met
   % builders and the metvariables contract).
   met.Properties.VariableUnits = icemodel.forcing.helpers.variableUnits( ...
      string(met.Properties.VariableNames));

   metadata = source_meta;
   metadata.checks = checks;
   metadata.albedo_policy = sprintf( ...
      "fillPromiceAlbedo(fillwinter=%d)", kwargs.fillwinter);
   metadata.precip_policy = ...
      "ppt = 0 (PROMICE has no precipitation sensor; swap in a gridded source)";
   metadata.humidity_policy = ...
      "rh from the PROMICE v3 rh_wrtwater channel (percent, w.r.t. water)";
   metadata.lwd_estimated = lwd_estimated;
   if lwd_estimated
      metadata.lwd_policy = ...
         "lwd ESTIMATED from tair + vapor pressure (no longwave sensor; " + ...
         "icemodel.surface.empirical_incoming_longwave_radiation, fill_lwd=true)";
   else
      metadata.lwd_policy = "lwd from the observed dlr channel";
   end
end

function data = channelOrNan(tt, name)
   %CHANNELORNAN Return a source channel or an all-NaN placeholder.
   if ismember(name, string(tt.Properties.VariableNames))
      data = tt.(name);
   else
      data = nan(height(tt), 1);
   end
end
