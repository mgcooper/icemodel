function [met, metadata] = buildPromiceMet(site, kwargs)
   %BUILDPROMICEMET Build an icemodel met timetable from PROMICE AWS data.
   %
   %  [met, metadata] = icemodel.forcing.buildPromiceMet(site)
   %  [met, metadata] = ... buildPromiceMet(site, source_dir=..., ...
   %     startdate=..., enddate=..., fillwinter=true)
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
   end

   [aws, source_meta] = icemodel.forcing.readPromiceAws(site, ...
      source_dir=kwargs.source_dir, timescale="hourly", ...
      startdate=kwargs.startdate, enddate=kwargs.enddate);

   % Assemble the met variables. PROMICE has no precipitation channel;
   % ppt is explicitly zero (see header note).
   Time = aws.Time;
   met = timetable(Time, ...
      aws.tair, aws.swd, aws.lwd, ...
      icemodel.forcing.fillPromiceAlbedo(aws.albedo, Time, ...
      fillwinter=kwargs.fillwinter), ...
      aws.wspd, aws.rh, aws.psfc, zeros(height(aws), 1), ...
      aws.wdir, aws.tsfc, aws.cfrac, aws.shf, aws.lhf, ...
      'VariableNames', {'tair', 'swd', 'lwd', 'albedo', 'wspd', 'rh', ...
      'psfc', 'ppt', 'wdir', 'tsfc', 'cfrac', 'shf', 'lhf'});

   % Standard QA/QC: gap fill + physical clamps.
   [met, checks] = icemodel.forcing.helpers.metchecks(met);

   met.Properties.VariableUnits = {'K', 'W/m2', 'W/m2', '-', 'm/s', ...
      '%', 'Pa', 'm', 'degrees', 'K', '-', 'W/m2', 'W/m2'};

   metadata = source_meta;
   metadata.checks = checks;
   metadata.albedo_policy = sprintf( ...
      "fillPromiceAlbedo(fillwinter=%d)", kwargs.fillwinter);
   metadata.precip_policy = ...
      "ppt = 0 (PROMICE has no precipitation sensor; swap in a gridded source)";
   metadata.humidity_policy = ...
      "rh from the PROMICE v3 rh_wrtwater channel (percent, w.r.t. water)";
end
