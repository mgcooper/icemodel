function out = stageRcmForcing(points, kwargs)
   %STAGERCMFORCING Stage RCM forcing + Data for a list of points.
   %
   %  colocation = icemodel.verification.setup.stageRcmForcing(points, ...
   %     legspec=L, forcing_sources=["mar","merra","racmo"], ...
   %     met_outdir=..., ...
   %     userdata_outdir=..., mar_dir=..., merra_dir=..., racmo_dir=...)
   %
   %  manifest = icemodel.verification.setup.stageRcmForcing( ...
   %     obs_manifest=<path-or-struct>, manifest_file=..., ...
   %     met_outdir=..., userdata_outdir=..., mar_dir=..., racmo_dir=...)
   %
   %  Single owner of RCM forcing/Data generation for the firn-evaluation
   %  staging. Observation import (importPromiceSites / importSumup) DELEGATES
   %  here instead of containing the RCM logic, and this builder is ALSO callable
   %  independently after observations are imported, so RCM products can be
   %  (re)built at any time - including when a new model/forcing collection is
   %  added - without re-importing observations.
   %
   %  WHAT IS WRITTEN (per source)
   %    * MAR / MERRA : the complete Data timetable via writeuserdata
   %      (userdata/<src>/) and the met-contract timetable via data2met+writemet
   %      (met/<src>/). The Data file carries every channel; the met file is the
   %      forcing subset.
   %    * RACMO       : the Data timetable via writeuserdata only (no met). The
   %      available RACMO 2.3p3 subsurface files lack the near-surface met STATE
   %      channels (tair/wspd/rh/psfc), so a RACMO met would fail validatemet;
   %      RACMO is staged as eval/reference Data, not a met source.
   %
   %  EXECUTION + ERROR HANDLING
   %    Each source is staged only for the files it still needs. Existing window
   %    files are reused only when they cover the full requested leg. A partial
   %    overlap triggers a wider rebuild when raw coverage is available; if that
   %    rebuild fails, the clipped cached fallback is retained with a warning and
   %    manifest note rather than being accepted silently. Points that
   %    still need source reads are grouped by identical year set, so one long
   %    station record does not force every shorter station to read the long union
   %    span. Sources are processed in order and each source's files are written
   %    IMMEDIATELY, so a later source's failure never rolls back an earlier
   %    source already on disk (MAR survives a RACMO failure). Manifest mode also
   %    checkpoints each completed source before starting the next one. A
   %    source-group build failure degrades only THAT source's participating legs
   %    to skip-with-reason; the other sources stand. The empty-coverage case is
   %    caught EARLIER (resolveLegWindows marks the leg staged=false before any
   %    build), so a model/site/window with no data is skipped without entering
   %    an hours-long build.
   %
   %  MODES
   %    Explicit (drivers): pass `points` (Nx2 [lat lon]) and `legspec` (a 1xN
   %      struct array; legspec(k).alias plus legspec(k).<src> = the
   %      resolveLegWindows leg for each requested source). Returns a 1xN cell of
   %      per-point colocation structs (colocation{k}.<src> = the leg record).
   %    Manifest convenience (standalone): pass `obs_manifest` (a family-manifest
   %      struct or its manifest.json path). Each case's point/window/alias is
   %      read from site_location/period/case_id, the legs are resolved from a
   %      fresh coverage probe, the explicit path stages the files, the colocation
   %      legs are merged back into every case, and the manifest is rewritten
   %      (manifest_file, or the supplied path). Returns the updated manifest.
   %
   %  Name-value
   %    legspec : struct array aligned to points (explicit mode). Each element has
   %        .alias (string, for file naming) and one field per requested source
   %        holding that source's leg (.staged/.years/.start/.end/.reason).
   %        Optional .discovery_start/.discovery_end broaden only cached-file
   %        search/ranking; .start/.end remain required build/reuse coverage.
   %    forcing_sources : sources requested by this call; subset of
   %        ["mar","merra","racmo"] (default all three). Repeated selectors are
   %        staged once in first-occurrence order. Manifest mode merges those
   %        results into each case and preserves omitted existing legs.
   %    met_outdir, userdata_outdir : output dirs (met_<...> / userdata files).
   %    mar_dir, merra_dir, racmo_dir, modis_dir : raw-source dirs per model.
   %    method : point sampling, "nearest" (default) or "natural".
   %    dt_out : model-met output timestep (default "15m"); pass "" to keep
   %        native cadence. Data/userdata is hourly at its shared writer.
   %    obs_manifest : path or struct selecting manifest-convenience mode.
   %    manifest_file : where to persist in manifest mode (when obs_manifest is a
   %        struct). Defaults to the obs_manifest path when that is a path.
   %    overwrite : logical (default false). Attempt to rebuild requested outputs
   %        even when broader artifacts exist; a compatible cache remains the
   %        fallback if the raw source build fails.
   %    overwrite_family : passed through to writeFamilyManifestMerge (manifest
   %        mode); default false (MERGE).
   %
   % See also: icemodel.verification.setup.resolveLegWindows,
   %  icemodel.verification.setup.rcmSourceCoverage,
   %  icemodel.verification.setup.importPromiceSites,
   %  icemodel.verification.setup.importSumup,
   %  icemodel.forcing.buildMarData, icemodel.forcing.data2met

   arguments
      points double = []
      kwargs.legspec struct = struct([])
      kwargs.forcing_sources (1, :) string ...
         {icemodel.verification.validators.mustBeRcmSourceSelection( ...
         kwargs.forcing_sources)} = ...
         icemodel.verification.namelists.rcmsources()
      kwargs.met_outdir (1, 1) string = ""
      kwargs.userdata_outdir (1, 1) string = ""
      kwargs.mar_dir (1, 1) string = ""
      kwargs.merra_dir (1, 1) string = ""
      kwargs.racmo_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.method (1, 1) string ...
         {mustBeMember(kwargs.method, ["nearest", "natural"])} = "nearest"
      kwargs.dt_out (1, 1) string ...
         {mustBeMember(kwargs.dt_out, ["", "15m"])} = "15m"
      kwargs.obs_manifest = ""
      kwargs.manifest_file (1, 1) string = ""
      kwargs.overwrite (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
   end

   % Remove blank optional selectors before dynamic field dispatch. Stable
   % deduplication then prevents repeated artifact and manifest work.
   sources = reshape(kwargs.forcing_sources, 1, []);
   sources = sources(strlength(strtrim(sources)) > 0);
   sources = unique(sources, 'stable');
   kwargs.modis_dir = defaultModisDir(kwargs.modis_dir);

   % Mode select: a struct or a non-empty path in obs_manifest -> manifest mode.
   om = kwargs.obs_manifest;
   manifest_mode = isstruct(om) || (~isstruct(om) && strlength(string(om)) > 0);

   if manifest_mode
      out = stageFromManifest(sources, kwargs);
   else
      out = stageExplicit(points, kwargs.legspec, sources, kwargs);
   end
end

%% Explicit mode
function modis_dir = defaultModisDir(modis_dir)
   %DEFAULTMODISDIR Use the standard GEUS MODIS cache when the caller is silent.
   if modis_dir ~= ""
      return
   end
   candidates = [
      string(getenv("ICEMODEL_MODIS_DIR"))
      "/Volumes/S03/DATA/greenland/geus/albedo/gris"
      fullfile(icemodel.internal.fullpath("data"), ...
      "verification", "geus", "albedo", "gris")
      fullfile(icemodel.internal.fullpath("data"), ...
      "test", "forcing", "geus", "albedo", "gris")];
   candidates = candidates(strlength(candidates) > 0);
   found = candidates(isfolder(candidates));
   if isempty(found)
      modis_dir = "";
   else
      modis_dir = found(1);
   end
end

function colocation = stageExplicit(points, legspec, sources, kwargs)
   %STAGEEXPLICIT Stage each requested source for the point list; return legs.
   n = size(points, 1);
   if numel(legspec) ~= n
      error('icemodel:verification:stageRcmForcing:legspecMismatch', ...
         'legspec (%d) must align with points (%d)', numel(legspec), n)
   end

   % colocation as a 1xN cell of scalar structs: per-source leg records have
   % heterogeneous fields (staged met/data legs vs skipped legs), so a cell
   % avoids forcing a single struct-array field union.
   colocation = repmat({struct()}, 1, n);

   % One source at a time, written immediately, so completed sources persist if
   % a later source fails.
   for src = sources
      colocation = stageOneSource(src, points, legspec, colocation, kwargs);
   end
end

function colocation = stageOneSource(src, points, legspec, colocation, kwargs)
   %STAGEONESOURCE Build one source by needed year groups, write per point.
   srcc = char(src);
   kind = rcmLegKind(src);
   n = size(points, 1);

   % Participants: points whose leg for this source resolved staged=true and
   % whose required output files are not already fully covered by an existing
   % staged window file. Cached-file discovery runs before the raw-coverage
   % gate so metadata-only RCM updates can reattach expensive artifacts even when
   % the raw archive probe is unavailable or outside the current source window.
   part = false(1, n);
   existing = repmat(emptyExistingFiles, 1, n);
   for k = 1:n
      L = legspec(k).(srcc);
      if hasDiscoveryWindow(L)
         existing(k) = existingRcmFiles(src, legspec(k).alias, L, ...
            points(k, :), kwargs);
      end
      if kwargs.overwrite && L.staged
         % Explicit overwrite forces a raw rebuild attempt, but a compatible
         % cache remains available to the source-failure fallback. Conflict
         % records have no reusable files; retain their reason so a concrete
         % mismatch cannot be hidden by a later raw-source failure.
         part(k) = true;
      elseif existing(k).conflict && L.staged
         % Any incompatible cache must yield to a requested raw build. If the
         % source is unavailable, the existing failure boundary records the
         % conflict without reusing or overwriting the incompatible artifact.
         part(k) = true;
      elseif existing(k).conflict
         colocation{k}.(srcc) = skippedLeg(kind, existing(k).reason);
      elseif completeExistingLegFromData(src, existing(k)) ...
            && (~L.staged || existingDataCoversLeg(existing(k), L))
         try
            existing(k) = deriveExistingMetFromData(src, legspec(k).alias, ...
               points(k, :), kwargs, existing(k));
            colocation{k}.(srcc) = existingRcmLeg(src, kind, L, ...
               existing(k), kwargs);
            warnExistingWindowReuse(src, legspec(k).alias, L, existing(k));
         catch derive_err
            if L.staged
               existing(k) = clearExistingData(existing(k));
               part(k) = true;
            else
               colocation{k}.(srcc) = skippedLeg(kind, derive_err.message);
            end
         end
      elseif existingLegComplete(src, existing(k)) ...
            && (~L.staged || (existingLegMatchesFullStage(src, existing(k)) ...
            && existingLegCoversRequestedWindow(src, existing(k), L)))
         colocation{k}.(srcc) = existingRcmLeg(src, kind, L, existing(k), kwargs);
         warnExistingWindowReuse(src, legspec(k).alias, L, existing(k));
      elseif L.staged
         part(k) = true;
      else
         colocation{k}.(srcc) = skippedLeg(kind, L.reason);
      end
   end
   idx = find(part);
   if isempty(idx)
      return
   end

   % Build each distinct year-set separately: this keeps the batch-reader
   % benefit for matching points without over-reading short-window points.
   groups = sameYearGroups(idx, legspec, srcc);
   for g = 1:numel(groups)
      gidx = groups(g).idx;
      pts = zeros(numel(gidx), 2);
      for j = 1:numel(gidx)
         pts(j, :) = points(gidx(j), :);
      end

      try
         Data = buildSourceData(src, pts, groups(g).years, kwargs);
      catch build_err
         % The source builder is the isolation boundary: missing, corrupt, or
         % otherwise unreadable raw inputs degrade this source/year group only.
         % Preserve any compatible met-only cache as a selectable forcing leg;
         % sources already written by this call remain untouched.
         for j = 1:numel(gidx)
            k = gidx(j);
            if endsWith(string(build_err.identifier), ...
                  ":pointOutsideValidDomain")
               % Native-mask rejection proves that the prior spatial payload is
               % invalid. Make the skipped result destructive so additive
               % manifest merging cannot resurrect an old off-mask artifact.
               leg = skippedLeg(kind, build_err.message);
               leg.replace_prior_artifacts = true;
               colocation{k}.(srcc) = leg;
            else
               colocation{k}.(srcc) = fallbackExistingOrSkipped( ...
                  src, kind, legspec(k).(srcc), existing(k), kwargs, ...
                  build_err.message);
            end
         end
         continue
      end

      % Per point: clip to its own window and write only the outputs that were
      % not already covered. A write failure degrades only that point's leg.
      for j = 1:numel(gidx)
          k = gidx(j);
          L = legspec(k).(srcc);
          try
            % MAR builders stamp per-day provenance on the whole source-year
            % axis. Preserve that exact axis before the case window removes rows.
            if src == "mar"
               source_days = unique(dateshift(Data{j}.Time(:), ...
                  'start', 'day'), 'stable');
            end
            d = windowSubset(Data{j}, L.start, L.end);
            if src == "mar"
               % Align once on hourly Data; the existing Data-to-met/writer path
               % then carries the same ledger into the derived 15-minute artifact.
               d.Properties.UserData = ...
                  icemodel.forcing.helpers.alignMarDailyMetadata( ...
                  d.Properties.UserData, source_days, d.Time);
            elseif src == "merra"
               % The builder proof describes its full source-year axis. Narrow it
               % to this case window before Data-to-met carries the same exact
               % native tavg3 inventory onto the 15-minute artifact.
               d.Properties.UserData = alignMerraTavg3Metadata( ...
                  d.Properties.UserData, d.Time);
            end
            write_existing = existing(k);
            if kwargs.overwrite
               % The discovered cache is failure fallback only. A successful
               % explicit overwrite must write the freshly built artifacts.
               write_existing = emptyExistingFiles;
            end
            colocation{k}.(srcc) = writeRcmLeg( ...
               src, d, legspec(k).alias, kind, L, points(k, :), ...
               kwargs, write_existing);
         catch write_err
            if ~isSkippableRcmStagingError(write_err)
               rethrow(write_err)
            end
            if endsWith(string(write_err.identifier), ...
                  ":pointOutsideValidDomain")
               % A non-local ice cell is a proven colocation conflict, not a
               % transient source failure. Replace any cached off-mask leg with
               % an explicit unavailable result instead of resurrecting it.
               leg = skippedLeg(kind, write_err.message);
               leg.replace_prior_artifacts = true;
               colocation{k}.(srcc) = leg;
            else
               colocation{k}.(srcc) = fallbackExistingOrSkipped( ...
                  src, kind, L, existing(k), kwargs, write_err.message);
            end
         end
      end
   end
end

function tf = isSkippableRcmStagingError(err)
   %ISSKIPPABLERCMSTAGINGERROR True only for absent source/window gaps.
   id = string(err.identifier);
   tf = endsWith(id, ":sourceNotFound") ...
      || endsWith(id, ":fileNotFound") ...
      || endsWith(id, ":emptyWindow") ...
      || endsWith(id, ":yearNotInArchive") ...
      || endsWith(id, ":pointOutsideValidDomain");
end

function existing = clearExistingData(existing)
   %CLEAREXISTINGDATA Force raw rebuild after cached Data proved unusable.
   existing.data_files = strings(1, 0);
   existing.data_start = NaT;
   existing.data_end = NaT;
end

function leg = fallbackExistingOrSkipped(src, kind, L, existing, kwargs, reason)
   %FALLBACKEXISTINGORSKIPPED Keep reusable cache when a raw refresh fails.
   if existingLegComplete(src, existing)
      leg = existingRcmLeg(src, kind, L, existing, kwargs);
      leg.note = sprintf('%s Raw refresh failed; preserved cached artifact: %s', ...
         char(string(leg.note)), char(string(reason)));
      [required_start, required_end] = requiredWindow(L);
      warning('icemodel:verification:stageRcmForcing:existingWindowFile', ...
         ['%s refresh failed for requested window %s to %s; preserving ' ...
         'the existing partial cached artifact: %s'], upper(char(src)), ...
         string(required_start), string(required_end), string(reason));
   elseif existing.conflict && contains(string(existing.reason), ...
         "artifact metadata does not match")
      % A proven point/method mismatch must remain destructive even when the
      % attempted overwrite later fails to read the raw source.
      leg = skippedLeg(kind, sprintf('%s; raw refresh failed: %s', ...
         char(string(existing.reason)), char(string(reason))));
   else
      leg = skippedLeg(kind, reason);
   end
end

function Data = buildSourceData(src, pts, years, kwargs)
   %BUILDSOURCEDATA Build the FULL Data timetable(s) for one source's points.
   % All three builders accept an Nx2 point list and return a 1xN cell of Data
   % timetables (one file open per source-year serving all points).
   switch char(src)
      case "mar"
         Data = asCell(icemodel.forcing.buildMarData(pts, years, ...
            source_dir=kwargs.mar_dir, modis_dir=kwargs.modis_dir, ...
            method=kwargs.method, fillgaps=false));
      case "merra"
         Data = asCell(icemodel.forcing.buildMerraData(pts, years, ...
            source_dir=kwargs.merra_dir, modis_dir=kwargs.modis_dir, ...
            method=kwargs.method, fillgaps=false));
      case "racmo"
         Data = asCell(icemodel.forcing.buildRacmoData(pts, years, ...
            source_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir, ...
            method=kwargs.method, dt="1hr"));
   end
end

function co = writeRcmLeg(src, d, alias, kind, L, point, kwargs, existing)
   %WRITERCMLEG Write one point's staged leg + build its colocation record.
   % MAR/MERRA write BOTH the full Data (userdata) AND the met (data2met); RACMO
   % writes Data only. Existing files are reused only when they cover the full
   % requested leg; partial files remain on disk while a wider artifact is built.
   % Every staged leg carries staged==true.
   if nargin < 8
      existing = emptyExistingFiles;
   end
   met_outdir = metOutdir(kwargs);
   userdata_outdir = userdataOutdir(kwargs);
   storage = rcmStorageToken(src);
   product_id = rcmProductId(src);
   if src == "racmo"
      d = tagSampleMethod(d, kwargs.method, point);
      if existingDataCoversLeg(existing, L)
         data_files = existing.data_files;
      else
         data_files = icemodel.forcing.helpers.writeuserdata( ...
            d, alias, storage, outdir=userdata_outdir, naming="window", ...
            overwrite=true);
         existing.met_files = strings(1, 0);
      end
      co = struct('kind', kind, 'staged', true, ...
         'source', char(src), ...
         'source_id', char(product_id), ...
         'data_files', ...
         icemodel.verification.setup.relpaths(data_files, userdata_outdir), ...
         'sample_method', char(kwargs.method), ...
         'window', icemodel.verification.setup.manifestWindow( ...
         L.start, L.end), ...
         'note', 'SMB/eval Data only; RACMO is not a met source.');
   else
      % MAR/MERRA: the COMPLETE Data timetable AND the met forcing it derives.
      d = tagSampleMethod(d, kwargs.method, point);
      if existingDataCoversLeg(existing, L)
         data_files = existing.data_files;
      else
         data_files = icemodel.forcing.helpers.writeuserdata( ...
            d, alias, storage, outdir=userdata_outdir, naming="window", ...
            overwrite=true);
      end
      if existingMetCoversLeg(existing, L)
         met_files = existing.met_files;
      else
         met = toMet(d);
         met = tagSampleMethod(met, kwargs.method, point);
         met_files = icemodel.forcing.helpers.writemet( ...
            met, alias, storage, outdir=met_outdir, naming="window", ...
            dt_out=kwargs.dt_out, overwrite=true);
      end
      co = struct('kind', kind, 'staged', true, ...
         'source', char(src), ...
         'source_id', char(product_id), ...
         'met_files', ...
         icemodel.verification.setup.relpaths(met_files, met_outdir), ...
         'data_files', ...
         icemodel.verification.setup.relpaths(data_files, userdata_outdir), ...
         'sample_method', char(kwargs.method), ...
         'window', icemodel.verification.setup.manifestWindow( ...
         L.start, L.end));
      % Readiness describes the exact saved path selected by the writer, which
      % may be a fresh file, an exact reuse, or a broader enclosing reuse.
      [forcing_ready, forcing_ready_reason, forcing_complete_windows] = ...
         icemodel.verification.setup.metArtifactReadiness(string(met_files));
      co.forcing_ready = logical(forcing_ready);
      co.forcing_ready_reason = char(forcing_ready_reason);
      co.forcing_complete_windows = forcing_complete_windows;
   end
end

function tf = existingMetCoversLeg(existing, L)
   %EXISTINGMETCOVERSLEG True when cached met spans the new raw-staged leg.
   [requested_start, requested_end] = requiredWindow(L);
   tf = ~isempty(existing.met_files) ...
      && ~isnat(existing.met_start) && ~isnat(existing.met_end) ...
      && existing.met_start <= requested_start ...
      && existing.met_end >= requested_end;
end

function tf = existingDataCoversLeg(existing, L)
   %EXISTINGDATACOVERSLEG True when cached Data spans the requested raw leg.
   [requested_start, requested_end] = requiredWindow(L);
   tf = ~isempty(existing.data_files) ...
      && ~isnat(existing.data_start) && ~isnat(existing.data_end) ...
      && existing.data_start <= requested_start ...
      && existing.data_end >= requested_end;
end

function tf = existingLegCoversRequestedWindow(src, existing, L)
   %EXISTINGLEGCOVERSREQUESTEDWINDOW Require full coverage for requested reuse.
   if src == "racmo"
      tf = existingDataCoversLeg(existing, L);
   else
      tf = existingMetCoversLeg(existing, L) ...
         && existingDataCoversLeg(existing, L);
   end
end

function warnExistingWindowReuse(src, alias, L, existing)
   %WARNEXISTINGWINDOWREUSE Report complete reuse or visible partial fallback.
   [required_start, required_end] = requiredWindow(L);
   if existingLegCoversRequestedWindow(src, existing, L)
      detail = "fully cover";
   else
      detail = "only partially cover";
   end
   warning('icemodel:verification:stageRcmForcing:existingWindowFile', ...
      ['Skipping %s build for %s because existing staged file(s) %s ' ...
      'the requested window %s to %s.'], upper(char(src)), alias, detail, ...
      string(required_start), string(required_end));
end

function met = toMet(d)
   %TOMET Convert a Data timetable to the native-cadence met contract.
   met = icemodel.forcing.data2met(d);
end

function metadata = alignMerraTavg3Metadata(metadata, times)
   %ALIGNMERRATAVG3METADATA Restrict the proven native grid to a saved window.
   expected = times(minute(times) == 0 & second(times) == 0 ...
      & mod(hour(times), 3) == 0);
   missing = metadata.merra_tavg3_missing_source_times(:);
   missing = missing(ismember(missing, expected));
   metadata.merra_tavg3_expected_source_row_count = numel(expected);
   metadata.merra_tavg3_source_time_gap_count = numel(missing);
   metadata.merra_tavg3_source_row_count = numel(expected) - numel(missing);
   metadata.merra_tavg3_missing_source_times = missing;
end

function T = tagSampleMethod(T, method, point)
   %TAGSAMPLEMETHOD Record reusable RCM sampling metadata in saved artifacts.
   metadata = T.Properties.UserData;
   if isempty(metadata) || ~isstruct(metadata)
      metadata = struct();
   end
   metadata.sample_method = char(method);
   metadata.lat_wgs84 = point(1);
   metadata.lon_wgs84 = point(2);
   custom = T.Properties.CustomProperties;
   if isprop(custom, 'Lat') && isprop(custom, 'Lon')
      metadata.source_lat_wgs84 = custom.Lat;
      metadata.source_lon_wgs84 = custom.Lon;
   end
   T.Properties.UserData = icemodel.forcing.helpers.columnizeMetadata(metadata);
end

%% Manifest-convenience mode
function manifest = stageFromManifest(sources, kwargs)
   %STAGEFROMMANIFEST Re-stage RCM forcing for every case in a staged manifest.
   % Reads each case's point (site_location), window (period) and alias
   % (case_id), resolves the legs from a fresh coverage probe, then stages and
   % checkpoints one source at a time through the shared explicit path.
   loaded_from_file = ~isstruct(kwargs.obs_manifest);
   if ~loaded_from_file
      manifest = kwargs.obs_manifest;
      manifest_file = kwargs.manifest_file;
   else
      manifest_file = string(kwargs.obs_manifest);
      manifest = jsondecode(fileread(manifest_file));
   end
   if strlength(string(manifest_file)) == 0
      error('icemodel:verification:stageRcmForcing:noManifestFile', ...
         'manifest_file is required when obs_manifest is a struct')
   end
   kwargs = inferManifestModeOutputDirs(kwargs, manifest_file);

   cases = manifest.cases;
   if isempty(cases)
      return
   end
   n = numel(cases);

   % Probe on-disk coverage once for the whole set (cheap, fail-early gate).
   coverage = icemodel.verification.setup.rcmSourceCoverage(sources, ...
      struct('mar', kwargs.mar_dir, 'merra', kwargs.merra_dir, ...
      'racmo', kwargs.racmo_dir));

   % Build points + legspec from the manifest cases.
   points = zeros(n, 2);
   legspec = repmat(emptyLeg(sources), 1, n);
   for k = 1:n
      c = cases(k);
      points(k, :) = [c.site_location.lat_wgs84, c.site_location.lon_wgs84];
      [t1, t2] = casePeriod(c);
      leg = icemodel.verification.setup.resolveLegWindows( ...
         sources, coverage, t1, t2);
      legspec(k).alias = string(c.case_id);
      for src = sources
         legspec(k).(char(src)) = ...
            legWithDiscoveryWindow(leg.(char(src)), t1, t2);
      end
   end

   % Persist one completed source at a time. This makes the manifest match the
   % already-written artifact state even when a later source raises unexpectedly.
   % Keep one no-source pass because callers also use this mode to persist edits
   % while preserving every omitted existing RCM leg.
   if ~isfield(manifest, 'skipped') || isempty(manifest.skipped)
      manifest.skipped = struct('site', {}, 'reason', {});
   end
   n_passes = max(1, numel(sources));
   ids = string({cases.case_id});
   for pass = 1:n_passes
      selected_sources = strings(1, 0);
      colocation = repmat({struct()}, 1, n);
      if ~isempty(sources)
         selected_sources = sources(pass);
         colocation = stageExplicit( ...
            points, legspec, selected_sources, kwargs);
      end

      % Merge only this source into the in-memory requested cases before its
      % checkpoint. Later passes therefore include every earlier completed leg.
      for k = 1:n
         colocation{k} = icemodel.verification.setup.preserveRcmLegs( ...
            cases(k).colocation, colocation{k}, selected_sources, ...
            cases(k).period, met_outdir=metOutdir(kwargs), ...
            userdata_outdir=userdataOutdir(kwargs), ...
            method=kwargs.method, point=points(k, :), ...
            overwrite=kwargs.overwrite);
         cases(k).colocation = icemodel.verification.setup.mergeColocation( ...
            cases(k).colocation, colocation{k});
         [forcing_sources, eval_sources] = ...
            icemodel.verification.setup.colocationSourceLists( ...
            cases(k).colocation);
         cases(k).forcing_sources = cellstr(forcing_sources);
         cases(k).eval_sources = cellstr(eval_sources);
      end
      manifest.cases = cases;

      % Path mode loaded these unrelated skips from the destination; later
      % checkpoints likewise received them from the prior merge. Pass them only
      % when they are genuinely incoming so ordinary refreshes do not duplicate.
      incoming = manifest;
      if ~kwargs.overwrite_family && (loaded_from_file || pass > 1)
         incoming.skipped = struct('site', {}, 'reason', {});
      end
      manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
         manifest_file, incoming, requested_ids=ids, ...
         overwrite_family=kwargs.overwrite_family);
   end

end

function [t1, t2] = casePeriod(c)
   %CASEPERIOD Parse a case's period{start,end} to UTC datetimes (NaT if blank).
   t1 = parseStamp(c.period.start);
   t2 = parseStamp(c.period.end);
end

function t = parseStamp(s)
   %PARSESTAMP "" / missing -> NaT (unbounded); else a UTC datetime.
   if isempty(s) || strlength(string(s)) == 0
      t = NaT;
   else
      t = icemodel.verification.setup.ensureUtc(string(s));
   end
end

function L = emptyLeg(sources)
   %EMPTYLEG Prototype legspec element (alias + one leg field per source).
   L = struct('alias', "");
   proto = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
      'reason', "");
   for src = sources
      L.(char(src)) = proto;
   end
end

%% Shared local helpers
function results = asCell(out)
   %ASCELL Normalize a builder output to a cell array (single point -> 1x1).
   if iscell(out)
      results = out;
   else
      results = {out};
   end
end

function groups = sameYearGroups(idx, legspec, srcc)
   %SAMEYEARGROUPS Partition participant indices by identical source years.
   groups = repmat(struct('idx', [], 'years', []), 1, numel(idx));
   ngroups = 0;
   for j = 1:numel(idx)
      k = idx(j);
      years = legspec(k).(srcc).years;
      found = false;
      for g = 1:ngroups
         if isequal(groups(g).years, years)
            groups(g).idx(end + 1) = k;
            found = true;
            break
         end
      end
      if ~found
         ngroups = ngroups + 1;
         groups(ngroups).idx = k;
         groups(ngroups).years = years;
      end
   end
   groups = groups(1:ngroups);
end

function existing = emptyExistingFiles
   %EMPTYEXISTINGFILES Prototype for optional pre-staged output paths.
   existing = struct('met_files', strings(1, 0), ...
      'met_start', NaT, 'met_end', NaT, ...
      'data_files', strings(1, 0), ...
      'data_start', NaT, 'data_end', NaT, ...
      'conflict', false, 'reason', "");
end

function existing = existingRcmFiles(src, alias, L, point, kwargs)
   %EXISTINGRCMFILES Locate staged files overlapping the discovery window.
   existing = emptyExistingFiles;
   storage = rcmStorageToken(src);
   [discovery_start, discovery_end] = discoveryWindow(L);
   [required_start, required_end] = requiredWindow(L);
   data_candidates = findExistingWindowRecords(userdataOutdir(kwargs), ...
      storage, sprintf('%s_%s', alias, storage), ".mat", ...
      discovery_start, discovery_end);
   data_records = reusableWindowRecords( ...
      data_candidates, src, kwargs.method, point);
   met_candidates = emptyWindowRecord;
   met_candidates = met_candidates([]);
   met_records = met_candidates;
   if src ~= "racmo"
      met_candidates = findExistingWindowRecords(metOutdir(kwargs), ...
         storage, sprintf('met_%s_%s', alias, storage), ...
         metWindowSuffix(kwargs.dt_out), discovery_start, discovery_end);
      met_records = reusableWindowRecords( ...
         met_candidates, src, kwargs.method, point);
   end

   if isempty(data_records) && isempty(met_records) ...
         && (~isempty(data_candidates) || ~isempty(met_candidates))
      existing.conflict = true;
      if allMissingArtifactMetadata([data_candidates; met_candidates])
         existing.reason = sprintf( ...
            '%s cache exists for %s but artifact metadata is missing', ...
            upper(char(src)), alias);
      else
         existing.reason = sprintf( ...
            '%s cache exists for %s but artifact metadata does not match requested point/method', ...
            upper(char(src)), alias);
      end
      return
   end

   if src == "racmo"
      existing = setExistingData(existing, ...
         bestWindowRecord(data_records, discovery_start, discovery_end, ...
         required_start, required_end));
   else
      [data_record, met_record] = bestMetDataWindowPair( ...
         data_records, met_records, discovery_start, discovery_end, ...
         required_start, required_end);
      existing = setExistingData(existing, data_record);
      existing = setExistingMet(existing, met_record);
   end
end

function tf = existingLegComplete(src, existing)
   %EXISTINGLEGCOMPLETE True when all outputs required for SRC already exist.
   if src == "racmo"
      tf = ~isempty(existing.data_files);
   else
      tf = ~isempty(existing.met_files);
   end
end

function tf = existingLegMatchesFullStage(src, existing)
   %EXISTINGLEGMATCHESFULLSTAGE True when cache equals a clean full RCM stage.
   if src == "racmo"
      tf = ~isempty(existing.data_files);
   else
      tf = ~isempty(existing.met_files) && ~isempty(existing.data_files);
   end
end

function tf = completeExistingLegFromData(src, existing)
   %COMPLETEEXISTINGLEGFROMDATA True when met can be derived from cached Data.
   tf = src ~= "racmo" && ~isempty(existing.data_files) ...
      && isempty(existing.met_files);
end

function existing = deriveExistingMetFromData(src, alias, point, kwargs, existing)
   %DERIVEEXISTINGMETFROMDATA Rebuild a missing MAR/MERRA met from cached Data.
   S = load(existing.data_files(1), 'Data');
   d = tagSampleMethod(S.Data, kwargs.method, point);
   met = toMet(d);
   met = tagSampleMethod(met, kwargs.method, point);
   existing.met_files = icemodel.forcing.helpers.writemet( ...
      met, alias, rcmStorageToken(src), outdir=metOutdir(kwargs), ...
      naming="window", dt_out=kwargs.dt_out, overwrite=true);
   existing.met_start = existing.data_start;
   existing.met_end = existing.data_end;
end

function co = existingRcmLeg(src, kind, L, existing, kwargs)
   %EXISTINGRCMLEG Manifest leg pointing at already-staged overlapping files.
   product_id = rcmProductId(src);
   [window_start, window_end] = existingOverlapWindow(src, L, existing);
   if existingLegCoversRequestedWindow(src, existing, L)
      data_note = 'Existing staged Data file fully covers requested window.';
      met_data_note = ...
         'Existing staged met/Data files fully cover requested window.';
   else
      data_note = ['Existing staged Data file overlaps but does not fully ' ...
         'cover requested window/output requirements.'];
      met_data_note = ['Existing staged met/Data files overlap but do not ' ...
         'fully cover requested window/output requirements.'];
   end
   if src == "racmo"
      co = struct('kind', kind, 'staged', true, ...
         'source', char(src), ...
         'source_id', char(product_id), ...
         'data_files', ...
         icemodel.verification.setup.relpaths( ...
         existing.data_files, userdataOutdir(kwargs)), ...
         'sample_method', char(kwargs.method), ...
         'window', icemodel.verification.setup.manifestWindow( ...
         window_start, window_end), ...
         'note', data_note);
   else
      % The reusable met file is the forcing artifact selected by this leg;
      % diagnose its bytes rather than inferring readiness from the request.
      [forcing_ready, forcing_ready_reason, forcing_complete_windows] = ...
         icemodel.verification.setup.metArtifactReadiness( ...
         string(existing.met_files));
      co = struct('kind', kind, 'staged', true, ...
         'source', char(src), ...
         'source_id', char(product_id), ...
         'sample_method', char(kwargs.method), ...
         'window', icemodel.verification.setup.manifestWindow( ...
         window_start, window_end), ...
         'note', met_data_note);
      co.forcing_ready = logical(forcing_ready);
      co.forcing_ready_reason = char(forcing_ready_reason);
      co.forcing_complete_windows = forcing_complete_windows;
      if ~isempty(existing.met_files)
         co.met_files = icemodel.verification.setup.relpaths( ...
            existing.met_files, metOutdir(kwargs));
      end
      if ~isempty(existing.data_files)
         co.data_files = icemodel.verification.setup.relpaths( ...
            existing.data_files, userdataOutdir(kwargs));
      else
         co.note = ['Existing staged met file overlaps requested window; ' ...
            'no matching Data file was found.'];
      end
   end
end

function records = findExistingWindowRecords( ...
      base, source, prefix, suffix, t1, t2)
   %FINDEXISTINGWINDOWRECORDS Staged files whose encoded windows overlap T1-T2.
   records = emptyWindowRecord;
   records = records([]);
   for d = icemodel.forcing.helpers.sourceSearchDirs(base, source)
      folder = string(d{1});
      if ~isfolder(folder)
         continue
      end
      listing = dir(fullfile(folder, char(string(prefix) + "_*_*" ...
         + string(suffix))));
      for k = 1:numel(listing)
         [ok, candidate_start, candidate_end] = parseWindowFilename( ...
            listing(k).name, prefix, suffix);
         if ~ok || ~windowsOverlap(candidate_start, candidate_end, t1, t2)
            continue
         end
         records(end + 1, 1) = struct( ...
            'filename', string(fullfile(listing(k).folder, listing(k).name)), ...
            'start', candidate_start, 'end', candidate_end); %#ok<AGROW>
      end
   end
end

function record = emptyWindowRecord()
   %EMPTYWINDOWRECORD Prototype for one cached window file.
   record = struct('filename', "", 'start', NaT, 'end', NaT);
end

function records = reusableWindowRecords(records, src, method, point)
   %REUSABLEWINDOWRECORDS Keep records with matching artifact metadata.
   keep = false(numel(records), 1);
   for k = 1:numel(records)
      keep(k) = isReusableRcmFile(records(k).filename, src, method, point);
   end
   records = records(keep);
end

function tf = allMissingArtifactMetadata(records)
   %ALLMISSINGARTIFACTMETADATA True when candidates predate metadata stamps.
   if isempty(records)
      tf = false;
      return
   end
   tf = true;
   for k = 1:numel(records)
      tf = tf && ~hasArtifactMetadata(records(k).filename);
   end
end

function tf = hasArtifactMetadata(filename)
   %HASARTIFACTMETADATA True when a MAT file carries artifact_metadata.
   tf = false;
   try
      info = whos('-file', filename);
   catch
      return
   end
   tf = any(string({info.name}) == "artifact_metadata");
end

function [data_record, met_record] = bestMetDataWindowPair( ...
      data_records, met_records, discovery_start, discovery_end, ...
      required_start, required_end)
   %BESTMETDATAWINDOWPAIR Prefer required-covering, then discovery-ranked pairs.
   data_record = emptyWindowRecord;
   met_record = emptyWindowRecord;
   best_duration = -Inf;
   best_width = -Inf;
   best_end = NaT;
   best_name = "";
   best_covers_required = false;
   for i = 1:numel(data_records)
      for j = 1:numel(met_records)
         [t1, t2] = commonOverlap(data_records(i), met_records(j), ...
            discovery_start, discovery_end);
         if t1 > t2
            continue
         end
         duration = seconds(t2 - t1);
         pair_start = max(data_records(i).start, met_records(j).start);
         pair_end = min(data_records(i).end, met_records(j).end);
         width = seconds(pair_end - pair_start);
         pair_name = data_records(i).filename + "|" ...
            + met_records(j).filename;
         covers_required = pair_start <= required_start ...
            && pair_end >= required_end;
         if betterCandidateRank(covers_required, duration, width, ...
               pair_end, pair_name, best_covers_required, best_duration, ...
               best_width, best_end, best_name)
            data_record = data_records(i);
            met_record = met_records(j);
            best_covers_required = covers_required;
            best_duration = duration;
            best_width = width;
            best_end = pair_end;
            best_name = pair_name;
         end
      end
   end
   if best_duration > -Inf
      return
   end
   if ~isempty(met_records)
      met_record = bestWindowRecord( ...
         met_records, discovery_start, discovery_end, ...
         required_start, required_end);
   else
      data_record = bestWindowRecord( ...
         data_records, discovery_start, discovery_end, ...
         required_start, required_end);
   end
end

function record = bestWindowRecord( ...
      records, discovery_start, discovery_end, required_start, required_end)
   %BESTWINDOWRECORD Prefer required coverage before discovery-window rank.
   record = emptyWindowRecord;
   best_duration = -Inf;
   best_width = -Inf;
   best_end = NaT;
   best_name = "";
   best_covers_required = false;
   for k = 1:numel(records)
      duration = overlapDuration(records(k).start, records(k).end, ...
         discovery_start, discovery_end);
      width = seconds(records(k).end - records(k).start);
      covers_required = records(k).start <= required_start ...
         && records(k).end >= required_end;
      if betterCandidateRank(covers_required, duration, width, ...
            records(k).end, records(k).filename, best_covers_required, ...
            best_duration, best_width, best_end, best_name)
         record = records(k);
         best_covers_required = covers_required;
         best_duration = duration;
         best_width = width;
         best_end = records(k).end;
         best_name = records(k).filename;
      end
   end
end

function tf = betterCandidateRank(covers_required, overlap, width, ...
      window_end, filename, best_covers_required, best_overlap, ...
      best_width, best_end, best_filename)
   %BETTERCANDIDATERANK Prefer sufficiency before broad discovery ranking.
   if covers_required ~= best_covers_required
      tf = covers_required;
      return
   end
   tf = betterWindowRank(overlap, width, window_end, filename, ...
      best_overlap, best_width, best_end, best_filename);
end

function tf = betterWindowRank(overlap, width, window_end, filename, ...
      best_overlap, best_width, best_end, best_filename)
   %BETTERWINDOWRANK Prefer overlap, widest coverage, latest end, then name.
   tf = overlap > best_overlap;
   if overlap ~= best_overlap
      return
   end
   tf = width > best_width;
   if width ~= best_width
      return
   end
   tf = isnat(best_end) || window_end > best_end;
   if isnat(best_end) || window_end ~= best_end
      return
   end
   ordered = sort([string(filename), string(best_filename)]);
   tf = best_filename == "" || ordered(1) == string(filename);
end

function [t1, t2] = commonOverlap( ...
      record_a, record_b, discovery_start, discovery_end)
   %COMMONOVERLAP Intersection of two records and the discovery window.
   t1 = max([record_a.start, record_b.start, discovery_start]);
   t2 = min([record_a.end, record_b.end, discovery_end]);
end

function existing = setExistingData(existing, record)
   %SETEXISTINGDATA Copy a selected data record into the existing-files struct.
   if strlength(record.filename) == 0
      return
   end
   existing.data_files = record.filename;
   existing.data_start = record.start;
   existing.data_end = record.end;
end

function existing = setExistingMet(existing, record)
   %SETEXISTINGMET Copy a selected met record into the existing-files struct.
   if strlength(record.filename) == 0
      return
   end
   existing.met_files = record.filename;
   existing.met_start = record.start;
   existing.met_end = record.end;
end

function [ok, t1, t2] = parseWindowFilename(name, prefix, suffix)
   %PARSEWINDOWFILENAME Parse <prefix>_<YYYYMMDD>_<YYYYMMDD><suffix>.
   ok = false;
   [t1, t2] = deal(NaT);
   pat = ['^' regexptranslate('escape', char(prefix)) ...
      '_(\d{8})_(\d{8})' regexptranslate('escape', char(suffix)) '$'];
   tok = regexp(name, pat, 'tokens', 'once');
   if isempty(tok)
      return
   end
   t1 = datetime(tok{1}, 'InputFormat', 'yyyyMMdd', 'TimeZone', 'UTC');
   t2 = datetime(tok{2}, 'InputFormat', 'yyyyMMdd', 'TimeZone', 'UTC') ...
      + hours(23) + minutes(59) + seconds(59);
   ok = true;
end

function tf = windowsOverlap(a_start, a_end, b_start, b_end)
   %WINDOWSOVERLAP True when two bounded windows share at least one instant.
   tf = ~isnat(a_start) && ~isnat(a_end) ...
      && ~isnat(b_start) && ~isnat(b_end) ...
      && a_start <= b_end && a_end >= b_start;
end

function duration = overlapDuration(a_start, a_end, b_start, b_end)
   %OVERLAPDURATION Duration of the common portion of two overlapping windows.
   duration = seconds(min(a_end, b_end) - max(a_start, b_start));
end

function [window_start, window_end] = existingOverlapWindow(src, L, existing)
   %EXISTINGOVERLAPWINDOW Common usable span for reused RCM artifacts.
   [window_start, window_end] = requiredWindow(L);
   if ~isempty(existing.data_files)
      window_start = max(window_start, existing.data_start);
      window_end = min(window_end, existing.data_end);
   end
   if src ~= "racmo" && ~isempty(existing.met_files)
      window_start = max(window_start, existing.met_start);
      window_end = min(window_end, existing.met_end);
   end
end

function L = legWithDiscoveryWindow(L, t1, t2)
   %LEGWITHDISCOVERYWINDOW Preserve the case period only for cache discovery.
   if ~isnat(t1) && ~isnat(t2)
      L.discovery_start = t1;
      L.discovery_end = t2;
   end
end

function tf = hasDiscoveryWindow(L)
   %HASDISCOVERYWINDOW True when L has a finite cached-artifact search window.
   [t1, t2] = discoveryWindow(L);
   tf = ~isnat(t1) && ~isnat(t2);
end

function [t1, t2] = discoveryWindow(L)
   %DISCOVERYWINDOW Return the broad window used to rank cached artifacts.
   if isfield(L, 'discovery_start') && isfield(L, 'discovery_end') ...
         && ~isnat(L.discovery_start) && ~isnat(L.discovery_end)
      t1 = L.discovery_start;
      t2 = L.discovery_end;
   else
      t1 = L.start;
      t2 = L.end;
   end
end

function [t1, t2] = requiredWindow(L)
   %REQUIREDWINDOW Return attainable source coverage required for cache reuse.
   if isfield(L, 'staged') && logical(L.staged)
      t1 = L.start;
      t2 = L.end;
   else
      [t1, t2] = discoveryWindow(L);
   end
end

function tf = isReusableRcmFile(filename, src, requested_method, point)
   %ISREUSABLERCMFILE True when staged metadata matches method and point.
   if isempty(filename) || strlength(string(filename)) == 0
      tf = false;
      return
   end

   % Reusable RCM files must carry the current artifact metadata contract.
   saved_method = savedSampleMethod(filename);
   tf = saved_method == requested_method;
   if tf
      tf = savedPointMatches(filename, point);
   end
   if tf && src == "racmo"
      % A pre-mask RACMO cache cannot prove which grid cell supplied its
      % payload. Force one canonical rebuild so off-mask coastal artifacts do
      % not bypass the current native-mask point-selection contract.
      tf = savedRacmoMaskContract(filename);
   end
end

function tf = savedRacmoMaskContract(filename)
   %SAVEDRACMOMASKCONTRACT True for artifacts built by the mask-aware selector.
   tf = false;
   metadata = ...
      icemodel.verification.setup.readRcmArtifactMetadata(string(filename));
   if ~isstruct(metadata) || ~isfield(metadata, 'racmo_ice_mask_applied') ...
         || ~isfield(metadata, 'racmo_point_max_distance_m')
      return
   end
   applied = metadata.racmo_ice_mask_applied;
   max_distance = double(metadata.racmo_point_max_distance_m);
   tf = isscalar(applied) && (islogical(applied) || isnumeric(applied)) ...
      && isfinite(double(applied)) && logical(applied) ...
      && isscalar(max_distance) && isfinite(max_distance) ...
      && max_distance > 0;
end

function method = savedSampleMethod(filename)
   %SAVEDSAMPLEMETHOD Read sample_method without loading full staged payloads.
   method = "";
   metadata = ...
      icemodel.verification.setup.readRcmArtifactMetadata(string(filename));
   if isstruct(metadata) && isfield(metadata, 'sample_method')
      method = string(metadata.sample_method);
   end
end

function tf = savedPointMatches(filename, point)
   %SAVEDPOINTMATCHES True when artifact metadata was written for this point.
   tf = false;
   metadata = ...
      icemodel.verification.setup.readRcmArtifactMetadata(string(filename));
   if ~isstruct(metadata) || ~isfield(metadata, 'lat_wgs84') ...
         || ~isfield(metadata, 'lon_wgs84')
      return
   end
   tf = abs(double(metadata.lat_wgs84) - point(1)) <= 1e-8 ...
      && abs(double(metadata.lon_wgs84) - point(2)) <= 1e-8;
end

function outdir = metOutdir(kwargs)
   %METOUTDIR Match writemet's default output root when callers omit it.
   [outdir, ~] = icemodel.verification.setup.rcmArtifactOutputDirs( ...
      kwargs.met_outdir, kwargs.userdata_outdir);
end

function kwargs = inferManifestModeOutputDirs(kwargs, manifest_file)
   %INFERMANIFESTMODEOUTPUTDIRS Pair <root>/eval manifests with <root>/input.
   [family_dir, ~, ~] = fileparts(char(manifest_file));
   [eval_dir, ~, ~] = fileparts(family_dir);
   [root_dir, eval_name, ~] = fileparts(eval_dir);
   if string(eval_name) ~= "eval"
      return
   end

   if kwargs.met_outdir == ""
      kwargs.met_outdir = string(fullfile(root_dir, 'input', 'met'));
   end
   if kwargs.userdata_outdir == ""
      kwargs.userdata_outdir = string(fullfile(root_dir, 'input', 'userdata'));
   end
end

function outdir = userdataOutdir(kwargs)
   %USERDATAOUTDIR Match writeuserdata's default output root when callers omit it.
   [~, outdir] = icemodel.verification.setup.rcmArtifactOutputDirs( ...
      kwargs.met_outdir, kwargs.userdata_outdir);
end

function suffix = metWindowSuffix(dt_out)
   %METWINDOWSUFFIX File suffix used by writemet for the requested met cadence.
   if dt_out == ""
      suffix = "_1hr.mat";
   else
      suffix = "_" + dt_out + ".mat";
   end
end

function kind = rcmLegKind(src)
   %RCMLEGKIND Manifest 'kind' label for a source's colocation leg.
   if src == "racmo"
      kind = 'point_data_smb_eval';
   else
      kind = 'point_met';
   end
end

function id = rcmProductId(src)
   %RCMPRODUCTID Explicit product identifier for one RCM source label.
   id = icemodel.verification.namelists.rcmProductIds(src);
end

function token = rcmStorageToken(src)
   %RCMSTORAGETOKEN Versioned token used in staged RCM paths and filenames.
   token = rcmProductId(src);
end

function leg = skippedLeg(kind, reason)
   %SKIPPEDLEG Manifest entry for a leg with no on-disk coverage or a failure.
   leg = struct('kind', kind, 'staged', false, ...
      'reason', char(string(reason)));
end

function tt = windowSubset(tt, t1, t2)
   %WINDOWSUBSET Clamp a timetable to [t1, t2] on a UTC-aware axis (no-op blank).
   if isnat(t1) || isnat(t2)
      return
   end
   t = tt.Time;
   if isempty(t.TimeZone)
      t.TimeZone = 'UTC';
   end
   keep = t >= t1 & t <= t2;
   tt = tt(keep, :);
end
