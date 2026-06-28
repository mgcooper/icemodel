function out = stageRcmForcing(points, kwargs)
   %STAGERCMFORCING Stage RCM forcing + Data for a list of points.
   %
   %  colocation = icemodel.verification.setup.stageRcmForcing(points, ...
   %     legspec=L, models=["mar","merra","racmo"], met_outdir=..., ...
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
   %    * MAR / MERRA : the COMPLETE Data timetable via writeuserdata
   %      (userdata/<src>/) AND the met-contract timetable via data2met+writemet
   %      (met/<src>/). ALWAYS BOTH. The Data file carries every channel; the met
   %      file is the forcing subset. (Earlier code wrote only met and discarded
   %      the Data - the regression this builder fixes.)
   %    * RACMO       : the Data timetable via writeuserdata only (no met). The
   %      available RACMO 2.3p3 subsurface files lack the near-surface met STATE
   %      channels (tair/wspd/rh/psfc), so a RACMO met would fail validatemet;
   %      RACMO is staged as eval/reference Data, not a met source.
   %
   %  EXECUTION + ERROR HANDLING
   %    Each source is staged only for the files it still needs. Existing window
   %    files whose encoded date range encloses the requested leg are reused and
   %    recorded in the manifest instead of building a narrower duplicate. Points
   %    that still need source reads are grouped by identical year set, so one
   %    long station record does not force every shorter station to read the long
   %    union span. Sources are processed in order and each source's files are
   %    written IMMEDIATELY, so a later source's failure never rolls back an
   %    earlier source already on disk (MAR survives a RACMO failure). A
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
   %    models : subset of ["mar","merra","racmo"] (default all three).
   %    met_outdir, userdata_outdir : output dirs (met_<...> / userdata files).
   %    mar_dir, merra_dir, racmo_dir, modis_dir : raw-source dirs per model.
   %    method : point sampling, "nearest" (default) or "natural".
   %    dt_out : optional met output timestep ("15m"); Data stays hourly.
   %    obs_manifest : path or struct selecting manifest-convenience mode.
   %    manifest_file : where to persist in manifest mode (when obs_manifest is a
   %        struct). Defaults to the obs_manifest path when that is a path.
   %    overwrite_family : passed through to writeFamilyManifestMerge (manifest
   %        mode); default false (MERGE).
   %
   % See also: icemodel.verification.setup.resolveLegWindows,
   %  icemodel.verification.setup.promiceSourceCoverage,
   %  icemodel.verification.setup.importPromiceSites,
   %  icemodel.verification.setup.importSumup,
   %  icemodel.forcing.buildMarData, icemodel.forcing.data2met

   arguments
      points double = []
      kwargs.legspec struct = struct([])
      kwargs.models (1, :) string {mustBeMember(kwargs.models, ...
         ["mar", "merra", "racmo"])} = ["mar", "merra", "racmo"]
      kwargs.met_outdir (1, 1) string = ""
      kwargs.userdata_outdir (1, 1) string = ""
      kwargs.mar_dir (1, 1) string = ""
      kwargs.merra_dir (1, 1) string = ""
      kwargs.racmo_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.method (1, 1) string = "nearest"
      kwargs.dt_out (1, 1) string = ""
      kwargs.obs_manifest = ""
      kwargs.manifest_file (1, 1) string = ""
      kwargs.overwrite_family (1, 1) logical = false
   end

   models = reshape(kwargs.models, 1, []);

   % Mode select: a struct or a non-empty path in obs_manifest -> manifest mode.
   om = kwargs.obs_manifest;
   manifest_mode = isstruct(om) || (~isstruct(om) && strlength(string(om)) > 0);

   if manifest_mode
      out = stageFromManifest(models, kwargs);
   else
      out = stageExplicit(points, kwargs.legspec, models, kwargs);
   end
end

%% Explicit mode
function colocation = stageExplicit(points, legspec, models, kwargs)
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
   for src = models
      colocation = stageOneSource(src, points, legspec, colocation, kwargs);
   end
end

function colocation = stageOneSource(src, points, legspec, colocation, kwargs)
   %STAGEONESOURCE Build one source by needed year groups, write per point.
   srcc = char(src);
   kind = rcmLegKind(src);
   n = size(points, 1);

   % Participants: points whose leg for this source resolved staged=true and
   % whose required output files are not already bracketed by an existing staged
   % window file. A non-participant gets a skip-with-reason now (its coverage gap
   % was found cheaply by resolveLegWindows, before any build).
   part = false(1, n);
   existing = repmat(emptyExistingFiles, 1, n);
   for k = 1:n
      L = legspec(k).(srcc);
      if L.staged
         existing(k) = existingRcmFiles(src, legspec(k).alias, L, kwargs);
         if existingLegComplete(src, existing(k))
            colocation{k}.(srcc) = existingRcmLeg(src, kind, L, existing(k), kwargs);
            warning('icemodel:verification:stageRcmForcing:existingWindowFile', ...
               ['Skipping %s build for %s because existing staged file(s) ' ...
               'already bracket %s to %s.'], upper(srcc), ...
               legspec(k).alias, string(L.start), string(L.end));
         else
            part(k) = true;
         end
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
         % The source build threw (vanished dir, read error): degrade this
         % source/year group only. Sources already written stand.
         for j = 1:numel(gidx)
            colocation{gidx(j)}.(srcc) = skippedLeg(kind, build_err.message);
         end
         continue
      end

      % Per point: clip to its own window and write only the outputs that were
      % not already covered. A write failure degrades only that point's leg.
      for j = 1:numel(gidx)
         k = gidx(j);
         L = legspec(k).(srcc);
         try
            d = windowSubset(Data{j}, L.start, L.end);
            colocation{k}.(srcc) = writeRcmLeg( ...
               src, d, legspec(k).alias, kind, L, kwargs, existing(k));
         catch write_err
            colocation{k}.(srcc) = skippedLeg(kind, write_err.message);
         end
      end
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
            method=kwargs.method));
      case "merra"
         Data = asCell(icemodel.forcing.buildMerraData(pts, years, ...
            source_dir=kwargs.merra_dir, modis_dir=kwargs.modis_dir, ...
            method=kwargs.method));
      case "racmo"
         Data = asCell(icemodel.forcing.buildRacmoData(pts, years, ...
            source_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir, ...
            method=kwargs.method, dt="1hr"));
   end
end

function co = writeRcmLeg(src, d, alias, kind, L, kwargs, existing)
   %WRITERCMLEG Write one point's staged leg + build its colocation record.
   % MAR/MERRA write BOTH the full Data (userdata) AND the met (data2met); RACMO
   % writes Data only. Existing enclosing files are reused so staging a narrower
   % window cannot create duplicate files. Every staged leg carries staged==true.
   if nargin < 7
      existing = emptyExistingFiles;
   end
   met_outdir = metOutdir(kwargs);
   userdata_outdir = userdataOutdir(kwargs);
   if src == "racmo"
      data_files = existing.data_files;
      if isempty(data_files)
         data_files = icemodel.forcing.helpers.writeuserdata( ...
            d, alias, "racmo", outdir=userdata_outdir, naming="window");
      end
      co = struct('kind', kind, 'staged', true, ...
         'data_files', ...
         icemodel.verification.setup.relpaths(data_files, userdata_outdir), ...
         'sample_method', 'nearest', ...
         'window', windowStruct(L.start, L.end), ...
         'note', 'SMB/eval Data only; RACMO is not a met source.');
   else
      % MAR/MERRA: the COMPLETE Data timetable AND the met forcing it derives.
      data_files = existing.data_files;
      if isempty(data_files)
         data_files = icemodel.forcing.helpers.writeuserdata( ...
            d, alias, src, outdir=userdata_outdir, naming="window");
      end
      met_files = existing.met_files;
      if isempty(met_files)
         met = toMet(d, kwargs.dt_out);
         met_files = icemodel.forcing.helpers.writemet( ...
            met, alias, src, outdir=met_outdir, naming="window");
      end
      co = struct('kind', kind, 'staged', true, ...
         'met_files', ...
         icemodel.verification.setup.relpaths(met_files, met_outdir), ...
         'data_files', ...
         icemodel.verification.setup.relpaths(data_files, userdata_outdir), ...
         'sample_method', 'nearest', ...
         'window', windowStruct(L.start, L.end));
   end
end

function met = toMet(d, dt_out)
   %TOMET Convert a Data timetable to the met contract (optional resample).
   met = icemodel.forcing.data2met(d);
   if dt_out ~= ""
      % icemodel.interpmet resamples one calendar year at a time.
      yrs = unique(year(met.Time));
      parts = cell(numel(yrs), 1);
      for n = 1:numel(yrs)
         parts{n} = icemodel.interpmet( ...
            met(year(met.Time) == yrs(n), :), char(dt_out));
      end
      met = vertcat(parts{:});
   end
end

%% Manifest-convenience mode
function manifest = stageFromManifest(models, kwargs)
   %STAGEFROMMANIFEST Re-stage RCM forcing for every case in a staged manifest.
   % Reads each case's point (site_location), window (period) and alias
   % (case_id), resolves the legs from a fresh coverage probe, stages via the
   % explicit path, merges the legs into every case and rewrites manifest.json.
   if isstruct(kwargs.obs_manifest)
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

   cases = manifest.cases;
   if isempty(cases)
      return
   end
   n = numel(cases);

   % Probe on-disk coverage once for the whole set (cheap, fail-early gate).
   coverage = icemodel.verification.setup.promiceSourceCoverage(models, ...
      struct('mar', kwargs.mar_dir, 'merra', kwargs.merra_dir, ...
      'racmo', kwargs.racmo_dir));

   % Build points + legspec from the manifest cases.
   points = zeros(n, 2);
   legspec = repmat(emptyLeg(models), 1, n);
   for k = 1:n
      c = cases(k);
      points(k, :) = [c.site_location.lat_wgs84, c.site_location.lon_wgs84];
      [t1, t2] = casePeriod(c);
      leg = icemodel.verification.setup.resolveLegWindows(models, coverage, t1, t2);
      legspec(k).alias = string(c.case_id);
      for src = models
         legspec(k).(char(src)) = leg.(char(src));
      end
   end

   colocation = stageExplicit(points, legspec, models, kwargs);

   % Merge the new RCM legs into every case + refresh the informational source
   % lists, then rewrite the manifest (MERGE by default).
   ids = strings(1, n);
   for k = 1:n
      cases(k).colocation = mergeFields(cases(k).colocation, colocation{k});
      [forcing_sources, eval_sources] = ...
         icemodel.verification.setup.colocationSourceLists( ...
         cases(k).colocation);
      cases(k).forcing_sources = cellstr(forcing_sources);
      cases(k).eval_sources = cellstr(eval_sources);
      ids(k) = string(cases(k).case_id);
   end
   manifest.cases = cases;

   manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
      manifest_file, manifest, requested_ids=ids, ...
      overwrite_family=kwargs.overwrite_family);
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

function s = mergeFields(s, add)
   %MERGEFIELDS Copy every field of `add` onto struct `s` (overwriting).
   f = fieldnames(add);
   for i = 1:numel(f)
      s.(f{i}) = add.(f{i});
   end
end

function L = emptyLeg(models)
   %EMPTYLEG Prototype legspec element (alias + one leg field per source).
   L = struct('alias', "");
   proto = struct('staged', false, 'years', [], 'start', NaT, 'end', NaT, ...
      'reason', "");
   for src = models
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
      'data_files', strings(1, 0));
end

function existing = existingRcmFiles(src, alias, L, kwargs)
   %EXISTINGRCMFILES Locate staged window files enclosing the requested leg.
   existing = emptyExistingFiles;
   existing.data_files = findExistingWindowPath(userdataOutdir(kwargs), src, ...
      sprintf('%s_%s', alias, src), ".mat", L.start, L.end);
   if src ~= "racmo"
      existing.met_files = findExistingWindowPath(metOutdir(kwargs), src, ...
         sprintf('met_%s_%s', alias, src), metWindowSuffix(kwargs.dt_out), ...
         L.start, L.end);
   end
end

function tf = existingLegComplete(src, existing)
   %EXISTINGLEGCOMPLETE True when all outputs required for SRC already exist.
   if src == "racmo"
      tf = ~isempty(existing.data_files);
   else
      tf = ~isempty(existing.data_files) && ~isempty(existing.met_files);
   end
end

function co = existingRcmLeg(src, kind, L, existing, kwargs)
   %EXISTINGRCMLEG Manifest leg pointing at already-staged enclosing files.
   if src == "racmo"
      co = struct('kind', kind, 'staged', true, ...
         'data_files', ...
         icemodel.verification.setup.relpaths( ...
         existing.data_files, userdataOutdir(kwargs)), ...
         'sample_method', 'nearest', ...
         'window', windowStruct(L.start, L.end), ...
         'note', 'Existing staged Data file brackets requested window.');
   else
      co = struct('kind', kind, 'staged', true, ...
         'met_files', ...
         icemodel.verification.setup.relpaths( ...
         existing.met_files, metOutdir(kwargs)), ...
         'data_files', ...
         icemodel.verification.setup.relpaths( ...
         existing.data_files, userdataOutdir(kwargs)), ...
         'sample_method', 'nearest', ...
         'window', windowStruct(L.start, L.end), ...
         'note', 'Existing staged met/Data files bracket requested window.');
   end
end

function path = findExistingWindowPath(base, source, prefix, suffix, t1, t2)
   %FINDEXISTINGWINDOWPATH First staged file whose encoded window encloses T1-T2.
   path = strings(1, 0);
   for d = icemodel.forcing.helpers.sourceSearchDirs(base, source)
      name = icemodel.forcing.helpers.findEnclosingWindowFile( ...
         d{1}, prefix, suffix, t1, t2);
      if strlength(name) > 0
         path = string(fullfile(d{1}, char(name)));
         return
      end
   end
end

function outdir = metOutdir(kwargs)
   %METOUTDIR Match writemet's default output root when callers omit it.
   outdir = kwargs.met_outdir;
   if outdir == ""
      outdir = string(fullfile(icemodel.getpath('input'), 'met'));
   end
end

function outdir = userdataOutdir(kwargs)
   %USERDATAOUTDIR Match writeuserdata's default output root when callers omit it.
   outdir = kwargs.userdata_outdir;
   if outdir == ""
      outdir = string(icemodel.getpath('userdata'));
   end
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

function leg = skippedLeg(kind, reason)
   %SKIPPEDLEG Manifest entry for a leg with no on-disk coverage or a failure.
   leg = struct('kind', kind, 'staged', false, ...
      'reason', char(string(reason)));
end

function w = windowStruct(t1, t2)
   %WINDOWSTRUCT JSON-friendly window record for the manifest.
   w = struct('start', char(string(t1)), 'end', char(string(t2)));
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
