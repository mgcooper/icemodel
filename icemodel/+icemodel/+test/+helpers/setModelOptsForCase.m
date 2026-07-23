function opts = setModelOptsForCase(c, kwargs)
   %SETMODELOPTSFORCASE Build resolved model OPTS for one case.
   %
   %  opts = icemodel.test.helpers.setModelOptsForCase(c)
   %  opts = icemodel.test.helpers.setModelOptsForCase(c, include_spinup=false)
   %  opts = icemodel.test.helpers.setModelOptsForCase(case_manifest)
   %  opts = icemodel.test.helpers.setModelOptsForCase(case_manifest, ...
   %     startdate=..., enddate=..., dt=900)
   %
   % Accepts either shape and dispatches:
   %
   %   Formal case (regression / perf / baseline matrices)
   %     Row table or scalar struct with fields smbmodel, sitename, forcings,
   %     userdata, uservars, simyear or simyears, optional n_spinup_years
   %     and solver. Resolves simyears + spinup policy, then applies
   %     case-specific overrides on top of the public icemodel.setopts
   %     contract.
   %
   %   Verification manifest (verification suite)
   %     Scalar struct with case_id and period. Uses the case_id as sitename,
   %     derives a selectable forcing from declared staged sources, derives
   %     simyears from a bounded window, and carries the window through
   %     opts.startdate / opts.enddate so the standard configureRun +
   %     createMetFileNames + loadmet chain resolves staged met files with no
   %     special-case bypass.
   %
   % Inputs
   %   c   Formal case row/struct OR verification case manifest.
   %
   % Name-value
   %   include_spinup : logical, default true. Formal-case only. When false,
   %       drop the n_spinup_years leading years from simyears.
   %   smbmodel : string, default "icemodel". Verification-manifest only.
   %   testname : string, default "verification". Verification-manifest only.
   %   dt : double seconds. Verification-manifest only. Default derived from
   %       the authoritative selected met filename, then native_timestep when
   %       no staged met path is available ("hourly" -> 3600).
   %   startdate / enddate : datetime, default NaT. Verification-manifest only.
   %       Explicit runtime window; defaults to the manifest's bounded period.
   %
   % See also: icemodel.setopts, icemodel.createMetFileNames,
   %  icemodel.verification.runIcemodelSnowCandidate

   arguments
      c
      kwargs.include_spinup (1, 1) logical = true
      kwargs.smbmodel (1, 1) string = "icemodel"
      kwargs.testname (1, 1) string = "verification"
      kwargs.dt double = []
      kwargs.startdate = NaT('TimeZone', 'UTC')
      kwargs.enddate = NaT('TimeZone', 'UTC')
   end

   % Accept either a one-row case table or an equivalent scalar struct.
   if istable(c)
      assert(height(c) == 1, ...
         'setModelOptsForCase expects a single table row')
      c = table2struct(c);
   end

   % Dispatch on input shape. Verification manifests carry case_id plus the
   % canonical period field. Formal cases carry simyear / simyears plus the full
   % positional contract for icemodel.setopts.
   if isVerificationManifest(c)
      opts = optsFromVerificationManifest(c, kwargs);
   else
      opts = optsFromFormalCase(c, kwargs);
   end
end

function tf = isVerificationManifest(c)
   %ISVERIFICATIONMANIFEST True for resolved verification case manifests.
   tf = isfield(c, 'case_id') && isfield(c, 'period');
end

function opts = optsFromFormalCase(c, kwargs)
   %OPTSFROMFORMALCASE Build opts from a formal regression/perf case row.

   [simyears, n_spinup_years] = resolveSimulationYears(c, kwargs.include_spinup);

   opts = icemodel.setopts(c.smbmodel, c.sitename, simyears, ...
      c.forcings, c.userdata, c.uservars, string.empty(), false, false, ...
      'n_spinup_years', n_spinup_years);

   if isfield(c, 'solver') && ~isempty(c.solver)
      opts = icemodel.resetopts(opts, 'solver', c.solver);
   end
end

function opts = optsFromVerificationManifest(case_manifest, kwargs)
   %OPTSFROMVERIFICATIONMANIFEST Build opts from a verification case manifest.

   sitename = char(string(case_manifest.case_id));
   forcings = verificationForcing(case_manifest, sitename);
   [window_start, window_end] = resolveComparisonWindow(case_manifest, ...
      kwargs.startdate, kwargs.enddate, forcings);

   % Resolve manifest-authoritative artifact paths before cadence. A partial
   % list can select a different covering file than its first stale record, so
   % runtime timestep must follow the exact path that metfname will load.
   metfname = {};
   staged_dt = [];
   userdatafname = {};
   has_input_root = isfield(case_manifest, 'input_data_root') ...
      && strlength(string(case_manifest.input_data_root)) > 0;
   if has_input_root
      [metfname, staged_dt] = stagedManifestMetFiles(case_manifest, forcings, ...
         case_manifest.input_data_root, window_start, window_end);
      userdatafname = stagedManifestUserdataFiles( ...
         case_manifest, forcings, case_manifest.input_data_root);
   end
   dt_seconds = resolveTimestep( ...
      case_manifest, kwargs.dt, forcings, metfname, staged_dt);

   smbmodel = char(kwargs.smbmodel);
   simyears = year(window_start):year(window_end);
   userdata = [];
   uservars = [];
   testname = char(kwargs.testname);
   saveflag = false;
   backupflag = false;

   % Assemble all runtime overrides before setopts performs its first
   % configureRun validation. A manifest-selected input root must not depend on
   % an unrelated active ICEMODEL_INPUT_PATH merely to reach the later reset.
   overrides = { ...
      'dt', dt_seconds, ...
      'startdate', window_start, ...
      'enddate', window_end};

   if has_input_root
      overrides = [overrides, { ...
         'pathinput', ...
         char(case_manifest.input_data_root), ...
         'pathuserdata', char(fullfile(case_manifest.input_data_root, ...
         'userdata')), ...
         'metfname', metfname, ...
         'userdatafname', userdatafname}];
   end

   % Build the finalized runtime contract once with the selected scoped paths.
   opts = icemodel.setopts(smbmodel, sitename, simyears, forcings, ...
      userdata, uservars, testname, saveflag, backupflag, overrides{:});
end

function forcings = verificationForcing(case_manifest, sitename)
   %VERIFICATIONFORCING Choose a selectable met source for a manifest case.
   valid = icemodel.namelists.forcings();
   if isfield(case_manifest, 'forcing_sources')
      declared = string(case_manifest.forcing_sources);
      declared = declared(ismember(declared, valid));
      if ~isempty(declared)
         % `forcing_ready` is advisory metadata for downstream fill/repair. A
         % declared source means a met file exists and can be selected here.
         forcings = char(declared(1));
         return
      end
      if ~allowsLegacyForcingFallback(case_manifest)
         error('icemodel:test:setModelOptsForCase:noRunnableForcing', ...
            'verification manifest %s does not declare a selectable forcing source', ...
            char(string(case_manifest.case_id)))
      end
   end

   candidates = string.empty();
   if isfield(case_manifest, 'dataset_family') ...
         && strlength(string(case_manifest.dataset_family)) > 0
      candidates(end+1) = string(case_manifest.dataset_family);
   end
   candidates(end+1) = string(sitename);

   runnable = candidates(ismember(candidates, valid));
   if ~isempty(runnable)
      forcings = char(runnable(1));
      return
   end

   error('icemodel:test:setModelOptsForCase:noRunnableForcing', ...
      'verification manifest %s does not declare a selectable forcing source', ...
      char(string(case_manifest.case_id)))
end

function tf = allowsLegacyForcingFallback(case_manifest)
   %ALLOWSLEGACYFORCINGFALLBACK True for legacy snow manifests without sources.
   tf = isfield(case_manifest, 'dataset_family') ...
      && string(case_manifest.dataset_family) == "esm_snowmip";
end

function [simyears, n_spinup_years] = resolveSimulationYears(c, include_spinup)
   %RESOLVESIMULATIONYEARS Resolve retained and spinup years for one formal case.

   if isfield(c, 'simyears') && ~isempty(c.simyears)
      simyears = c.simyears;
   elseif isfield(c, 'simyear') && ~isempty(c.simyear)
      simyears = [c.simyear - 1, c.simyear];
   else
      error('formal case is missing simyear/simyears')
   end

   if isfield(c, 'n_spinup_years') && ~isempty(c.n_spinup_years)
      n_spinup_years = c.n_spinup_years;
   elseif isfield(c, 'simyear') && ~isempty(c.simyear)
      n_spinup_years = 1;
   else
      n_spinup_years = 0;
   end

   if ~include_spinup
      if n_spinup_years >= numel(simyears)
         error('n_spinup_years must be smaller than numel(simyears)')
      end
      if n_spinup_years > 0
         simyears = simyears(n_spinup_years+1:end);
      end
      n_spinup_years = 0;
   end
end

function [window_start, window_end] = resolveComparisonWindow(manifest, ...
      startdate, enddate, forcing_source)
   %RESOLVECOMPARISONWINDOW Resolve the canonical window for a manifest case.

   if xor(isnat(startdate), isnat(enddate))
      error('icemodel:test:setModelOptsForCase:halfWindow', ...
         'startdate and enddate must be provided together')
   end

   has_explicit_window = ~isnat(startdate) && ~isnat(enddate);
   if has_explicit_window
      explicit_start = icemodel.verification.setup.ensureUtc(startdate);
      explicit_end = icemodel.verification.setup.ensureUtc(enddate);
   end

   colocation_window = colocationWindow(manifest, forcing_source);
   [colocation_start, colocation_end] = parseWindow(colocation_window);
   manifest_window = manifestWindow(manifest);
   [manifest_start, manifest_end] = parseWindow(manifest_window);

   if ~isnat(colocation_start) && ~isnat(colocation_end)
      if ~isnat(manifest_start) && ~isnat(manifest_end)
         available_start = max([colocation_start, manifest_start]);
         available_end = min([colocation_end, manifest_end]);
         if available_start > available_end
            error('icemodel:test:setModelOptsForCase:emptyWindow', ...
               ['verification manifest %s has no overlap between %s ' ...
               'forcing and case period'], char(string(manifest.case_id)), ...
               char(string(forcing_source)))
         end
      else
         available_start = colocation_start;
         available_end = colocation_end;
      end
   else
      available_start = manifest_start;
      available_end = manifest_end;
   end

   if has_explicit_window
      if isnat(available_start) || isnat(available_end)
         window_start = explicit_start;
         window_end = explicit_end;
      else
         window_start = max([explicit_start, available_start]);
         window_end = min([explicit_end, available_end]);
         if window_start > window_end
            error('icemodel:test:setModelOptsForCase:emptyWindow', ...
               ['explicit verification window for %s has no overlap with ' ...
               '%s forcing and case period'], char(string(manifest.case_id)), ...
               char(string(forcing_source)))
         end
      end
      return
   end

   window_start = available_start;
   window_end = available_end;
   if ~isnat(window_start) && ~isnat(window_end)
      return
   end
   error('icemodel:test:setModelOptsForCase:unboundedWindow', ...
      ['verification manifest %s needs a bounded period, a bounded ' ...
      'colocation window for %s, or an explicit window'], ...
      char(string(manifest.case_id)), char(string(forcing_source)))
end

function window = manifestWindow(manifest)
   %MANIFESTWINDOW Return the manifest runtime window.
   window = manifest.period;
end

function window = colocationWindow(manifest, forcing_source)
   %COLOCATIONWINDOW Return the selected forcing leg's staged window if present.
   window = struct('start', '', 'end', '');
   if ~isfield(manifest, 'colocation')
      return
   end
   for source = colocationSourceKeys(forcing_source)
      key = char(source);
      if isfield(manifest.colocation, key)
         leg = manifest.colocation.(key);
         if isfield(leg, 'window')
            window = leg.window;
            return
         end
      end
   end
end

function keys = colocationSourceKeys(forcing_source)
   %COLOCATIONSOURCEKEYS Return product-id and short RCM keys for lookup.
   source = string(forcing_source);
   labels = icemodel.verification.namelists.rcmsources();
   product_ids = icemodel.verification.namelists.rcmProductIds(labels);
   match = product_ids == source;
   if any(match)
      keys = unique([source, labels(match)], 'stable');
   else
      keys = source;
   end
end

function [paths, dt_seconds] = stagedManifestMetFiles(case_manifest, forcing_source, ...
      input_root, window_start, window_end)
   %STAGEDMANIFESTMETFILES Resolve explicitly recorded staged met artifacts.
   % Recorded paths are authoritative. Select a recorded file whose saved Time
   % axis encloses the run, or prove that an ordered subset of recorded files
   % supplies one continuous, non-overlapping, single-cadence runtime axis.
   dt_seconds = [];
   paths = stagedManifestArtifactFiles(case_manifest, forcing_source, ...
      input_root, "met", "met_files");
   if isempty(paths)
      return
   end

   % A manifest-selected list must be provable from its recorded artifacts. Fail
   % here instead of letting configureRun discover an unrelated filename sibling.
   exists = cellfun(@isfile, paths);
   if ~any(exists)
      error('icemodel:test:setModelOptsForCase:invalidManifestMetFiles', ...
         ['recorded met_files for %s do not exist and cannot cover the ' ...
         'requested runtime: %s'], char(string(forcing_source)), ...
         char(strjoin(string(paths), ', ')))
   end

   % Inspect the saved timetable rather than trusting filename dates or cadence.
   % Each candidate must have a strictly increasing, uniform Time axis before it
   % can participate in either single-file or split-file coverage.
   candidates = reshape(string(paths(exists)), [], 1);
   starts = NaT(numel(candidates), 1, 'TimeZone', window_start.TimeZone);
   ends = starts;
   cadences = nan(numel(candidates), 1);
   for n = 1:numel(candidates)
      payload = load(candidates(n), 'met');
      if ~isfield(payload, 'met') || ~istimetable(payload.met) ...
            || height(payload.met) < 2
         error('icemodel:test:setModelOptsForCase:invalidManifestMetFiles', ...
            'recorded met_file does not contain a provable timetable: %s', ...
            char(candidates(n)))
      end
      times = icemodel.verification.setup.ensureUtc(payload.met.Time);
      steps = seconds(diff(times));
      cadence = median(steps);
      tolerance = max(1e-6, abs(cadence) * 1e-9);
      if ~isfinite(cadence) || cadence <= 0 ...
            || any(abs(steps - cadence) > tolerance)
         error('icemodel:test:setModelOptsForCase:invalidManifestMetFiles', ...
            'recorded met_file has a nonuniform Time axis: %s', ...
            char(candidates(n)))
      end
      starts(n) = times(1);
      ends(n) = times(end);
      cadences(n) = cadence;
   end

   % Prefer one saved artifact that encloses the complete runtime. This retains
   % the deterministic widest/latest/path policy without loading list siblings.
   covering = starts <= window_start & ends >= window_end;
   if any(covering)
      selected_paths = candidates(covering);
      selected_starts = starts(covering);
      selected_ends = ends(covering);
      selected_cadences = cadences(covering);
      durations = seconds(selected_ends - selected_starts);
      end_key = year(selected_ends) .* 10000 + ...
         month(selected_ends) .* 100 + day(selected_ends);
      rank = table(-durations, -end_key, selected_paths, ...
         'VariableNames', {'neg_duration', 'neg_end', 'path'});
      [~, order] = sortrows(rank, {'neg_duration', 'neg_end', 'path'});
      paths = {char(selected_paths(order(1)))};
      dt_seconds = selected_cadences(order(1));
      return
   end

   % A legitimate split list contributes only files whose saved support meets
   % the requested window. Sort that subset by actual support before validating
   % its cadence, overlap, adjacency, and endpoint coverage.
   selected = starts <= window_end & ends >= window_start;
   if ~any(selected)
      error('icemodel:test:setModelOptsForCase:invalidManifestMetFiles', ...
         ['recorded met_files do not cover the requested runtime for %s: ' ...
         '%s'], char(string(forcing_source)), ...
         char(strjoin(string(paths), ', ')))
   end
   selected_paths = candidates(selected);
   selected_starts = starts(selected);
   selected_ends = ends(selected);
   selected_cadences = cadences(selected);
   rank = table(selected_starts, selected_ends, selected_paths, ...
      'VariableNames', {'start_time', 'end_time', 'path'});
   [~, order] = sortrows(rank, {'start_time', 'end_time', 'path'});
   selected_paths = selected_paths(order);
   selected_starts = selected_starts(order);
   selected_ends = selected_ends(order);
   selected_cadences = selected_cadences(order);

   cadence = selected_cadences(1);
   tolerance = max(1e-6, abs(cadence) * 1e-9);
   if any(abs(selected_cadences - cadence) > tolerance)
      error('icemodel:test:setModelOptsForCase:invalidManifestMetFiles', ...
         'recorded met_files have inconsistent saved cadences: %s', ...
         char(strjoin(selected_paths, ', ')))
   end
   for n = 2:numel(selected_paths)
      separation = seconds(selected_starts(n) - selected_ends(n - 1));
      if separation <= 0
         error('icemodel:test:setModelOptsForCase:invalidManifestMetFiles', ...
            'recorded met_files overlap or duplicate Time rows: %s', ...
            char(strjoin(selected_paths, ', ')))
      end
      if abs(separation - cadence) > tolerance
         error('icemodel:test:setModelOptsForCase:invalidManifestMetFiles', ...
            'recorded met_files leave a gap in the requested runtime: %s', ...
            char(strjoin(selected_paths, ', ')))
      end
   end
   if selected_starts(1) > window_start || selected_ends(end) < window_end
      error('icemodel:test:setModelOptsForCase:invalidManifestMetFiles', ...
         ['recorded met_files do not cover the requested runtime for %s: ' ...
         '%s'], char(string(forcing_source)), ...
         char(strjoin(string(paths), ', ')))
   end
   paths = cellstr(reshape(selected_paths, 1, []));
   dt_seconds = cadence;
end

function paths = stagedManifestUserdataFiles(case_manifest, forcing_source, ...
      input_root)
   %STAGEDMANIFESTUSERDATAFILES Resolve exact staged Data artifact variants.
   paths = stagedManifestArtifactFiles(case_manifest, forcing_source, ...
      input_root, "userdata", "data_files");
end

function paths = stagedManifestArtifactFiles(case_manifest, forcing_source, ...
      input_root, folder, field)
   %STAGEDMANIFESTARTIFACTFILES Resolve one selected leg's recorded file field.
   paths = {};
   files = recordedManifestArtifactFiles( ...
      case_manifest, forcing_source, field);
   if isempty(files)
      return
   end
   candidates = arrayfun(@(file) fullfile(string(input_root), folder, file), ...
      files, UniformOutput=false);
   % Preserve every recorded path; family-specific callers may select a proven
   % covering member but never replace this manifest whitelist with discovery.
   paths = cellfun(@char, candidates, UniformOutput=false);
end

function [window_start, window_end] = parseWindow(window)
   %PARSEWINDOW Convert a manifest window struct into UTC datetimes.
   window_start = icemodel.verification.setup.ensureUtc(window.start);
   window_end = icemodel.verification.setup.ensureUtc(window.end);
end

function dt_seconds = resolveTimestep(case_manifest, dt_override, ...
      forcing_source, selected_met_files, staged_dt)
   %RESOLVETIMESTEP Map the selected met artifact cadence to dt seconds.

   if ~isempty(staged_dt)
      % Explicit manifest met is authoritative. A caller may restate its saved
      % cadence, but cannot silently configure a different model timestep.
      if ~isempty(dt_override) && (~isscalar(dt_override) ...
            || ~isfinite(dt_override) || dt_override ~= staged_dt)
         error('icemodel:test:setModelOptsForCase:metCadenceConflict', ...
            ['requested dt=%s conflicts with authoritative saved met ' ...
            'cadence %g seconds for %s'], mat2str(dt_override), staged_dt, ...
            char(strjoin(string(selected_met_files), ', ')))
      end
      dt_seconds = staged_dt;
      return
   end
   if ~isempty(dt_override)
      % Legacy manifests without explicit met artifacts retain caller control.
      dt_seconds = dt_override;
      return
   end
   label = metFileTimestepLabel( ...
      case_manifest, forcing_source, selected_met_files);
   if label == "" && isfield(case_manifest, 'native_timestep')
      label = lower(string(case_manifest.native_timestep));
   end
   switch label
      case {"hourly", "1hr", "1 hour", "1-hour"}
         dt_seconds = 3600;
      case {"15-min", "15 min", "15m", "15-minute", "15 minute"}
         dt_seconds = 900;
      case {"30-min", "30 min", "30m", "30-minute", "30 minute", ...
            "00:30:00"}
         dt_seconds = 1800;
      otherwise
         dt_seconds = 3600;
   end
end

function label = metFileTimestepLabel( ...
      case_manifest, forcing_source, selected_met_files)
   %METFILETIMESTEPLABEL Read the authoritative met filename cadence suffix.
   label = "";
   files = reshape(string(selected_met_files), 1, []);
   if isempty(files)
      % Legacy/no-root manifests have no resolved absolute path, so retain the
      % declared-leg fallback used before exact-path selection was available.
      files = recordedManifestArtifactFiles( ...
         case_manifest, forcing_source, "met_files");
   end
   if isempty(files)
      return
   end
   token = regexp(files(1), "_(15m|30m|1hr)\.mat$", "tokens", "once");
   if ~isempty(token)
      label = string(token{1});
   end
end

function files = recordedManifestArtifactFiles(case_manifest, forcing_source, field)
   %RECORDEDMANIFESTARTIFACTFILES Return one selected leg's nonblank records.
   files = strings(1, 0);
   if ~isfield(case_manifest, 'colocation')
      return
   end
   for source = colocationSourceKeys(forcing_source)
      key = char(source);
      if ~isfield(case_manifest.colocation, key)
         continue
      end
      leg = case_manifest.colocation.(key);
      name = char(field);
      if ~isfield(leg, name) || isempty(leg.(name))
         continue
      end
      files = string(leg.(name));
      files = files(strlength(files) > 0);
      if ~isempty(files)
         return
      end
   end
end
