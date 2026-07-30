function opts = verifyPromiceFilledReadiness(opts, fileiter)
   %VERIFYPROMICEFILLEDREADINESS Gate derived PROMICE forcing by coverage.
   %
   %  opts = icemodel.verifyPromiceFilledReadiness(opts) verifies the
   %  producer-manifest identity of the configured promice_filled artifacts
   %  and proves complete required-channel coverage of every requested
   %  timestep between opts.startdate and opts.enddate (or the opts.simyears
   %  span when no dates are set) directly against the filled met files'
   %  samples (POLICY A4; water years and arbitrary windows are
   %  first-class). The canonical product and runtime timestep are both
   %  fixed at 15 minutes. Calendar-year ledger verdicts are producer
   %  bookkeeping and never the runtime gate. The returned options are
   %  marked for code-generation-safe loading.
   %
   % See also: icemodel.loadmet, icemodel.setopts

   if ~strcmpi(string(opts.forcings), "promice_filled")
      return
   end
   site = lower(string(opts.sitename));
   icemodel.forcing.reconstruct.mustBeStationToken(site)
   if nargin < 2 || isempty(fileiter)
      fileiter = 1:numel(opts.metfname);
   end
   met_files = string(opts.metfname(fileiter));
   if numel(fileiter) < numel(opts.metfname)
      opts = narrowRequestToSelectedFiles(opts, met_files);
   end
   if ~isfield(opts, 'dt') || ~isscalar(opts.dt) || opts.dt ~= 900
      error('icemodel:loadmet:promiceFilledCadenceMismatch', ...
         'promice_filled is published only at the 900-second cadence');
   end

   % Resolve the producer ledger without placing table I/O in loadmet's
   % code-generation surface. The ledger no longer gates the run (A4) but
   % remains a required, manifest-pinned bookkeeping artifact.
   readiness_file = "";
   if isfield(opts, 'readiness_file')
      readiness_file = string(opts.readiness_file);
   end
   if readiness_file == ""
      readiness_file = fullfile(icemodel.internal.fullpath, 'data', ...
         'preview', 'qa', 'gapfill', 'ledger', ...
         site + "-readiness.csv");
   end
   if ~isfile(readiness_file)
      error('icemodel:loadmet:promiceFilledNotReady', ...
         'PROMICE readiness ledger is unavailable: %s', readiness_file);
   end

   % The producer manifest binds both the verdict and every configured filled
   % met file to the bytes published by the reconstruction transaction.
   report_inputs_file = "";
   if isfield(opts, 'report_inputs_file')
      report_inputs_file = string(opts.report_inputs_file);
   end
   if report_inputs_file == ""
      qa_dir = fileparts(fileparts(readiness_file));
      report_inputs_file = fullfile(qa_dir, 'plans', ...
         site + "-report-inputs.json");
   end
   verifyProducerManifest(report_inputs_file, site, readiness_file, ...
      met_files(:));

   % A4 runtime gate: the filled samples themselves must cover every
   % requested timestep with the required channels. Ledger verdict strings
   % are never consulted here.
   verifyRequestedWindowCoverage(opts, site, met_files);

   opts.promice_filled_readiness_verified = true;
   opts.promice_filled_manifest_verified = true;
   opts.promice_filled_verified_forcing = char(string(opts.forcings));
   opts.promice_filled_verified_site = char(site);
   opts.promice_filled_verified_simyears = double(opts.simyears(:)).';
   opts.promice_filled_verified_metfname = opts.metfname(fileiter);
   opts.promice_filled_verified_startdate = opts.startdate;
   opts.promice_filled_verified_enddate = opts.enddate;
   opts.promice_filled_verified_calendar_type = char(string(opts.calendar_type));
   opts.promice_filled_verified_smbmodel = char(string(opts.smbmodel));
   opts.promice_filled_verified_dt = double(opts.dt);
   opts.report_inputs_file = report_inputs_file;
end

function opts = narrowRequestToSelectedFiles(opts, met_files)
   %NARROWREQUESTTOSELECTEDFILES Make fileiter the verified run contract.
   times_by_file = cell(numel(met_files), 1);
   for k = 1:numel(met_files)
      payload = load(met_files(k), 'met');
      if ~isfield(payload, 'met') || ~istimetable(payload.met) ...
            || isempty(payload.met)
         error('icemodel:loadmet:promiceFilledIdentityMismatch', ...
            'file does not contain a nonempty met timetable: %s', ...
            met_files(k));
      end
      times_by_file{k} = payload.met.Time;
   end
   times = vertcat(times_by_file{:});
   selected_years = intersect(double(opts.simyears(:)).', ...
      unique(year(times)).', 'stable');
   if isempty(selected_years)
      error('icemodel:loadmet:promiceFilledWindowUncovered', ...
         'selected met files do not cover any requested simulation year');
   end
   opts.simyears = selected_years;
   if isfield(opts, 'numyears')
      opts.numyears = numel(selected_years);
   end
   selected_start = min(times);
   selected_end = max(times);
   if isfield(opts, 'startdate') && isdatetime(opts.startdate) ...
         && ~isempty(opts.startdate) && ~isnat(opts.startdate)
      opts.startdate = max([opts.startdate; selected_start]);
   end
   if isfield(opts, 'enddate') && isdatetime(opts.enddate) ...
         && ~isempty(opts.enddate) && ~isnat(opts.enddate)
      opts.enddate = min([opts.enddate; selected_end]);
   end
end

function verifyRequestedWindowCoverage(opts, site, met_files)
   %VERIFYREQUESTEDWINDOWCOVERAGE Prove filled samples cover the request.
   %
   % POLICY A4: the runtime gate is complete required-channel coverage of
   % every requested timestep, checked against the actual samples in the
   % configured (manifest-verified) filled met files.

   % The gate channel set plus the conditional snow-model precipitation
   % requirement (A5: finite total ppt OR snowf per timestep).
   required = icemodelRequiredChannels();
   need_precip = requiresSnowfallForcing(opts.smbmodel);

   % Collect sample times and per-channel validity from every configured
   % file; per-file cells avoid growing arrays inside the loop.
   met_files = string(met_files(:));
   n_files = numel(met_files);
   times_by_file = cell(n_files, 1);
   flags_by_file = cell(n_files, 1);
   precip_by_file = cell(n_files, 1);
   for k = 1:n_files
      payload = load(met_files(k), 'met');
      % Manifest hashing proves bytes, not payload shape; refuse files
      % that do not carry the met timetable contract.
      if ~isfield(payload, 'met') || ~istimetable(payload.met)
         error('icemodel:loadmet:promiceFilledIdentityMismatch', ...
            'file does not contain a met timetable: %s', met_files(k));
      end
      met = payload.met;
      % Match the loader's timezone normalization so window comparisons
      % against opts.startdate/enddate (UTC) are exact.
      met.Time.TimeZone = 'UTC';
      if height(met) < 2 || any(diff(met.Time) ~= minutes(15))
         error('icemodel:loadmet:promiceFilledCadenceMismatch', ...
            'promice_filled file is not a contiguous 15-minute product: %s', ...
            met_files(k));
      end
      if any(mod(minute(met.Time), 15) ~= 0 | second(met.Time) ~= 0)
         error('icemodel:loadmet:promiceFilledCadenceMismatch', ...
            ['promice_filled file is not aligned to the fixed UTC ' ...
            'quarter-hour grid: %s'], met_files(k));
      end
      variables = string(met.Properties.VariableNames);
      flags = false(height(met), numel(required));
      for c = 1:numel(required)
         % An absent required channel counts as wholly missing (A5);
         % finite values count only inside the A15 scalar registry.
         if ismember(required(c), variables)
            flags(:, c) = ...
               icemodel.forcing.reconstruct.scalarValidity( ...
               required(c), met.(required(c)));
         end
      end
      precip_ok = false(height(met), 1);
      if need_precip
         % A5 snow-model requirement: scalar-valid total ppt OR snowf.
         % Snowfall shares ppt's nonnegative accumulation bounds. Any
         % finite source phases must also satisfy A10: nonnegative values,
         % no phase above a finite total, and exact mass balance for a
         % complete split. Missing phases stay admissible because the
         % runtime threshold option can resolve them.
         ppt = nan(height(met), 1);
         rainf = nan(height(met), 1);
         snowf = nan(height(met), 1);
         if ismember("ppt", variables)
            ppt = met.ppt;
            precip_ok = precip_ok ...
               | icemodel.forcing.reconstruct.scalarValidity( ...
               "ppt", met.ppt);
         end
         if ismember("rainf", variables)
            rainf = met.rainf;
         end
         if ismember("snowf", variables)
            snowf = met.snowf;
            precip_ok = precip_ok ...
               | icemodel.forcing.reconstruct.scalarValidity( ...
               "ppt", met.snowf);
         end
         phase_ok = icemodel.forcing.helpers.precipitationValidity( ...
            ppt, rainf, snowf);
         precip_ok = precip_ok & phase_ok;
      end
      times_by_file{k} = met.Time;
      flags_by_file{k} = flags;
      precip_by_file{k} = precip_ok;
   end
   times = vertcat(times_by_file{:});
   flags = vertcat(flags_by_file{:});
   precip_ok = vertcat(precip_by_file{:});

   % Merge duplicate timestamps across files: a timestep is covered when
   % any configured file covers it.
   [times, ~, index] = unique(times);
   merged = false(numel(times), size(flags, 2));
   for c = 1:size(flags, 2)
      merged(:, c) = accumarray(index, double(flags(:, c)), ...
         [numel(times), 1], @max) > 0;
   end
   flags = merged;
   precip_ok = accumarray(index, double(precip_ok), ...
      [numel(times), 1], @max) > 0;

   [window_start, window_end] = requestedWindow(opts);
   if numel(times) < 2
      % A posting cadence cannot be derived from fewer than two samples,
      % and such a product cannot cover any real request.
      error('icemodel:loadmet:promiceFilledWindowUncovered', ...
         ['promice_filled files provide %d samples; cannot cover the ' ...
         'requested window %s..%s for %s'], numel(times), ...
         string(window_start), string(window_end), site);
   end

   % The published product and requested runtime cadence are both fixed
   % at 15 minutes; per-file validation above prevents a filename suffix
   % or missing rows from weakening that contract.
   cadence = 900;

   % Requested timesteps: the fixed UTC quarter-hour grid across the
   % requested window, independent of the first product sample so a missing
   % leading row cannot redefine the coverage contract.
   anchor = dateshift(window_start, 'start', 'day');
   offset_start = ceil(seconds(window_start - anchor) / cadence);
   offset_end = floor(seconds(window_end - anchor) / cadence);
   requested = anchor + seconds(cadence * (offset_start:offset_end).');
   requested = requested(ismember(year(requested), ...
      double(opts.simyears(:))));
   if strcmp('noleap', opts.calendar_type)
      requested = requested(~(month(requested) == 2 ...
         & day(requested) == 29));
   end
   if isempty(requested)
      error('icemodel:loadmet:promiceFilledWindowUncovered', ...
         ['requested window %s..%s contains no requested timesteps for ' ...
         '%s (check startdate/enddate against simyears)'], ...
         string(window_start), string(window_end), site);
   end

   % Coverage: a requested timestep is covered per channel when a
   % configured file posts a scalar-valid sample exactly at that time.
   [present, where] = ismember(requested, times);
   labels = required;
   covered = false(numel(requested), numel(required));
   covered(present, :) = flags(where(present), :);
   if need_precip
      labels = [required, "ppt|snowf"];
      covered_precip = false(numel(requested), 1);
      covered_precip(present) = precip_ok(where(present));
      covered = [covered, covered_precip];
   end
   if all(covered(:))
      return
   end

   % Refuse with one compact per-channel account of the uncovered request:
   % absent rows first, then each channel's count and first..last times.
   problems = strings(1 + numel(labels), 1);
   n_problems = 0;
   if any(~present)
      n_problems = n_problems + 1;
      problems(n_problems) = describeUncovered("timesteps absent", ...
         requested(~present));
   end
   for c = 1:numel(labels)
      bad = ~covered(:, c);
      if any(bad)
         n_problems = n_problems + 1;
         problems(n_problems) = describeUncovered(labels(c), ...
            requested(bad));
      end
   end
   error('icemodel:loadmet:promiceFilledWindowUncovered', ...
      'promice_filled does not cover the requested window %s..%s for %s: %s', ...
      string(window_start), string(window_end), site, ...
      strjoin(problems(1:n_problems), '; '));
end

function [window_start, window_end] = requestedWindow(opts)
   %REQUESTEDWINDOW Resolve the requested coverage window (POLICY A4).
   %
   % startdate/enddate are first-class arbitrary bounds (water years);
   % each unset bound defaults to the corresponding simyears calendar
   % edge, so a dateless run requests its complete simulation years.
   years = double(opts.simyears(:));
   window_start = datetime(min(years), 1, 1, 0, 0, 0, 'TimeZone', 'UTC');
   % Inclusive cadence-agnostic end bound: the last instant of the final
   % requested year (grid construction floors onto the posting cadence).
   window_end = datetime(max(years) + 1, 1, 1, 0, 0, 0, ...
      'TimeZone', 'UTC') - seconds(1);
   if isfield(opts, 'startdate') && ~isnat(opts.startdate)
      window_start = opts.startdate;
   end
   if isfield(opts, 'enddate') && ~isnat(opts.enddate)
      window_end = opts.enddate;
   end
end

function description = describeUncovered(label, uncovered_times)
   %DESCRIBEUNCOVERED Compact account: channel, count, first..last times.
   description = string(label) + ": " + numel(uncovered_times) ...
      + " uncovered (" + string(min(uncovered_times)) + " .. " ...
      + string(max(uncovered_times)) + ")";
end

function channels = icemodelRequiredChannels()
   %ICEMODELREQUIREDCHANNELS POLICY A5 runtime forcing channel set.
   %
   % The set itself lives once in the reconstruct namespace SSOT so the
   % runtime gate and the ledger default can never silently diverge; the
   % gate keeps this thin accessor because callers must not tune it the
   % way reconstruct.setopts required_channels can be tuned per product.
   channels = icemodel.forcing.reconstruct.icemodelRequiredChannels();
end

function tf = requiresSnowfallForcing(smbmodel)
   %REQUIRESSNOWFALLFORCING True when the smbmodel consumes precip mass.
   %
   % POLICY A5: ready_icemodel requires the seven forcing channels; a snow
   % model additionally requires finite total ppt OR snowf. The two
   % historical models keep the D-0b zero-rain contract and consume no
   % precipitation mass, so any other smbmodel is treated as
   % snowfall-consuming and gates on the wider set.
   tf = ~ismember(lower(string(smbmodel)), ["icemodel", "skinmodel"]);
end

function verifyProducerManifest(filename, site, readiness_file, met_files)
   %VERIFYPRODUCERMANIFEST Match runtime inputs to transaction-stamped hashes.
   if ~isfile(filename)
      error('icemodel:loadmet:promiceFilledIdentityMismatch', ...
         'PROMICE producer manifest is unavailable: %s', filename);
   end
   try
      manifest = jsondecode(fileread(filename));
   catch
      error('icemodel:loadmet:promiceFilledIdentityMismatch', ...
         'PROMICE producer manifest is unreadable: %s', filename);
   end
   valid_manifest = isstruct(manifest) && isscalar(manifest) ...
       && isfield(manifest, 'site') && isfield(manifest, 'artifacts') ...
       && isfield(manifest, 'path_base') ...
       && string(manifest.path_base) == "selected_data_root" ...
       && strcmpi(string(manifest.site), site) ...
       && isstruct(manifest.artifacts);
   if ~valid_manifest
      error('icemodel:loadmet:promiceFilledIdentityMismatch', ...
         'PROMICE producer manifest has invalid identity: %s', filename);
   end

   artifacts = manifest.artifacts;
   required = ["role", "path", "bytes", "sha256"];
   if isempty(artifacts) || ~all(isfield(artifacts, required))
      error('icemodel:loadmet:promiceFilledIdentityMismatch', ...
         'PROMICE producer manifest lacks artifact hashes: %s', filename);
   end
   roles = string({artifacts.role});
   relative_paths = string({artifacts.path});
   if isempty(met_files)
      error('icemodel:loadmet:promiceFilledIdentityMismatch', ...
         'PROMICE producer manifest cannot resolve an empty met selection');
   end
   [data_root, ~] = icemodel.forcing.reconstruct.selectedDataRoot( ...
      string(fileparts(met_files(1))));
   paths = strings(size(relative_paths));
   for k = 1:numel(relative_paths)
      relative_file = java.io.File(char(relative_paths(k)));
      candidate = string(fullfile(data_root, relative_paths(k)));
      if relative_file.isAbsolute() || ~icemodel.internal.isPathInside(candidate, data_root)
         error('icemodel:loadmet:promiceFilledIdentityMismatch', ...
            'PROMICE producer manifest path escapes its selected root: %s', ...
            relative_paths(k));
      end
      paths(k) = canonicalPath(candidate);
   end
   expected_paths = [string(readiness_file); met_files];
   for k = 1:numel(expected_paths)
      expected_paths(k) = canonicalPath(expected_paths(k));
   end
   expected_roles = ["readiness"; repmat("filled", numel(met_files), 1)];
   for k = 1:numel(expected_paths)
      match = roles(:) == expected_roles(k) & paths(:) == expected_paths(k);
      if nnz(match) ~= 1 || ~isfile(expected_paths(k))
         error('icemodel:loadmet:promiceFilledIdentityMismatch', ...
            'runtime artifact is not uniquely pinned: %s', expected_paths(k));
      end
      record = artifacts(match);
      info = dir(expected_paths(k));
      digest = icemodel.verification.setup.fileSha256(expected_paths(k));
      if info.bytes ~= double(record.bytes) ...
            || digest ~= lower(string(record.sha256))
         error('icemodel:loadmet:promiceFilledIdentityMismatch', ...
            'runtime artifact differs from producer manifest: %s', ...
            expected_paths(k));
      end
   end

function pathname = canonicalPath(pathname)
   %CANONICALPATH Resolve one artifact path for identity comparison.
   pathname = string(java.io.File(char(pathname)).getCanonicalPath());
end

end
