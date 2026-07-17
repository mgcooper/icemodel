function report = writePromiceSnowModelReadyYears(kwargs)
   %WRITEPROMICESNOWMODELREADYYEARS Write the annual PROMICE snow handoff.
   %
   %  report = icemodel.verification.setup.writePromiceSnowModelReadyYears( ...
   %     evaluation_data_root=eval_root, input_data_root=input_root, ...
   %     output_dir=qa_root)
   %
   % Audits complete Gregorian years declared by the PROMICE family manifest.
   % Immediately usable forcing requires a complete exact-grid MAR 3.11 met
   % year, including precipitation. Practical evaluation requires an exact
   % hourly PROMICE snow-depth grid with at least 95 percent finite samples and
   % one finite sample every day; the strict subset requires every sample.
   % PROMICE met plus MAR precipitation is retained as a separate diagnostic.
   %
   % Inputs
   %  evaluation_data_root  Root containing promice/manifest.json and targets.
   %  input_data_root       Root containing met/<source> staged forcing files.
   %  output_dir            Directory receiving CSV, Markdown, and JSON output.
   %
   % Outputs
   %  report  Struct containing coverage/ready tables, the self-check, and paths.

   arguments
      kwargs.evaluation_data_root (1, 1) string
      kwargs.input_data_root (1, 1) string
      kwargs.output_dir (1, 1) string
   end

   % Keep every path explicit so candidate and canonical stages use the same
   % implementation without embedding a machine-specific workspace location.
   manifest_file = fullfile(kwargs.evaluation_data_root, ...
      "promice", "manifest.json");
   if ~isfile(manifest_file)
      error('icemodel:verification:promiceReadyYears:manifestMissing', ...
         'PROMICE manifest does not exist: %s', manifest_file)
   end
   if ~isfolder(kwargs.output_dir)
      mkdir(kwargs.output_dir)
   end
   manifest = icemodel.verification.helpers.readFamilyManifest(manifest_file);

   % Process one case at a time so multi-decade met artifacts are released
   % before the next station is loaded.
   rows = repmat(rowTemplate(), 0, 1);
   sites_without_rcm = strings(0, 1);
   for icase = 1:numel(manifest.cases)
      c = manifest.cases(icase);
      [promice_met, promice_exists] = loadMet(c, "promice", ...
         kwargs.input_data_root);
      [mar_met, mar_exists] = loadMet(c, "mar", kwargs.input_data_root);
      observations = loadObservations(c, manifest.family_root);
      if ~mar_exists
         sites_without_rcm(end + 1, 1) = string(c.site_id); %#ok<AGROW>
      end

      % Audit only complete years inside the declared MAR leg: precipitation is
      % mandatory, so years outside that source window are not candidates.
      years = candidateYears(c, mar_exists);
      for y = reshape(years, 1, [])
         row = annualRow(c, y, promice_met, mar_met, observations, ...
            promice_exists, mar_exists);
         rows(end + 1, 1) = row; %#ok<AGROW>
      end
      clear promice_met mar_met observations
   end

   % Stable station/year ordering makes CSV bytes and grouped Markdown
   % deterministic even if a future manifest writer reorders case records.
   coverage = struct2table(rows);
   coverage = sortrows(coverage, {'site_id', 'year'});
   ready = coverage(coverage.practical_ready, :);
   sites_without_rcm = sort(unique(sites_without_rcm));

   coverage_file = fullfile(kwargs.output_dir, ...
      "promice_snow_model_site_year_coverage.csv");
   ready_file = fullfile(kwargs.output_dir, ...
      "promice_snow_model_ready_site_years.csv");
   summary_file = fullfile(kwargs.output_dir, ...
      "promice_snow_model_ready_site_years.md");
   check_file = fullfile(kwargs.output_dir, ...
      "promice_snow_model_ready_site_years_check.json");

   % CSV is the machine-readable contract; Markdown and JSON summarize exactly
   % those written rows instead of maintaining separate scientific logic.
   writetable(coverage, coverage_file)
   writetable(ready, ready_file)
   writeText(summary_file, summaryMarkdown(coverage, ready, ...
      sites_without_rcm));
   check = buildCheck(coverage, ready, sites_without_rcm, coverage_file, ...
      ready_file, summary_file);
   writeText(check_file, string(jsonencode(check, PrettyPrint=true)));

   report = struct('coverage', coverage, 'ready', ready, 'check', check, ...
      'files', struct('coverage_csv', string(coverage_file), ...
      'ready_csv', string(ready_file), 'summary_md', string(summary_file), ...
      'check_json', string(check_file)));
end

function years = candidateYears(c, mar_exists)
   %CANDIDATEYEARS Intersect the case period with its mandatory MAR leg.
   if ~mar_exists || ~isfield(c.colocation.mar, 'window')
      years = zeros(1, 0);
      return
   end
   case_start = parseUtc(c.period.start);
   case_end = parseUtc(c.period.end);
   mar_start = parseUtc(c.colocation.mar.window.start);
   mar_end = parseUtc(c.colocation.mar.window.end);
   period = struct('start', formatUtc(max(case_start, mar_start)), ...
      'end', formatUtc(min(case_end, mar_end)));
   years = fullCalendarYears(period);
end

function text = formatUtc(value)
   %FORMATUTC Encode one UTC datetime for reuse by the period helper.
   text = string(value, 'yyyy-MM-dd HH:mm:ss');
end

%% Local functions
function row = annualRow(c, y, promice_met, mar_met, observations, ...
      promice_exists, mar_exists)
   %ANNUALROW Derive one annual readiness record from staged artifacts.
   met_time = calendarGrid(y, minutes(15));
   snow_time = calendarGrid(y, hours(1));
   required = requiredPromiceChannels();

   % Required-channel missing counts include absent coordinate samples, not
   % merely explicit NaNs in the rows that happen to remain in a file.
   promice_grid = exactGrid(promice_met, met_time);
   promice_missing = zeros(numel(required), 1);
   for k = 1:numel(required)
      promice_missing(k) = missingCount(promice_met, required(k), met_time);
   end
   promice_complete = promice_grid && all(promice_missing == 0);

   mar_required = [required; "ppt"];
   mar_grid = exactGrid(mar_met, met_time);
   mar_missing = zeros(numel(mar_required), 1);
   for k = 1:numel(mar_required)
      mar_missing(k) = missingCount(mar_met, mar_required(k), met_time);
   end
   ppt_missing = mar_missing(end);
   mar_precip_complete = mar_exists && mar_grid && ppt_missing == 0;
   mar_complete = mar_exists && mar_grid && all(mar_missing == 0);

   % Snow depth is evaluated independently of met coverage so exclusions remain
   % diagnostic: callers can distinguish unavailable evaluation from forcing.
   snow_present = ismember("snow_depth", ...
      string(observations.Properties.VariableNames));
   snow_grid = exactGrid(observations, snow_time);
   snow_finite = finiteMask(observations, "snow_depth", snow_time);
   snow_finite_samples = nnz(snow_finite);
   snow_expected_samples = numel(snow_time);
   snow_fraction = snow_finite_samples / snow_expected_samples;
   snow_expected_days = daysInYear(y);
   snow_finite_days = nnz(any(reshape(snow_finite, 24, []), 1));
   snow_daily_fraction = snow_finite_days / snow_expected_days;
   snow_practical = snow_present && snow_grid && snow_fraction >= 0.95 ...
      && snow_finite_days == snow_expected_days;
   snow_strict = snow_present && snow_grid ...
      && snow_finite_samples == snow_expected_samples;

   hybrid_forcing = promice_exists && promice_complete ...
      && mar_precip_complete;
   forcing_ready = mar_complete;
   strict_ready = forcing_ready && snow_strict;
   practical_ready = forcing_ready && snow_practical;

   row = rowTemplate();
   row.site_id = string(c.site_id);
   row.case_id = string(c.case_id);
   row.surface_zone = string(c.surface_zone);
   row.year = y;
   row.expected_met_samples = numel(met_time);
   row.promice_time_grid_complete = promice_grid;
   row.promice_required_complete = promice_complete;
   row.rcm_time_grid_complete = mar_grid;
   row.rcm_ppt_missing = ppt_missing;
   row.rcm_precip_complete = mar_precip_complete;
   row.rcm_required_missing_total = sum(mar_missing);
   row.rcm_full_forcing_complete = mar_complete;
   row.promice_hybrid_forcing_ready = hybrid_forcing;
   row.promice_hybrid_practical_ready = hybrid_forcing && snow_practical;
   row.forcing_ready = forcing_ready;
   row.snow_depth_present = snow_present;
   row.snow_time_grid_complete = snow_grid;
   row.snow_finite_samples = snow_finite_samples;
   row.snow_expected_samples = snow_expected_samples;
   row.snow_finite_fraction = snow_fraction;
   row.snow_finite_days = snow_finite_days;
   row.snow_expected_days = snow_expected_days;
   row.snow_daily_presence_fraction = snow_daily_fraction;
   row.strict_ready = strict_ready;
   row.practical_ready = practical_ready;
   row.promice_hybrid_exclusion_reason = hybridReason(promice_exists, ...
      promice_grid, promice_missing, mar_exists, mar_grid, ppt_missing);
   row.exclusion_reason = readinessReason(mar_exists, mar_grid, ...
      mar_missing, snow_present, snow_grid, snow_fraction, ...
      snow_finite_days, snow_expected_days);

   % Preserve the seven per-channel PROMICE counts used to diagnose whether the
   % optional hybrid lane can ever become useful.
   for k = 1:numel(required)
      row.(char(required(k) + "_missing")) = promice_missing(k);
   end
end

function row = rowTemplate()
   %ROWTEMPLATE Define the stable public CSV schema and column order.
   row = struct( ...
      'site_id', "", 'case_id', "", 'surface_zone', "", 'year', 0, ...
      'rcm_precip_source', "mar3.11", 'expected_met_samples', 0, ...
      'promice_time_grid_complete', false, ...
      'promice_required_complete', false, ...
      'rcm_time_grid_complete', false, 'rcm_ppt_missing', 0, ...
      'rcm_precip_complete', false, 'rcm_required_missing_total', 0, ...
      'rcm_full_forcing_complete', false, ...
      'promice_hybrid_forcing_ready', false, ...
      'promice_hybrid_practical_ready', false, 'forcing_ready', false, ...
      'snow_depth_present', false, 'snow_time_grid_complete', false, ...
      'snow_finite_samples', 0, 'snow_expected_samples', 0, ...
      'snow_finite_fraction', 0, 'snow_finite_days', 0, ...
      'snow_expected_days', 0, 'snow_daily_presence_fraction', 0, ...
      'strict_ready', false, 'practical_ready', false, ...
      'promice_hybrid_exclusion_reason', "", 'exclusion_reason', "", ...
      'tair_missing', 0, 'swd_missing', 0, 'lwd_missing', 0, ...
      'albedo_missing', 0, 'wspd_missing', 0, 'rh_missing', 0, ...
      'psfc_missing', 0);
end

function channels = requiredPromiceChannels()
   %REQUIREDPROMICECHANNELS Return non-precipitation model forcing channels.
   channels = ["tair"; "swd"; "lwd"; "albedo"; "wspd"; "rh"; "psfc"];
end

function [met, exists] = loadMet(c, source, input_root)
   %LOADMET Load one manifest-referenced met timetable, if staged.
   met = timetable();
   exists = false;
   source = char(source);
   if ~isfield(c, 'colocation') || ~isfield(c.colocation, source)
      return
   end
   leg = c.colocation.(source);
   if ~isfield(leg, 'met_files') || isempty(leg.met_files)
      return
   end
   met_files = reshape(string(leg.met_files), [], 1);
   met_file = fullfile(input_root, "met", met_files(1));
   if ~isfile(met_file)
      return
   end

   % Load only the canonical payload and reject malformed referenced bytes.
   saved = load(met_file, 'met');
   if ~isfield(saved, 'met') || ~istimetable(saved.met)
      error('icemodel:verification:promiceReadyYears:badMet', ...
         'Referenced met artifact lacks timetable named met: %s', met_file)
   end
   met = ensureUtcTime(saved.met);
   exists = true;
end

function observations = loadObservations(c, family_root)
   %LOADOBSERVATIONS Load the PROMICE target timetable when it exists.
   observations = timetable();
   if ~isfield(c, 'evaluation_file') || strlength(c.evaluation_file) == 0
      return
   end
   evaluation_file = fullfile(family_root, c.evaluation_file);
   if ~isfile(evaluation_file)
      return
   end
   targets = icemodel.verification.helpers.loadArtifact( ...
      evaluation_file, "targets");
   if ~isstruct(targets) || ~isfield(targets, 'data') ...
         || ~istimetable(targets.data)
      error('icemodel:verification:promiceReadyYears:badObservations', ...
         'PROMICE evaluation artifact lacks targets.data timetable: %s', ...
         evaluation_file)
   end
   observations = ensureUtcTime(targets.data);
end

function tt = ensureUtcTime(tt)
   %ENSUREUTCTIME Normalize staged row times without shifting naive values.
   tt.Time.TimeZone = 'UTC';
end

function years = fullCalendarYears(period)
   %FULLCALENDARYEARS Return complete hourly calendar years in a case period.
   start_time = parseUtc(period.start);
   end_time = parseUtc(period.end);
   first_year = year(start_time);
   if start_time > datetime(first_year, 1, 1, 'TimeZone', 'UTC')
      first_year = first_year + 1;
   end
   last_year = year(end_time);
   last_required = datetime(last_year, 12, 31, 23, 0, 0, ...
      'TimeZone', 'UTC');
   if end_time < last_required
      last_year = last_year - 1;
   end
   years = first_year:last_year;
end

function value = parseUtc(text)
   %PARSEUTC Decode the manifest's portable second-resolution datetime text.
   value = datetime(string(text), 'InputFormat', 'yyyy-MM-dd HH:mm:ss', ...
      'TimeZone', 'UTC');
end

function time = calendarGrid(y, step)
   %CALENDARGRID Build one exact UTC calendar-year coordinate.
   start_time = datetime(y, 1, 1, 'TimeZone', 'UTC');
   end_time = datetime(y + 1, 1, 1, 'TimeZone', 'UTC') - step;
   time = (start_time:step:end_time)';
end

function tf = exactGrid(tt, expected_time)
   %EXACTGRID Test ordered equality for one year without interpolation.
   if height(tt) == 0
      tf = false;
      return
   end
   y = year(expected_time(1));
   actual_time = tt.Time(year(tt.Time) == y);
   tf = isequal(actual_time(:), expected_time(:));
end

function count = missingCount(tt, variable, expected_time)
   %MISSINGCOUNT Count absent postings and nonfinite values on the exact grid.
   count = numel(expected_time) - nnz(finiteMask(tt, variable, expected_time));
end

function mask = finiteMask(tt, variable, expected_time)
   %FINITEMASK Align finite support to an expected coordinate without filling.
   mask = false(numel(expected_time), 1);
   if height(tt) == 0 ...
         || ~ismember(variable, string(tt.Properties.VariableNames))
      return
   end
   [present, location] = ismember(expected_time, tt.Time);
   values = tt.(char(variable));
   row_finite = all(isfinite(values), 2);
   mask(present) = row_finite(location(present));
end

function n = daysInYear(y)
   %DAYSINYEAR Return the Gregorian calendar-day count.
   n = day(datetime(y, 12, 31), 'dayofyear');
end

function reason = hybridReason(promice_exists, promice_grid, missing, ...
      mar_exists, mar_grid, ppt_missing)
   %HYBRIDREASON Explain why PROMICE met plus MAR precipitation is unusable.
   parts = strings(0, 1);
   if ~promice_exists
      parts(end + 1) = "PROMICE met leg absent";
   elseif ~promice_grid
      parts(end + 1) = "PROMICE time grid incomplete";
   end
   channels = requiredPromiceChannels();
   for k = 1:numel(channels)
      if missing(k) > 0
         parts(end + 1) = sprintf('%s missing=%d', ...
            channels(k), missing(k)); %#ok<AGROW>
      end
   end
   if ~mar_exists
      parts(end + 1) = "MAR precipitation leg absent";
   elseif ~mar_grid
      parts(end + 1) = "MAR time grid incomplete";
   elseif ppt_missing > 0
      parts(end + 1) = sprintf('MAR ppt missing=%d', ...
         ppt_missing);
   end
   reason = strjoin(parts, '; ');
end

function reason = readinessReason(mar_exists, mar_grid, mar_missing, ...
      snow_present, snow_grid, snow_fraction, snow_days, expected_days)
   %READINESSREASON Explain each failed practical-readiness criterion.
   parts = strings(0, 1);
   if ~mar_exists
      parts(end + 1) = "MAR 3.11 met leg absent";
   elseif ~mar_grid
      parts(end + 1) = "MAR 3.11 time grid incomplete";
   end
   if any(mar_missing > 0)
      parts(end + 1) = sprintf('MAR 3.11 required missing=%d', ...
         sum(mar_missing));
   end
   if ~snow_present
      parts(end + 1) = "PROMICE snow_depth absent";
   elseif ~snow_grid
      parts(end + 1) = "PROMICE snow_depth time grid incomplete";
   elseif snow_fraction < 0.95
      parts(end + 1) = "snow_depth finite fraction <95%";
   elseif snow_days < expected_days
      parts(end + 1) = "snow_depth lacks daily presence";
   end
   reason = strjoin(parts, '; ');
end

function lines = summaryMarkdown(coverage, ready, sites_without_rcm)
   %SUMMARYMARKDOWN Render the concise human-readable annual handoff.
   strict = ready(ready.strict_ready, :);
   ready_sites = unique(ready.site_id, 'sorted');
   strict_sites = unique(strict.site_id, 'sorted');
   hybrid_count = nnz(coverage.promice_hybrid_practical_ready);
   lines = [ ...
      "# PROMICE snow-model-ready site-years"; ""; ...
      "Generated from `promice/manifest.json` under the supplied " + ...
      "evaluation-data root."; ""; "## Criteria"; ""; ...
      "- **Immediately usable forcing:** complete exact-grid MAR 3.11 " + ...
      "15-minute met, including precipitation and all required channels."; ...
      "- **PROMICE-hybrid audit:** complete PROMICE non-precipitation met " + ...
      "plus complete MAR precipitation; reported separately."; ...
      "- **Practical evaluation:** exact hourly PROMICE `snow_depth`, at " + ...
      "least 95% finite samples, and one finite sample every day."; ...
      "- **Strict subset:** the same forcing criteria and 100% finite " + ...
      "hourly `snow_depth`."; ""; "## Result"; ""; ...
      sprintf('- Candidate full site-years audited: **%d**', height(coverage)); ...
      sprintf('- Sites without a staged RCM precipitation leg: **%d** (%s)', ...
      numel(sites_without_rcm), strjoin(sites_without_rcm, ', ')); ...
      sprintf('- Recommended practical site-years: **%d** across **%d** sites', ...
      height(ready), numel(ready_sites)); ...
      sprintf('- Strict 100%% snow-depth site-years: **%d** across **%d** sites', ...
      height(strict), numel(strict_sites)); ...
      sprintf('- PROMICE-hybrid practical site-years: **%d**', hybrid_count); ...
      ""; "## Recommended ranges"; ""; ...
      "| Site | Practical years | Strict 100% snow-depth years |"; ...
      "|---|---|---|"];

   % Render one stable station row with compact contiguous ranges.
   for site = reshape(ready_sites, 1, [])
      practical_years = ready.year(ready.site_id == site);
      strict_years = strict.year(strict.site_id == site);
      lines(end + 1) = sprintf('| %s | %s | %s |', site, ...
         compactYears(practical_years), compactYears(strict_years)); %#ok<AGROW>
   end
   lines(end + 1:end + 3) = [""; ...
      "Detailed exclusions and channel counts: " + ...
      "`promice_snow_model_site_year_coverage.csv`."; ""];
end

function text = compactYears(years)
   %COMPACTYEARS Render sorted years as deterministic contiguous ranges.
   years = sort(unique(reshape(years, 1, [])));
   if isempty(years)
      text = "none";
      return
   end
   parts = strings(0, 1);
   first = years(1);
   last = first;
   for k = 2:numel(years)
      if years(k) == last + 1
         last = years(k);
      else
         parts(end + 1) = yearRange(first, last); %#ok<AGROW>
         first = years(k);
         last = first;
      end
   end
   parts(end + 1) = yearRange(first, last);
   text = strjoin(parts, ', ');
end

function text = yearRange(first, last)
   %YEARRANGE Render one year or one inclusive year range.
   if first == last
      text = string(first);
   else
      text = sprintf('%d-%d', first, last);
   end
end

function check = buildCheck(coverage, ready, sites_without_rcm, ...
      coverage_file, ready_file, summary_file)
   %BUILDCHECK Build portable provenance and verify the written deliverables.
   strict = ready(ready.strict_ready, :);
   check = struct();
   check.manifest = "promice/manifest.json";
   check.required_promice_channels = requiredPromiceChannels();
   check.rcm_precip_source = "mar3.11";
   check.recommended_forcing_source = "mar3.11 full met";
   check.practical_snow_minimum_fraction = 0.95;
   check.practical_snow_requires_every_day = true;
   check.candidate_site_year_count = height(coverage);
   check.ready_site_year_count = height(ready);
   check.ready_site_count = numel(unique(ready.site_id));
   check.strict_site_year_count = height(strict);
   check.strict_site_count = numel(unique(strict.site_id));
   check.promice_hybrid_ready_site_year_count = ...
      nnz(coverage.promice_hybrid_practical_ready);
   check.sites_without_rcm_precip_leg = sites_without_rcm;
   check.coverage_csv = "promice_snow_model_site_year_coverage.csv";
   check.ready_csv = "promice_snow_model_ready_site_years.csv";
   check.summary_md = "promice_snow_model_ready_site_years.md";
   check.coverage_csv_sha256 = ...
      icemodel.verification.setup.fileSha256(coverage_file);
   check.ready_csv_sha256 = ...
      icemodel.verification.setup.fileSha256(ready_file);
   check.summary_md_sha256 = ...
      icemodel.verification.setup.fileSha256(summary_file);

   % The self-check proves the ready CSV is exactly the practical subset and
   % every advertised hash names one written artifact.
   check.all_checks_passed = isequal(ready, ...
      coverage(coverage.practical_ready, :)) ...
      && all(strlength([check.coverage_csv_sha256; ...
      check.ready_csv_sha256; check.summary_md_sha256]) == 64);
end

function writeText(pathname, lines)
   %WRITETEXT Write UTF-8 text with stable Unix newlines.
   fid = fopen(pathname, 'w', 'n', 'UTF-8');
   if fid < 0
      error('icemodel:verification:promiceReadyYears:writeFailed', ...
         'Could not open output file: %s', pathname)
   end
   cleanup = onCleanup(@() fclose(fid));
   for line = reshape(string(lines), 1, [])
      fprintf(fid, '%s\n', line);
   end
   clear cleanup
end
