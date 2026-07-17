function filenames = writemet(met, site, forcings, kwargs)
   %WRITEMET Validate and save an icemodel met file.
   %
   %  filenames = icemodel.forcing.helpers.writemet(met, site, forcings)
   %  filenames = ... writemet(_, outdir=..., naming="yearly", ...
   %     dt_out="15m", overwrite=true, validate=false)
   %
   % Saves MET (a timetable satisfying the icemodel met contract) under
   % the standard naming convention so icemodel.createMetFileNames and
   % icemodel.loadmet resolve it without special cases:
   %
   %  naming="window" (default): one file spanning the full time axis,
   %     met_<site>_<forcings>_<YYYYMMDD>_<YYYYMMDD>_<dt>.mat
   %
   %  naming="yearly": one file per calendar year (legacy form),
   %     met_<site>_<forcings>_<YYYY>_<dt>.mat
   %
   % OUTDIR defaults to fullfile(icemodel.getpath('input'), 'met'), the
   % met directory of the active icemodel workspace (demo/data when the
   % demo or test config is active). The directory is created when it
   % does not exist. The saved .mat file holds one variable named met. Model met
   % defaults to 15-minute output across repository writers. Pass dt_out="" to
   % preserve an explicit native cadence. Existing targets are additive no-ops
   % unless overwrite=true; explicit replacement emits a warning. A successful
   % wider window write removes only strictly contained files for the same
   % site/source/cadence class and warns with the exact paths removed.
   % Yearly naming validates/resamples the full source before calendar slicing,
   % then applies compact exact per-year support/count summaries. A guarded
   % 15-minute input must carry those summaries and constant cadence blocks.
   %
   % Inputs
   %  met      - timetable of forcing variables (see
   %             icemodel.forcing.helpers.metvariables)
   %  site     - site name encoded in the filename (e.g. "kanm")
   %  forcings - forcing source encoded in the filename (e.g. "mar")
   %
   % Outputs
   %  filenames - string column of selected full paths (including reused files)
   %
   % See also: icemodel.forcing.helpers.metfilename,
   %  icemodel.forcing.helpers.validatemet, icemodel.loadmet

   arguments
      met timetable
      site (1, 1) string
      forcings (1, 1) string
      kwargs.outdir (1, 1) string = ""
      kwargs.naming (1, 1) string {mustBeMember(kwargs.naming, ...
         ["window", "yearly"])} = "window"
      kwargs.dt_out (1, 1) string {mustBeMember(kwargs.dt_out, ...
         ["", "15m"])} = "15m"
      kwargs.overwrite (1, 1) logical = false
      kwargs.validate (1, 1) logical = true
   end

   % Validate and derive the complete output before any directory mutation. This
   % keeps rejected input side-effect free for both window and yearly naming.
   met = prepareMet(met, kwargs.dt_out, kwargs.validate);
   dt = seconds(met.Time(2) - met.Time(1));
   guarded = false;
   if kwargs.naming == "yearly"
      % Guarded yearly provenance is another validation boundary and must fail
      % before the writer creates its per-source directory.
      metadata = met.Properties.UserData;
      guarded = isstruct(metadata) ...
         && isfield(metadata, 'met_resample_policy');
      if guarded
         validateGuardedYearlyMet(met)
      end
   end

   outdir = kwargs.outdir;
   if outdir == ""
      outdir = string(fullfile(icemodel.getpath('input'), 'met'));
   end
   % Stage into the per-source subfolder met/<forcings>/ so the flat met/ folder
   % does not sprawl as verification staging grows; the runtime resolves this
   % subfolder first (icemodel.configureRun / createMetFileNames).
   outdir = fullfile(outdir, char(forcings));
   if ~isfolder(outdir)
      mkdir(outdir)
   end

   switch kwargs.naming
      case "window"
         % Centralize model-met cadence at the final shared writer. Userdata
         % uses its separate hourly writer boundary, and outages remain missing.
         identity_matches = @(filename) ...
            icemodel.forcing.helpers.artifactIdentityMatches( ...
            filename, met, "met");
         cadence_matches = @(filename) ...
            icemodel.forcing.helpers.artifactCadenceMatches( ...
            filename, "met", dt);
         candidate_matches = @(filename) cadence_matches(filename) ...
            && identity_matches(filename);
         name = icemodel.forcing.helpers.metfilename(site, forcings, ...
            met.Time(1), met.Time(end), dt);
         filenames = fullfile(outdir, name);

         % Runtime gives an exact window name precedence over broader files, so
         % validate and return that same artifact before considering enclosure.
         % A conflicting exact file requires explicit overwrite even when a
         % compatible broad artifact also exists.
         if isfile(filenames) && ~kwargs.overwrite
            if ~cadence_matches(filenames)
               cadenceConflict(filenames, "window")
            end
            if ~identity_matches(filenames)
               identityConflict(filenames, "window")
            end
            return
         end

         % With no exact file, a broader current window already satisfies a
         % narrower ordinary request without creating a duplicate artifact.
         if ~kwargs.overwrite
            prefix = "met_" + site + "_" + forcings;
            suffix = metSuffix(dt);
            enclosing = ...
               icemodel.forcing.helpers.findEnclosingWindowFile( ...
               outdir, prefix, suffix, met.Time(1), met.Time(end), ...
               accept_candidate=candidate_matches);
            % Return the same path normal runtime resolution would select.
            if enclosing ~= ""
               filenames = fullfile(outdir, enclosing);
               return
            end
         end
         wrote = savemet(filenames, met, kwargs.overwrite);
         if wrote
            % A successful wider refresh supersedes only contained shorter
            % windows for this exact site/source/cadence naming class.
             icemodel.forcing.helpers.pruneSupersededWindowFiles( ...
                filenames, "met_" + site + "_" + forcings, metSuffix(dt), ...
                accept_candidate=candidate_matches);
         end

      case "yearly"
         % The full source was validated above before slicing, so a gap or
         % irregular step crossing Jan 1 remains visible in the yearly files.
         identity_matches = @(filename) ...
            icemodel.forcing.helpers.artifactIdentityMatches( ...
            filename, met, "met");
         cadence_matches = @(filename) ...
            icemodel.forcing.helpers.artifactCadenceMatches( ...
            filename, "met", dt);
         years_present = unique(year(met.Time));
         filenames = strings(numel(years_present), 1);
         for n = 1:numel(years_present)
            yyyy = years_present(n);
            year_met = met(year(met.Time) == yyyy, :);
            if guarded
               year_met = localizeYearlyProvenance(year_met, yyyy);
            end
            name = icemodel.forcing.helpers.metfilename(site, forcings, ...
               yyyy, [], dt);
            filenames(n) = fullfile(outdir, name);
            if isfile(filenames(n)) && ~kwargs.overwrite
               if ~cadence_matches(filenames(n))
                  cadenceConflict(filenames(n), "year")
               end
               if ~identity_matches(filenames(n))
                  identityConflict(filenames(n), "year")
               end
            end
            savemet(filenames(n), year_met, kwargs.overwrite)
         end
   end
end

function cadenceConflict(filename, span)
   %CADENCECONFLICT Reject an exact filename whose saved time axis disagrees.
   error('icemodel:forcing:writemet:cadenceConflict', ...
      ['Existing model-met artifact %s has a different or unprovable saved ' ...
      'cadence. Pass overwrite=true to replace that exact %s.'], filename, span)
end

function identityConflict(filename, span)
   %IDENTITYCONFLICT Surface an exact-name collision with other provenance.
   error('icemodel:forcing:writemet:identityConflict', ...
      ['Existing model-met artifact %s has conflicting source, product, ' ...
      'schema, sampling-method, or point metadata. Pass overwrite=true ' ...
      'to replace that exact %s.'], filename, span)
end

%% Local functions
function met = prepareMet(met, dt_out, validate)
   %PREPAREMET Resample, validate, and stamp one independently saved timetable.
   met = icemodel.forcing.helpers.resampleMetTimestep(met, dt_out);
   if validate
      icemodel.forcing.helpers.validatemet(met)
   end

   % Stamp after validation so caller-provided unit mismatches stay visible.
   met = icemodel.forcing.helpers.stampMetadata(met, strict=false);
end

function validateGuardedYearlyMet(met)
   %VALIDATEGUARDEDYEARLYMET Prove guarded support before calendar slicing.
   metadata = met.Properties.UserData;
   required = ["met_resample_policy", ...
      "met_resample_source_cadence_seconds", ...
      "met_resample_time_semantics", ...
      "met_resample_support_start_inclusive", ...
      "met_resample_support_end_exclusive", ...
      "met_resample_source_gap_intervals", ...
      "met_resample_yearly_summaries"];
   if ~all(isfield(metadata, required))
      error('icemodel:forcing:writemet:guardedYearlyProvenanceMissing', ...
         'guarded 15-minute met lacks exact yearly resample summaries')
   end

   cadence_s = double(metadata.met_resample_source_cadence_seconds);
   rows_per_interval = cadence_s / 900;
   allowed = ["interval_start_zero_order_hold", "native_15m_unchanged"];
   valid_layout = ismember(string(metadata.met_resample_policy), allowed) ...
      && string(metadata.met_resample_time_semantics) == "interval_start" ...
      && isscalar(cadence_s) && isfinite(cadence_s) && cadence_s >= 900 ...
      && abs(rows_per_interval - round(rows_per_interval)) <= 1e-9 ...
      && height(met) > 1 && all(diff(met.Time) == minutes(15)) ...
      && mod(height(met), round(rows_per_interval)) == 0 ...
      && isequal(metadata.met_resample_support_start_inclusive, met.Time(1)) ...
      && isequal(metadata.met_resample_support_end_exclusive, ...
      met.Time(end) + minutes(15));
   if ~valid_layout
      error('icemodel:forcing:writemet:guardedYearlyProvenanceInvalid', ...
         'guarded met has invalid support or source-cadence provenance')
   end
   rows_per_interval = round(rows_per_interval);

   % Every recorded source interval (and every omitted all-NaN interval) must
   % remain constant across its full 15-minute support block.
   first_rows = 1:rows_per_interval:height(met);
   for name = string(met.Properties.VariableNames)
      values = met.(char(name));
      first_values = values(first_rows, :);
      for offset = 2:rows_per_interval
         if ~isequaln(values(first_rows + offset - 1, :), first_values)
            error('icemodel:forcing:writemet:guardedSourceNotConstant', ...
               ['guarded met variable %s is not constant within recorded ' ...
               'source-cadence blocks'], name)
         end
      end
   end

   % Gap intervals are the compact distinction between omitted all-NaN blocks
   % and real all-NaN source rows. Validate them and every source-side aggregate
   % before any per-year summary is copied into an artifact.
   validateGuardedSourceSummaries(met, metadata, cadence_s, first_rows)
end

function validateGuardedSourceSummaries(met, metadata, cadence_s, first_rows)
   %VALIDATEGUARDEDSOURCESUMMARIES Prove source/gap counts from guarded blocks.
   gaps = metadata.met_resample_source_gap_intervals;
   if ~isstruct(gaps) || (~isempty(gaps) ...
         && ~all(isfield(gaps, ["start", "end"])))
      error('icemodel:forcing:writemet:guardedYearlyProvenanceInvalid', ...
         'guarded met source-gap intervals are malformed')
   end

   block_times = met.Time(first_rows);
   block_support_end = block_times + seconds(cadence_s);
   gap_blocks = false(numel(block_times), 1);
   for n = 1:numel(gaps)
      start_time = gaps(n).start;
      end_time = gaps(n).end;
      aligned = isdatetime(start_time) && isscalar(start_time) ...
         && ~isnat(start_time) && isdatetime(end_time) ...
         && isscalar(end_time) && ~isnat(end_time) ...
         && start_time < end_time ...
         && abs(seconds(start_time - met.Time(1)) / cadence_s ...
         - round(seconds(start_time - met.Time(1)) / cadence_s)) <= 1e-9 ...
         && abs(seconds(end_time - met.Time(1)) / cadence_s ...
         - round(seconds(end_time - met.Time(1)) / cadence_s)) <= 1e-9;
      if ~aligned
         error('icemodel:forcing:writemet:guardedYearlyProvenanceInvalid', ...
            'guarded met source-gap intervals are not cadence aligned')
      end
      in_gap = block_times >= start_time & block_times < end_time;
      output_gap = met.Time >= start_time & met.Time < end_time;
      if ~any(in_gap) || ~allNumericMissing(met(output_gap, :))
         error('icemodel:forcing:writemet:guardedYearlyProvenanceInvalid', ...
            'guarded met source-gap interval is not explicit missing support')
      end
      gap_blocks = gap_blocks | in_gap;
   end

   source_blocks = ~gap_blocks;
   top_valid = metadata.met_resample_source_row_count == nnz(source_blocks) ...
      && metadata.met_resample_source_time_gap_count == numel(gaps) ...
      && isequaln(metadata.met_resample_source_missing_counts, ...
      missingCounts(met(first_rows(source_blocks), :))) ...
      && isequaln(metadata.met_resample_expected_missing_counts, ...
      missingCounts(met));
   if ~top_valid
      error('icemodel:forcing:writemet:guardedYearlyProvenanceInvalid', ...
         'guarded met full-window source/gap counts are inconsistent')
   end

   summaries = metadata.met_resample_yearly_summaries;
   required = ["year", "source_row_count", "source_time_gap_count", ...
      "source_missing_counts", "source_gap_intervals", ...
      "expected_missing_counts", "support_start_inclusive", ...
      "support_end_exclusive"];
   if ~isstruct(summaries) || ~all(isfield(summaries, required))
      error('icemodel:forcing:writemet:guardedYearlyProvenanceInvalid', ...
         'guarded met yearly summaries are malformed')
   end
   for n = 1:numel(summaries)
      summary = summaries(n);
      output = met(year(met.Time) == summary.year, :);
      if isempty(output)
         error('icemodel:forcing:writemet:guardedYearlyProvenanceInvalid', ...
            'guarded met summary year is absent from the timetable')
      end
      support_start = output.Time(1);
      support_end = output.Time(end) + minutes(15);
      overlaps = block_times < support_end ...
         & block_support_end > support_start;
      source_overlaps = overlaps & source_blocks;
      gap_overlaps = false(numel(gaps), 1);
      for k = 1:numel(gaps)
         gap_overlaps(k) = gaps(k).start < support_end ...
            && gaps(k).end > support_start;
      end
      valid = isequal(summary.support_start_inclusive, support_start) ...
         && isequal(summary.support_end_exclusive, support_end) ...
         && summary.source_row_count == nnz(source_overlaps) ...
         && summary.source_time_gap_count == nnz(gap_overlaps) ...
         && isequaln(summary.source_missing_counts, ...
         missingCounts(met(first_rows(source_overlaps), :))) ...
         && isequaln(summary.source_gap_intervals, gaps(gap_overlaps)) ...
         && isequaln(summary.expected_missing_counts, missingCounts(output));
      if ~valid
         error('icemodel:forcing:writemet:guardedYearlyProvenanceInvalid', ...
            'guarded met summary source/gap counts are inconsistent')
      end
   end
end

function tf = allNumericMissing(T)
   %ALLNUMERICMISSING True when every numeric channel is non-finite.
   tf = true;
   for name = string(T.Properties.VariableNames)
      values = T.(char(name));
      if isnumeric(values) && any(isfinite(values), 'all')
         tf = false;
         return
      end
   end
end

function year_met = localizeYearlyProvenance(year_met, yyyy)
   %LOCALIZEYEARLYPROVENANCE Replace full-window facts with one exact summary.
   metadata = year_met.Properties.UserData;
   summaries = metadata.met_resample_yearly_summaries;
   required = ["year", "source_row_count", "source_time_gap_count", ...
      "source_missing_counts", "source_gap_intervals", ...
      "expected_missing_counts", ...
      "support_start_inclusive", "support_end_exclusive"];
   if ~isstruct(summaries) || ~all(isfield(summaries, required))
      error('icemodel:forcing:writemet:guardedYearlyProvenanceInvalid', ...
         'guarded met yearly summaries are malformed')
   end
   hit = find([summaries.year] == yyyy);
   if numel(hit) ~= 1
      error('icemodel:forcing:writemet:guardedYearlyProvenanceInvalid', ...
         'guarded met must carry exactly one summary for year %d', yyyy)
   end
   summary = summaries(hit);
   actual_missing = missingCounts(year_met);
   valid_slice = isequal(summary.support_start_inclusive, year_met.Time(1)) ...
      && isequal(summary.support_end_exclusive, ...
      year_met.Time(end) + minutes(15)) ...
      && isequaln(summary.expected_missing_counts, actual_missing);
   if ~valid_slice
      error('icemodel:forcing:writemet:guardedYearlyProvenanceInvalid', ...
         'guarded met summary does not match the saved year %d slice', yyyy)
   end

   metadata.met_resample_source_row_count = summary.source_row_count;
   metadata.met_resample_source_time_gap_count = ...
      summary.source_time_gap_count;
   metadata.met_resample_source_missing_counts = summary.source_missing_counts;
   metadata.met_resample_source_gap_intervals = summary.source_gap_intervals;
   metadata.met_resample_expected_missing_counts = ...
      summary.expected_missing_counts;
   metadata.met_resample_support_start_inclusive = ...
      summary.support_start_inclusive;
   metadata.met_resample_support_end_exclusive = ...
      summary.support_end_exclusive;
   metadata.met_resample_yearly_summaries = summary;
   year_met.Properties.UserData = metadata;
end

function counts = missingCounts(T)
   %MISSINGCOUNTS Count unavailable values in numeric timetable variables.
   counts = struct();
   for name = string(T.Properties.VariableNames)
      values = T.(char(name));
      if isnumeric(values)
         counts.(char(name)) = nnz(~isfinite(values));
      end
   end
end

function suffix = metSuffix(dt)
   %METSUFFIX Filename suffix paired with a validated model-met cadence.
   tag = icemodel.forcing.helpers.metTimestepSuffix(dt);
   suffix = "_" + tag + ".mat";
end

function wrote = savemet(filename, met, overwrite)
   %SAVEMET Save MET to FILENAME as a variable named met.
   exists = isfile(filename);
   wrote = false;
   % Ordinary repeated writes reuse the exact current artifact bytes.
   if exists && ~overwrite
      return
   end
   % Explicit replacement is intentionally visible to setup callers.
   if exists
      warning('icemodel:forcing:writemet:overwrite', ...
         'Replacing existing model-met artifact %s.', filename);
   end
   % Persist one exact provenance record in both supported read locations. The
   % adapter may derive cadence/location facts absent from incoming UserData.
   S.artifact_metadata = icemodel.forcing.helpers.artifactMetadata(met);
   met.Properties.UserData = S.artifact_metadata;
   S.met = met;
   save(filename, '-struct', 'S')
   wrote = true;
end
