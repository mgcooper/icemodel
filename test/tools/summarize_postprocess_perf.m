function report = summarize_postprocess_perf(kwargs)
   %SUMMARIZE_POSTPROCESS_PERF Compare hourly RETIME implementations.
   %
   %  report = summarize_postprocess_perf()
   %  report = summarize_postprocess_perf(output_file="/tmp/postprocess_perf.mat")

   arguments
      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      kwargs.output_file (1, :) string = ""
   end

   % Install the formal test config once so the fixture helpers resolve the
   % same environment used by the accepted suite.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Build one representative quarter-hour surface-output timetable.
   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(kwargs.simyear, ...
      configure=true, nsteps=96, dt_seconds=900);
   cleanup = onCleanup(@() icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      workspace));
   opts = icemodel.test.helpers.buildSyntheticOpts(workspace, ...
      'skinmodel', kwargs.simyear, dt=900, solver=1, ...
      testname='postprocess_perf_summary');
   [ice1_raw, ~, opts] = icemodel.test.helpers.runSmbModel(opts);
   met = icemodel.loadmet(opts);
   ice1_tt = rawIce1ToTimetable(ice1_raw, met.Time);
   ice1_tt = repeatQuarterHourTimetable(ice1_tt, 35040);

   % Measure the legacy timetable RETIME branch against the fixed-step path.
   legacy_s = timeit(@() legacyHourlyMean(ice1_tt));
   fixed_s = timeit(@() icemodel.retimeHourlyFixedStep(ice1_tt));

   % Confirm the fixed-step helper reproduces the legacy output exactly on
   % the aligned 15-minute timetable used by postprocess.
   legacy = legacyHourlyMean(ice1_tt);
   fixed = icemodel.retimeHourlyFixedStep(ice1_tt);

   report = struct();
   report.timing = table(["retime_legacy"; "retime_fixed_step"], ...
      [legacy_s; fixed_s], ...
      'VariableNames', {'variant', 'seconds_per_call'});
   report.timing.speedup_vs_legacy = legacy_s ./ report.timing.seconds_per_call;
   report.accuracy = table("retime_fixed_step", ...
      max(abs(double(fixed{:,:}) - double(legacy{:,:})), [], 'all'), ...
      'VariableNames', {'variant', 'max_abs_err'});
   report.timestamp_utc = datetime('now', 'TimeZone', 'UTC');

   if ~isblanktext(kwargs.output_file)
      outdir = fileparts(char(kwargs.output_file));
      if ~isempty(outdir) && exist(outdir, 'dir') ~= 7
         mkdir(outdir);
      end
      save(char(kwargs.output_file), 'report');
   end

   disp('Postprocess timing summary:')
   disp(report.timing)
   disp('Postprocess accuracy summary:')
   disp(report.accuracy)
   clear cleanup
end

function ice1_tt = rawIce1ToTimetable(ice1_raw, time)
   %RAWICE1TOTIMETABLE Convert raw model output to POSTPROCESS timetable form.

   if isfield(ice1_raw, 'Tice_converged')
      ice1_raw.Tice_converged = single(ice1_raw.Tice_converged);
   end
   if isfield(ice1_raw, 'Tsfc_converged')
      ice1_raw.Tsfc_converged = single(ice1_raw.Tsfc_converged);
   end

   time.TimeZone = 'UTC';
   ice1_tt = struct2table(ice1_raw);
   ice1_tt = table2timetable(ice1_tt, 'RowTimes', time);
end

function TT = repeatQuarterHourTimetable(TT, n_rows)
   %REPEATQUARTERHOURTIMETABLE Tile a short timetable to a target row count.

   vars = TT.Properties.VariableNames;
   data = TT{:, vars};
   n_repeat = ceil(n_rows / size(data, 1));
   data = repmat(data, n_repeat, 1);
   data = data(1:n_rows, :);

   time = datetime(2016, 1, 1, 0, 0, 0, 'TimeZone', 'UTC') + ...
      minutes(15) * transpose(0:n_rows-1);
   TT = array2timetable(data, 'RowTimes', time, 'VariableNames', vars);
end

function ice1_hourly = legacyHourlyMean(ice1_tt)
   %LEGACYHOURLYMEAN Reproduce the timetable RETIME branch from POSTPROCESS.

   ice1_hourly = retime(ice1_tt, 'hourly', 'mean');
   ice1_hourly = ice1_hourly(~(month(ice1_hourly.Time) == 2 ...
      & day(ice1_hourly.Time) == 29), :);
end
