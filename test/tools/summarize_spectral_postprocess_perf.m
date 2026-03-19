function report = summarize_spectral_postprocess_perf(kwargs)
   %SUMMARIZE_SPECTRAL_POSTPROCESS_PERF Compare spectral and retime variants.
   %
   %  report = summarize_spectral_postprocess_perf()
   %  report = summarize_spectral_postprocess_perf(output_file="/tmp/perf.mat")
   %
   % This diagnostic tool uses TIMEIT to compare the current inline spectral
   % path, an exact helper-call decomposition, an approximate lookup-table
   % variant, and the hourly postprocess RETIME hotspot.

   arguments
      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      kwargs.output_file (1, :) string = ""
   end

   % Install the canonical test config once so the synthetic fixtures and
   % namespaced helpers resolve the same environment as the formal suite.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Time the spectral alternatives on one shared daytime synthetic state.
   spectral = summarizeSpectralPerf(kwargs.simyear);

   % Time the hourly postprocess hotspot on a representative full-year-sized
   % surface-output timetable.
   postprocess = summarizePostprocessPerf(kwargs.simyear);

   % Return the combined report and optionally save it for later inspection.
   report = struct();
   report.spectral = spectral;
   report.postprocess = postprocess;
   report.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   report.matlab_version = string(version);
   report.host = string(computer);

   if ~isblanktext(kwargs.output_file)
      outdir = fileparts(char(kwargs.output_file));
      if ~isempty(outdir) && exist(outdir, 'dir') ~= 7
         mkdir(outdir);
      end
      save(char(kwargs.output_file), 'report');
   end

   % Print the compact timing tables so interactive use mirrors the report.
   disp('Spectral timing summary:')
   disp(report.spectral.timing)
   disp('Spectral accuracy summary:')
   disp(report.spectral.accuracy)
   disp('Postprocess timing summary:')
   disp(report.postprocess.timing)
   disp('Postprocess accuracy summary:')
   disp(report.postprocess.accuracy)
end

function report = summarizeSpectralPerf(simyear)
   %SUMMARIZESPECTRALPERF Measure exact and approximate spectral variants.

   % Build one daytime spectral state so all spectral variants see the same
   % geometry, forcing, and optical coefficients.
   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(simyear, ...
      configure=true, nsteps=96, dt_seconds=900);
   cleanup = onCleanup(@() icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      workspace));
   state = icemodel.test.fixtures.makeSyntheticColumnState(workspace, ...
      'icemodel', solver=3, metstep=49, include_spectral=true, ...
      testname='spectral_postprocess_summary');

   % Precompute the shared spectral transforms and lookup table once.
   ro_sno = state.ro_ice * state.f_ice + state.ro_liq * state.f_liq + ...
      state.ro_air * (1 - state.f_ice - state.f_liq);
   grid_thermal = cumsum(state.dz) - state.dz / 2;
   grid_spectral = cumsum(state.dz_spect) - state.dz_spect / 2;
   ro_sno_spect = max(GRIDFORWARD(ro_sno, grid_thermal, grid_spectral), 300);
   [~, ~, ~, z_walls] = CVMESH(state.opts.z0_spectral, state.opts.dz_spectral);
   bulk_lookup = MAKEBULKEXTCOEFSLOOKUP(state.dz_spect, state.spect_N, ...
      state.spect_S, state.solardwavl);

   % Measure the exact bulk transform, the lookup approximation, and the
   % full UPDATEEXTCOEFS variants as per-call timings.
   bulk_exact_s = timeit(@() repeatCall(@() BULKEXTCOEFS(state.dz_spect, ...
      ro_sno_spect, state.spect_N, state.spect_S, state.solardwavl), 32)) / 32;
   bulk_lookup_s = timeit(@() repeatCall(@() BULKEXTCOEFSLOOKUP(ro_sno_spect, ...
      bulk_lookup), 256)) / 256;
   inline_s = timeit(@() repeatCall(@() UPDATEEXTCOEFSINLINELEGACY( ...
      state.swd, ...
      state.albedo, state.Q0, state.dz_spect, state.spect_N, ...
      state.spect_S, state.solardwavl, state.Sc, state.dz, ro_sno), 8)) / 8;
   function_exact_s = timeit(@() repeatCall(@() UPDATEEXTCOEFSDECOMPOSED( ...
      state.swd, state.albedo, state.Q0, state.dz_spect, state.spect_N, ...
      state.spect_S, state.solardwavl, state.Sc, state.dz, ro_sno), ...
      8)) / 8;
   function_exact_cached_s = timeit(@() repeatCall(@() ...
      UPDATEEXTCOEFSDECOMPOSEDCACHED(state.swd, state.albedo, state.Q0, ...
      state.dz_spect, state.spect_N, state.spect_S, state.solardwavl, ...
      state.Sc, state.dz, ro_sno, grid_thermal, grid_spectral, z_walls), ...
      8)) / 8;
   function_lookup_s = timeit(@() repeatCall(@() UPDATEEXTCOEFSLOOKUP( ...
      state.swd, state.albedo, state.Q0, state.dz_spect, state.Sc, ...
      state.dz, ro_sno, grid_thermal, grid_spectral, z_walls, ...
      bulk_lookup), 8)) / 8;

   % Compute exact and approximate agreement metrics alongside the timings.
   bulk_exact = BULKEXTCOEFS(state.dz_spect, ro_sno_spect, state.spect_N, ...
      state.spect_S, state.solardwavl);
   bulk_approx = BULKEXTCOEFSLOOKUP(ro_sno_spect, bulk_lookup);
   [Sc_inline, chi_inline] = UPDATEEXTCOEFSINLINELEGACY(state.swd, ...
      state.albedo, ...
      state.Q0, state.dz_spect, state.spect_N, state.spect_S, ...
      state.solardwavl, state.Sc, state.dz, ro_sno);
   [Sc_exact, chi_exact] = UPDATEEXTCOEFSDECOMPOSED(state.swd, ...
      state.albedo, state.Q0, state.dz_spect, state.spect_N, ...
      state.spect_S, state.solardwavl, state.Sc, state.dz, ro_sno);
   [Sc_exact_cached, chi_exact_cached] = UPDATEEXTCOEFSDECOMPOSEDCACHED( ...
      state.swd, state.albedo, state.Q0, state.dz_spect, state.spect_N, ...
      state.spect_S, state.solardwavl, state.Sc, state.dz, ro_sno, ...
      grid_thermal, grid_spectral, z_walls);
   [Sc_lookup, chi_lookup] = UPDATEEXTCOEFSLOOKUP(state.swd, ...
      state.albedo, state.Q0, state.dz_spect, state.Sc, state.dz, ro_sno, ...
      grid_thermal, grid_spectral, z_walls, bulk_lookup);

   timing = table( ...
      ["bulk_exact"; "bulk_lookup"; "update_inline"; ...
      "update_function_exact"; "update_function_exact_cached"; ...
      "update_function_lookup"], ...
      [bulk_exact_s; bulk_lookup_s; inline_s; function_exact_s; ...
      function_exact_cached_s; function_lookup_s], ...
      'VariableNames', {'variant', 'seconds_per_call'});
   timing.speedup_vs_inline = inline_s ./ timing.seconds_per_call;

   accuracy = table( ...
      ["bulk_lookup"; "update_function_exact"; ...
      "update_function_exact_cached"; "update_function_lookup"], ...
      [max(abs(bulk_approx - bulk_exact) ./ max(abs(bulk_exact), 1e-12)); ...
      max(abs(Sc_exact - Sc_inline) ./ max(abs(Sc_inline), 1e-12)); ...
      max(abs(Sc_exact_cached - Sc_inline) ./ max(abs(Sc_inline), 1e-12)); ...
      max(abs(Sc_lookup - Sc_inline) ./ max(abs(Sc_inline), 1e-12))], ...
      [NaN; abs(chi_exact - chi_inline); ...
      abs(chi_exact_cached - chi_inline); abs(chi_lookup - chi_inline)], ...
      'VariableNames', {'variant', 'max_rel_err', 'chi_abs_err'});

   report = struct();
   report.timing = timing;
   report.accuracy = accuracy;
   clear cleanup
end

function report = summarizePostprocessPerf(simyear)
   %SUMMARIZEPOSTPROCESSPERF Measure the hourly postprocess retime hotspot.

   % Build one representative quarter-hour surface-output timetable.
   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(simyear, ...
      configure=true, nsteps=96, dt_seconds=900);
   cleanup = onCleanup(@() icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      workspace));
   opts = icemodel.test.helpers.buildSyntheticOpts(workspace, ...
      'skinmodel', simyear, dt=900, solver=1, ...
      testname='postprocess_perf_summary');
   [ice1_raw, ~, opts] = icemodel.test.helpers.runSmbModel(opts);
   met = icemodel.loadmet(opts);
   ice1_tt = rawIce1ToTimetable(ice1_raw, met.Time);
   ice1_tt = repeatQuarterHourTimetable(ice1_tt, 35040);

   % Measure the legacy timetable RETIME branch against the fixed-step path.
   legacy_s = timeit(@() legacyHourlyMean(ice1_tt));
   fixed_s = timeit(@() icemodel.internal.retimeHourlyFixedStep(ice1_tt));

   % Confirm the fixed-step helper reproduces the legacy output within the
   % single-precision rounding noise introduced by the raw model output.
   legacy = legacyHourlyMean(ice1_tt);
   fixed = icemodel.internal.retimeHourlyFixedStep(ice1_tt);
   timing = table(["retime_legacy"; "retime_fixed_step"], ...
      [legacy_s; fixed_s], 'VariableNames', {'variant', 'seconds_per_call'});
   timing.speedup_vs_legacy = legacy_s ./ timing.seconds_per_call;

   accuracy = table("retime_fixed_step", ...
      max(abs(double(fixed{:,:}) - double(legacy{:,:})), [], 'all'), ...
      'VariableNames', {'variant', 'max_abs_err'});

   report = struct();
   report.timing = timing;
   report.accuracy = accuracy;
   clear cleanup
end

function repeatCall(func, n_reps)
   %REPEATCALL Run one timing target enough times for a stable TIMEIT sample.

   % Repeat a function handle locally so TIMEIT can compare per-call cost on
   % small kernels without the measurement being dominated by timer overhead.
   for i = 1:n_reps
      func();
   end
end

function ice1_tt = rawIce1ToTimetable(ice1_raw, time)
   %RAWICE1TOTIMETABLE Convert raw model output to POSTPROCESS timetable form.

   % Match POSTPROCESS by converting logical convergence flags to single
   % before the timetable conversion.
   if isfield(ice1_raw, 'Tice_converged')
      ice1_raw.Tice_converged = single(ice1_raw.Tice_converged);
   end
   if isfield(ice1_raw, 'Tsfc_converged')
      ice1_raw.Tsfc_converged = single(ice1_raw.Tsfc_converged);
   end

   % Build the timetable shape timed by the hourly aggregation hotspot.
   time.TimeZone = 'UTC';
   ice1_tt = struct2table(ice1_raw);
   ice1_tt = table2timetable(ice1_tt, 'RowTimes', time);
end

function TT = repeatQuarterHourTimetable(TT, n_rows)
   %REPEATQUARTERHOURTIMETABLE Tile a short timetable to a target row count.

   % Repeat the source data enough times to hit the requested representative
   % benchmark length, then assign one continuous quarter-hour time axis.
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

   % Match the legacy path exactly, including the leap-day removal.
   ice1_hourly = retime(ice1_tt, 'hourly', 'mean');
   ice1_hourly = ice1_hourly(~(month(ice1_hourly.Time) == 2 ...
      & day(ice1_hourly.Time) == 29), :);
end
