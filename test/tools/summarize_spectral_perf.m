function report = summarize_spectral_perf(kwargs)
   %SUMMARIZE_SPECTRAL_PERF Compare the spectral source-term variants.
   %
   %  report = summarize_spectral_perf()
   %  report = summarize_spectral_perf(output_file="/tmp/spectral_perf.mat")
   %
   % This diagnostic reports:
   %  1. kernel timings for the inlined, exact, and lookup paths
   %  2. direct whole-model timings for exact vs lookup
   %  3. agreement metrics against the historical inlined path using the same
   %     scalar summary semantics as the formal regression suite
   %
   % The kernel section preserves all three variants (inlined, exact, lookup)
   % because those functions are called directly. The direct-model section
   % compares only exact vs lookup via opts.lookup_k_bulk.

   arguments
      kwargs.simyear (1, 1) double {mustBeInteger, mustBePositive} ...
         = 2016

      kwargs.output_file (1, :) string ...
         = ""

      kwargs.smoke_site (1, :) string ...
         = "kanm"

      kwargs.solver (1, 1) double {mustBeMember(kwargs.solver, [1 2 3])} ...
         = 2

      kwargs.include_full_model (1, 1) logical ...
         = true

      kwargs.include_direct_model (1, 1) logical ...
         = true

      kwargs.n_direct_runs (1, 1) double {mustBeInteger, mustBePositive} = 1
   end

   % Install the formal test config once so the synthetic fixtures and perf
   % runner share the same environment as the accepted test suite.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Measure the narrow kernel path first because it isolates the spectral
   % transforms from the rest of the model runtime.
   kernels = summarizeSpectralKernelPerf(kwargs.simyear);

   % Direct TIC/TOC measurement wraps only the model call.
   direct_model = struct();
   if kwargs.include_full_model && kwargs.include_direct_model
      direct_model = summarizeSpectralDirectModelPerf(kwargs.simyear, ...
         kwargs.smoke_site, kwargs.solver, kwargs.n_direct_runs);
   end

   % Return the report and optionally save it for later comparison.
   report = struct();
   report.kernels = kernels;
   report.direct_model = direct_model;
   report.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   report.matlab_version = string(version);
   report.host = string(computer);

   % Save the file.
   if ~isblanktext(kwargs.output_file)
      outdir = fileparts(char(kwargs.output_file));
      if ~isempty(outdir) && exist(outdir, 'dir') ~= 7
         mkdir(outdir);
      end
      save(char(kwargs.output_file), 'report');
   end

   % Print compact tables so interactive use mirrors the saved report without
   % flooding the command window with fields that are only useful in the MAT
   % artifact.
   disp('Spectral kernel timing summary:')
   disp(report.kernels.timing(:, {'variant', 'seconds_per_call', ...
      'ref_variant', 'speedup_vs_ref'}))
   disp('Spectral kernel accuracy summary:')
   disp(report.kernels.accuracy)
   if isfield(report.direct_model, 'timing')
      disp('Spectral direct-model timing summary:')
      disp(report.direct_model.timing(:, {'variant', 'median_wall_s', ...
         'speedup_vs_exact'}))
      disp('Spectral direct-model output agreement:')
      disp(report.direct_model.accuracy)
   end
end

function report = summarizeSpectralKernelPerf(simyear)
   %SUMMARIZESPECTRALKERNELPERF Measure exact and approximate kernel variants.

   % Configure the synthetic workspace
   workspace = icemodel.test.fixtures.makeSyntheticWorkspace(simyear, ...
      configure=true, nsteps=96, dt_seconds=900);
   cleanup = onCleanup(@() icemodel.test.fixtures.cleanupSyntheticWorkspace( ...
      workspace));

   % Build one daytime spectral state so all variants see the same geometry,
   % forcing, and optical coefficients.
   state = icemodel.test.fixtures.makeSyntheticColumnState(workspace, ...
      'icemodel', solver=3, metstep=49, include_spectral=true, ...
      testname='spectral_summary');

   % Precompute the remapped spectral density and lookup table once.
   ro_sno = state.ro_ice * state.f_ice + state.ro_liq * state.f_liq + ...
      state.ro_air * (1 - state.f_ice - state.f_liq);

   ro_sno_spect = max(interp1(state.z_nodes, ro_sno, ...
      state.z_nodes_spect, 'nearest', 'extrap'), 300);

   k_lookup = icemodel.radiation.make_bulk_extinction_lookup(state.dz_spect, ...
      state.tau_N, state.tau_S, state.solar_dwavel);

   % Empty lookup table for the exact path
   k_lookup_empty = struct([]);

   % Measure the exact bulk transform, the lookup approximation, and the full
   % spectral source-term variants as per-call timings.
   bulk_exact_s = timeit(@() repeatCall(@() ...
      icemodel.radiation.bulk_extinction_coefficients(state.dz_spect, ...
      ro_sno_spect, state.tau_N, state.tau_S, state.solar_dwavel), 32)) / 32;

   bulk_lookup_s = timeit(@() repeatCall(@() ...
      icemodel.radiation.bulk_extinction_coefficients( ...
      state.dz_spect, ro_sno_spect, state.tau_N, state.tau_S, ...
      state.solar_dwavel, k_lookup), 256)) / 256;

   inlined_s = timeit(@() repeatCall(@() SPECTRALSOURCETERM_INLINE( ...
      state.swd, state.albedo, state.I0, state.dz_spect, state.tau_N, ...
      state.tau_S, state.solar_dwavel, state.dz, ro_sno, state.z_nodes, ...
      state.z_nodes_spect), 8)) / 8;

   functions_s = timeit(@() repeatCall(@() ...
      icemodel.column.shortwave_source_term(state.swd, state.albedo, ...
      state.I0, state.dz_spect, state.tau_N, state.tau_S, ...
      state.solar_dwavel, state.dz, state.z_nodes, state.z_nodes_spect, ...
      state.z_edges_spect, ro_sno, k_lookup_empty), 8)) / 8;

   lookup_s = timeit(@() repeatCall(@() ...
      icemodel.column.shortwave_source_term(state.swd, state.albedo, ...
      state.I0, state.dz_spect, state.tau_N, state.tau_S, ...
      state.solar_dwavel, state.dz, state.z_nodes, state.z_nodes_spect, ...
      state.z_edges_spect, ro_sno, k_lookup), 8)) / 8;

   % Compute exact and approximate agreement metrics alongside the timings.
   bulk_exact = icemodel.radiation.bulk_extinction_coefficients(state.dz_spect, ...
      ro_sno_spect, state.tau_N, state.tau_S, state.solar_dwavel);

   bulk_lookup = icemodel.radiation.bulk_extinction_coefficients(state.dz_spect, ...
      ro_sno_spect, state.tau_N, state.tau_S, state.solar_dwavel, k_lookup);

   [Sc_inlined, chi_inlined] = SPECTRALSOURCETERM_INLINE( ...
      state.swd, state.albedo, state.I0, state.dz_spect, state.tau_N, ...
      state.tau_S, state.solar_dwavel, state.dz, ro_sno, state.z_nodes, ...
      state.z_nodes_spect);

   [Sc_functions, chi_functions] = icemodel.column.shortwave_source_term( ...
      state.swd, state.albedo, state.I0, state.dz_spect, state.tau_N, ...
      state.tau_S, state.solar_dwavel, state.dz, state.z_nodes, ...
      state.z_nodes_spect, state.z_edges_spect, ro_sno, k_lookup_empty);

   [Sc_lookup, chi_lookup] = icemodel.column.shortwave_source_term( ...
      state.swd, state.albedo, state.I0, state.dz_spect, state.tau_N, ...
      state.tau_S, state.solar_dwavel, state.dz, state.z_nodes, ...
      state.z_nodes_spect, state.z_edges_spect, ro_sno, k_lookup);

   timing = table( ...
      ["bulk_exact"; "bulk_lookup"; "source_inlined"; ...
      "source_functions"; "source_lookup"], ...
      [bulk_exact_s; bulk_lookup_s; inlined_s; functions_s; lookup_s], ...
      ["bulk_exact"; "bulk_exact"; "source_inlined"; ...
      "source_inlined"; "source_inlined"], ...
      'VariableNames', {'variant', 'seconds_per_call', 'ref_variant'});

   timing.speedup_vs_ref = [ ...
      bulk_exact_s / bulk_exact_s; ...
      bulk_exact_s / bulk_lookup_s; ...
      inlined_s / inlined_s; ...
      inlined_s / functions_s; ...
      inlined_s / lookup_s];

   accuracy = table( ...
      ["bulk_lookup"; "source_functions"; "source_lookup"], ...
      [max(abs(bulk_lookup - bulk_exact) ./ max(abs(bulk_exact), 1e-12)); ...
      max(abs(Sc_functions - Sc_inlined) ./ max(abs(Sc_inlined), 1e-12)); ...
      max(abs(Sc_lookup - Sc_inlined) ./ max(abs(Sc_inlined), 1e-12))], ...
      [NaN; abs(chi_functions - chi_inlined); abs(chi_lookup - chi_inlined)], ...
      'VariableNames', {'variant', 'max_rel_err', 'chi_abs_err'});

   report = struct();
   report.timing = timing;
   report.accuracy = accuracy;
   report.measurement_note = [ ...
      "bulk_* rows measure only the bulk-extinction transform. ", ...
      "source_* rows measure the full spectral source-term path."];
   clear cleanup
end

function report = summarizeSpectralDirectModelPerf(simyear, smoke_site, ...
      solver, n_runs)
   %SUMMARIZESPECTRALDIRECTMODELPERF Time exact vs lookup whole-model runs.

   % Build one formal smoke case and reuse it across both variants.
   cases = icemodel.test.helpers.getPerfCaseMatrix(tier="smoke", ...
      smbmodel="icemodel", solver=solver, simyear=simyear, ...
      smoke_sites=smoke_site);
   if height(cases) ~= 1
      error('spectral summary expected one smoke perf case, found %d', ...
         height(cases))
   end
   c = cases(1, :);
   variants = ["exact"; "lookup"];
   outputs = cell(numel(variants), 1);
   wall_s = nan(numel(variants), n_runs);

   % Time each variant directly around RUNSMBMODEL so the summary exposes the
   % wall-clock behavior without any perf-framework sampling wrapper.
   for i = 1:numel(variants)
      opts = icemodel.test.helpers.setModelOptsForCase(c);
      opts.lookup_k_bulk = (variants(i) == "lookup");
      opts.testname = "direct_spectral_" + variants(i);

      icemodel.test.helpers.runSmbModel(opts);

      for j = 1:n_runs
         t0 = tic;
         [ice1, ice2] = icemodel.test.helpers.runSmbModel(opts);
         wall_s(i, j) = toc(t0);
         if j == 1
            outputs{i} = struct('ice1', ice1, 'ice2', ice2);
         end
      end
   end

   timing = table(variants, median(wall_s, 2), mean(wall_s, 2), ...
      min(wall_s, [], 2), max(wall_s, [], 2), ...
      'VariableNames', {'variant', 'median_wall_s', 'mean_wall_s', ...
      'min_wall_s', 'max_wall_s'});
   timing.speedup_vs_exact = timing.median_wall_s(1) ./ timing.median_wall_s;

   accuracy = compareVariantOutputs(outputs, variants);

   report = struct();
   report.timing = timing;
   report.accuracy = accuracy;
   report.case_id = string(c.case_id);
   report.n_runs = n_runs;
   report.measurement_note = [ ...
      "Direct TIC/TOC around RUNSMBMODEL for one formal smoke ICEMODEL ", ...
      "case. This includes the full model runtime but no perf-framework ", ...
      "sampling."];
end

function T = compareVariantOutputs(outputs, variants)
   %COMPAREVARIANTOUTPUTS Summarize full-model output differences by variant.

   % Summarize the direct outputs with the same scalar metric helper used by
   % the formal regression suite so the study stays aligned with accepted
   % report semantics.
   summaries = cellfun(@(out) icemodel.test.helpers.summarizeIce1Metrics( ...
      out.ice1), outputs, 'UniformOutput', false);
   ref = summaries{1};

   pct_mean_T = nan(numel(variants), 1);
   pct_seb_rmse = nan(numel(variants), 1);
   max_Tice_out = nan(numel(variants), 1);
   n_not_converged = nan(numel(variants), 1);
   for i = 1:numel(variants)
      S = summaries{i};
      pct_mean_T(i) = pctDelta(S.mean_Tice_numiter, ref.mean_Tice_numiter);
      pct_seb_rmse(i) = pctDelta(S.closure_seb_rmse, ...
         ref.closure_seb_rmse);
      max_Tice_out(i) = S.max_Tice_numiter;
      n_not_converged(i) = S.n_not_converged;
   end

   T = table(variants, pct_mean_T, pct_seb_rmse, max_Tice_out, ...
      n_not_converged, ...
      'VariableNames', {'variant', 'pct_dif_mean_Tice_numiter', ...
      'pct_dif_seb_rmse', 'max_Tice_numiter', ...
      'n_not_converged'});
end

function y = pctDelta(x, ref)
   %PCTDELTA Return percent change relative to one reference value.

   if abs(ref) <= 1e-12
      y = NaN;
   else
      y = 100 * (x - ref) / ref;
   end
end

function repeatCall(fcn, n_repeat)
   %REPEATCALL Run one function handle repeatedly for TIMEIT batching.

   for i = 1:n_repeat
      fcn();
   end
end
