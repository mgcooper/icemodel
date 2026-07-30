function tests = test_reconstruct_harness
   %TEST_RECONSTRUCT_HARNESS Verify the held-out validation harness contracts.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install the verification path for namespace resolution.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   testCase.TestData.tmp = tempname;
   mkdir(testCase.TestData.tmp);
end

function teardown(testCase)
   % Remove temporary split manifests.
   if isfolder(testCase.TestData.tmp)
      rmdir(testCase.TestData.tmp, 's')
   end
   clear testCase.TestData.cleanup
end

%% gapCensus

function test_census_counts_runs_and_buckets(testCase)
   % Known inserted gaps produce exact run counts, buckets, and totals.
   series = icemodel.test.fixtures.makeReconstructSeries();
   series.tair(10:15) = NaN;          % 6 h run  -> bucket 1 (<=6 h)
   series.tair(100:123) = NaN;        % 24 h run -> bucket 2 (6-24 h)
   series.lwd(200:399) = NaN;         % 200 h run -> bucket 5 (>168 h)

   returned = icemodel.forcing.reconstruct.gapCensus(series, ...
      channels=["tair", "lwd"]);

   tair_row = returned.summary(returned.summary.channel == "tair", :);
   testCase.verifyEqual(tair_row.n_runs, 2);
   testCase.verifyEqual(tair_row.runs_bucket1, 1);
   testCase.verifyEqual(tair_row.runs_bucket2, 1);
   testCase.verifyEqual(tair_row.n_missing, 30);
   testCase.verifyEqual(tair_row.samples_fixable_within_cap, 6);
   bucket_columns = startsWith( ...
      string(tair_row.Properties.VariableNames), "runs_bucket");
   testCase.verifyEqual(sum(tair_row{1, bucket_columns}), tair_row.n_runs);
   lwd_row = returned.summary(returned.summary.channel == "lwd", :);
   testCase.verifyEqual(lwd_row.runs_bucket5, 1);
   testCase.verifyEqual(lwd_row.longest_run_hours, 200);
   testCase.verifyEqual(height(returned.runs), 3);
   testCase.verifyTrue(all(ismember(returned.runs.season, ...
      ["DJF", "MAM", "JJA", "SON"])));
end

function test_census_uses_configured_interpolation_cap(testCase)
   % Census tier sizing must follow the same central cap as reconstruction.
   series = icemodel.test.fixtures.makeReconstructSeries();
   series.lwd(10:12) = NaN;

   returned = icemodel.forcing.reconstruct.gapCensus(series, ...
      channels="lwd", cap_hours=2);

   testCase.verifyEqual(returned.cap_hours, 2);
   testCase.verifyEqual(returned.summary.samples_fixable_within_cap, 0);
   testCase.verifyEqual( ...
      returned.summary.samples_needing_donor_or_proxy, 3);
end

function test_census_uses_channel_specific_interpolation_cap(testCase)
   % Planning evidence must size albedo with the same D-41 cap as execution.
   series = icemodel.test.fixtures.makeReconstructSeries();
   series.albedo = 0.7 * ones(height(series), 1);
   series.albedo(10:35) = NaN;

   returned = icemodel.forcing.reconstruct.gapCensus(series, ...
      channels="albedo", cap_hours=6, ...
      cap_hours_by_channel=struct('albedo', 30));

   testCase.verifyEqual(returned.cap_hours_by_channel.albedo, 30);
   testCase.verifyEqual( ...
      returned.summary.samples_fixable_within_cap, 26);
   testCase.verifyEqual( ...
      returned.summary.samples_needing_donor_or_proxy, 0);
end

function test_census_bounds_to_observed_record(testCase)
   % Leading spans remain pre-deployment until every core channel is finite.
   series = icemodel.test.fixtures.makeReconstructSeries();
   series.tair(1:50) = NaN;
   series.rh(1:60) = NaN;
   series.psfc(1:70) = NaN;
   series.lwd(1:70) = NaN;

   returned = icemodel.forcing.reconstruct.gapCensus(series, ...
      channels="lwd");

   testCase.verifyEqual(returned.record_start, ...
      series.Properties.RowTimes(71));
   testCase.verifyEqual(returned.summary.n_missing, 0);
end

function test_census_rejects_unknown_channel(testCase)
   % Unknown channels fail loudly instead of censusing nothing.
   testCase.verifyError(@() icemodel.forcing.reconstruct.gapCensus( ...
      icemodel.test.fixtures.makeReconstructSeries(), channels="nope"), ...
      'icemodel:reconstruct:gapCensus:unknownChannel');
end

function test_census_rejects_absent_configured_core_channel(testCase)
   % Joint-core support is impossible when any configured channel is absent.
   series = icemodel.test.fixtures.makeReconstructSeries();
   series.psfc = [];
   testCase.verifyError(@() icemodel.forcing.reconstruct.gapCensus( ...
      series, channels="lwd", core_channels=["tair", "rh", "psfc"]), ...
      'icemodel:reconstruct:gapCensus:noFiniteCore');
end

%% validationSplit

function test_split_is_deterministic_and_disjoint(testCase)
   % Same seed reproduces the split; the partition is disjoint and complete.
   years = 2009:2023;
   a = icemodel.forcing.reconstruct.validationSplit(years, ...
      station="kanm", seed=42);
   b = icemodel.forcing.reconstruct.validationSplit(years, ...
      station="kanm", seed=42);

   testCase.verifyEqual(a, b);
   testCase.verifyEmpty(intersect(a.years_selection, a.years_evaluation));
   testCase.verifyEqual(sort([a.years_selection, a.years_evaluation]), years);
   testCase.verifyEqual(numel(a.years_selection), round(0.7 * numel(years)));
end

function test_split_manifest_persists_and_wins(testCase)
   % A persisted manifest is returned as-is on replay, even with a new seed.
   manifest_file = fullfile(testCase.TestData.tmp, 'kanm-split.json');
   first = icemodel.forcing.reconstruct.validationSplit(2009:2023, ...
      station="kanm", seed=1, manifest_file=manifest_file);

   replay = icemodel.forcing.reconstruct.validationSplit(2009:2023, ...
      station="kanm", seed=999, manifest_file=manifest_file);

   testCase.verifyEqual(replay.years_selection, first.years_selection);
   testCase.verifyEqual(replay.seed, 1);
   testCase.verifyError(@() icemodel.forcing.reconstruct.validationSplit( ...
      2009:2023, station="kanl", seed=1, manifest_file=manifest_file), ...
      'icemodel:reconstruct:validationSplit:stationMismatch');
end

function test_split_manifest_rejects_overlap_and_stale_years(testCase)
   % A persisted replay may win over a new seed, but it cannot overlap its
   % protocol sets or reference a year outside the current station record.
   manifest_file = fullfile(testCase.TestData.tmp, 'strict-split.json');
   original = icemodel.forcing.reconstruct.validationSplit(2009:2023, ...
      station="kanm", seed=1, manifest_file=manifest_file);

   changed = original;
   changed.years_evaluation(1) = changed.years_selection(1);
   writeJson(manifest_file, changed);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.validationSplit(2009:2023, ...
      station="kanm", seed=1, manifest_file=manifest_file), ...
      'icemodel:reconstruct:validationSplit:invalidManifest');

   changed = original;
   changed.years_selection(1) = 1999;
   writeJson(manifest_file, changed);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.validationSplit(2009:2023, ...
      station="kanm", seed=1, manifest_file=manifest_file), ...
      'icemodel:reconstruct:validationSplit:invalidManifest');
end

%% syntheticMissingness

function test_draws_are_blocked_and_deterministic(testCase)
   % Draws are reproducible, land only on finite samples in allowed years,
   % and keep observed context on both sides of every gap.
   series = icemodel.test.fixtures.makeReconstructSeries();
   series_years = year(series.Properties.RowTimes);

   a = icemodel.forcing.reconstruct.syntheticMissingness(series, "tair", ...
      seededRuns(), years=2020, seed=7, n_gaps=5);
   b = icemodel.forcing.reconstruct.syntheticMissingness(series, "tair", ...
      seededRuns(), years=2020, seed=7, n_gaps=5);

   testCase.verifyEqual(a.mask, b.mask);
   testCase.verifyEqual(a.inserted, 5);
   testCase.verifyTrue(all(isfinite(series.tair(a.mask))));
   testCase.verifyTrue(all(series_years(a.mask) == 2020));
   % Context margin: the samples immediately outside each gap stay unmasked.
   for g = 1:height(a.gaps)
      idx = find(series.Properties.RowTimes == a.gaps.start_time(g));
      testCase.verifyFalse(a.mask(idx - 1));
   end
   testCase.verifyEqual(nnz(a.mask), ...
      sum(round(a.gaps.duration_hours)));
   % Gaps arrive in time order: validationMetrics maps time-ordered masked
   % samples to gap rows sequentially, so any other order would misalign.
   testCase.verifyTrue(issorted(a.gaps.start_time));
end

function test_draws_reject_empty_duration_pool(testCase)
   % A stratum with no census runs cannot be drawn from.
   series = icemodel.test.fixtures.makeReconstructSeries();
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.syntheticMissingness(series, "tair", ...
      seededRuns(), years=2020, seed=7, season="XXX"), ...
      'icemodel:reconstruct:syntheticMissingness:emptyDurationPool');
end

function test_draws_stay_in_requested_season(testCase)
   % A season-stratified duration pool must also constrain insertion times.
   series = icemodel.test.fixtures.makeReconstructSeries();
   draws = icemodel.forcing.reconstruct.syntheticMissingness(series, ...
      "tair", seededRuns(), years=2020, seed=9, n_gaps=3, ...
      bucket=1, season="JJA");

   testCase.verifyGreaterThan(draws.inserted, 0);
   testCase.verifyTrue(all(draws.gaps.season == "JJA"));
   testCase.verifyTrue(all( ...
      icemodel.forcing.reconstruct.seasonOf( ...
       series.Properties.RowTimes(draws.mask)) == "JJA"));
end

function test_persistence_masks_every_heldout_draw(testCase)
   % A sample held out by another stratum cannot act as the active gap's
   % persistence anchor; without that second draw the native anchor is valid.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:4)).';
   values = (10:10:50).';
   active_mask = false(5, 1);
   active_mask(3:4) = true;
   gaps = table(times(3), times(4), ...
      'VariableNames', {'start_time', 'end_time'});
   active = struct('mask', active_mask, 'gaps', gaps);
   other_mask = false(5, 1);
   other_mask(2) = true;
   other = struct('mask', other_mask);

   blocked = icemodel.forcing.reconstruct.persistenceEstimate( ...
      values, times, active, {active, [], other});
   native = icemodel.forcing.reconstruct.persistenceEstimate( ...
      values, times, active, {active});

   testCase.verifyTrue(all(isnan(blocked)));
   testCase.verifyEqual(native, [20; 20]);
end

%% validationMetrics and admissionGate

function test_options_reject_inverted_variability_range(testCase)
   % The central options contract rejects an impossible admission interval.
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.setopts( ...
      min_variability_ratio=2, max_variability_ratio=1), ...
      'icemodel:reconstruct:setopts:variabilityRange');
end

function test_options_reject_out_of_range_coverage(testCase)
   % Admission and native-triage fractions are probabilities, not arbitrary
   % finite scalars.
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.setopts(min_coverage=-0.1), ...
      'icemodel:reconstruct:setopts:coverageRange');
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.setopts(min_native_core_coverage=1.1), ...
      'icemodel:reconstruct:setopts:coverageRange');
end

function test_options_reject_cap_above_policy_ceiling(testCase)
   % The ordinary cap stays six hours while evidenced channel overrides remain
   % narrow and independently bounded.
   edges = icemodel.forcing.reconstruct.bucketEdges();
   opts = icemodel.forcing.reconstruct.setopts();
   testCase.verifyEqual(opts.cap_hours, edges(2));
   testCase.verifyEqual(opts.cap_hours_by_channel.swd, 9);
   testCase.verifyEqual(opts.cap_hours_by_channel.rh, 9);
   testCase.verifyEqual(opts.cap_hours_by_channel.albedo, 30);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.setopts(cap_hours=8), ...
        'icemodel:reconstruct:setopts:capHours');
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.setopts( ...
      cap_hours_by_channel=struct('swd', 10)), ...
      'icemodel:reconstruct:setopts:capHours');
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.setopts( ...
      cap_hours_by_channel=struct('rh', 10)), ...
      'icemodel:reconstruct:setopts:capHours');
end

function test_options_reject_required_unplanned_channel(testCase)
   % Every required reconstructed output needs a planned provenance column;
   % adopted total precipitation is the sole policy exception.
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.setopts( ...
      required_channels=["tair", "swu"], plan_channels="tair"), ...
      'icemodel:reconstruct:setopts:requiredChannelNotPlanned');
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels=["tair", "ppt"], plan_channels="tair");
   testCase.verifyEqual(opts.required_channels, ["tair", "ppt"]);
end

function test_options_reject_unsupported_method_schema(testCase)
   % Dedicated precipitation, interpolation, and proxy policies must reject
   % configurations that their execution layers cannot honor.
   testCase.verifyEqual( ...
      icemodel.forcing.helpers.precipitationVariables(), ...
      ["ppt", "rainf", "snowf"]);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.setopts(plan_channels=["tair", "ppt"]), ...
      'icemodel:reconstruct:setopts:plannedPrecipitation');
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.setopts(interp_channels="swu"), ...
      'icemodel:reconstruct:setopts:unsupportedInterpolationChannel');
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.setopts(proxy_sources="racmo"), ...
      'icemodel:reconstruct:setopts:unsupportedProxySource');
   opts = icemodel.forcing.reconstruct.setopts( ...
      proxy_sources=["MAR", "MERRA"]);
   testCase.verifyEqual(opts.proxy_sources, ["mar", "merra"]);
end

function test_duration_buckets_include_exact_upper_boundary(testCase)
   % Exact policy boundaries stay in the shorter right-closed stratum.
   returned = icemodel.forcing.reconstruct.gapDurationBucket( ...
      [0.25, 6, 6 + eps(6), 24, 24 + eps(24)]);

   testCase.verifyEqual(returned, [1, 1, 2, 2, 3]);
end

function test_metrics_score_known_errors(testCase)
   % A constructed reconstruction yields exact bias/RMSE/coverage and flags
   % the deliberate bound violation and boundary jump.
   series = icemodel.test.fixtures.makeReconstructSeries();
   draws = icemodel.forcing.reconstruct.syntheticMissingness(series, ...
      "tair", seededRuns(), years=2020, seed=7, n_gaps=3);
   truth = series.tair(draws.mask);
   filled = truth + 1;                          % constant +1 K error
   filled(1) = 500;                             % out of tair bounds

   metrics = icemodel.forcing.reconstruct.validationMetrics( ...
      truth, filled, draws.gaps, series, "tair");

   testCase.verifyEqual(metrics.overall.n, numel(truth));
   testCase.verifyEqual(metrics.overall.coverage, 1);
   testCase.verifyGreaterThan(metrics.overall.bias, 1 - 1e-6);
   testCase.verifyTrue(isfinite(metrics.overall.correlation));
   testCase.verifyGreaterThan(metrics.overall.variability_ratio, 0);
   testCase.verifyGreaterThan(metrics.overall.typical_magnitude, 0);
   testCase.verifyEqual(metrics.overall.bound_violations, 1);
   % The 500 K first sample is a boundary jump at its gap's leading edge.
   testCase.verifyGreaterThan(metrics.overall.boundary_jump_rate, 0);
    testCase.verifyEqual(sum(metrics.by_stratum.n), metrics.overall.n);
end

function test_metrics_and_gate_require_filled_sample_provenance(testCase)
   % A candidate that leaves any supplied value with the missing code fails
   % the policy's 100% provenance-accounting invariant.
   series = icemodel.test.fixtures.makeReconstructSeries();
   draws = icemodel.forcing.reconstruct.syntheticMissingness(series, ...
      "tair", seededRuns(), years=2020, seed=17, n_gaps=2);
   truth = series.tair(draws.mask);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   provenance = repmat(double(codes.station_transfer), numel(truth), 1);
   provenance(1) = double(codes.missing);

   metrics = icemodel.forcing.reconstruct.validationMetrics( ...
      truth, truth, draws.gaps, series, "tair", ...
      provenance=provenance);
   testCase.verifyLessThan(metrics.overall.provenance_accounting, 1);
   gate = icemodel.forcing.reconstruct.admissionGate( ...
      "tair", metrics.overall, NaN);
   testCase.verifyFalse(gate.admit);
   testCase.verifyTrue(any(contains(gate.reasons, ...
      "provenance accounting")));
end

function test_metrics_fail_closed_without_boundary_anchors(testCase)
   % A record-spanning held-out gap has no measurable native seam and must
   % not be scored as a clean boundary.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   gaps = table(times(1), times(end), hours(times(end) - times(1)), 4, ...
      icemodel.forcing.reconstruct.seasonOf(times(1)), ...
      'VariableNames', {'start_time', 'end_time', 'duration_hours', ...
      'bucket', 'season'});

   metrics = icemodel.forcing.reconstruct.validationMetrics( ...
      series.tair, series.tair, gaps, series, "tair");

   testCase.verifyEqual(metrics.overall.boundary_jump_rate, 1);
end

function test_metrics_context_must_exclude_withheld_truth(testCase)
   % Held-out values can inflate the native step scale enough to hide a bad
   % seam, so planner callers must mask them before requesting metrics.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:9)).';
   native = (250:259).';
   gap = (3:8).';
   gaps = table(times(gap(1)), times(gap(end)), 6, 1, "DJF", ...
      'VariableNames', {'start_time', 'end_time', 'duration_hours', ...
      'bucket', 'season'});
   truth = native(gap);
   filled = truth + 10;
   leaked = timetable(times, native, 'VariableNames', {'tair'});
   leaked.tair(gap) = repmat([200; 300], 3, 1);
   masked = leaked;
   masked.tair(gap) = NaN;

   leaked_metrics = icemodel.forcing.reconstruct.validationMetrics( ...
      truth, filled, gaps, leaked, "tair");
   masked_metrics = icemodel.forcing.reconstruct.validationMetrics( ...
      truth, filled, gaps, masked, "tair");

   testCase.verifyEqual(leaked_metrics.overall.boundary_jump_rate, 0);
   testCase.verifyEqual(masked_metrics.overall.boundary_jump_rate, 1);
end

function test_metrics_count_relational_shortwave_violation(testCase)
    % Held-out scoring counts swu above its paired swd as a hard violation.
    series = icemodel.test.fixtures.makeReconstructSeries();
    n = height(series);
    series.swd = 100 * ones(n, 1);
    series.swu = 40 * ones(n, 1);
    times = series.Properties.RowTimes;
    gaps = table(times(5), times(6), 0.5, 2, ...
       icemodel.forcing.reconstruct.seasonOf(times(5)), ...
       'VariableNames', {'start_time', 'end_time', 'duration_hours', ...
       'bucket', 'season'});
    metrics = icemodel.forcing.reconstruct.validationMetrics( ...
       [40; 40], [50; 120], gaps, series, "swu");
    testCase.verifyEqual(metrics.overall.bound_violations, 1);
end

function test_metrics_score_uncertainty_calibration(testCase)
   % Perfectly sized 1-sigma errors give ~68% nominal coverage bounds; a
   % constant error equal to sigma sits exactly at the 1-sigma edge.
   series = icemodel.test.fixtures.makeReconstructSeries();
   draws = icemodel.forcing.reconstruct.syntheticMissingness(series, ...
      "tair", seededRuns(), years=2020, seed=11, n_gaps=2);
   truth = series.tair(draws.mask);
   filled = truth + 0.5;
   sigma = ones(size(truth));

   metrics = icemodel.forcing.reconstruct.validationMetrics( ...
      truth, filled, draws.gaps, series, "tair", sigma=sigma);

   testCase.verifyEqual(metrics.overall.sigma1_coverage, 1);
   testCase.verifyEqual(metrics.overall.sigma2_coverage, 1);
   testCase.verifyEqual(metrics.overall.correlation, 1, 'AbsTol', 1e-12);
   testCase.verifyEqual(metrics.overall.variability_ratio, 1, ...
      'AbsTol', 1e-12);
end

function test_metrics_variability_removes_between_gap_offsets(testCase)
   % Large differences between gap means must not disguise a candidate that
   % is constant within every outage.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:7)).';
   native = (250:257).';
   series = timetable(times, native, 'VariableNames', {'tair'});
   gaps = table(times([2; 6]), times([3; 7]), [1; 1], [1; 1], ...
      ["DJF"; "DJF"], 'VariableNames', {'start_time', 'end_time', ...
      'duration_hours', 'bucket', 'season'});
   truth = [250; 252; 270; 272];
   filled = [251; 251; 271; 271];

   metrics = icemodel.forcing.reconstruct.validationMetrics( ...
      truth, filled, gaps, series, "tair");

   testCase.verifyEqual(metrics.overall.variability_ratio, 0);
   testCase.verifyEqual(metrics.by_stratum.variability_ratio, 0);
end

function test_metrics_duration_uses_only_complete_gaps(testCase)
   % A finite fragment in a longer draw cannot authorize that draw length.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   first = (10:11).';
   second = (20:23).';
   gaps = table(times([first(1); second(1)]), ...
      times([first(end); second(end)]), [2; 4], [1; 1], ["DJF"; "DJF"], ...
      'VariableNames', {'start_time', 'end_time', 'duration_hours', ...
      'bucket', 'season'});
   truth = series.lwd([first; second]);
   filled = truth;
   filled(end) = NaN;

   returned = icemodel.forcing.reconstruct.validationMetrics( ...
      truth, filled, gaps, series, "lwd");

   testCase.verifyEqual(returned.max_complete_gap_hours, 2);
end

function test_metrics_handle_constant_spread(testCase)
   % Variability diagnostics stay explicit for constant truth/predictions
   % instead of dividing by zero or inventing a correlation.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   gaps = table(times(5), times(6), 2, 2, ...
      icemodel.forcing.reconstruct.seasonOf(times(5)), ...
      'VariableNames', {'start_time', 'end_time', 'duration_hours', ...
      'bucket', 'season'});

    constant_truth = icemodel.forcing.reconstruct.validationMetrics( ...
       [1; 1], [1; 1], gaps, series, "tair");
    testCase.verifyTrue(isnan(constant_truth.overall.correlation));
    testCase.verifyTrue(isnan( ...
       constant_truth.overall.variability_ratio));
    testCase.verifyEqual( ...
       constant_truth.overall.within_gap_observed_spread, 0);

   compressed = icemodel.forcing.reconstruct.validationMetrics( ...
      [1; 2], [1; 1], gaps, series, "tair");
   testCase.verifyTrue(isnan(compressed.overall.correlation));
    testCase.verifyEqual(compressed.overall.variability_ratio, 0);
end

function test_gate_rejects_unmeasurable_within_gap_variability(testCase)
   % Multiple singleton gaps cannot prove a constant held-out stratum.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:4)).';
   series = timetable(times, (259:263).', 'VariableNames', {'tair'});
   gaps = table(times([2; 4]), times([2; 4]), [0.25; 0.25], [1; 1], ...
      ["DJF"; "DJF"], 'VariableNames', {'start_time', 'end_time', ...
      'duration_hours', 'bucket', 'season'});
   metrics = icemodel.forcing.reconstruct.validationMetrics( ...
      [260; 262], [260; 262], gaps, series, "tair", provenance=[1; 1]);

   testCase.verifyTrue(isnan( ...
      metrics.overall.within_gap_observed_spread));
   gate = icemodel.forcing.reconstruct.admissionGate( ...
      "tair", metrics.overall, 1);
   testCase.verifyFalse(gate.admit);
   testCase.verifyTrue(any(contains(gate.reasons, ...
      "within-gap variability evidence unavailable")));

   % One measurable gap cannot hide a singleton gap in the same stratum.
   mixed_gaps = table(times([2; 4]), times([2; 5]), [1; 2], [1; 1], ...
      ["DJF"; "DJF"], 'VariableNames', {'start_time', 'end_time', ...
      'duration_hours', 'bucket', 'season'});
   mixed = icemodel.forcing.reconstruct.validationMetrics( ...
      [260; 262; 263], [260; 262; 263], mixed_gaps, series, "tair", ...
      provenance=[1; 1; 1]);
   testCase.verifyTrue(isnan(mixed.overall.within_gap_observed_spread));
   mixed_gate = icemodel.forcing.reconstruct.admissionGate( ...
      "tair", mixed.overall, 1);
   testCase.verifyFalse(mixed_gate.admit);
   testCase.verifyTrue(any(contains(mixed_gate.reasons, ...
      "within-gap variability evidence unavailable")));
end

function test_gate_admits_and_denies_per_policy(testCase)
   % The gate admits a clean row and records one reason per failed criterion.
   clean = table(100, 1.0, 0.1, 1.0, 1.0, 1.0, 1.0, 0, 0.0, ...
       NaN, NaN, 1.0, ...
       'VariableNames', {'n', 'coverage', 'bias', 'rmse', 'correlation', ...
       'variability_ratio', 'within_gap_observed_spread', ...
       'bound_violations', 'boundary_jump_rate', ...
       'sigma1_coverage', 'sigma2_coverage', 'provenance_accounting'});

   admitted = icemodel.forcing.reconstruct.admissionGate("tair", clean, 2.0);
   testCase.verifyTrue(admitted.admit);
   testCase.verifyEmpty(admitted.reasons);

   dirty = clean;
   dirty.bias = 2;                 % exceeds 0.5 K cap
   dirty.rmse = 1.95;              % not 10% better than baseline 2.0
   dirty.variability_ratio = 0.2;  % suppresses held-out variability
   dirty.bound_violations = 3;
   dirty.coverage = 0.05;          % below the 0.10 usefulness floor
   denied = icemodel.forcing.reconstruct.admissionGate("tair", dirty, 2.0);
   testCase.verifyFalse(denied.admit);
   testCase.verifyEqual(numel(denied.reasons), 5);

   inflated = clean;
   inflated.variability_ratio = 2;
   high_spread = icemodel.forcing.reconstruct.admissionGate( ...
      "tair", inflated, 2.0);
   testCase.verifyFalse(high_spread.admit);
   testCase.verifyTrue(any(contains(high_spread.reasons, ...
      "variability ratio")));

   % A NaN baseline disables only the improvement criterion and records it.
   no_baseline = icemodel.forcing.reconstruct.admissionGate( ...
      "tair", clean, NaN);
    testCase.verifyTrue(no_baseline.admit);
    testCase.verifyFalse(no_baseline.baseline_available);

    % Missing or nonfinite provenance accounting is a failed invariant, not
    % permission to skip the provenance gate.
    missing_provenance = removevars(clean, 'provenance_accounting');
    missing_gate = icemodel.forcing.reconstruct.admissionGate( ...
       "tair", missing_provenance, 2.0);
    testCase.verifyFalse(missing_gate.admit);
    testCase.verifyTrue(any(contains(missing_gate.reasons, ...
       "provenance accounting unavailable")));
    nonfinite_provenance = clean;
    nonfinite_provenance.provenance_accounting = NaN;
    nonfinite_gate = icemodel.forcing.reconstruct.admissionGate( ...
       "tair", nonfinite_provenance, 2.0);
    testCase.verifyFalse(nonfinite_gate.admit);
end

function test_common_support_skill_uses_identical_samples(testCase)
   % Unsupported hard samples cannot inflate the baseline used to admit a
   % candidate that supplies only one easy sample.
   skill = icemodel.forcing.reconstruct.commonSupportSkill( ...
      [0; 10; 20; 30], [1; NaN; NaN; NaN], [0.5; 100; 100; 100]);
   testCase.verifyEqual(skill.n, 1);
   testCase.verifyEqual(skill.candidate_rmse, 1);
   testCase.verifyEqual(skill.baseline_rmse, 0.5);
   testCase.verifyEqual(skill.fractional_improvement, -1);

   unavailable = icemodel.forcing.reconstruct.commonSupportSkill( ...
      [0; 1], [NaN; NaN], [0; 1]);
   testCase.verifyEqual(unavailable.n, 0);
   testCase.verifyTrue(isnan(unavailable.rmse_ratio));

   tied = icemodel.forcing.reconstruct.commonSupportSkill(0, 0, 0);
   worse = icemodel.forcing.reconstruct.commonSupportSkill(0, 1, 0);
   testCase.verifyEqual(tied.fractional_improvement, 0);
   testCase.verifyEqual(worse.rmse_ratio, Inf);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.commonSupportSkill([0; 1], 0, 0), ...
      'icemodel:reconstruct:commonSupportSkill:sizeMismatch');
end

function test_census_defaults_to_all_numeric_channels(testCase)
   % Omitting channels censuses every numeric variable in the series.
   returned = icemodel.forcing.reconstruct.gapCensus(icemodel.test.fixtures.makeReconstructSeries());
   testCase.verifyEqual(sort(string(returned.summary.channel)), ...
      sort(["tair"; "rh"; "psfc"; "lwd"]));
end

function test_census_rejects_all_missing_core(testCase)
   % A series with no finite core samples cannot define an observed record.
   series = icemodel.test.fixtures.makeReconstructSeries();
   series.tair(:) = NaN;
   series.rh(:) = NaN;
   series.psfc(:) = NaN;
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.gapCensus(series, channels="lwd"), ...
      'icemodel:reconstruct:gapCensus:noFiniteCore');
end

function test_census_daylight_cut_ignores_night_shortwave(testCase)
   % With the site point supplied, nighttime-missing swd is not outage;
   % only the deliberate daytime gap counts.
   series = icemodel.test.fixtures.makeReconstructSeries();
   elevation = icemodel.forcing.helpers.solarElevation( ...
      series.Properties.RowTimes, 67.0, -48.8);
   swd = max(0, 600 * sind(max(elevation(:), 0)));
   swd(elevation(:) <= 0) = NaN;                 % night screening pattern
   series.swd = swd;
   % Mask six CONTIGUOUS daylight hours — midsummer midday at 67 N is
   % guaranteed daylight, so the gap cannot straddle a night (which would
   % correctly census as two runs).
   times = series.Properties.RowTimes;
   gap_idx = find(times >= datetime(2020, 6, 21, 9, 0, 0, ...
      'TimeZone', 'UTC') & times <= datetime(2020, 6, 21, 14, 0, 0, ...
      'TimeZone', 'UTC'));
   series.swd(gap_idx) = NaN;

   returned = icemodel.forcing.reconstruct.gapCensus(series, ...
      channels="swd", latitude=67.0, longitude=-48.8);

   testCase.verifyEqual(returned.summary.n_missing, 6);
   testCase.verifyEqual(returned.summary.n_runs, 1);
end

function test_shortwave_draws_exclude_dark_samples(testCase)
   % SWD validation must never mask finite polar-night zeros.
   times = (datetime(2020, 3, 20, 0, 0, 0, 'TimeZone', 'UTC'): ...
      hours(1):datetime(2020, 3, 21, 23, 0, 0, 'TimeZone', 'UTC')).';
   series = timetable(times, zeros(numel(times), 1), ...
      'VariableNames', {'swd'});
   runs = table("swd", times(1), times(1), 1, 1, "MAM", ...
      'VariableNames', {'channel', 'start_time', 'end_time', ...
      'duration_hours', 'bucket', 'season'});

   draws = icemodel.forcing.reconstruct.syntheticMissingness( ...
      series, "swd", runs, years=2020, seed=3, n_gaps=20, ...
      bucket=1, season="MAM", context_hours=0, ...
      latitude=67.0, longitude=-48.8);

   toa = icemodel.forcing.reconstruct.toaIrradiance(times, 67.0, -48.8);
   testCase.verifyGreaterThan(draws.inserted, 0);
   testCase.verifyTrue(all(toa(draws.mask) >= ...
      icemodel.forcing.reconstruct.setopts().toa_dark_wm2));
end

function test_draws_honor_bucket_filter(testCase)
   % A bucket-filtered draw resamples only that bucket's durations.
   series = icemodel.test.fixtures.makeReconstructSeries();
   draws = icemodel.forcing.reconstruct.syntheticMissingness(series, ...
      "tair", seededRuns(), years=2020, seed=3, n_gaps=6, bucket=2);
   testCase.verifyTrue(all(ismember(draws.gaps.duration_hours, [12, 24])));
end

function test_draws_report_shortfall_when_pool_exhausts(testCase)
   % When observed spans cannot host every requested gap, the shortfall is
   % reported through inserted < requested instead of erroring.
   series = icemodel.test.fixtures.makeReconstructSeries();
   series.tair(200:end) = NaN;                   % tiny observed span
   draws = icemodel.forcing.reconstruct.syntheticMissingness(series, ...
      "tair", seededRuns(), years=2020, seed=5, n_gaps=50, ...
      context_hours=12);
   testCase.verifyLessThan(draws.inserted, draws.requested);
end

function test_metrics_reject_malformed_inputs(testCase)
   % Size mismatches, unsorted gap tables, and gap/sample disagreements all
   % fail loudly rather than mis-scoring.
   series = icemodel.test.fixtures.makeReconstructSeries();
   draws = icemodel.forcing.reconstruct.syntheticMissingness(series, ...
      "tair", seededRuns(), years=2020, seed=13, n_gaps=3);
   truth = series.tair(draws.mask);

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.validationMetrics( ...
      truth, truth(1:end - 1), draws.gaps, series, "tair"), ...
      'icemodel:reconstruct:validationMetrics:sizeMismatch');
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.validationMetrics( ...
      truth, truth, draws.gaps, series, "tair", provenance=1), ...
      'icemodel:reconstruct:validationMetrics:provenanceSizeMismatch');
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.validationMetrics( ...
      truth, truth, flipud(draws.gaps), series, "tair"), ...
      'icemodel:reconstruct:validationMetrics:unsortedGaps');
   short_gaps = draws.gaps(1:end - 1, :);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.validationMetrics( ...
      truth, truth, short_gaps, series, "tair"), ...
      'icemodel:reconstruct:validationMetrics:gapSampleMismatch');
end

function test_gate_relative_wspd_cap_and_unknown_channel(testCase)
   % The wspd cap widens to 10% of the typical magnitude when supplied, and
   % unknown channels fail loudly.
   row = table(100, 1.0, 1.4, 1.0, 1.0, 1.0, 1.0, 0, 0.0, ...
       NaN, NaN, 1.0, ...
       'VariableNames', {'n', 'coverage', 'bias', 'rmse', 'correlation', ...
       'variability_ratio', 'within_gap_observed_spread', ...
       'bound_violations', 'boundary_jump_rate', ...
       'sigma1_coverage', 'sigma2_coverage', 'provenance_accounting'});

   fixed = icemodel.forcing.reconstruct.admissionGate("wspd", row, 2.0);
   testCase.verifyFalse(fixed.admit);            % 1.4 > 1 m/s fixed cap
   relative = icemodel.forcing.reconstruct.admissionGate("wspd", row, ...
      2.0, typical_magnitude=15);
   testCase.verifyTrue(relative.admit);          % cap max(1, 1.5) = 1.5
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.admissionGate("nope", row, 1.0), ...
      'icemodel:reconstruct:admissionGate:unknownChannel');
end

function test_physical_bounds_registry(testCase)
   % Known channels return [lower upper]; unknown channels fail loudly.
   returned = icemodel.forcing.reconstruct.physicalBounds("tair");
   testCase.verifyEqual(returned, [193, 300]);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.physicalBounds("nope"), ...
      'icemodel:reconstruct:physicalBounds:unknownChannel');
end

function test_scalar_validity_uses_physical_bounds(testCase)
   % Readiness-grade validity uses the scalar registry without applying
   % relational shortwave rules.
   returned = icemodel.forcing.reconstruct.scalarValidity( ...
      "tair", [193; 250; 300; 301; NaN]);
   testCase.verifyEqual(returned, [true; true; true; false; false]);
   swd = icemodel.forcing.reconstruct.scalarValidity( ...
      "swd", [0; 500; Inf; -1]);
   testCase.verifyEqual(swd, [true; true; false; false]);
end

function test_solar_elevation_bands_contract(testCase)
   % The swd solar-geometry constants are single-sourced (D-28): strictly
   % increasing calibration bin edges spanning the full elevation range,
   % the civil-twilight boundary as an interior edge shared with the
   % darkness pre-pass, and a positive per-band support floor.
   returned = icemodel.forcing.reconstruct.solarElevationBands();
   edges = returned.calibration_bin_edges_deg;
   testCase.verifyEqual(edges(1), -90);
   testCase.verifyEqual(edges(end), 90);
   testCase.verifyTrue(all(diff(edges) > 0));
   testCase.verifyTrue(ismember(returned.civil_twilight_deg, ...
      edges(2:end - 1)));
   testCase.verifyGreaterThan(returned.min_bin_samples, 0);
   testCase.verifyGreaterThan(returned.toa_ceiling_floor_wm2, 0);
   testCase.verifyGreaterThan(returned.toa_ceiling_multiplier, 1);
   testCase.verifyEqual(returned.solar_constant_wm2, 1361);
   testCase.verifyGreaterThan(returned.twilight_ceiling_wm2, ...
      returned.toa_ceiling_floor_wm2);
end

function test_plan_prefers_calibrated_proxy_over_climatology_for_swd(testCase)
   % D-29: when a stratum admits both, the calibrated RCM proxy must be
   % walked before climatology for swd, even when climatology's held-out
   % skill ranks higher (its day-of-year median structurally inserts
   % clearer-than-context days into cloudy weeks). Constructed so the
   % periodic truth makes climatology near-perfect while the proxy's
   % aperiodic wobble leaves residual error: both beat the bucket-1
   % persistence baseline, climatology out-skills the proxy, and only the
   % deterministic post-ranking swap can list the proxy first.
   times = (datetime(2019, 1, 1, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2021, 12, 31, 23, 0, 0, 'TimeZone', 'UTC')).';
   toa = icemodel.forcing.reconstruct.toaIrradiance(times, 45.0, -48.8);
   doy = day(times, 'dayofyear');
   % Finite core channels support the year split; only swd is planned.
   tair = 260 + 10 * sin(2 * pi * (doy - 30) / 365);
   rh = 80 + zeros(numel(times), 1);
   psfc = 90000 + zeros(numel(times), 1);
   % Perfectly periodic truth: the day-of-year x hour climatology
   % reproduces it almost exactly, guaranteeing top-ranked climatology.
   swd = 0.7 * toa;
   % The proxy tracks the truth through a fixed ratio plus a small
   % 11-day multiplicative wobble the calendar-locked climatology cannot
   % see: calibration recovers the ratio, the wobble leaves error.
   wobble = 1 + 0.10 * sin(2 * pi * (1:numel(times)).' / (24 * 11));
   model = (swd / 0.8) .* wobble;
   % Real JJA daylight gaps seed the census duration pool for bucket 1
   % in every year, so both selection and evaluation draws exist.
   durations = [2, 3, 4, 3, 2];
   for y = [2019, 2020, 2021]
      for k = 1:numel(durations)
         t0 = datetime(y, 7, 2 + 3 * k, 13, 0, 0, 'TimeZone', 'UTC');
         swd(times >= t0 & times < t0 + hours(durations(k))) = NaN;
      end
   end
   series = timetable(times, tair, rh, psfc, swd, ...
      'VariableNames', {'tair', 'rh', 'psfc', 'swd'});
   target = struct('series', series, 'station', "swdtest", 'location', ...
      struct('lat_wgs84', 45.0, 'lon_wgs84', -48.8, 'elev_m', 1000));
   donors = struct('series', {}, 'station', {}, 'family', {}, ...
      'location', {});
   proxies = struct('series', timetable(times, model, ...
      'VariableNames', {'swd'}), 'name', "mar", 'code_name', "mar");

   % swd draws demand CONTINUOUS daylight across the gap plus both
   % context margins; a mid-latitude day cannot host the default 24 h
   % margin, and the margin length is not what this test exercises.
   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      donors, proxies, channels="swd", seed=7, synthetic_context_hours=2);

   methods = plan.channels(1).methods;
   testCase.assertNotEmpty(methods);
   names = string({methods.name}).';
   seasons = string({methods.seasons}).';
   buckets = [methods.buckets].';
   proxy_index = find(names == "proxy:mar" & seasons == "JJA" ...
      & buckets == 1, 1);
   clim_index = find(names == "climatology" & seasons == "JJA" ...
      & buckets == 1, 1);
   testCase.assertNotEmpty(proxy_index);
   testCase.assertNotEmpty(clim_index);
   % Climatology genuinely out-skills the proxy on the periodic truth,
   % so pure skill ranking would list it first...
   testCase.verifyGreaterThan( ...
      methods(clim_index).selection.fractional_improvement, ...
      methods(proxy_index).selection.fractional_improvement);
   % ...but the D-29 swap lists the calibrated proxy first in the walk
   % order the engine consumes.
   testCase.verifyLessThan(proxy_index, clim_index);
   % The persisted calibration registry carries the version-2 binned
   % schema (D-28) for downstream consumers.
   calibrations = plan.channels(1).proxy_calibrations;
   testCase.verifyEqual(string(calibrations(1).source), "mar");
   testCase.verifyEqual(calibrations(1).parameters.version, 2);
   testCase.verifyTrue(isfield(calibrations(1).parameters, ...
      'binned_corrections'));
end

%% Fixtures


function runs = seededRuns()
   %SEEDEDRUNS Small census-run pool with short and daily durations.
   runs = table( ...
      repmat("tair", 4, 1), ...
      repmat(datetime(2020, 6, 1, 'TimeZone', 'UTC'), 4, 1), ...
      repmat(datetime(2020, 6, 2, 'TimeZone', 'UTC'), 4, 1), ...
      [3; 5; 12; 24], [1; 1; 2; 2], ...
      ["JJA"; "JJA"; "DJF"; "SON"], ...
      'VariableNames', {'channel', 'start_time', 'end_time', ...
      'duration_hours', 'bucket', 'season'});
end

function writeJson(filename, value)
   %WRITEJSON Persist one compact fixture manifest.
   fid = fopen(filename, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', jsonencode(value));
   clear cleaner
end
