function tests = test_reconstruct_qa_screen
   %TEST_RECONSTRUCT_QA_SCREEN Verify the buried/rime-encased sensor screen.
   %
   % Exercises icemodel.forcing.reconstruct.flatRunScreen against a clean
   % synthetic year and the same year with the verified FRE 2018 buried-run
   % signature stamped over a spring window (bead icemodel-g1n.45).
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Install the verification path for namespace resolution.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
end

function teardown(testCase)
   % Dropping the stored handle destroys the onCleanup object, which
   % restores the caller's configuration deterministically per test.
   testCase.TestData.cleanup = [];
end

%% Screen behavior

function test_clean_series_not_flagged(testCase)
   % A clean diurnal year produces no flags and the documented empty
   % findings schema, so per-station results vertcat cleanly.
   met = makeCleanMet();

   [returned, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met, 67, -50);

   testCase.verifyEqual(returned, false(height(met), 1));
   testCase.verifyEqual(height(findings), 0);
   expected = {'start_time', 'end_time', 'n_days', 'channels', ...
      'tair_range_max_k', 'swd_max_wm2', 'toa_max_wm2', ...
      'lwd_dev_wm2', 'rh_range_pct'};
   testCase.verifyEqual(findings.Properties.VariableNames, expected);
end

function test_buried_run_flagged_with_full_evidence(testCase)
   % The FRE signature (flat tair, dead swd under daylight, blackbody
   % lwd, pinned rh) is flagged over exactly its 22-day window, and the
   % findings row carries the evidence statistics.
   [met, window] = makeBuriedMet();

   [returned, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met, 67, -50);

   % Every sample in the stamped window is flagged, nothing outside it.
   testCase.verifyEqual(returned, window);
   testCase.verifyEqual(height(findings), 1);
   % The run spans the stamped days and implicates all four channels.
   times = met.Properties.RowTimes;
   testCase.verifyEqual(findings.start_time, min(times(window)));
   testCase.verifyEqual(findings.end_time, max(times(window)));
   testCase.verifyEqual(findings.n_days, 22);
   testCase.verifyEqual(findings.channels, "tair,swd,lwd,rh");
   % Evidence statistics match the stamped synthetic signature: 0.1 K
   % tair wiggle, 0.5 W/m2 swd, ~1 W/m2 blackbody deviation, rh at 85.
   testCase.verifyEqual(findings.tair_range_max_k, 0.1, 'AbsTol', 1e-2);
   testCase.verifyEqual(findings.swd_max_wm2, 0.5, 'AbsTol', 1e-12);
   testCase.verifyGreaterThan(findings.toa_max_wm2, 100);
   testCase.verifyEqual(findings.lwd_dev_wm2, 1, 'AbsTol', 1e-6);
   testCase.verifyEqual(findings.rh_range_pct, 0, 'AbsTol', 1e-12);
end

function test_short_run_not_flagged(testCase)
   % A two-day flat run stays under the three-day floor: persistent
   % synoptic overcast must not trip the screen.
   [met, ~] = makeBuriedMet(datetime(2020, 3, 25, 'TimeZone', 'UTC'), ...
      datetime(2020, 3, 27, 'TimeZone', 'UTC'));

   [returned, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met, 67, -50);

   testCase.verifyFalse(any(returned));
   testCase.verifyEqual(height(findings), 0);
end

function test_lwd_only_corroboration(testCase)
   % With a working pyranometer (realistic swd) the blackbody lwd alone
   % corroborates the flat tair, and swd is not implicated.
   [met, window] = makeBuriedMet();
   % Restore clean shortwave inside the window: the sensor still sees sun.
   clean = makeCleanMet();
   met.swd(window) = clean.swd(window);

   [returned, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met, 67, -50);

   testCase.verifyEqual(returned, window);
   testCase.verifyEqual(findings.channels, "tair,lwd,rh");
end

function test_swd_only_corroboration_without_rh(testCase)
   % With lwd well below blackbody the dead daylight pyranometer alone
   % corroborates; a met table without rh reports NaN evidence for it.
   [met, window] = makeBuriedMet();
   sigma = icemodel.physicalConstant('SB');
   % Clear-sky longwave sits far below the blackbody envelope.
   met.lwd(window) = sigma * met.tair(window) .^ 4 - 60;
   met = removevars(met, 'rh');

   [returned, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met, 67, -50);

   testCase.verifyEqual(returned, window);
   testCase.verifyEqual(findings.channels, "tair,swd");
   testCase.verifyTrue(isnan(findings.rh_range_pct));
end

function test_tair_only_requires_corroboration_toggle(testCase)
   % A flat tair with no radiation channels only flags when the policy
   % corroboration requirement is explicitly relaxed (diagnostic mode).
   [met, window] = makeBuriedMet();
   met = removevars(met, {'swd', 'lwd', 'rh'});

   [returned, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met, 67, -50);
   % Default policy: no corroboration available, so nothing is flagged.
   testCase.verifyFalse(any(returned));
   testCase.verifyEqual(height(findings), 0);

   [returned, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met, 67, -50, require_corroboration=false);
   % Diagnostic mode screens on the core condition alone.
   testCase.verifyEqual(returned, window);
   testCase.verifyEqual(findings.channels, "tair");
   testCase.verifyTrue(isnan(findings.swd_max_wm2));
   testCase.verifyTrue(isnan(findings.lwd_dev_wm2));
end

function test_polar_night_darkness_not_corroborating(testCase)
   % Dark swd during polar night is instrument reality, not evidence:
   % the TOA daylight gate must block condition (b) in December.
   met = makeCleanMet();
   times = met.Properties.RowTimes;
   window = times >= datetime(2020, 12, 10, 'TimeZone', 'UTC') ...
      & times < datetime(2020, 12, 20, 'TimeZone', 'UTC');
   sigma = icemodel.physicalConstant('SB');
   % Flat tair and dead swd, but lwd stays clear-sky (no corroboration
   % path): polar night TOA never reaches the daylight gate at 67 N.
   met.tair(window) = 250;
   met.swd(window) = 0;
   met.lwd(window) = sigma * met.tair(window) .^ 4 - 60;

   [returned, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met, 67, -50);

   testCase.verifyFalse(any(returned));
   testCase.verifyEqual(height(findings), 0);
end

function test_low_coverage_day_splits_run(testCase)
   % A day with mostly-missing tair is not evaluable and splits the run
   % into two findings; the unevaluable day itself is never flagged.
   [met, window] = makeBuriedMet();
   times = met.Properties.RowTimes;
   split_day = dateshift(times, 'start', 'day') ...
      == datetime(2020, 4, 1, 'TimeZone', 'UTC');
   met.tair(split_day) = NaN;

   [returned, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met, 67, -50);

   % Mar 25-31 (7 days) and Apr 2-15 (14 days) survive as separate runs.
   expected = window & ~split_day;
   testCase.verifyEqual(returned, expected);
   testCase.verifyEqual(height(findings), 2);
   testCase.verifyEqual(findings.n_days, [7; 14]);
end

function test_axis_hole_splits_run(testCase)
   % Removing a day of rows entirely breaks calendar consecutiveness:
   % the screen must not bridge a hole in the time axis.
   [met, ~] = makeBuriedMet();
   times = met.Properties.RowTimes;
   hole = dateshift(times, 'start', 'day') ...
      == datetime(2020, 4, 1, 'TimeZone', 'UTC');
   met = met(~hole, :);

   [~, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met, 67, -50);

   testCase.verifyEqual(height(findings), 2);
   testCase.verifyEqual(findings.n_days, [7; 14]);
end

function test_threshold_overrides_respected(testCase)
   % Tightening either the flatness threshold or the run-length floor
   % below/above the stamped signature disables the flag, proving the
   % name-value knobs actually plumb into the conditions.
   [met, ~] = makeBuriedMet();

   [returned, ~] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met, 67, -50, max_tair_range_k=0.05);
   % The stamped 0.1 K daily range exceeds a 0.05 K flatness threshold.
   testCase.verifyFalse(any(returned));

   [returned, ~] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met, 67, -50, min_run_days=30);
   % A 22-day run stays under a 30-day floor.
   testCase.verifyFalse(any(returned));
end

function test_empty_and_scalar_met(testCase)
   % Zero rows return empty outputs; a single sample can never form a
   % multi-day run (and exercises the scalar-cadence branch).
   met = makeCleanMet();

   [returned, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met([], :), 67, -50);
   testCase.verifyEqual(returned, false(0, 1));
   testCase.verifyEqual(height(findings), 0);

   [returned, findings] = icemodel.forcing.reconstruct.flatRunScreen( ...
      met(1, :), 67, -50);
   testCase.verifyEqual(returned, false);
   testCase.verifyEqual(height(findings), 0);
end

function test_missing_tair_errors(testCase)
   % The screen is meaningless without air temperature and must refuse.
   met = removevars(makeCleanMet(), 'tair');
   testCase.verifyError( ...
      @() icemodel.forcing.reconstruct.flatRunScreen(met, 67, -50), ...
      'icemodel:reconstruct:flatRunScreen:missingTair');
end

%% Shortwave seam quality

function test_shortwave_seam_repairs_only_synthetic_sample(testCase)
   % A single implausible synthetic posting between smooth native anchors
   % is reconnected without changing any native observation.
   [times, swd, provenance, spike, codes] = makeSeamFixture();
   native = provenance == codes.observed;
   original = swd;

   [filled, repaired_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=1);

   testCase.verifyEqual(filled(spike), 100);
   testCase.verifyEqual(filled(native), original(native));
   testCase.verifyEqual(repaired_provenance(spike), ...
      codes.bounded_interp);
   testCase.verifyEqual(height(audit), 1);
   testCase.verifyEqual(quality.status, "pass");
   testCase.verifyGreaterThan(quality.initial_outliers, 0);
   testCase.verifyEqual(quality.final_outliers, 0);
   testCase.verifyEqual(quality.changed_samples, 1);
   testCase.verifyEqual(quality.passes, 1);
end

function test_shortwave_seam_preserves_protected_source_codes(testCase)
   % Native SWD and the separately validated twilight-climatology tier may
   % diagnose a seam, but post-final repair must never replace either.
   [times, swd, provenance, spike, codes] = makeSeamFixture();
   original = swd;
   for native_code = [codes.observed, codes.raw_shortwave, ...
         codes.clamped_shortwave, codes.twilight_climatology, ...
         codes.darkness]
      native_provenance = provenance;
      native_provenance(spike) = native_code;
      [returned, returned_provenance] = ...
         icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
         times, swd, native_provenance, 0, 0, min_reference_steps=1);
      testCase.verifyEqual(returned, original);
      testCase.verifyEqual(returned_provenance, native_provenance);
   end
end

function test_shortwave_seam_ignores_solar_regime_transitions(testCase)
   % Dawn/dusk changes across dark, twilight, and daylight supports are
   % solar geometry, not method-seam discontinuities.
   times = (datetime(2020, 3, 20, 0, 0, 0, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2020, 3, 22, 23, 0, 0, 'TimeZone', 'UTC')).';
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, 0, 0, hours(1));
   daylight = elevation >= bands.horizon_deg;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   swd = 100 * double(daylight);
   provenance = repmat(codes.darkness, size(times));
   provenance(daylight) = codes.observed;
   dawn = daylight & [true; ~daylight(1:end - 1)];
   provenance(dawn) = codes.mar;

   [returned, returned_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=1);

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(returned_provenance, provenance);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyEqual(quality.initial_outliers, 0);
   testCase.verifyEqual(quality.status, "pass");
end

function test_shortwave_seam_ignores_elevation_bin_transitions(testCase)
   % A method run occupying one complete calibration band may jump at its
   % solar-band edges without borrowing evidence from either neighboring
   % band.
   times = (datetime(2020, 3, 20, 0, 0, 0, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2020, 3, 22, 23, 0, 0, 'TimeZone', 'UTC')).';
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, 0, 0, hours(1));
   elevation_bin = discretize(elevation, ...
      bands.calibration_bin_edges_deg);
   target_bin = find(bands.calibration_bin_edges_deg == 5);
   synthetic = elevation_bin == target_bin;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   swd = repmat(100, size(times));
   swd(synthetic) = 400;
   provenance = repmat(codes.observed, size(times));
   provenance(synthetic) = codes.mar;

   [returned, returned_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=1);

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(returned_provenance, provenance);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyEqual(quality.initial_outliers, 0);
   testCase.verifyEqual(quality.status, "pass");
end

function test_shortwave_seam_ignores_regime_edge_neighborhood(testCase)
   % The first same-daylight-regime seam at dawn is still part of the
   % physical solar transition and must not be diagnosed or repaired.
   times = (datetime(2020, 3, 20, 0, 0, 0, 'TimeZone', 'UTC'):minutes(5): ...
      datetime(2020, 3, 22, 23, 0, 0, 'TimeZone', 'UTC')).';
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, 0, 0, minutes(5));
   daylight = elevation >= bands.horizon_deg;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   swd = 100 * double(daylight);
   provenance = repmat(codes.darkness, size(times));
   provenance(daylight) = codes.observed;
   dawn = find(daylight & [true; ~daylight(1:end - 1)], 1);
   provenance(dawn) = codes.mar;
   swd(dawn) = 400;

   [returned, returned_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=1);

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(returned_provenance, provenance);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyEqual(quality.initial_outliers, 0);
   testCase.verifyEqual(quality.changed_samples, 0);
   testCase.verifyEqual(quality.final_outliers, 0);
   testCase.verifyEqual(quality.status, "pass");
end

function test_shortwave_seam_accepts_smooth_method_ramp(testCase)
   % A method transition embedded in a locally linear solar ramp is not a
   % discontinuity even when each step exceeds the seasonal percentile
   % and its narrow stratum does not meet the reference-count minimum.
   [times, swd, provenance, synthetic] = makeSeamFixture();
   ramp = synthetic + (-2:2);
   swd(ramp) = (0:100:400).';

   [returned, returned_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, percentile=90, ...
      min_reference_steps=numel(times));

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(returned_provenance, provenance);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyEqual(quality.initial_outliers, 0);
   testCase.verifyEqual(quality.final_unsupported_boundaries, 0);
   testCase.verifyEqual(quality.changed_samples, 0);
   testCase.verifyEqual(quality.status, "pass");
end

function test_shortwave_seam_accepts_local_ramp_across_fine_bin(testCase)
   % Local mathematical continuity may cross a fine calibration-bin edge
   % when the diagnosed method boundary itself remains in one band.
   times = (datetime(2020, 3, 20, 0, 0, 0, 'TimeZone', 'UTC'):minutes(5): ...
      datetime(2020, 3, 22, 23, 0, 0, 'TimeZone', 'UTC')).';
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, 0, 0, minutes(5));
   elevation_bin = discretize(elevation, ...
      bands.calibration_bin_edges_deg);
   regime = 1 + (elevation > bands.civil_twilight_deg) ...
      + (elevation >= bands.horizon_deg);
   candidate = find( ...
      elevation_bin(1:end - 2) ~= elevation_bin(2:end - 1) ...
      & elevation_bin(2:end - 1) == elevation_bin(3:end) ...
      & regime(1:end - 2) == regime(2:end - 1) ...
      & regime(2:end - 1) == regime(3:end), 1) + 1;
   testCase.assertNotEmpty(candidate);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   swd = repmat(100, size(times));
   swd(candidate + (-1:2)) = (0:100:300).';
   provenance = repmat(codes.observed, size(times));
   provenance(candidate) = codes.mar;

   [returned, returned_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, percentile=90, ...
      min_reference_steps=1);

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(returned_provenance, provenance);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyEqual(quality.initial_outliers, 0);
   testCase.verifyEqual(quality.status, "pass");
end

function test_shortwave_seam_repairs_one_side_of_method_transition(testCase)
   % A seam between two synthetic methods still changes only one posting
   % during one narrow repair pass.
   [times, swd, provenance, synthetic, codes] = makeSeamFixture();
   provenance(synthetic) = codes.mar;
   provenance(synthetic + 1) = codes.merra2;
   swd(synthetic) = 100;
   swd(synthetic + 1) = 400;

   [returned, ~, ~, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, percentile=90, ...
      min_reference_steps=1, max_passes=1);

   testCase.verifyEqual(returned(synthetic), 100);
   testCase.verifyEqual(returned(synthetic + 1), 100);
   testCase.verifyEqual(quality.changed_samples, 1);
end

function test_shortwave_seam_requires_reference_in_matching_stratum(testCase)
   % Native evidence in another season must not authorize a DJF repair.
   times = (datetime(2020, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2020, 6, 30, 23, 0, 0, 'TimeZone', 'UTC')).';
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   swd = repmat(100, size(times));
   provenance = repmat(codes.observed, size(times));
   winter = icemodel.forcing.reconstruct.seasonOf(times) == "DJF";
   provenance(winter) = codes.mar;
   synthetic = find(winter, 1) + 24;
   provenance(synthetic) = codes.merra2;
   swd(synthetic) = 400;

   [returned, returned_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=1);

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(returned_provenance, provenance);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyGreaterThan(quality.final_unsupported_boundaries, 0);
   testCase.verifyEqual(quality.status, "review");

   % Native anchors and their bridge floor must not convert an
   % empirically unsupported stratum into supported evidence.
   bridge_provenance = repmat(codes.observed, size(times));
   bridge_provenance(synthetic) = codes.mar;
   bridge_swd = repmat(100, size(times));
   bridge_swd(synthetic) = 400;
   [~, ~, ~, bridge_quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, bridge_swd, bridge_provenance, 0, 0, ...
      min_reference_steps=numel(times));
   testCase.verifyGreaterThan( ...
      bridge_quality.final_unsupported_boundaries, 0);
   testCase.verifyEqual(bridge_quality.status, "review");

   % The same unsupported method transition is harmless when it is exactly
   % continuous; other seasons still provide the global QA reference.
   swd(synthetic) = 100;
   [~, ~, ~, continuous_quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=1);
   testCase.verifyEqual( ...
      continuous_quality.final_unsupported_boundaries, 0);
   testCase.verifyEqual(continuous_quality.status, "pass");
end

function test_shortwave_seam_rejects_synthetic_continuity_floor(testCase)
   % A same-direction synthetic step is not independent evidence for an
   % otherwise unsupported method offset.
   [times, swd, provenance, synthetic, codes] = makeSeamFixture();
   swd(synthetic:synthetic + 1) = [200; 400];
   provenance(synthetic:synthetic + 1) = codes.mar;
   swd(synthetic + 2) = NaN;
   provenance(synthetic + 2) = codes.missing;

   [returned, returned_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, ...
      min_reference_steps=numel(times));

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(returned_provenance, provenance);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyGreaterThan(quality.final_unsupported_boundaries, 0);
   testCase.verifyEqual(quality.status, "review");
end

function test_shortwave_seam_accepts_exact_finite_anchor_interpolation(testCase)
   % An exact bounded interpolation between two finite postings is
   % continuous even when their source labels differ and the narrow
   % stratum lacks native support.
   [times, swd, provenance, center, codes] = makeSeamFixture();
   support = center + (-2:2);
   swd(support) = [100; 100; 200; 300; 300];
   provenance(support(2)) = codes.mar;
   provenance(support(3)) = codes.bounded_interp;
   provenance(support(4)) = codes.observed;

   [returned, returned_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, ...
      min_reference_steps=numel(times));

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(returned_provenance, provenance);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyEqual(quality.final_unsupported_boundaries, 0);
   testCase.verifyEqual(quality.status, "pass");

   % A mislabeled or nonlinear synthetic point is not exempt.
   swd(center) = 250;
   [~, ~, ~, nonlinear_quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, ...
      min_reference_steps=numel(times));
   testCase.verifyGreaterThan( ...
      nonlinear_quality.final_unsupported_boundaries, 0);
   testCase.verifyEqual(nonlinear_quality.status, "review");
end

function test_shortwave_seam_exact_bridge_support_is_boundary_local(testCase)
   % An exact bridge may begin in twilight; a later endpoint boundary is
   % continuous when that boundary itself is safely inside daylight.
   times = (datetime(2020, 3, 20, 0, 0, 0, 'TimeZone', 'UTC'):minutes(5): ...
      datetime(2020, 3, 22, 23, 55, 0, 'TimeZone', 'UTC')).';
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, 0, 0, minutes(5));
   regime = 1 + (elevation > bands.civil_twilight_deg) ...
      + (elevation >= bands.horizon_deg);
   elevation_bin = discretize(elevation, ...
      bands.calibration_bin_edges_deg);
   bridge_start = find(regime(1:end - 4) == 2 ...
      & regime(2:end - 3) == 3 ...
      & regime(3:end - 2) == 3 ...
      & regime(4:end - 1) == 3 ...
      & elevation_bin(3:end - 2) == elevation_bin(4:end - 1), 1);
   testCase.assertNotEmpty(bridge_start);
   support = bridge_start + (0:4);
   swd = repmat(100, size(times));
   swd(support) = [0; 10; 20; 30; 30];
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   provenance = repmat(codes.observed, size(times));
   provenance(support(2:3)) = codes.bounded_interp;
   provenance(support(4)) = codes.mar;
   same_regime = regime(2:end) == regime(1:end - 1);
   near_regime_edge = ~same_regime;
   near_regime_edge(2:end) = near_regime_edge(2:end) ...
      | ~same_regime(1:end - 1);
   near_regime_edge(1:end - 1) = near_regime_edge(1:end - 1) ...
      | ~same_regime(2:end);
   boundary_step = support(3);
   testCase.verifyNotEqual( ...
      provenance(boundary_step), provenance(boundary_step + 1));
   testCase.verifyEqual( ...
      elevation_bin(boundary_step), elevation_bin(boundary_step + 1));
   testCase.verifyTrue(same_regime(boundary_step));
   testCase.verifyFalse(near_regime_edge(boundary_step));
   testCase.verifyEqual( ...
      icemodel.forcing.reconstruct.seasonOf(times(boundary_step)), ...
      icemodel.forcing.reconstruct.seasonOf(times(boundary_step + 1)));

   [returned, returned_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, ...
      min_reference_steps=numel(times));

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(returned_provenance, provenance);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyGreaterThan(quality.initial_boundaries, 0);
   testCase.verifyEqual(quality.final_unsupported_boundaries, 0);
   testCase.verifyEqual(quality.status, "pass");
end

function test_shortwave_seam_uses_same_regime_reference_fallback(testCase)
   % A sparse fine elevation bin may borrow only same-season native steps
   % from its coarse solar regime.
   times = (datetime(2020, 3, 20, 0, 0, 0, 'TimeZone', 'UTC'):minutes(5): ...
      datetime(2020, 3, 22, 23, 55, 0, 'TimeZone', 'UTC')).';
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, 0, 0, minutes(5));
   bin = discretize(elevation, bands.calibration_bin_edges_deg);
   regime = 1 + (elevation > bands.civil_twilight_deg) ...
      + (elevation >= bands.horizon_deg);
   day_three = day(times) == 22;
   candidate = find(day_three(2:end - 1) ...
      & regime(2:end - 1) == 3 ...
      & bin(1:end - 2) == bin(2:end - 1) ...
      & bin(2:end - 1) == bin(3:end), 1) + 1;
   testCase.assertNotEmpty(candidate);
   target_bin = bin(candidate);
   step = repmat(20, numel(times) - 1, 1);
   step(bin(2:end) == target_bin) = 10;
   swd = cumsum([100; step]);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   provenance = repmat(codes.observed, size(times));
   provenance(candidate:end) = codes.mar;
   swd(candidate:end) = swd(candidate:end) + 5;

   [returned, returned_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=100);

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(returned_provenance, provenance);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyEqual(quality.final_unsupported_boundaries, 0);
   testCase.verifyEqual(quality.status, "pass");
end

function test_shortwave_seam_no_reference_requires_supported_boundaries(testCase)
   % Missing station-native seam references are not themselves a defect;
   % exact/no-op boundaries pass, while a nonlinear boundary still reviews.
   [times, swd, provenance, center, codes] = makeSeamFixture();
   provenance(:) = codes.mar;

   [returned, returned_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=1);

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(returned_provenance, provenance);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyEqual(quality.observed_reference_steps, 0);
   testCase.verifyEqual(quality.final_boundaries, 0);
   testCase.verifyEqual(quality.status, "pass");

   support = center + (-1:1);
   swd(support) = [100; 200; 300];
   provenance(center) = codes.bounded_interp;
   [~, ~, ~, exact_quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=1);
   testCase.verifyEqual(exact_quality.observed_reference_steps, 0);
   testCase.verifyGreaterThan(exact_quality.final_boundaries, 0);
   testCase.verifyEqual(exact_quality.final_unsupported_boundaries, 0);
   testCase.verifyEqual(exact_quality.status, "pass");

   swd(center) = 250;
   [~, ~, ~, nonlinear_quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=1);
   testCase.verifyGreaterThan( ...
      nonlinear_quality.final_unsupported_boundaries, 0);
   testCase.verifyEqual(nonlinear_quality.observed_reference_steps, 0);
   testCase.verifyEqual(nonlinear_quality.status, "review");
end

function test_shortwave_seam_zero_pass_is_diagnostic_only(testCase)
   % max_passes=0 exposes the empirical verdict without mutating data.
   [times, swd, provenance] = makeSeamFixture();

   [returned, returned_provenance, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=1, ...
      max_passes=0);

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(returned_provenance, provenance);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyEqual(quality.status, "review");
   testCase.verifyEqual(quality.changed_samples, 0);
   testCase.verifyEqual(quality.passes, 0);
end

function test_shortwave_seam_does_not_bridge_deep_darkness(testCase)
   % A synthetic midnight spike is diagnosed but not smoothed across the
   % policy's civil-twilight boundary.
   [times, swd, provenance, ~, codes] = makeSeamFixture();
   provenance(:) = codes.observed;
   midnight = find(day(times) == 2 & hour(times) == 0, 1);
   provenance(midnight) = codes.mar;
   swd(midnight) = 400;

   [returned, ~, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=1);

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyGreaterThan(quality.final_outliers, ...
      quality.expected_tail_outliers);
   testCase.verifyEqual(quality.status, "review");
end

function test_shortwave_seam_respects_gap_cap(testCase)
   % A two-posting synthetic run is not broadened into a repair when the
   % caller's one-hour short-gap cap cannot admit it.
   [times, swd, provenance, spike, codes] = makeSeamFixture();
   second = spike + 1;
   provenance(second) = codes.mar;
   swd(second) = swd(spike);

   [returned, ~, audit, quality] = ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times, swd, provenance, 0, 0, min_reference_steps=1, ...
      cap_hours=1);

   testCase.verifyEqual(returned, swd);
   testCase.verifyEqual(audit, cell(0, 1));
   testCase.verifyEqual(quality.status, "review");
end

function test_shortwave_seam_rejects_invalid_axes(testCase)
   % Time, values, and provenance must share an axis with a cadence.
   [times, swd, provenance] = makeSeamFixture();
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times(1:end - 1), swd, provenance, 0, 0), ...
      'icemodel:reconstruct:smoothShortwaveSeams:sizeMismatch');
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.smoothShortwaveSeams( ...
      times(1), swd(1), provenance(1), 0, 0), ...
      'icemodel:reconstruct:smoothShortwaveSeams:shortAxis');
end

%% Local fixtures

function [times, swd, provenance, spike, codes] = makeSeamFixture()
   %MAKESEAMFIXTURE Smooth native flux with one midday synthetic spike.
   times = (datetime(2020, 6, 1, 0, 0, 0, 'TimeZone', 'UTC'):hours(1): ...
      datetime(2020, 6, 3, 23, 0, 0, 'TimeZone', 'UTC')).';
   swd = repmat(100, size(times));
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   provenance = repmat(codes.observed, size(times));
   spike = find(day(times) == 2 & hour(times) == 12, 1);
   swd(spike) = 400;
   provenance(spike) = codes.mar;
end

function met = makeCleanMet()
   %MAKECLEANMET One clean synthetic year: diurnal tair, sunny swd.
   %
   % Extends the shared reconstruction fixture with a 0.7 clear-sky-index
   % shortwave channel so the daylight-darkness corroboration has a
   % realistic working-pyranometer baseline to contrast against.
   met = icemodel.test.fixtures.makeReconstructSeries();
   toa = icemodel.forcing.reconstruct.toaIrradiance( ...
      met.Properties.RowTimes, 67, -50);
   met.swd = 0.7 * toa;
end

function [met, window] = makeBuriedMet(window_start, window_end)
   %MAKEBURIEDMET Clean year with the FRE buried-run signature stamped in.
   %
   % Stamps the verified FRE 2018-03-25..04-15 signature over a spring
   % window of the clean fixture: 0.1 K daily tair range, 0.5 W/m2 swd
   % under high spring TOA, lwd one W/m2 above blackbody, rh pinned 85.
   % Defaults reproduce the 22-day reference window on the fixture year.
   arguments
      window_start (1, 1) datetime = ...
         datetime(2020, 3, 25, 'TimeZone', 'UTC')
      window_end (1, 1) datetime = ...
         datetime(2020, 4, 16, 'TimeZone', 'UTC')
   end
   met = makeCleanMet();
   times = met.Properties.RowTimes;
   window = times >= window_start & times < window_end;
   % A 0.05 K diurnal wiggle keeps the range near the verified < 0.15 K
   % while avoiding a degenerate exactly-constant series.
   met.tair(window) = 250 + 0.05 * sin(2 * pi * hour(times(window)) / 24);
   met.swd(window) = 0.5;
   sigma = icemodel.physicalConstant('SB');
   met.lwd(window) = sigma * met.tair(window) .^ 4 + 1;
   met.rh(window) = 85;
end
