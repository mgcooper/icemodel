function tests = test_forcing_destep
   %TEST_FORCING_DESTEP Verify the surface step-detection / de-stepping transform.
   %
   % Exercises icemodel.forcing.destepSurface and the companion
   % icemodel.forcing.helpers.surfaceFlags on synthetic series (no data
   % dependency) so the detection/classification/correction contract is pinned
   % independent of the staged PROMICE product.
   tests = functiontests(localfunctions);
end

function test_destep_detects_and_levels_unambiguous_step(testCase)
   % A gross single-timestep jump (>> max_step) is detected, classified
   % UNAMBIGUOUS (magnitude + implausible), and leveled away by the default
   % mode, leaving the surrounding gradual trend intact.

   t = (datetime(2015, 1, 1):hours(1):datetime(2015, 1, 1, 9, 0, 0))';
   % A flat pre-step segment, a single clean +5 m jump at sample 6, then a flat
   % post-step segment, so the leveled series is exactly the pre-step level
   % (no ramp increment conflated into the step).
   surf = [zeros(5, 1); 5 * ones(5, 1)];

   [corrected, record, flags] = icemodel.forcing.destepSurface(t, surf);

   % Exactly one unambiguous, applied step at the jump.
   unamb = record([record.classification] == "unambiguous");
   testCase.verifyEqual(numel(unamb), 1);
   testCase.verifyTrue(unamb.applied);
   testCase.verifyEqual(unamb.magnitude, 5, 'AbsTol', 1e-9);

   % The leveled series recovers the flat baseline (the jump is removed).
   testCase.verifyEqual(corrected, zeros(10, 1), 'AbsTol', 1e-9);

   % The per-sample masks mark the post-step sample correctable.
   testCase.verifyEqual(flags.step_correctable(6), 1);
   testCase.verifyEqual(flags.step_magnitude(6), 5, 'AbsTol', 1e-9);
end

function test_destep_flags_ambiguous_step_without_correcting(testCase)
   % A jump that clears the RATE gate but is NOT grossly implausible, not near a
   % transition, and seasonally plausible is AMBIGUOUS: flagged, never
   % corrected. The series is returned unchanged.

   t = (datetime(2015, 7, 1):hours(1):datetime(2015, 7, 1, 5, 0, 0))';
   % A 0.5 m melt-season jump: above the hourly rate bound but below max_step.
   surf = [0; 0.01; 0.02; 0.52; 0.53; 0.54];

   [corrected, record, flags] = icemodel.forcing.destepSurface(t, surf, ...
      season="ablation");

   testCase.verifyEqual(numel(record), 1);
   testCase.verifyEqual(string(record.classification), "ambiguous");
   testCase.verifyFalse(record.applied);
   testCase.verifyEqual(corrected, surf, 'AbsTol', 1e-12);   % unchanged
   testCase.verifyEqual(flags.step_correctable(4), 0);
   testCase.verifyEqual(flags.step_detected(4), 1);
end

function test_destep_transition_corroborates_unambiguous(testCase)
   % A sub-max_step jump that coincides with a known station transition is
   % corroborated -> UNAMBIGUOUS, and corrected by the default mode.

   t = (datetime(2015, 7, 1):hours(1):datetime(2015, 7, 1, 5, 0, 0))';
   surf = [0; 0.01; 0.02; 0.52; 0.53; 0.54];   % 0.5 m jump at sample 4

   [corrected, record, ~] = icemodel.forcing.destepSurface(t, surf, ...
      transition_times=t(4), season="ablation");

   testCase.verifyEqual(string(record.classification), "unambiguous");
   testCase.verifyTrue(record.applied);
   testCase.verifyTrue(any(record.evidence == "transition"));
   % Leveled: the post-step segment drops by the 0.5 m offset.
   testCase.verifyEqual(corrected(end), 0.04, 'AbsTol', 1e-9);
end

function test_destep_detect_mode_leaves_series_unaltered(testCase)
   % mode="detect" never changes the series (staging-faithful contract): the
   % flags/record are produced but corrected == surf.

   t = (datetime(2015, 1, 1):hours(1):datetime(2015, 1, 1, 9, 0, 0))';
   surf = (0:9)' * 0.01;
   surf(6:end) = surf(6:end) + 5;

   [corrected, record, flags] = icemodel.forcing.destepSurface(t, surf, ...
      mode="detect");

   testCase.verifyEqual(corrected, surf, 'AbsTol', 1e-12);
   testCase.verifyFalse(any([record.applied]));
   testCase.verifyEqual(flags.step_detected(6), 1);   % still detected
end

function test_destep_gap_interior_is_not_a_step(testCase)
   % A large change accumulated ACROSS a NaN gap (a long slope-bridged interval)
   % is NOT a single-timestep step: between the two finite endpoints the rate is
   % plausible, so nothing is flagged.

   t = (datetime(2015, 6, 1):days(1):datetime(2015, 8, 9, 0, 0, 0))';
   surf = nan(numel(t), 1);
   surf(1) = 0;
   surf(end) = 2;   % +2 m over ~69 days -> ~0.03 m/day, plausible

   [corrected, record, ~] = icemodel.forcing.destepSurface(t, surf);

   testCase.verifyEmpty(record);
   testCase.verifyEqual(corrected, surf, 'AbsTol', 1e-12);   % untouched
end

function test_destep_gap_flagged_sample_is_excluded(testCase)
   % A sharp jump whose endpoint is gap-bridged (gap_flag==1) is an
   % interpolation artifact, not a direct-observation step. With gap_flag
   % supplied the bridged sample is excluded from detection; without it the
   % same spike would be detected.
   t = (datetime(2015, 1, 1):days(1):datetime(2015, 1, 10))';
   surf = (0:9)' * 0.01;          % smooth, plausible daily change
   surf(6) = surf(6) + 3;         % a +3 m spike landing on one sample
   gap = zeros(numel(t), 1);
   gap(6) = 1;                    % that sample is gap-bridged

   [~, rec_no_gap] = icemodel.forcing.destepSurface(t, surf, mode="detect");
   testCase.verifyNotEmpty(rec_no_gap);  % spike detected without the flag

   [~, rec_gap] = icemodel.forcing.destepSurface(t, surf, mode="detect", ...
      gap_flag=gap);
   testCase.verifyEmpty(rec_gap);        % excluded once flagged as bridged
end

function test_surfaceFlags_gap_from_sensors_not_just_z_nan(testCase)
   % surfaceFlags marks a sample gap-bridged when z is finite but every sensor
   % is NaN (slope-bridged), which the bare z-NaN heuristic misses.

   z = [0; 0.1; 0.2; 0.3; 0.4];        % all finite (slope-bridged interior)
   sensors = [0.0 NaN; 0.1 0.5; NaN NaN; 0.3 0.6; 0.4 NaN];
   t = (datetime(2015, 1, 1):hours(1):datetime(2015, 1, 1, 4, 0, 0))';

   flags = icemodel.forcing.helpers.surfaceFlags(z, sensors, t);

   % Sample 3 (both sensors NaN, z finite) is gap-bridged; samples with any
   % finite sensor are direct observations.
   testCase.verifyEqual(flags.gap, [0; 0; 1; 0; 0]);
   % z-NaN-only would have flagged none here, missing sample 3.
   testCase.verifyGreaterThan(nnz(flags.gap), nnz(~isfinite(z)));
end

function test_surfaceFlags_station_transition_window(testCase)
   % With a known handover time, samples within tol_days are flagged a station
   % transition; the gap flag is independent of it.

   t = (datetime(2015, 1, 1):days(1):datetime(2015, 1, 31))';
   z = (0:numel(t) - 1)' * 0.01;
   sensors = z;   % all observed (no gaps)

   flags = icemodel.forcing.helpers.surfaceFlags(z, sensors, t, ...
      transition_times=datetime(2015, 1, 15), tol_days=3);

   inwin = abs(days(t - datetime(2015, 1, 15))) <= 3;
   testCase.verifyEqual(flags.station_transition, double(inwin));
   testCase.verifyEqual(nnz(flags.gap), 0);   % all observed -> no gap
end
