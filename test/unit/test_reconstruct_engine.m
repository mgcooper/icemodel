function tests = test_reconstruct_engine
   %TEST_RECONSTRUCT_ENGINE Verify the reconstruction engine contracts.
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

%% fillShortGaps

function test_shortgaps_fills_interior_3h_gap(testCase)
   % A 3 h interior gap is bridged linearly between its finite anchors,
   % flagged in the fill mask, and audited as bounded_interp.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   truth = series.tair;
   x = truth;
   gap = (4000:4002).';                          % 3 h interior run
   x(gap) = NaN;

   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "tair");

   % The fill is the straight line between the two anchor samples.
   expected = interp1([gap(1) - 1; gap(end) + 1], ...
      [truth(gap(1) - 1); truth(gap(end) + 1)], gap);
   testCase.verifyEqual(returned(gap), expected, 'AbsTol', 1e-12);
   testCase.verifyEqual(find(filled), gap);
   % Native finite samples outside the gap are never modified.
   keep = true(size(x));
   keep(gap) = false;
   testCase.verifyEqual(returned(keep), truth(keep));
   % One audit row records the channel, duration, and method label.
   testCase.verifyEqual(numel(audit), 1);
   testCase.verifyEqual(audit{1}{1}, 'tair');
   testCase.verifyEqual(audit{1}{4}, 3);
   testCase.verifyEqual(audit{1}{5}, 'bounded_interp');
end

function test_shortgaps_cap_rejects_12h_gap(testCase)
   % A 12 h run exceeds the tier-1 wall-clock cap: the gap stays missing
   % for later tiers and nothing is audited.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   x = series.tair;
   gap = (4000:4011).';                          % 12 h run, over the cap
   x(gap) = NaN;

   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "tair");

   testCase.verifyTrue(all(isnan(returned(gap))));
   testCase.verifyFalse(any(filled));
   testCase.verifyEmpty(audit);
end

function test_shortgaps_season_boundary_rules(testCase)
   % State channels remain season-contained except the explicitly
   % evidence-backed RH and SWD calendar-label rules.
   times = (datetime(2020, 2, 29, 22, 0, 0, 'TimeZone', 'UTC'): ...
      hours(1):datetime(2020, 3, 1, 4, 0, 0, ...
      'TimeZone', 'UTC')).';
   x = (260:266).';
   gap = (3:5).';
   x(gap) = NaN;

   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "tair");

   testCase.verifyTrue(all(isnan(returned(gap))));
   testCase.verifyFalse(any(filled));
   testCase.verifyEmpty(audit);

   % D-48/D-50 apply RH's validated nine-hour bridge cap across the
   % nonphysical calendar-season label while other state channels remain
   % season-contained.
   rh_times = (datetime(2020, 2, 29, 19, 0, 0, 'TimeZone', 'UTC'): ...
      hours(1):datetime(2020, 3, 1, 5, 0, 0, ...
      'TimeZone', 'UTC')).';
   rh_gap = (2:10).';
   rh = [90; nan(9, 1); 99];
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(rh_times, rh, "rh", ...
      cap_hours=9);
   expected = linspace(90, 99, 11).';
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   testCase.verifyTrue(all(filled(rh_gap)));
   testCase.verifyEqual(audit{1}{5}, 'bounded_interp');

   % The primary CSI path admits the evidenced two-hour SWD exception.
   swd = 100 * ones(size(times));
   short_gap = (3:4).';
   swd(short_gap) = NaN;
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, swd, "swd", ...
      latitude=0, longitude=180, allow_swd_flux_fallback=true);
   testCase.verifyTrue(all(isfinite(returned(short_gap))));
   testCase.verifyTrue(all(filled(short_gap)));
   testCase.verifyTrue(contains(audit{1}{6}, 'csi'));

   % Force CSI undefined so the production flux fallback exercises the
   % same two-hour exception.
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, swd, "swd", ...
      latitude=0, longitude=180, toa_dark_wm2=2000, ...
      allow_swd_flux_fallback=true);
   testCase.verifyEqual(returned(short_gap), 100 * ones(2, 1));
   testCase.verifyTrue(all(filled(short_gap)));
   testCase.verifyTrue(contains(audit{1}{6}, 'flux-linear fallback'));

   % The exact three-hour edge class is also evidence-supported.
   swd = 100 * ones(size(times));
   swd(gap) = NaN;
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, swd, "swd", ...
      latitude=0, longitude=180, toa_dark_wm2=2000, ...
      allow_swd_flux_fallback=true);
   testCase.verifyEqual(returned(gap), 100 * ones(3, 1));
   testCase.verifyTrue(all(filled(gap)));
   testCase.verifyTrue(contains(audit{1}{6}, 'flux-linear fallback'));

   % The same SWD ceiling covers a gap whose missing samples straddle the
   % bookkeeping boundary itself.
   swd = 100 * ones(size(times));
   straddling_gap = (2:3).';
   swd(straddling_gap) = NaN;
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, swd, "swd", ...
      latitude=0, longitude=180, toa_dark_wm2=2000, ...
      allow_swd_flux_fallback=true);
   testCase.verifyEqual(returned(straddling_gap), 100 * ones(2, 1));
   testCase.verifyTrue(all(filled(straddling_gap)));
   testCase.verifyTrue(contains(audit{1}{6}, 'flux-linear fallback'));

   % A boundary CSI bridge cannot borrow an extraterrestrially impossible
   % anchor even when interpolation would bring every candidate in bounds.
   guard_times = (datetime(2020, 2, 29, 23, 0, 0, ...
      'TimeZone', 'UTC'):hours(1):datetime(2020, 3, 1, 2, 0, 0, ...
      'TimeZone', 'UTC')).';
   guarded = [1500; NaN; NaN; 0];
   guarded_gap = (2:3).';
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps( ...
      guard_times, guarded, "swd", latitude=0, longitude=180, ...
      allow_swd_flux_fallback=true);
   testCase.verifyTrue(all(isnan(returned(guarded_gap))));
   testCase.verifyFalse(any(filled));
   testCase.verifyEmpty(audit);

   % D-49 uses the full existing SWD cap at the boundary.
   long_gap = (3:6).';
   swd = 100 * ones(size(times));
   swd(long_gap) = NaN;
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, swd, "swd", ...
      latitude=0, longitude=180, toa_dark_wm2=2000, ...
      allow_swd_flux_fallback=true);
   testCase.verifyEqual(returned(long_gap), 100 * ones(4, 1));
   testCase.verifyTrue(all(filled(long_gap)));
   testCase.verifyTrue(contains(audit{1}{6}, 'flux-linear fallback'));

   % The boundary exception still cannot exceed SWD's nine-hour ceiling.
   long_times = (datetime(2020, 2, 29, 18, 0, 0, ...
      'TimeZone', 'UTC'):hours(1):datetime(2020, 3, 1, 5, 0, 0, ...
      'TimeZone', 'UTC')).';
   over_cap = 100 * ones(size(long_times));
   ten_hour_gap = (2:11).';
   over_cap(ten_hour_gap) = NaN;
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps( ...
      long_times, over_cap, "swd", latitude=0, longitude=180, ...
      toa_dark_wm2=2000, allow_swd_flux_fallback=true, cap_hours=9);
   testCase.verifyTrue(all(isnan(returned(ten_hour_gap))));
   testCase.verifyFalse(any(filled));
   testCase.verifyEmpty(audit);
end

function test_shortgaps_flux_fallback_closes_blocked_twilight_csi(testCase)
   % D-32 uses the scalar-valid local flux bridge when a low-sun CSI
   % candidate is unavailable or rejected by the diagnostic solar ceiling.
   times = (datetime(2019, 10, 12, 19, 45, 0, 'TimeZone', 'UTC'): ...
      minutes(15):datetime(2019, 10, 12, 21, 0, 0, ...
      'TimeZone', 'UTC')).';
   x = [206; nan(4, 1); 65];
   gap = (2:5).';

   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=64.4828, longitude=-49.5391, ...
      allow_swd_flux_fallback=true);

   expected = interp1([1; 6], [206; 65], gap);
   testCase.verifyEqual(returned(gap), expected, 'AbsTol', 1e-12);
   testCase.verifyTrue(all(filled(gap)));
   testCase.verifyTrue(contains(audit{1}{6}, 'flux-linear fallback'));
   testCase.verifyEqual(returned([1, 6]), x([1, 6]));

   % A scalar-finite but extraterrestrially impossible pair of daylight
   % anchors cannot bypass the same final guard.
   noon_times = (datetime(2019, 7, 1, 12, 0, 0, ...
      'TimeZone', 'UTC'):minutes(15):datetime(2019, 7, 1, 13, 15, 0, ...
      'TimeZone', 'UTC')).';
   impossible = [2000; nan(4, 1); 2000];
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps( ...
      noon_times, impossible, "swd", latitude=64.4828, ...
      longitude=-49.5391, allow_swd_flux_fallback=true);
   testCase.verifyTrue(all(isnan(returned(gap))));
   testCase.verifyFalse(any(filled(gap)));
   testCase.verifyEmpty(audit);
end

function test_twilight_climatology_fills_only_supported_edge_postings(testCase)
   % D-44 reuses the native day-of-year/posting climatology only for one
   % residual civil-twilight posting adjacent to known darkness.
   times = (datetime(2017, 1, 1, 0, 0, 0, 'TimeZone', 'UTC'): ...
      hours(1):datetime(2020, 12, 31, 23, 0, 0, ...
      'TimeZone', 'UTC')).';
   latitude = 66.48;
   longitude = -42.50;
   maximum_elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, latitude, longitude, hours(1));
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   dark = maximum_elevation <= bands.civil_twilight_deg;
   twilight = maximum_elevation > bands.civil_twilight_deg ...
      & maximum_elevation <= bands.horizon_deg;
   edge = false(size(times));
   edge(2:end - 1) = twilight(2:end - 1) ...
      & xor(dark(1:end - 2), dark(3:end));
   dawn = find(edge & [false; dark(1:end - 1)] ...
      & year(times) == 2020, 1);
   dusk = find(edge & [dark(2:end); false] ...
      & year(times) == 2020, 1);
   testCase.assertNotEmpty(dawn);
   testCase.assertNotEmpty(dusk);

   native = max(0, min(40, 4 * (maximum_elevation + 6)));
   targets = [dawn; dusk];
   x = native;
   x(targets) = NaN;
   native(targets) = NaN;
   expected = icemodel.forcing.reconstruct.climatologyFill( ...
      times, native, times(targets), min_support=30);

   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillTwilightClimatology( ...
      times, x, native, latitude, longitude);

   testCase.verifyEqual(returned(targets), expected, 'AbsTol', 1e-12);
   testCase.verifyTrue(all(filled(targets)));
   testCase.verifyEqual(numel(audit), 2);
   testCase.verifyTrue(all(cellfun( ...
      @(row) strcmp(row{5}, 'twilight_climatology'), audit)));
   keep = true(size(x));
   keep(targets) = false;
   testCase.verifyEqual(returned(keep), x(keep));

   % A true-dark posting and an under-supported twilight posting decline.
   refused = x;
   dark_target = find(dark, 1);
   refused(dark_target) = NaN;
   [returned, filled] = ...
      icemodel.forcing.reconstruct.fillTwilightClimatology( ...
      times, refused, native, latitude, longitude, min_support=10000);
   testCase.verifyTrue(isnan(returned(dark_target)));
   testCase.verifyTrue(all(isnan(returned(targets))));
   testCase.verifyFalse(any(filled));

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillTwilightClimatology( ...
      times(1:end - 1), x, native, latitude, longitude), ...
      'icemodel:reconstruct:fillTwilightClimatology:sizeMismatch');
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillTwilightClimatology( ...
      times(1), x(1), native(1), latitude, longitude), ...
      'icemodel:reconstruct:fillTwilightClimatology:missingCadence');
end

function test_shortgaps_csi_fill_preserves_diurnal_shape(testCase)
   % Shortwave built as a constant 0.7 clear-sky index of top-of-atmosphere
   % irradiance makes the CSI-space fill exactly predictable: interpolating
   % the constant index reproduces the diurnal truth a raw line would cut.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   elevation = icemodel.forcing.helpers.solarElevation(times, 67.0, -48.8);
   toa = 1361 * max(0, sind(elevation(:)));
   truth = 0.7 * toa;
   x = truth;
   % Mask five contiguous high-sun hours; solar noon at 48.8 W is ~15:15
   % UTC, so this window sits inside guaranteed midsummer daylight.
   gap = find(times >= datetime(2020, 6, 21, 13, 0, 0, ...
      'TimeZone', 'UTC') & times <= datetime(2020, 6, 21, 17, 0, 0, ...
      'TimeZone', 'UTC'));
   x(gap) = NaN;

   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8);

   % The CSI bridge reconverts to the exact diurnal truth.
   testCase.verifyEqual(returned(gap), truth(gap), 'AbsTol', 1e-6);
   testCase.verifyTrue(all(filled(gap)));
   % A raw line between the anchors undercuts the midday hump; the CSI
   % fill sits above it at the gap center, preserving the cycle shape.
   rawline = interp1([gap(1) - 1; gap(end) + 1], ...
      [truth(gap(1) - 1); truth(gap(end) + 1)], gap);
   testCase.verifyGreaterThan(returned(gap(3)), rawline(3));
   % Dark samples are native zeros: work-space NaN runs with no native
   % gap must pass through untouched rather than be "filled".
   keep = true(size(x));
   keep(gap) = false;
   testCase.verifyEqual(returned(keep), truth(keep));
   % The audit detail records the CSI variant of the method.
   testCase.verifyTrue(contains(audit{1}{6}, 'csi'));
end

function test_shortgaps_csi_does_not_bridge_across_night(testCase)
   % Darkness is a hard interpolation boundary: a daylight gap beginning at
   % dawn has no preceding daylight anchor and must remain for later tiers.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   elevation = icemodel.forcing.helpers.solarElevation(times, 67.0, -48.8);
   toa = 1361 * max(0, sind(elevation(:)));
   truth = 0.7 * toa;
   x = truth;
   % Equinox morning: night exists at 67 N in late March. The gap starts
   % at the day's first meaningful-sun sample (toa >= 10, the CSI mask
   % threshold), so the preceding masked night joins the work-space run.
   day_start = find(times >= datetime(2020, 3, 21, 0, 0, 0, ...
      'TimeZone', 'UTC') & toa >= 10, 1, 'first');
   gap = (day_start:day_start + 2).';
   x(gap) = NaN;

   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8);

   testCase.verifyFalse(any(filled(gap)));
   testCase.verifyTrue(all(isnan(returned(gap))));
   testCase.verifyEmpty(audit);
end

function test_shortgaps_post_final_flux_fallback_closes_dawn_sliver(testCase)
   % D-32: after later tiers supply finite neighbors, a short dawn sliver
   % for which CSI has no preceding daylight denominator uses a local
   % flux-linear bridge. Deep darkness remains a hard boundary.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   elevation = icemodel.forcing.helpers.solarElevation(times, 67.0, -48.8);
   toa = icemodel.forcing.reconstruct.toaIrradiance( ...
      times, 67.0, -48.8);
   truth = 0.7 * toa;
   day_start = find(times >= datetime(2020, 3, 21, 0, 0, 0, ...
      'TimeZone', 'UTC') & toa >= 10, 1, 'first');
   gap = (day_start:day_start + 2).';
   x = truth;
   x(gap) = NaN;

   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8, ...
      allow_swd_flux_fallback=true);

   expected = interp1([gap(1) - 1; gap(end) + 1], ...
      truth([gap(1) - 1; gap(end) + 1]), gap);
   testCase.verifyEqual(returned(gap), expected, 'AbsTol', 1e-12);
   testCase.verifyTrue(all(filled(gap)));
   testCase.verifyTrue(contains(audit{1}{6}, 'flux-linear fallback'));

   % Enabling the fallback must not turn a polar-night hole into light.
   dark = find(elevation < ...
      icemodel.forcing.reconstruct.solarElevationBands().civil_twilight_deg, ...
      1);
   x = truth;
   x(dark) = NaN;
   [returned, filled] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8, ...
      allow_swd_flux_fallback=true);
   testCase.verifyTrue(isnan(returned(dark)));
   testCase.verifyFalse(filled(dark));

   % The one-posting darkness edge remains reserved for the separately
   % validated twilight-climatology tier.
   interval_maximum = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, 67.0, -48.8, hours(1));
   civil = ...
      icemodel.forcing.reconstruct.solarElevationBands().civil_twilight_deg;
   twilight_start = find(interval_maximum(2:end) >= civil ...
      & interval_maximum(1:end - 1) < civil, 1) + 1;
   x = truth;
   x(twilight_start) = NaN;
   [returned, filled] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8, ...
      allow_swd_flux_fallback=true);
   testCase.verifyTrue(isnan(returned(twilight_start)));
   testCase.verifyFalse(filled(twilight_start));

   % Two or more still-sunlit postings may use one known darkness zero and
   % the opposite daylight anchor's CSI. This preserves the solar curve
   % without filling any deep-darkness sample.
   location = struct('lat_wgs84', 67.0, 'lon_wgs84', -48.8);
   [truth_csi, ~] = ...
      icemodel.forcing.reconstruct.clearSkyIndex(times, truth, location);
   darkness_anchor = find( ...
      interval_maximum(1:end - 3) <= civil ...
      & interval_maximum(2:end - 2) > civil ...
      & interval_maximum(3:end - 1) > civil ...
      & isfinite(truth_csi(4:end)), 1);
   testCase.assertNotEmpty(darkness_anchor);
   twilight_gap = darkness_anchor + (1:2).';
   x = truth;
   x(twilight_gap) = NaN;
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8, ...
      allow_swd_flux_fallback=true);
   testCase.verifyEqual( ...
      returned(twilight_gap), truth(twilight_gap), 'AbsTol', 1e-12);
   testCase.verifyTrue(all(filled(twilight_gap)));
   testCase.verifyTrue(contains(audit{1}{6}, 'one-sided CSI'));

   % Exercise the reverse (sunset) direction explicitly.
   sunset_gap_start = find( ...
      isfinite(truth_csi(1:end - 3)) ...
      & interval_maximum(2:end - 2) > civil ...
      & interval_maximum(3:end - 1) > civil ...
      & interval_maximum(4:end) <= civil, 1) + 1;
   testCase.assertNotEmpty(sunset_gap_start);
   sunset_gap = sunset_gap_start + (0:1).';
   x = truth;
   x(sunset_gap) = NaN;
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8, ...
      allow_swd_flux_fallback=true);
   testCase.verifyEqual( ...
      returned(sunset_gap), truth(sunset_gap), 'AbsTol', 1e-12);
   testCase.verifyTrue(all(filled(sunset_gap)));
   testCase.verifyTrue(contains(audit{1}{6}, 'one-sided CSI'));

   % The caller's CSI-darkness threshold remains authoritative.
   before_sunset = sunset_gap(1) - 1;
   x = truth;
   x(sunset_gap) = NaN;
   [returned, filled] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8, ...
      toa_dark_wm2=toa(before_sunset) + 1, ...
      allow_swd_flux_fallback=true);
   testCase.verifyTrue(all(isnan(returned(sunset_gap))));
   testCase.verifyFalse(any(filled(sunset_gap)));

   % A finite but over-ceiling daylight anchor must not transfer an
   % impossible CSI into missing dawn postings.
   ceiling_gap = darkness_anchor + (1:3).';
   after_dawn = ceiling_gap(end) + 1;
   x = truth;
   x(after_dawn) = 1000 * toa(after_dawn);
   x(ceiling_gap) = NaN;
   [returned, filled] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8, ...
      allow_swd_flux_fallback=true);
   testCase.verifyTrue(all(isnan(returned(ceiling_gap))));
   testCase.verifyFalse(any(filled(ceiling_gap)));
end

function test_shortgaps_swd_override_closes_nine_hour_sliver(testCase)
   % D-39: observed-only holdouts support a narrow nine-hour SWD
   % extension; the ordinary six-hour cap must still leave this run alone.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   toa = icemodel.forcing.reconstruct.toaIrradiance( ...
      times, 67.0, -48.8);
   truth = 0.7 * toa;
   day_start = find(times >= datetime(2020, 3, 21, 0, 0, 0, ...
      'TimeZone', 'UTC') & toa >= 10, 1, 'first');
   gap = (day_start:day_start + 8).';
   x = truth;
   x(gap) = NaN;
   expected = interp1([gap(1) - 1; gap(end) + 1], ...
      truth([gap(1) - 1; gap(end) + 1]), gap);

   [six_hour, filled_six] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8, cap_hours=6, ...
      allow_swd_flux_fallback=true);
   testCase.verifyTrue(all(isnan(six_hour(gap))));
   testCase.verifyFalse(any(filled_six(gap)));

   [eight_hour, filled_eight] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8, cap_hours=8, ...
      allow_swd_flux_fallback=true);
   testCase.verifyTrue(all(isnan(eight_hour(gap))));
   testCase.verifyFalse(any(filled_eight(gap)));

   [nine_hour, filled_nine, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8, cap_hours=9, ...
      allow_swd_flux_fallback=true);
   testCase.verifyEqual(nine_hour(gap), expected, 'AbsTol', 1e-12);
   testCase.verifyTrue(all(filled_nine(gap)));
   testCase.verifyTrue(contains(audit{1}{6}, 'cap 9 h'));

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillShortGaps( ...
      times, x, "tair", cap_hours=9), ...
      'icemodel:reconstruct:fillShortGaps:capHours');
end

function test_shortgaps_albedo_override_closes_day_scale_sliver(testCase)
   % D-41: albedo seam remnants use the observed-only-supported
   % 30-hour linear bridge and cannot silently widen beyond it.
   times = (datetime(2020, 7, 1, 0, 0, 0, 'TimeZone', 'UTC'): ...
      hours(1):datetime(2020, 7, 2, 3, 0, 0, ...
      'TimeZone', 'UTC')).';
   x = [0.50; nan(26, 1); 0.56];
   expected = linspace(0.50, 0.56, 28).';
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps( ...
      times, x, "albedo", cap_hours=30);
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   testCase.verifyTrue(all(filled(2:27)));
   testCase.verifyTrue(contains(audit{1}{6}, 'cap 30 h'));
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillShortGaps( ...
      times, x, "albedo", cap_hours=31), ...
      'MATLAB:validators:mustBeLessThanOrEqual');
end

function test_shortgaps_rh_override_closes_nine_hour_sliver(testCase)
   % D-42/D-50: observed-only holdouts support a nine-hour RH bridge, while
   % longer outages remain outside tier 1.
   times = (datetime(2020, 7, 1, 0, 0, 0, 'TimeZone', 'UTC'): ...
      hours(1):datetime(2020, 7, 1, 10, 0, 0, ...
      'TimeZone', 'UTC')).';
   x = [60; nan(9, 1); 78];
   expected = linspace(60, 78, 11).';
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps( ...
      times, x, "rh", cap_hours=9);
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   testCase.verifyTrue(all(filled(2:10)));
   testCase.verifyTrue(contains(audit{1}{6}, 'cap 9 h'));
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillShortGaps( ...
      times, x, "rh", cap_hours=10), ...
      'icemodel:reconstruct:fillShortGaps:capHours');

   % The approved nine-hour cap itself, not just option validation,
   % leaves a ten-hour outage for later tiers.
   long_times = (datetime(2020, 7, 1, 0, 0, 0, 'TimeZone', 'UTC'): ...
      hours(1):datetime(2020, 7, 1, 11, 0, 0, ...
      'TimeZone', 'UTC')).';
   long_x = [60; nan(10, 1); 78];
   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps( ...
      long_times, long_x, "rh", cap_hours=9);
   testCase.verifyTrue(all(isnan(returned(2:11))));
   testCase.verifyFalse(any(filled));
   testCase.verifyEmpty(audit);
end

function test_shortgaps_csi_caps_separate_daylight_gaps(testCase)
   % Two short native gaps joined only by finite night values are separate
   % cap decisions even though CSI masks the whole night between them.
   times = (datetime(2020, 3, 20, 0, 0, 0, 'TimeZone', 'UTC'): ...
      hours(1):datetime(2020, 3, 23, 0, 0, 0, ...
      'TimeZone', 'UTC')).';
   toa = icemodel.forcing.reconstruct.toaIrradiance( ...
      times, 67.0, -48.8);
   truth = 0.7 * toa;
   daylight = toa >= 10;
   edge = diff([false; daylight; false]);
   day_starts = find(edge == 1);
   day_stops = find(edge == -1) - 1;
   gap = [day_stops(1) - 3; day_stops(1) - 2; ...
      day_starts(2) + 2; day_starts(2) + 3];
   x = truth;
   x(gap) = NaN;

   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "swd", ...
      latitude=67.0, longitude=-48.8);

   testCase.verifyTrue(all(filled(gap)));
   testCase.verifyTrue(all(isfinite(returned(gap))));
   testCase.verifyTrue(all(returned(gap) >= 0 ...
      & returned(gap) <= 1.05 * toa(gap)));
   testCase.verifyEqual(numel(audit), 2);
   testCase.verifyEqual([audit{1}{4}; audit{2}{4}], [2; 2]);
end

function test_shortgaps_swd_requires_solar_geometry(testCase)
   % The CSI variant cannot run without the site point: omitting latitude
   % and longitude for swd fails loudly instead of raw-interpolating.
   series = icemodel.test.fixtures.makeReconstructSeries();
   testCase.verifyError(@() icemodel.forcing.reconstruct.fillShortGaps( ...
      series.Properties.RowTimes, zeros(height(series), 1), "swd"), ...
      'icemodel:reconstruct:fillShortGaps:missingSolarGeometry');
end

function test_shortgaps_rejects_out_of_bounds_fill(testCase)
   % Anchors above the tair upper bound make the bridge unphysical: the
   % jump check passes (tiny anchor-to-anchor step) but the bounds check
   % leaves the gap missing for the next tier.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   x = series.tair;
   gap = (4000:4002).';
   x(gap(1) - 1) = 302;                          % suspect native anchors
   x(gap(end) + 1) = 305;                        % above the 300 K bound
   x(gap) = NaN;

   [returned, filled] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "tair");

   testCase.verifyTrue(all(isnan(returned(gap))));
   testCase.verifyFalse(any(filled));
   % The suspect native anchors themselves are never modified.
   testCase.verifyEqual(returned(gap(1) - 1), 302);
end

function test_shortgaps_bridge_floor_handles_native_jump(testCase)
   % A large but in-bounds native anchor step raises the feasible seam
   % floor; the linear bridge is retained instead of refused.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   truth = series.tair;
   x = truth;
   gap = (500:502).';                            % January, ~252 K locally
   x(gap(end) + 1) = truth(gap(end) + 1) + 40;   % in-bounds but jumpy
   x(gap) = NaN;

   [returned, filled] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "tair");

   expected = interp1([gap(1) - 1; gap(end) + 1], ...
      x([gap(1) - 1; gap(end) + 1]), gap);
   testCase.verifyEqual(returned(gap), expected, 'AbsTol', 1e-12);
   testCase.verifyTrue(all(filled(gap)));
end

function test_step_scale_is_zero_for_constant_channel(testCase)
   % Finite zero steps define an exact-continuity scale of zero, not NaN.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:99)).';

   scale = icemodel.forcing.reconstruct.stepScale(times, ...
      7 * ones(numel(times), 1));

   testCase.verifyEqual([scale.DJF, scale.MAM, scale.JJA, scale.SON], ...
      zeros(1, 4));
end

function test_step_scale_is_zero_without_adjacent_finite_samples(testCase)
   % Sparse observations with no adjacent pair must produce a finite
   % fail-closed seam scale rather than allowing comparisons against NaN.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:99)).';
   x = nan(numel(times), 1);
   x([1, 3]) = [1, 9];

   scale = icemodel.forcing.reconstruct.stepScale(times, x);

   testCase.verifyEqual([scale.DJF, scale.MAM, scale.JJA, scale.SON], ...
      zeros(1, 4));
end

function test_step_scale_retains_zero_observed_steps(testCase)
   % A mostly flat quantized record must have the documented median step of
   % zero; one rare change cannot inflate the seam threshold.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:99)).';
   x = zeros(numel(times), 1);
   x(end) = 10;

   scale = icemodel.forcing.reconstruct.stepScale(times, x);

   testCase.verifyEqual(scale.DJF, 0);
end

function test_shortgaps_skips_leading_and_trailing_gaps(testCase)
   % Edge gaps have no interior anchors: both stay missing for later
   % tiers and the fill mask and audit stay empty.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   x = series.tair;
   x(1:3) = NaN;                                 % leading run
   x(end - 2:end) = NaN;                         % trailing run

   [returned, filled, audit] = ...
      icemodel.forcing.reconstruct.fillShortGaps(times, x, "tair");

   testCase.verifyTrue(all(isnan(returned([1:3, ...
      numel(x) - 2:numel(x)]))));
   testCase.verifyFalse(any(filled));
   testCase.verifyEmpty(audit);
end

%% climatologyFill

function test_climatology_median_honors_doy_window(testCase)
   % A day-of-year ramp makes the pooled median analytic: the default
   % 7-day half-window pools 15 same-hour days centered on the query, and
   % a narrower window shrinks the pool while keeping the centered median.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   x = day(times, 'dayofyear');
   query = datetime(2020, 4, 9, 6, 0, 0, 'TimeZone', 'UTC');   % doy 100

   [returned, n_support] = ...
      icemodel.forcing.reconstruct.climatologyFill(times, x, query);
   testCase.verifyEqual(returned, 100);
   testCase.verifyEqual(n_support, 15);

   [narrow, n_narrow] = icemodel.forcing.reconstruct.climatologyFill( ...
      times, x, query, window_days=2);
   testCase.verifyEqual(narrow, 100);
   testCase.verifyEqual(n_narrow, 5);
end

function test_climatology_diurnal_pooling_toggle(testCase)
   % An hour-of-day ramp separates the pooling modes: diurnal pooling
   % returns the query hour, pooling all hours returns the daily median.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   x = hour(times);
   query = datetime(2020, 4, 9, 6, 0, 0, 'TimeZone', 'UTC');

   [returned, n_support] = ...
      icemodel.forcing.reconstruct.climatologyFill(times, x, query);
   testCase.verifyEqual(returned, 6);
   testCase.verifyEqual(n_support, 15);

   [pooled, n_pooled] = icemodel.forcing.reconstruct.climatologyFill( ...
      times, x, query, diurnal=false);
   testCase.verifyEqual(pooled, 11.5);           % median of hours 0..23
   testCase.verifyEqual(n_pooled, 360);          % 15 days x 24 hours
end

function test_climatology_preserves_subhour_postings(testCase)
   % Half-hour K-transect postings retain distinct diurnal bins rather than
   % collapsing both samples into the same whole-hour climatology.
   days_in = (1:15).';
   times = sort([ ...
      datetime(2020, 1, days_in, 0, 0, 0, 'TimeZone', 'UTC'); ...
      datetime(2020, 1, days_in, 0, 30, 0, 'TimeZone', 'UTC')]);
   x = minute(times);
   query = [ ...
      datetime(2020, 1, 8, 0, 0, 0, 'TimeZone', 'UTC'); ...
      datetime(2020, 1, 8, 0, 30, 0, 'TimeZone', 'UTC')];

   [returned, n_support] = ...
      icemodel.forcing.reconstruct.climatologyFill(times, x, query);

   testCase.verifyEqual(returned, [0; 30]);
   testCase.verifyEqual(n_support, [15; 15]);
end

function test_climatology_normalizes_post_leap_calendar_days(testCase)
   % The same March calendar window must pool the same dates in leap and
   % non-leap years; raw day-of-year bins would return 10.5 and 9.5 here.
   days_in = (8:12).';
   times = [datetime(2020, 3, days_in, 6, 0, 0, 'TimeZone', 'UTC'); ...
      datetime(2021, 3, days_in, 6, 0, 0, 'TimeZone', 'UTC')];
   x = [days_in; days_in];
   query = [datetime(2020, 3, 10, 6, 0, 0, 'TimeZone', 'UTC'); ...
      datetime(2021, 3, 10, 6, 0, 0, 'TimeZone', 'UTC')];

   [returned, n_support] = ...
      icemodel.forcing.reconstruct.climatologyFill(times, x, query, ...
      window_days=1);

   testCase.verifyEqual(returned, [10; 10]);
   testCase.verifyEqual(n_support, [6; 6]);
end

function test_climatology_min_support_returns_nan(testCase)
   % Insufficient pooled support yields NaN while still reporting the
   % count so callers can diagnose why the estimate declined.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   x = day(times, 'dayofyear');
   query = datetime(2020, 4, 9, 6, 0, 0, 'TimeZone', 'UTC');

   [returned, n_support] = icemodel.forcing.reconstruct.climatologyFill( ...
      times, x, query, min_support=16);          % pool only holds 15

   testCase.verifyTrue(isnan(returned));
   testCase.verifyEqual(n_support, 15);
end

%% elevationAdjust

function test_elevation_tair_lapse(testCase)
   % A 500 m rise cools the donor by the default lapse rate.
   returned = icemodel.forcing.reconstruct.elevationAdjust( ...
      "tair", [250; 260], 500);
   expected = [250; 260] - 3;                    % -0.0060 K/m x 500 m
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
end

function test_elevation_psfc_barometric(testCase)
   % Pressure follows the barometric factor, with the scale height taken
   % from the fallback temperature by default and a supplied coincident
   % temperature when the caller has one.
   % Constants come from the canonical physical-constant source so the
   % expectation can never drift from the implementation.
   [Rd, g] = icemodel.physicalConstant('Rd', 'gravity');
   x = [80000; 82000];
   returned = icemodel.forcing.reconstruct.elevationAdjust("psfc", x, 1000);
   expected = x * exp(-1000 * g / (Rd * 255));
   testCase.verifyEqual(returned, expected, 'RelTol', 1e-12);

   warm = icemodel.forcing.reconstruct.elevationAdjust("psfc", x, 1000, ...
      tair_for_pressure=270);
   expected_warm = x * exp(-1000 * g / (Rd * 270));
   testCase.verifyEqual(warm, expected_warm, 'RelTol', 1e-12);
end

function test_elevation_threshold_noop(testCase)
   % |dz| at or below the threshold applies no adjustment at all.
   returned = icemodel.forcing.reconstruct.elevationAdjust( ...
      "tair", [250; 260], 100);
   testCase.verifyEqual(returned, [250; 260]);
end

function test_elevation_passthrough_channel(testCase)
   % Channels without approved elevation dependence pass through even for
   % large elevation differences.
   returned = icemodel.forcing.reconstruct.elevationAdjust( ...
      "wspd", [5; 8], 1000);
   testCase.verifyEqual(returned, [5; 8]);
end

%% fitDonorTransfer / applyDonorTransfer

function test_donor_linear_recovers_known_transfer(testCase)
   % An exactly affine target recovers the slope and intercept, and the
   % applied transfer reproduces the target over the fitted range.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   donor = series.tair;
   target = 1.2 * donor + 5;

   transfer = icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, target, donor, "tair", fit_years=2020, lag_search=false);

   testCase.verifyEqual(transfer.lag_hours, 0);
   testCase.verifyEqual(transfer.models.JJA.slope, 1.2, 'AbsTol', 1e-9);
   testCase.verifyEqual(transfer.models.JJA.intercept, 5, 'AbsTol', 1e-6);
   testCase.verifyEqual(transfer.n_overlap, numel(donor));

   returned = icemodel.forcing.reconstruct.applyDonorTransfer( ...
      transfer, times, donor);
   testCase.verifyEqual(returned, target, 'AbsTol', 1e-6);
end

function test_donor_fits_per_season_models(testCase)
   % A summer-only offset appears in the JJA fit while MAM stays neutral;
   % seasons under a quarter of the overlap requirement (DJF and SON in a
   % leap year, 2184 h each) inherit the identical all-season fit.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   donor = series.tair;
   season = icemodel.forcing.reconstruct.seasonOf(times);
   target = donor + 3 * (season == "JJA");

   transfer = icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, target, donor, "tair", fit_years=2020, lag_search=false);

   testCase.verifyEqual(transfer.models.JJA.slope, 1, 'AbsTol', 1e-9);
   testCase.verifyEqual(transfer.models.JJA.intercept, 3, 'AbsTol', 1e-6);
   testCase.verifyEqual(transfer.models.MAM.intercept, 0, 'AbsTol', 1e-6);
   testCase.verifyEqual(transfer.models.DJF, transfer.models.SON);
end

function test_donor_lag_search_detects_shift(testCase)
   % A donor leading the target by 3 h is detected by the lag search; a
   % stricter minimum gain rejects the same lag and records zero.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   n = numel(times);
   % Diurnal-dominant signal so a 3 h misalignment costs ~0.29 of the
   % overlap correlation (cos 45 deg), well above the default gain gate.
   target = 250 + 20 * sin(2 * pi * (0:n - 1).' / 24);
   donor = [target(4:end); nan(3, 1)];           % donor leads by 3 h

   detected = icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, target, donor, "tair", fit_years=2020);
   testCase.verifyEqual(detected.lag_hours, 3);
   testCase.verifyGreaterThan(detected.overlap_r_lagged, 0.999);
   testCase.verifyEqual(detected.models.JJA.slope, 1, 'AbsTol', 1e-6);

   % A gain requirement above the achievable ~0.29 keeps the lag at zero.
   rejected = icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, target, donor, "tair", fit_years=2020, min_lag_gain=0.5);
   testCase.verifyEqual(rejected.lag_hours, 0);
   testCase.verifyEqual(rejected.overlap_r_lagged, rejected.overlap_r_base);
end

function test_donor_insufficient_overlap_errors(testCase)
   % Overlap below the one-year policy requirement fails loudly instead
   % of fitting an under-constrained transfer.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   donor = series.tair;
   target = donor;
   target(101:end) = NaN;                        % only 100 h of overlap

   testCase.verifyError(@() icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, target, donor, "tair", fit_years=2020, lag_search=false), ...
      'icemodel:reconstruct:fitDonorTransfer:insufficientOverlap');
end

function test_donor_coarse_short_support_counts_elapsed_hours(testCase)
   % Support-holding ninety days of hourly donor postings onto a 15-minute
   % target creates four rows per posting but still represents only 2160 h.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC'):minutes(15): ...
      datetime(2020, 12, 31, 23, 45, 0, 'TimeZone', 'UTC')).';
   target = 260 + sin((1:numel(times)).' / 100);
   donor = nan(size(target));
   n_support = 90 * 24 * 4;
   donor(1:n_support) = repelem(target(1:4:n_support), 4);

   transfer = icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, target, donor, "tair", fit_years=2020, ...
      lag_search=false, min_overlap_hours=2000);
   testCase.verifyEqual(transfer.n_overlap_hours, 90 * 24);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, target, donor, "tair", fit_years=2020, ...
      lag_search=false, min_overlap_hours=8000), ...
      'icemodel:reconstruct:fitDonorTransfer:insufficientOverlap');
end

function test_donor_short_record_lag_search_fails_as_overlap(testCase)
   % A lag wider than the record must decline for overlap, not index past it.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:3)).';
   target = (1:4).';
   donor = target;

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, target, donor, "tair", fit_years=2020, ...
      max_lag_hours=18, min_overlap_hours=5), ...
      'icemodel:reconstruct:fitDonorTransfer:insufficientOverlap');
end

function test_donor_single_sample_fails_as_overlap(testCase)
   % One timestamp cannot define cadence and must use the documented fit error.
   times = datetime(2020, 1, 1, 'TimeZone', 'UTC');

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, 1, 1, "tair", fit_years=2020), ...
      'icemodel:reconstruct:fitDonorTransfer:insufficientOverlap');
end

function test_donor_apply_short_axis_returns_missing(testCase)
   % A fitted lag wider than a later application axis retains no samples.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:3)).';
   model = struct('kind', "linear", 'slope', 1, 'intercept', 0);
   models = struct('DJF', model, 'MAM', model, 'JJA', model, 'SON', model);
   transfer = struct('lag_hours', 5, 'donor_range', [0, 10], ...
      'models', models);

   returned = icemodel.forcing.reconstruct.applyDonorTransfer( ...
      transfer, times, (1:4).');

   testCase.verifySize(returned, [4, 1]);
   testCase.verifyTrue(all(isnan(returned)));
end

function test_swd_donor_transfer_uses_station_specific_csi(testCase)
   % Identical cloud transmissivity at distinct stations must transfer
   % exactly despite their different raw TOA flux cycles.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   target_location = struct('lat_wgs84', 67.0, 'lon_wgs84', -48.8);
   donor_location = struct('lat_wgs84', 64.0, 'lon_wgs84', -42.0);
   target_toa = icemodel.forcing.reconstruct.toaIrradiance(times, ...
      target_location.lat_wgs84, target_location.lon_wgs84);
   donor_toa = icemodel.forcing.reconstruct.toaIrradiance(times, ...
      donor_location.lat_wgs84, donor_location.lon_wgs84);
   transmissivity = 0.55 + 0.1 * sin((1:numel(times)).' / 200);
   target = transmissivity .* target_toa;
   donor = transmissivity .* donor_toa;

   transfer = icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, target, donor, "swd", fit_years=2020, lag_search=false, ...
      min_overlap_hours=1000, target_location=target_location, ...
      donor_location=donor_location, toa_dark_wm2=20);
   returned = icemodel.forcing.reconstruct.applyDonorTransfer( ...
      transfer, times, donor);

   lit = isfinite(returned);
   testCase.verifyEqual(string(transfer.transfer_space), "clear_sky_index");
   testCase.verifyGreaterThan(mean(abs(target(lit) - donor(lit))), 1);
   testCase.verifyEqual(returned(lit), target(lit), 'AbsTol', 1e-8);
end

function test_donor_apply_one_sample_respects_lag(testCase)
   % A singleton application axis supports a zero-lag transfer but cannot
   % silently treat a fitted nonzero lag as zero.
   time = datetime(2020, 1, 1, 'TimeZone', 'UTC');
   model = struct('kind', "linear", 'slope', 2, 'intercept', 1);
   models = struct('DJF', model, 'MAM', model, 'JJA', model, 'SON', model);
   transfer = struct('lag_hours', 1, 'donor_range', [0, 10], ...
      'models', models);

   shifted = icemodel.forcing.reconstruct.applyDonorTransfer( ...
      transfer, time, 3);
   transfer.lag_hours = 0;
   unshifted = icemodel.forcing.reconstruct.applyDonorTransfer( ...
      transfer, time, 3);

   testCase.verifyTrue(isnan(shifted));
   testCase.verifyEqual(unshifted, 7);
end

function test_donor_piecewise_recovers_monotone_map(testCase)
   % A quadratic, monotone-over-range target is captured by the knotted
   % monotone map: node values are sorted and the applied estimate tracks
   % the curve.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   donor = series.tair;
   % Monotone because donor > 230 K everywhere in the fixture.
   target = 200 + 0.05 * (donor - 230).^2;

   transfer = icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, target, donor, "tair", fit_years=2020, lag_search=false, ...
      knots=6);

   testCase.verifyEqual(transfer.models.JJA.kind, "piecewise");
   testCase.verifyTrue(issorted(transfer.models.JJA.node_values));

   returned = icemodel.forcing.reconstruct.applyDonorTransfer( ...
      transfer, times, donor);
   testCase.verifyTrue(all(isfinite(returned)));
   testCase.verifyLessThan(mean(abs(returned - target)), 3);
end

function test_donor_piecewise_degenerate_falls_back_linear(testCase)
   % A two-valued donor collapses the quantile breakpoints to one bin;
   % the degenerate map falls back to the linear fit instead of crashing
   % or fabricating a flat transfer.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   n = numel(times);
   donor = repmat(270, n, 1);
   donor(1:50:end) = 260;                        % sparse second level
   target = donor + 1;

   transfer = icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, target, donor, "tair", fit_years=2020, lag_search=false, ...
      knots=6);

   testCase.verifyEqual(transfer.models.JJA.kind, "linear");
   testCase.verifyEqual(transfer.models.JJA.slope, 1, 'AbsTol', 1e-9);
   testCase.verifyEqual(transfer.models.JJA.intercept, 1, 'AbsTol', 1e-6);
end

function test_donor_apply_enforces_extrapolation_limit(testCase)
   % Donor samples beyond 110% of the fitted range return NaN; samples
   % just outside the fitted range but inside the widened window still
   % transfer, and missing donors stay missing.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   donor = series.tair;
   target = 1.2 * donor + 5;
   transfer = icemodel.forcing.reconstruct.fitDonorTransfer( ...
      times, target, donor, "tair", fit_years=2020, lag_search=false);

   span = diff(transfer.donor_range);
   inside = transfer.donor_range(2) + 0.05 * span;   % inside the +10% edge
   outside = transfer.donor_range(2) + 0.20 * span;  % beyond the +10% edge
   returned = icemodel.forcing.reconstruct.applyDonorTransfer( ...
      transfer, times(1:3), [inside; outside; NaN]);

   testCase.verifyEqual(returned(1), 1.2 * inside + 5, 'AbsTol', 1e-6);
   testCase.verifyTrue(all(isnan(returned(2:3))));
end

%% fitProxyCalibration / applyProxyCalibration

function test_proxy_recovers_additive_bias(testCase)
   % A constant -2.5 model bias is recovered as a +2.5 additive median
   % correction per season and removed exactly on application.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   model = series.lwd;
   obs = model + 2.5;

   calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
      times, obs, model, "lwd", fit_years=2020);

   testCase.verifyEqual(calibration.mode, "additive");
   testCase.verifyFalse(calibration.identity);
   testCase.verifyEqual(calibration.corrections.JJA, 2.5, 'AbsTol', 1e-12);

   returned = icemodel.forcing.reconstruct.applyProxyCalibration( ...
      calibration, times, model);
   testCase.verifyEqual(returned, obs, 'AbsTol', 1e-12);
end

function test_proxy_multiplicative_swd_ratio(testCase)
   % Shortwave calibrates by ratio only where target TOA is meaningful,
   % regardless of a bright proxy during target-station darkness.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   % Idealized diurnal shortwave: daylight sinusoid, zero at night.
   model = max(0, 600 * sin(2 * pi * (hour(times) - 6) / 24));
   obs = 0.8 * model;
   target_toa = icemodel.forcing.reconstruct.toaIrradiance( ...
      times, 0, 120);
   dark_proxy = target_toa < 10 & model > 10;
   testCase.verifyTrue(any(dark_proxy));
   obs(dark_proxy) = 2 * model(dark_proxy);

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fitProxyCalibration( ...
      times, obs, model, "swd", fit_years=2020), ...
      'icemodel:reconstruct:fitProxyCalibration:targetToaRequired');

   calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
      times, obs, model, "swd", fit_years=2020, target_toa=target_toa);

   testCase.verifyEqual(calibration.mode, "multiplicative");
   testCase.verifyEqual(calibration.corrections.DJF, 0.8, 'AbsTol', 1e-12);
   expected_overlap = year(times) == 2020 & target_toa >= 10 & model > 0;
   testCase.verifyEqual(calibration.n_overlap, nnz(expected_overlap));

   returned = icemodel.forcing.reconstruct.applyProxyCalibration( ...
      calibration, times, model);
   testCase.verifyEqual(returned(~dark_proxy), obs(~dark_proxy), ...
      'AbsTol', 1e-9);
   testCase.verifyEqual(returned(dark_proxy), ...
      0.8 * model(dark_proxy), 'AbsTol', 1e-9);
end

function test_proxy_multiplicative_wind_preserves_support(testCase)
   % Wind speed scales by a positive overlap ratio instead of subtracting
   % an additive bias that can turn physically valid low winds negative.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   model = linspace(0.1, 10, height(series)).';
   obs = 0.4 * model;

   calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
      times, obs, model, "wspd", fit_years=2020);
   returned = icemodel.forcing.reconstruct.applyProxyCalibration( ...
      calibration, times, model);

   testCase.verifyEqual(calibration.mode, "multiplicative");
   testCase.verifyEqual(returned, obs, 'AbsTol', 1e-12);
   testCase.verifyGreaterThanOrEqual(min(returned), 0);
end

function test_proxy_seasonal_fallback_to_annual(testCase)
   % A season with too few overlap samples inherits the annual correction
   % instead of trusting its thin local bias.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   model = series.lwd;
   season = icemodel.forcing.reconstruct.seasonOf(times);
   obs = nan(numel(model), 1);
   % JJA carries a +3 bias with full-season support; DJF carries a +1
   % bias on only 100 samples, below the 300-sample seasonal minimum.
   jja = find(season == "JJA");
   obs(jja) = model(jja) + 3;
   djf = find(season == "DJF");
   obs(djf(1:100)) = model(djf(1:100)) + 1;

   calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
      times, obs, model, "lwd", fit_years=2020);

   testCase.verifyEqual(calibration.corrections.JJA, 3, 'AbsTol', 1e-12);
   % The annual median is +3 (JJA dominates), so thin DJF inherits it.
   testCase.verifyEqual(calibration.corrections.DJF, 3, 'AbsTol', 1e-12);
   testCase.verifyEqual(calibration.corrections.MAM, 3, 'AbsTol', 1e-12);
end

function test_proxy_identity_without_overlap(testCase)
   % No usable overlap yields the recorded identity: the model passes
   % through unchanged and missing model samples stay missing.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   model = series.lwd;
   model(5) = NaN;
   obs = nan(numel(model), 1);

   calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
      times, obs, model, "lwd", fit_years=2020);

   testCase.verifyTrue(calibration.identity);
   testCase.verifyEqual(calibration.corrections.JJA, 0);

   returned = icemodel.forcing.reconstruct.applyProxyCalibration( ...
      calibration, times, model);
   testCase.verifyEqual(returned, model);
end

function test_proxy_swd_binned_calibration_round_trip(testCase)
   % D-28: swd fits one ratio per solar-elevation band per season, so
   % regime-dependent biases (twilight LOW, shoulder HIGH) are corrected
   % where a single seasonal ratio mixes them. 15-minute cadence at 45 N
   % keeps every band populated in every season.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') ...
      + minutes(0:15:(60 * 24 * 366 - 15))).';
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   edges = bands.calibration_bin_edges_deg;
   elevation = icemodel.forcing.helpers.solarElevation(times, 45, 0);
   toa = icemodel.forcing.reconstruct.toaIrradiance(times, 45, 0);
   band = discretize(elevation, edges);
   % A strictly positive proxy gives every band a usable denominator;
   % production twilight RCM output is near-zero, but the fit only needs
   % model > 0.
   model = max(2, 0.8 * toa);
   % One known ratio per band: strong twilight boost, shoulder inflation
   % decaying toward a mild high-sun deflation.
   ratios_true = [0.5, 6, 2, 1.5, 1.1, 0.9];
   obs = model .* ratios_true(band).';
   % Deep darkness never has observed overlap here (mirrors production,
   % where the darkness pre-pass owns those samples), forcing the
   % insufficient-band fallback path in band 1.
   obs(band == 1) = NaN;

   calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
      times, obs, model, "swd", fit_years=2020, target_toa=toa, ...
      target_elevation=elevation);

   % Version-2 schema: edges plus per-season binned ratios alongside the
   % legacy seasonal scalar that elevation-less consumers keep reading.
   testCase.verifyEqual(calibration.version, 2);
   testCase.verifyEqual(calibration.bin_edges_deg, edges);
   testCase.verifyTrue(isfield(calibration, 'corrections'));
   season = icemodel.forcing.reconstruct.seasonOf(times);
   recovered = false(numel(times), 1);
   for name = ["DJF", "MAM", "JJA", "SON"]
      counts = calibration.binned_counts.(char(name));
      returned = calibration.binned_corrections.(char(name));
      enough = counts >= bands.min_bin_samples;
      % Bands with enough overlap recover their exact regime ratio.
      testCase.verifyTrue(any(enough));
      testCase.verifyEqual(returned(enough), ratios_true(enough), ...
         'AbsTol', 1e-9);
      % Erased band 1 (and any season's naturally thin band) inherits the
      % season's single ratio instead of a noisy median.
      testCase.verifyEqual(counts(1), 0);
      testCase.verifyEqual(returned(~enough), repmat( ...
         calibration.corrections.(char(name)), 1, nnz(~enough)));
      recovered(season == name & ismember(band, find(enough))) = true;
   end

   % Applying with elevation reproduces the banded observations exactly
   % wherever the sample's season x band ratio was actually fit.
   applied = icemodel.forcing.reconstruct.applyProxyCalibration( ...
      calibration, times, model, target_elevation=elevation);
   testCase.verifyEqual(applied(recovered), obs(recovered), ...
      'AbsTol', 1e-9);

   % Elevation-less application falls back to the same record's
   % per-season scalar, so v2 records stay readable by callers without
   % station geometry.
   fallback = icemodel.forcing.reconstruct.applyProxyCalibration( ...
      calibration, times, model);
   for name = ["DJF", "MAM", "JJA", "SON"]
      in_season = season == name;
      testCase.verifyEqual(fallback(in_season), model(in_season) ...
         * calibration.corrections.(char(name)), 'AbsTol', 1e-9);
   end
end

function test_proxy_apply_legacy_single_ratio_record(testCase)
   % D-28 backward compatibility: an old persisted single-ratio swd
   % record — no version, edges, or binned fields — still applies its
   % per-season scalar exactly, with or without a supplied elevation.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   model = max(0, 600 * sin(2 * pi * (hour(times) - 6) / 24));
   legacy = struct('channel', "swd", 'mode', "multiplicative", ...
      'fit_years', 2020, 'corrections', struct('DJF', 0.8, ...
      'MAM', 0.9, 'JJA', 1.1, 'SON', 1.2), ...
      'n_overlap', 500, 'identity', false);
   elevation = icemodel.forcing.helpers.solarElevation( ...
      times, 67.0, -48.8);

   season = icemodel.forcing.reconstruct.seasonOf(times);
   expected = nan(numel(model), 1);
   for name = ["DJF", "MAM", "JJA", "SON"]
      in_season = season == name;
      expected(in_season) = model(in_season) ...
         * legacy.corrections.(char(name));
   end

   returned = icemodel.forcing.reconstruct.applyProxyCalibration( ...
      legacy, times, model);
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   % The field-presence guard means a supplied elevation cannot invent
   % bands a legacy record never fit.
   returned = icemodel.forcing.reconstruct.applyProxyCalibration( ...
      legacy, times, model, target_elevation=elevation);
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
end

%% lwdEstimator

function test_lwd_estimator_plausible_and_propagates_nan(testCase)
   % Sane polar inputs produce finite lwd inside the physical registry
   % bounds, and missing inputs propagate to NaN estimates.
   tair = linspace(250, 285, 50).';
   rh = linspace(50, 100, 50).';
   tair(3) = NaN;
   rh(7) = NaN;

   returned = icemodel.forcing.reconstruct.lwdEstimator(tair, rh);

   testCase.verifyTrue(all(isnan(returned([3, 7]))));
   good = setdiff(1:50, [3, 7]);
   bounds = icemodel.forcing.reconstruct.physicalBounds("lwd");
   testCase.verifyTrue(all(isfinite(returned(good))));
   testCase.verifyTrue(all(returned(good) > bounds(1) & ...
      returned(good) < bounds(2)));
end

%% deriveUpwardShortwave

function test_derive_upward_shortwave_uses_final_operands(testCase)
   % Native swu survives; missing samples derive only where both final
   % operands exist and receive the dedicated algebraic provenance.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:3)).';
   swd = [100; 200; NaN; 0];
   albedo = [0.5; 0.6; 0.7; 0.8];
   swu = [40; NaN; NaN; NaN];
   filled = timetable(times, swd, albedo, swu);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   swu_code = [codes.observed; codes.missing; codes.missing; codes.missing];
   provenance = timetable(times, swu_code, 'VariableNames', {'swu'});
   audit = table('Size', [0 7], 'VariableTypes', {'cellstr', ...
      'datetime', 'datetime', 'double', 'cellstr', 'cellstr', 'cellstr'}, ...
      'VariableNames', {'channel', 'start_time', 'end_time', ...
      'duration_hours', 'method', 'detail', 'context_id'});
   audit.start_time.TimeZone = 'UTC';
   audit.end_time.TimeZone = 'UTC';

   [filled, provenance, audit] = ...
      icemodel.forcing.reconstruct.deriveUpwardShortwave( ...
      filled, provenance, audit, codes);

   testCase.verifyEqual(filled.swu, [40; 120; NaN; 0]);
   testCase.verifyEqual(provenance.swu, [codes.observed; ...
      codes.derived_shortwave; codes.missing; codes.derived_shortwave]);
   testCase.verifyEqual(string(audit.method), ...
      repmat("derived_shortwave", 2, 1));
end

%% provenanceCodes

function test_provenance_registry_values(testCase)
   % The registry is the append-only published mapping; any renumbering
   % would silently re-label staged products, so pin every value.
   returned = icemodel.forcing.reconstruct.provenanceCodes();
   expected = struct( ...
      'observed', uint8(0), ...
      'bounded_interp', uint8(1), ...
      'station_transfer', uint8(2), ...
      'spline_adjustment', uint8(3), ...
      'mar', uint8(4), ...
      'merra2', uint8(5), ...
      'racmo', uint8(6), ...
      'modis', uint8(7), ...
      'climatology', uint8(8), ...
      'constant', uint8(9), ...
      'empirical_estimator', uint8(10), ...
       'darkness', uint8(11), ...
       'derived_shortwave', uint8(12), ...
      'raw_shortwave', uint8(13), ...
      'clamped_shortwave', uint8(14), ...
      'twilight_climatology', uint8(15), ...
       'missing', uint8(255));
   testCase.verifyEqual(returned, expected);
end

function test_reconstruct_preserves_native_source_provenance(testCase)
   % Explicit source-backed native codes override only the default observed
   % label, and malformed vectors fail at the engine boundary.
   series = icemodel.test.fixtures.makeReconstructSeries();
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   cm = struct('channel', "tair", 'methods', struct.empty);
   source = struct('tair', repmat(codes.raw_shortwave, height(series), 1));

   result = icemodel.forcing.reconstruct.reconstructSeries( ...
      series, cm, native_provenance=source);

   testCase.verifyTrue(all(result.provenance.tair == codes.raw_shortwave));
   source.tair = codes.raw_shortwave;
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.reconstructSeries( ...
      series, cm, native_provenance=source), ...
      'icemodel:reconstruct:reconstructSeries:invalidNativeProvenance');
end

%% reconstructSeries

function test_reconstruct_composes_tiers_with_provenance(testCase)
   % Three deliberate summer gaps compose through the tiers: a 3 h gap
   % fills by tier-1 interpolation, a 12 h gap by the first admitted
   % method (donor), and a 48 h gap by climatology after the donor
   % declines; stratum-restricted decoys never fire outside their
   % season/bucket despite offering finite estimates.
   series = icemodel.test.fixtures.makeReconstructSeries();
   truth = series.tair;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   n = numel(truth);
   gap_a = (3865:3867).';                        % 3 h  -> tier 1 (Jun 10)
   gap_b = (4000:4011).';                        % 12 h -> donor method
   gap_c = (4200:4247).';                        % 48 h -> climatology
   series.tair([gap_a; gap_b; gap_c]) = NaN;

   % Donor declines everywhere except the 12 h gap; climatology offers
   % the truth everywhere; the decoys are blocked by bucket and season.
   est_donor = nan(n, 1);
   est_donor(gap_b) = truth(gap_b);
   m_bucket = struct('name', "bucket1_decoy", 'code', codes.constant, ...
      'estimate', truth, 'seasons', "all", 'buckets', 1, ...
      'audit_context_id', "plan-bucket");
   m_season = struct('name', "winter_decoy", 'code', codes.mar, ...
      'estimate', truth, 'seasons', "DJF", 'buckets', [], ...
      'audit_context_id', "plan-season");
   m_donor = struct('name', "donor:test", ...
      'code', codes.station_transfer, 'estimate', est_donor, ...
      'seasons', "all", 'buckets', [], ...
      'audit_context_id', "plan-donor");
   m_clim = struct('name', "clim:test", 'code', codes.climatology, ...
      'estimate', truth, 'seasons', "all", 'buckets', [], ...
      'audit_context_id', "plan-climatology");
   cm = struct('channel', "tair", ...
      'methods', [m_bucket, m_season, m_donor, m_clim]);

   result = icemodel.forcing.reconstruct.reconstructSeries(series, cm);

   % Per-sample provenance is exact: native, tier-1, donor, climatology.
   expected = repmat(codes.observed, n, 1);
   expected(gap_a) = codes.bounded_interp;
   expected(gap_b) = codes.station_transfer;
   expected(gap_c) = codes.climatology;
   testCase.verifyEqual(result.provenance.tair, expected);
   % Native samples are untouched and every gap sample is now finite.
   native = expected == codes.observed;
   testCase.verifyEqual(result.series.tair(native), truth(native));
   testCase.verifyTrue(all(isfinite(result.series.tair)));
   % Method-filled samples carry the exact admitted estimates.
   testCase.verifyEqual(result.series.tair(gap_b), truth(gap_b));
   testCase.verifyEqual(result.series.tair(gap_c), truth(gap_c));
   % The audit records one row per fill with a plan-context join key.
   testCase.verifyEqual(result.audit.Properties.VariableNames, ...
      {'channel', 'start_time', 'end_time', 'duration_hours', ...
      'method', 'detail', 'context_id'});
   testCase.verifyEqual(string(result.audit.method), ...
      ["bounded_interp"; "donor:test"; "clim:test"]);
   testCase.verifyEqual(string(result.audit.context_id), ...
      ["bounded_interp"; "plan-donor"; "plan-climatology"]);
   % The registry ships with the product for later decoding.
   testCase.verifyEqual(result.registry, codes);
end

function test_reconstruct_blends_offset_method_seams(testCase)
   % Policy rule 7: a finite, in-bounds estimate with a large anchor
   % offset fills with seam blending — the excess offsets taper across
   % the blend window and the method keeps the fill.
   series = icemodel.test.fixtures.makeReconstructSeries();
   truth = series.tair;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   n = numel(truth);
   gap = (4000:4011).';                          % 12 h summer run
   series.tair(gap) = NaN;

   est_offset = nan(n, 1);
   est_offset(gap) = truth(gap) + 25;            % in-bounds, big offset
   m_offset = struct('name', "offset", 'code', codes.station_transfer, ...
      'estimate', est_offset, 'seasons', "all", 'buckets', []);
   m_clim = struct('name', "clim:test", 'code', codes.climatology, ...
      'estimate', truth, 'seasons', "all", 'buckets', []);
   cm = struct('channel', "tair", 'methods', [m_offset, m_clim]);

   result = icemodel.forcing.reconstruct.reconstructSeries(series, cm);

   % The offset method fills the run (no demotion to climatology).
   testCase.verifyEqual(result.provenance.tair(gap), ...
      repmat(codes.station_transfer, numel(gap), 1));
   % Both anchored seams land inside the jump limit (the blend corrects
   % the excess mismatch, targeting half the limit).
   scale = icemodel.forcing.reconstruct.stepScale( ...
      series.Properties.RowTimes, series.tair);
   limit = 3 * scale.JJA;
   returned = result.series.tair;
   testCase.verifyLessThanOrEqual( ...
      abs(returned(gap(1)) - truth(gap(1) - 1)), limit);
   testCase.verifyLessThanOrEqual( ...
      abs(returned(gap(end)) - truth(gap(end) + 1)), limit);
   % Each independent long-run taper reaches zero before the untouched
   % interior, so no correction jump is introduced inside the gap.
   testCase.verifyEqual(returned(gap(6:7)), est_offset(gap(6:7)), ...
      'AbsTol', 1e-12);
   % The audit names the fill and records the seam blend offsets.
   row = result.audit(strcmp(result.audit.method, 'offset'), :);
   testCase.verifyEqual(height(row), 1);
   testCase.verifySubstring(row.detail{1}, 'seam blend');
end

function test_blend_seams_ignores_segment_without_usable_samples(testCase)
   % A method usable elsewhere can leave one residual segment wholly
   % unusable; seam handling must leave it unchanged for the next method.
   x = [1; NaN; NaN; 4];
   candidate = [10; 11];
   [returned, note] = icemodel.forcing.reconstruct.blendSeams( ...
      x, candidate, [2; 3], [false; false], 2, 3, 1, 4, 2, 1);

   testCase.verifyEqual(returned, candidate);
   testCase.verifyEmpty(note);
end

function test_reconstruct_blend_preserves_interior_signal(testCase)
   % Beyond the taper window the method's interior signal survives
   % untouched: a 48 h run with a constant +25 offset keeps that offset
   % in the middle while both 6 h seams decay to the anchors.
   series = icemodel.test.fixtures.makeReconstructSeries();
   truth = series.tair;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   n = numel(truth);
   gap = (4200:4247).';                          % 48 h run
   series.tair(gap) = NaN;

   est_offset = nan(n, 1);
   est_offset(gap) = truth(gap) + 25;
   m_offset = struct('name', "offset", 'code', codes.station_transfer, ...
      'estimate', est_offset, 'seasons', "all", 'buckets', []);
   cm = struct('channel', "tair", 'methods', m_offset);

   result = icemodel.forcing.reconstruct.reconstructSeries(series, cm);

   % Interior samples (outside both 6-sample tapers) carry the method
   % value exactly; both seams land inside the jump limit.
   scale = icemodel.forcing.reconstruct.stepScale( ...
      series.Properties.RowTimes, series.tair);
   limit = 3 * scale.JJA;
   returned = result.series.tair;
   interior = gap(8:41);
   testCase.verifyEqual(returned(interior), truth(interior) + 25, ...
      'AbsTol', 1e-12);
   testCase.verifyLessThanOrEqual( ...
      abs(returned(gap(1)) - truth(gap(1) - 1)), limit);
   testCase.verifyLessThanOrEqual( ...
      abs(returned(gap(end)) - truth(gap(end) + 1)), limit);
end

function test_reconstruct_blends_run_shorter_than_taper(testCase)
   % A run shorter than two taper windows takes one straight correction
   % line between the two seam offsets, so a 3 h two-sided offset gap
   % fills with both seams inside the jump limit. Tier 1 is disabled so
   % the short gap reaches the method walk (mirroring tier-1 declines
   % and channels with no tier-1 path).
   series = icemodel.test.fixtures.makeReconstructSeries();
   truth = series.tair;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   n = numel(truth);
   gap = (4000:4002).';                          % 3 h run < 6 h taper
   series.tair(gap) = NaN;

   est_offset = nan(n, 1);
   est_offset(gap) = truth(gap) + 25;            % offsets both seams
   m_offset = struct('name', "offset", 'code', codes.station_transfer, ...
      'estimate', est_offset, 'seasons', "all", 'buckets', []);
   cm = struct('channel', "tair", 'methods', m_offset);

   result = icemodel.forcing.reconstruct.reconstructSeries(series, cm, ...
      interp_channels=string.empty(1, 0));

   % The offset method fills the whole short run.
   testCase.verifyEqual(result.provenance.tair(gap), ...
      repmat(codes.station_transfer, numel(gap), 1));
   % Both seams land inside the jump limit — the overlap-free correction
   % applies each seam offset exactly once at its own edge.
   scale = icemodel.forcing.reconstruct.stepScale( ...
      series.Properties.RowTimes, series.tair);
   limit = 3 * scale.JJA;
   returned = result.series.tair;
   testCase.verifyLessThanOrEqual( ...
      abs(returned(gap(1)) - truth(gap(1) - 1)), limit);
   testCase.verifyLessThanOrEqual( ...
      abs(returned(gap(end)) - truth(gap(end) + 1)), limit);
end

function test_reconstruct_bucket1_requires_own_admission(testCase)
   % A sub-6 h run reaching the method walk uses only methods that beat
   % persistence in its own held-out bucket; bucket-2 skill cannot rescue it.
   series = icemodel.test.fixtures.makeReconstructSeries();
   truth = series.tair;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   n = numel(truth);
   gap = (4000:4002).';                          % 3 h run -> bucket 1
   series.tair(gap) = NaN;

   est = nan(n, 1);
   est(gap) = truth(gap);
   m_b2 = struct('name', "donor:b2only", ...
      'code', codes.constant, 'estimate', est, ...
      'seasons', "all", 'buckets', 2);
   m_b1 = struct('name', "donor:b1", ...
      'code', codes.station_transfer, 'estimate', est, ...
      'seasons', "all", 'buckets', 1);
   cm = struct('channel', "tair", 'methods', [m_b2, m_b1]);

   result = icemodel.forcing.reconstruct.reconstructSeries(series, cm, ...
      interp_channels=string.empty(1, 0));

   testCase.verifyEqual(result.provenance.tair(gap), ...
      repmat(codes.station_transfer, numel(gap), 1));
   testCase.verifyEqual(result.series.tair(gap), truth(gap));
end

function test_reconstruct_splits_gap_at_season_boundary(testCase)
   % A short cross-season outage bypasses tier-1 bridging and each segment
   % uses only the method admitted for its own held-out season.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   seasons = icemodel.forcing.reconstruct.seasonOf(times);
   boundary = find(seasons(1:end - 1) == "DJF" ...
      & seasons(2:end) == "MAM", 1) + 1;
   testCase.assertNotEmpty(boundary);
   djf_gap = (boundary - 2:boundary - 1).';
   mam_gap = (boundary:boundary + 1).';
   truth = series.tair;
   series.tair([djf_gap; mam_gap]) = NaN;

   % Each estimate is finite only on its validated seasonal segment.
   djf_estimate = nan(height(series), 1);
   djf_estimate(djf_gap) = truth(djf_gap);
   mam_estimate = nan(height(series), 1);
   mam_estimate(mam_gap) = truth(mam_gap);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   djf_method = struct('name', "djf", 'code', codes.station_transfer, ...
      'estimate', djf_estimate, 'seasons', "DJF", 'buckets', 1);
   mam_method = struct('name', "mam", 'code', codes.climatology, ...
      'estimate', mam_estimate, 'seasons', "MAM", 'buckets', 1);
   request = struct('channel', "tair", ...
      'methods', [djf_method, mam_method]);

   result = icemodel.forcing.reconstruct.reconstructSeries(series, request);

   testCase.verifyEqual(result.provenance.tair(djf_gap), ...
      repmat(codes.station_transfer, numel(djf_gap), 1));
   testCase.verifyEqual(result.provenance.tair(mam_gap), ...
      repmat(codes.climatology, numel(mam_gap), 1));
   testCase.verifyEqual(result.series.tair([djf_gap; mam_gap]), ...
      truth([djf_gap; mam_gap]));
end

function test_reconstruct_later_seams_use_prefill_native_scale(testCase)
   % Tier-1 values must not tighten or inflate the observed step scale used
   % to constrain a later donor/proxy gap.
   series = icemodel.test.fixtures.makeReconstructSeries();
   series.tair = 250 + 10 * mod((1:height(series)).', 2);
   series.tair(99) = 240;
   series.tair(100:105) = NaN;                  % tier-1 6 h interpolation
   series.tair(106) = 280;
   later_gap = (200:211).';                     % bucket-2 method fill
   series.tair(later_gap) = NaN;
   estimate = nan(height(series), 1);
   estimate(later_gap) = 280;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   method = struct('name', "donor:native-scale", ...
      'code', codes.station_transfer, 'estimate', estimate, ...
      'seasons', "all", 'buckets', 2);
   channel_methods = struct('channel', "tair", 'methods', method);

   result = icemodel.forcing.reconstruct.reconstructSeries( ...
      series, channel_methods, interp_channels="tair");

   testCase.verifyEqual(result.series.tair(later_gap), ...
      repmat(280, numel(later_gap), 1), 'AbsTol', 1e-12);
   testCase.verifyEqual(result.provenance.tair(later_gap), ...
      repmat(codes.station_transfer, numel(later_gap), 1));
end

function test_reconstruct_tier1_uses_prefill_native_scale(testCase)
   % Darkness zeros inserted before Tier 1 must not inflate the native seam
   % scale used by the shortwave CSI bridge.
   full_times = (datetime(2020, 3, 21, 0, 0, 0, 'TimeZone', 'UTC'): ...
      hours(1):datetime(2020, 3, 21, 23, 0, 0, ...
      'TimeZone', 'UTC')).';
   full_toa = icemodel.forcing.reconstruct.toaIrradiance( ...
      full_times, 67.0, -48.8);
   full_dark = full_toa < 10;
   dawn = find(full_dark(1:end - 1) & ~full_dark(2:end), 1, 'first');
   keep = dawn:(dawn + 5);
   times = full_times(keep);
   toa = full_toa(keep);
   dark = toa < 10;
   anchor = 4;
   x = nan(size(toa));
   x(2) = 200;
   x([anchor, anchor + 2]) = 100;
   series = timetable(times, x, 'VariableNames', {'swd'});
   methods = struct('name', {}, 'code', {}, 'estimate', {}, ...
      'seasons', {}, 'buckets', {});
   channel_methods = struct('channel', "swd", 'methods', methods);

   result = icemodel.forcing.reconstruct.reconstructSeries( ...
      series, channel_methods, latitude=67.0, longitude=-48.8, ...
      interp_channels="swd");

   after_dark = x;
   after_dark(dark & ~isfinite(after_dark)) = 0;
   contaminated = icemodel.forcing.reconstruct.stepScale(times, after_dark);
   [contaminated_fill, filled] = ...
      icemodel.forcing.reconstruct.fillShortGaps( ...
      times, after_dark, "swd", latitude=67.0, longitude=-48.8, ...
      step_scale=contaminated);

   testCase.verifyTrue(filled(anchor + 1));
   testCase.verifyGreaterThan( ...
      abs(contaminated_fill(anchor + 1) - 100), 1e-9);
   testCase.verifyEqual(result.series.swd(anchor + 1), 100, ...
      'AbsTol', 1e-12);
end

function test_reconstruct_bridge_floor_rescues_single_sample_gap(testCase)
   % When the native record itself steps across a gap by more than the
   % seasonal seam limit, no fill value can sit inside the limit of both
   % anchors — the per-run limit is floored at the linear bridge's
   % implied per-step change so a bridge-smooth fill is never refused.
   series = icemodel.test.fixtures.makeReconstructSeries();
   truth = series.tair;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   n = numel(truth);
   gap = 4000;                                   % single-sample run
   series.tair(gap) = NaN;
   % Force an anchor step far beyond 3x the median step: a synthetic
   % frontal jump directly after the gap.
   series.tair(gap + 1) = truth(gap + 1) + 30;

   est = nan(n, 1);
   est(gap) = truth(gap);
   m_est = struct('name', "donor:test", ...
      'code', codes.station_transfer, 'estimate', est, ...
      'seasons', "all", 'buckets', []);
   cm = struct('channel', "tair", 'methods', m_est);

   result = icemodel.forcing.reconstruct.reconstructSeries(series, cm, ...
      interp_channels=string.empty(1, 0));

   % The single sample fills instead of being refused as jointly
   % infeasible; the value lands between the two anchors.
   testCase.verifyEqual(result.provenance.tair(gap), ...
      codes.station_transfer);
   returned = result.series.tair(gap);
   testCase.verifyTrue(isfinite(returned));
   lo = min(series.tair(gap - 1), series.tair(gap + 1));
   hi = max(series.tair(gap - 1), series.tair(gap + 1));
   testCase.verifyGreaterThanOrEqual(returned, lo - 1e-9);
   testCase.verifyLessThanOrEqual(returned, hi + 1e-9);
end

function test_reconstruct_audits_residual_missing_runs(testCase)
   % Runs no admitted method fills stay missing AND leave an audit row
   % naming every refusal reason.
   series = icemodel.test.fixtures.makeReconstructSeries();
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   n = height(series);
   gap = (4200:4247).';                          % 48 h run -> bucket 3
   series.tair(gap) = NaN;

   % Admitted only for bucket 1, so the run has no eligible method; the
   % second method offers no finite samples at all.
   m_bucket = struct('name', "bucket1_only", 'code', codes.constant, ...
      'estimate', series.tair, 'seasons', "all", 'buckets', 1);
   m_empty = struct('name', "empty", 'code', codes.climatology, ...
      'estimate', nan(n, 1), 'seasons', "all", 'buckets', []);
   cm = struct('channel', "tair", 'methods', [m_bucket, m_empty]);

   result = icemodel.forcing.reconstruct.reconstructSeries(series, cm);

   % The run is still missing and provenance says so.
   testCase.verifyTrue(all(~isfinite(result.series.tair(gap))));
   testCase.verifyEqual(result.provenance.tair(gap), ...
      repmat(codes.missing, numel(gap), 1));
   % One 'unfilled' audit row explains both refusals.
   row = result.audit(strcmp(result.audit.method, 'unfilled'), :);
   testCase.verifyEqual(height(row), 1);
   testCase.verifySubstring(row.detail{1}, 'not admitted');
   testCase.verifySubstring(row.detail{1}, 'no finite estimate');
end

function test_reconstruct_partial_fill_advances_methods(testCase)
   % A method that declines the tail of a run fills its half; the
   % leftover samples advance to the next method in the same walk.
   series = icemodel.test.fixtures.makeReconstructSeries();
   truth = series.tair;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   n = numel(truth);
   gap = (4400:4411).';                          % 12 h run (Jul 2)
   series.tair(gap) = NaN;

   est_half = nan(n, 1);
   est_half(gap(1:6)) = truth(gap(1:6));         % declines the tail
   m_half = struct('name', "donor:half", ...
      'code', codes.station_transfer, 'estimate', est_half, ...
      'seasons', "all", 'buckets', []);
   m_clim = struct('name', "clim:test", 'code', codes.climatology, ...
      'estimate', truth, 'seasons', "all", 'buckets', []);
   cm = struct('channel', "tair", 'methods', [m_half, m_clim]);

   result = icemodel.forcing.reconstruct.reconstructSeries(series, cm);

   testCase.verifyEqual(result.provenance.tair(gap(1:6)), ...
      repmat(codes.station_transfer, 6, 1));
   testCase.verifyEqual(result.provenance.tair(gap(7:12)), ...
      repmat(codes.climatology, 6, 1));
   testCase.verifyEqual(result.series.tair(gap), truth(gap));
end

function test_reconstruct_blends_disjoint_leftovers_separately(testCase)
   % A middle fill can leave two edge fragments; the next method must taper
   % each fragment in wall-clock order without crossing the filled middle.
   series = icemodel.test.fixtures.makeReconstructSeries();
   truth = series.tair;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   n = height(series);
   gap = (4400:4413).';
   series.tair(gap) = NaN;

   middle = nan(n, 1);
   middle(gap(5:10)) = truth(gap(5:10));
   edges = nan(n, 1);
   edges(gap([1:4, 11:14])) = truth(gap([1:4, 11:14])) + 25;
   m_middle = struct('name', "middle", 'code', codes.station_transfer, ...
      'estimate', middle, 'seasons', "all", 'buckets', []);
   m_edges = struct('name', "edges", 'code', codes.climatology, ...
      'estimate', edges, 'seasons', "all", 'buckets', []);
   request = struct('channel', "tair", ...
      'methods', [m_middle, m_edges]);

   result = icemodel.forcing.reconstruct.reconstructSeries( ...
      series, request, interp_channels=string.empty(1, 0));

   testCase.verifyEqual(result.provenance.tair(gap(5:10)), ...
      repmat(codes.station_transfer, 6, 1));
   testCase.verifyEqual(result.provenance.tair(gap([1:4, 11:14])), ...
      repmat(codes.climatology, 8, 1));
   testCase.verifyEqual(result.series.tair(gap([4, 11])), ...
      truth(gap([4, 11])) + 25, 'AbsTol', 1e-12);
end

function test_toa_irradiance_bounds_and_night(testCase)
   % The shared clear-sky reference is zero at night, positive at midday,
   % and never exceeds the solar constant.
   times = datetime(2020, 3, 21, 0, 0, 0, 'TimeZone', 'UTC') ...
      + hours(0:23).';
   returned = icemodel.forcing.reconstruct.toaIrradiance(times, ...
      67.0, -48.8);
   % Local solar midnight (~03:15 UTC at 48.8 W) is dark at the equinox;
   % local solar noon (~15:15 UTC) carries meaningful sun.
   testCase.verifyEqual(returned(4), 0);
   testCase.verifyGreaterThan(returned(16), 100);
   testCase.verifyLessThanOrEqual(max(returned), 1361);
end

function test_reconstruct_darkness_zero_fills_swd_nights(testCase)
   % Policy swd darkness rule (B2/D-28): missing samples with the sun
   % below civil twilight are KNOWN zeros. A multi-day March swd outage
   % with no admitted methods gets its deep-dark samples zero-filled with
   % the darkness provenance code, while twilight-band and daylight
   % samples stay honestly missing — decomposed into per-day fragments
   % for the bucketed methods (none here).
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   toa = icemodel.forcing.reconstruct.toaIrradiance(times, 67.0, -48.8);
   % Give the series a shortwave channel with real nights (March window).
   series.swd = 0.7 * toa;
   gap = find(times >= datetime(2020, 3, 10, 0, 0, 0, 'TimeZone', ...
      'UTC') & times < datetime(2020, 3, 13, 0, 0, 0, 'TimeZone', 'UTC'));
   series.swd(gap) = NaN;

   % No admitted methods and tier 1 disabled: only the darkness pre-pass
   % may act, isolating its behavior.
   methods = struct('name', {}, 'code', {}, 'estimate', {}, ...
      'seasons', {}, 'buckets', {});
   cm = struct('channel', "swd", 'methods', methods);
   result = icemodel.forcing.reconstruct.reconstructSeries(series, cm, ...
      latitude=67.0, longitude=-48.8, ...
      interp_channels=string.empty(1, 0));

   % The darkness boundary is the civil-twilight solar elevation, never
   % the TOA threshold: the twilight band carries real diffuse light the
   % pre-pass must not erase (D-28).
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   elevation = icemodel.forcing.helpers.solarElevation(times, 67.0, -48.8);
   maximum_elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, 67.0, -48.8, hours(1));
   dark = gap(maximum_elevation(gap) <= bands.civil_twilight_deg);
   twilight = gap(maximum_elevation(gap) > bands.civil_twilight_deg ...
      & elevation(gap) < 0);
   light = gap(maximum_elevation(gap) > bands.civil_twilight_deg);
   testCase.verifyNotEmpty(dark);
   testCase.verifyNotEmpty(twilight);
   testCase.verifyNotEmpty(light);
   returned = result.series.swd;
   % Deep-dark samples are exact zeros with the darkness code.
   testCase.verifyEqual(returned(dark), zeros(numel(dark), 1));
   testCase.verifyEqual(result.provenance.swd(dark), ...
      repmat(codes.darkness, numel(dark), 1));
   % Twilight and daylight samples stay missing — twilight is left to
   % the fill tiers (none admitted here), never hard-zeroed.
   testCase.verifyTrue(all(~isfinite(returned(light))));
   testCase.verifyEqual(result.provenance.swd(light), ...
      repmat(codes.missing, numel(light), 1));
    % The audit carries one row per contiguous night, never a false span
    % across intervening daylight.
    row = result.audit(strcmp(result.audit.method, 'darkness_zero'), :);
    testCase.verifyGreaterThan(height(row), 1);
    testCase.verifyEqual(sum(row.duration_hours), ...
       numel(dark) * hours(median(diff(times))));
end

function test_reconstruct_darkness_leaves_twilight_to_methods(testCase)
   % D-28 twilight handling: samples between civil twilight and sunrise
   % are NOT known zeros — stations measure real diffuse light there —
   % so an admitted method (not the darkness pre-pass) fills them, while
   % deep darkness still zero-fills first.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   toa = icemodel.forcing.reconstruct.toaIrradiance(times, 67.0, -48.8);
   series.swd = 0.7 * toa;
   % One full March night including both evening and morning twilight.
   gap = find(times >= datetime(2020, 3, 11, 12, 0, 0, 'TimeZone', ...
      'UTC') & times < datetime(2020, 3, 12, 12, 0, 0, 'TimeZone', 'UTC'));
   series.swd(gap) = NaN;

   % The method offers a small positive twilight-compatible value
   % everywhere (kept inside the darkness-era swd ceiling), so any
   % twilight sample it fills is visibly non-zero with its provenance.
   est = nan(height(series), 1);
   est(gap) = max(4, 0.7 * toa(gap));
   method = struct('name', "proxy:test", 'code', codes.mar, ...
      'estimate', est, 'seasons', "all", 'buckets', []);
   cm = struct('channel', "swd", 'methods', method);
   result = icemodel.forcing.reconstruct.reconstructSeries(series, cm, ...
      latitude=67.0, longitude=-48.8, ...
      interp_channels=string.empty(1, 0));

   bands = icemodel.forcing.reconstruct.solarElevationBands();
   elevation = icemodel.forcing.helpers.solarElevation(times, 67.0, -48.8);
   maximum_elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, 67.0, -48.8, hours(1));
   dark = gap(maximum_elevation(gap) <= bands.civil_twilight_deg);
   twilight = gap(maximum_elevation(gap) > bands.civil_twilight_deg ...
      & elevation(gap) < 0);
   testCase.verifyNotEmpty(dark);
   testCase.verifyNotEmpty(twilight);
   returned = result.series.swd;
   % Deep darkness is still the pre-pass's known zero.
   testCase.verifyEqual(returned(dark), zeros(numel(dark), 1));
   testCase.verifyEqual(result.provenance.swd(dark), ...
      repmat(codes.darkness, numel(dark), 1));
   % Twilight samples carry the METHOD's finite fill and provenance —
   % values may be seam-blended, so provenance is the honest witness.
   testCase.verifyTrue(all(isfinite(returned(twilight))));
   testCase.verifyEqual(result.provenance.swd(twilight), ...
      repmat(codes.mar, numel(twilight), 1));
end

function test_reconstruct_wholly_missing_shortwave_keeps_known_darkness(testCase)
   % B2 depends on solar geometry, not native sensor support: deep darkness
   % is known zero even when the SWD channel is otherwise wholly missing.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   series.swd = nan(height(series), 1);
   methods = struct('name', {}, 'code', {}, 'estimate', {}, ...
      'seasons', {}, 'buckets', {});
   request = struct('channel', "swd", 'methods', methods);

   result = icemodel.forcing.reconstruct.reconstructSeries( ...
      series, request, latitude=67.0, longitude=-48.8, ...
      interp_channels=string.empty(1, 0));

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   bands = icemodel.forcing.reconstruct.solarElevationBands();
   maximum_elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, 67.0, -48.8, hours(1));
   dark = maximum_elevation <= bands.civil_twilight_deg;
   testCase.verifyTrue(any(dark));
   testCase.verifyEqual(result.series.swd(dark), zeros(nnz(dark), 1));
   testCase.verifyEqual(result.provenance.swd(dark), ...
      repmat(codes.darkness, nnz(dark), 1));
   testCase.verifyTrue(all(isnan(result.series.swd(~dark))));
   testCase.verifyTrue(all(result.provenance.swd(~dark) == codes.missing));
   testCase.verifyTrue(any(string(result.audit.method) == "darkness_zero"));
end

function test_reconstruct_rejects_unknown_channel(testCase)
   % Channels missing from the series fail loudly instead of silently
   % composing nothing.
   series = icemodel.test.fixtures.makeReconstructSeries();
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   methods = struct('name', "m", 'code', codes.climatology, ...
      'estimate', nan(height(series), 1), 'seasons', "all", ...
      'buckets', []);
   cm = struct('channel', "nope", 'methods', methods);

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.reconstructSeries(series, cm), ...
      'icemodel:reconstruct:reconstructSeries:unknownChannel');
end

function test_reconstruct_refuses_duration_beyond_validation(testCase)
   % An admitted bucket does not authorize a longer real gap than the
   % longest held-out block actually tested; the explicit factor can widen
   % that horizon when policy changes.
   series = icemodel.test.fixtures.makeReconstructSeries();
   gap = (100:107).';
   truth = series.tair;
   series.tair(gap) = NaN;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   method = struct('name', "validated_shorter", ...
      'code', codes.climatology, 'estimate', truth, 'seasons', "all", ...
      'buckets', 2, 'max_validated_hours', 6);
   request = struct('channel', "tair", 'methods', method);

   refused = icemodel.forcing.reconstruct.reconstructSeries( ...
      series, request, interp_channels=string.empty(1, 0));
   testCase.verifyTrue(all(isnan(refused.series.tair(gap))));
   testCase.verifyTrue(any(contains(string(refused.audit.detail), ...
      "held-out duration limit")));

   admitted = icemodel.forcing.reconstruct.reconstructSeries( ...
      series, request, interp_channels=string.empty(1, 0), ...
      max_validation_duration_factor=2);
   testCase.verifyEqual(admitted.series.tair(gap), truth(gap));
end

function test_reconstruct_empty_audit_shape(testCase)
   % A gap-free channel with no methods produces the documented empty
   % audit shape, all-observed provenance, and an untouched series.
   series = icemodel.test.fixtures.makeReconstructSeries();
   truth = series.tair;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   % An empty method list is a legal composition request: nothing to walk.
   methods = struct('name', {}, 'code', {}, 'estimate', {}, ...
      'seasons', {}, 'buckets', {});
   cm = struct('channel', "tair", 'methods', methods);

   result = icemodel.forcing.reconstruct.reconstructSeries(series, cm);

   testCase.verifyEqual(height(result.audit), 0);
   testCase.verifyEqual(result.audit.Properties.VariableNames, ...
      {'channel', 'start_time', 'end_time', 'duration_hours', ...
      'method', 'detail', 'context_id'});
   testCase.verifyEqual(result.provenance.tair, ...
      repmat(codes.observed, numel(truth), 1));
   testCase.verifyEqual(result.series.tair, truth);
end


function test_donor_fit_handles_constant_donor(testCase)
   % A near-constant donor must fall back to the linear fit instead of
   % crashing discretize with collapsed quantile edges.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:8783)).';
   x = 260 + 5 * sin(2 * pi * (1:8784).' / 8784);
   d = 250 * ones(8784, 1);

   transfer = icemodel.forcing.reconstruct.fitDonorTransfer(times, x, ...
      d, "tair", fit_years=2020, knots=6, min_overlap_hours=4000);

    testCase.verifyEqual(string(transfer.models.DJF.kind), "linear");
end

function test_physical_validity_enforces_relational_bounds(testCase)
    % Shared validity rejects shortwave above TOA and upward flux above swd.
    times = datetime(2020, 3, 21, [3; 15], 0, 0, 'TimeZone', 'UTC');
    toa = icemodel.forcing.reconstruct.toaIrradiance(times, 67.0, -48.8);
    values = [0; 1.04 * toa(2)];
    testCase.verifyTrue(all( ...
       icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", values, times, latitude=67.0, longitude=-48.8)));
    values(2) = 1.06 * toa(2);
    valid = icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", values, times, latitude=67.0, longitude=-48.8);
    testCase.verifyFalse(valid(2));
    testCase.verifyEqual( ...
       icemodel.forcing.reconstruct.physicalValidity( ...
       "swu", [50; 110], times, swd=[100; 100]), [true; false]);
    testCase.verifyEqual( ...
       icemodel.forcing.reconstruct.physicalValidity( ...
       "tair", [193; 301], times), [true; false]);
end

function test_physical_validity_rejects_missing_context(testCase)
    % Relational checks fail loudly without geometry or a paired swd axis.
    times = datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:1).';
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", [0; 0], times), ...
       'icemodel:reconstruct:physicalValidity:missingSolarGeometry');
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.physicalValidity( ...
       "swu", [0; 0], times), ...
       'icemodel:reconstruct:physicalValidity:missingShortwaveReference');
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.physicalValidity( ...
       "tair", 260, times), ...
       'icemodel:reconstruct:physicalValidity:sizeMismatch');
end

function test_physical_validity_accepts_precomputed_geometry(testCase)
    % Hot-loop callers pass axis-sliced precomputed solar geometry; the
    % verdict and both per-sample limits must match the internally
    % computed path exactly, and wrong-length vectors must be refused
    % rather than silently recomputed.
    times = datetime(2020, 3, 21, 'TimeZone', 'UTC') + hours(0:23).';
    values = 120 * ones(numel(times), 1);
    toa = icemodel.forcing.reconstruct.toaIrradiance(times, 67.0, -48.8);
    elevation = icemodel.forcing.helpers.solarElevation( ...
       times, 67.0, -48.8);
    [expected, expected_lower, expected_upper] = ...
       icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", values, times, latitude=67.0, longitude=-48.8);
    [returned, returned_lower, returned_upper] = ...
       icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", values, times, latitude=67.0, longitude=-48.8, ...
       toa=toa, elevation=elevation);
    testCase.verifyEqual(returned, expected);
    testCase.verifyEqual(returned_lower, expected_lower);
    testCase.verifyEqual(returned_upper, expected_upper);
    % Supplying only one of the two vectors computes the other internally.
    returned = icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", values, times, latitude=67.0, longitude=-48.8, toa=toa);
    testCase.verifyEqual(returned, expected);
    returned = icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", values, times, latitude=67.0, longitude=-48.8, ...
       elevation=elevation);
    testCase.verifyEqual(returned, expected);
    % Length mismatches expose caller slicing bugs loudly on both fields.
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", values, times, latitude=67.0, longitude=-48.8, ...
       toa=toa(1:3)), ...
       'icemodel:reconstruct:physicalValidity:precomputedGeometrySize');
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", values, times, latitude=67.0, longitude=-48.8, ...
       elevation=elevation(1:3)), ...
       'icemodel:reconstruct:physicalValidity:precomputedGeometrySize');
end

function test_physical_validity_uses_complete_posting_twilight(testCase)
    % A posting that starts below civil twilight but reaches twilight within
    % its support receives the same diffuse-light allowance as staging.
    times = datetime(2020, 3, 20, 0, 0, 0, 'TimeZone', 'UTC') ...
       + minutes(0:15:24 * 60 - 15).';
    point_elevation = icemodel.forcing.helpers.solarElevation( ...
       times, 67.0, -48.8);
    maximum_elevation = ...
       icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
       times, 67.0, -48.8, hours(1));
    bands = icemodel.forcing.reconstruct.solarElevationBands();
    target = find(point_elevation < bands.civil_twilight_deg ...
       & maximum_elevation >= bands.civil_twilight_deg ...
       & maximum_elevation < 0, 1);
    testCase.assertNotEmpty(target);

    value = 20;
    point_valid = icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", value, times(target), latitude=67.0, longitude=-48.8);
    interval_valid = icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", value, times(target), latitude=67.0, longitude=-48.8, ...
       interval=hours(1));
    precomputed_valid = icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", value, times(target), latitude=67.0, longitude=-48.8, ...
       elevation=maximum_elevation(target));
    testCase.verifyFalse(point_valid);
    testCase.verifyTrue(interval_valid);
    testCase.verifyTrue(precomputed_valid);
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.physicalValidity( ...
       "swd", value, times(target), latitude=67.0, longitude=-48.8, ...
       interval=-hours(1)), ...
       'icemodel:reconstruct:physicalValidity:negativeInterval');
end

function test_reconstruct_rejects_upward_shortwave_above_downwelling(testCase)
    % The production engine applies the same swu<=swd relation as metrics.
    series = icemodel.test.fixtures.makeReconstructSeries();
    n = height(series);
    series.swd = 100 * ones(n, 1);
    series.swu = 40 * ones(n, 1);
    gap = 100:110;
    series.swu(gap) = NaN;
    codes = icemodel.forcing.reconstruct.provenanceCodes();
    methods = struct('name', "too_bright", 'code', codes.climatology, ...
       'estimate', 120 * ones(n, 1), 'seasons', "all", 'buckets', 2:5);
    cm = struct('channel', "swu", 'methods', methods);
    result = icemodel.forcing.reconstruct.reconstructSeries(series, cm, ...
       interp_channels=string.empty(1, 0));
    testCase.verifyTrue(all(isnan(result.series.swu(gap))));
    testCase.verifyTrue(all(result.provenance.swu(gap) == codes.missing));
end

function test_reconstruct_processes_swd_before_swu(testCase)
   % Reverse caller order must still let swu use newly reconstructed swd.
   series = icemodel.test.fixtures.makeReconstructSeries();
   times = series.Properties.RowTimes;
   start = find(month(times) == 6 & day(times) == 15 ...
      & hour(times) == 8, 1);
   gap = (start:start + 7).';
   series.swd(:) = 50;
   series.swu = 25 * ones(height(series), 1);
   swd_estimate = series.swd;
   swu_estimate = series.swu;
   series.swd(gap) = NaN;
   series.swu(gap) = NaN;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   swd_method = struct('name', "swd_method", 'code', codes.climatology, ...
      'estimate', swd_estimate, 'seasons', "all", 'buckets', 2);
   swu_method = struct('name', "swu_method", 'code', codes.climatology, ...
      'estimate', swu_estimate, 'seasons', "all", 'buckets', 2);
   requests = [ ...
      struct('channel', "swu", 'methods', swu_method)
      struct('channel', "swd", 'methods', swd_method)];

   returned = icemodel.forcing.reconstruct.reconstructSeries( ...
      series, requests, latitude=67, longitude=-48.8, ...
      interp_channels=string.empty(1, 0));

   testCase.verifyEqual(returned.series.swd(gap), swd_estimate(gap));
   testCase.verifyEqual(returned.series.swu(gap), swu_estimate(gap));
end

function test_partition_precipitation_preserves_total(testCase)
    % Warm/cold totals enter exactly one phase; invalid inputs stay missing.
    Tf = icemodel.physicalConstant('Tf');
    [rain, snow] = icemodel.forcing.reconstruct.partitionPrecipitation( ...
       [1; 2; 0; -1; 3], [Tf; Tf - 1; Tf + 1; Tf; NaN]);
    testCase.verifyEqual(rain(1:3), [1; 0; 0]);
    testCase.verifyEqual(snow(1:3), [0; 2; 0]);
    testCase.verifyTrue(all(isnan(rain(4:5))));
    testCase.verifyTrue(all(isnan(snow(4:5))));
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.partitionPrecipitation([1; 2], Tf), ...
       'icemodel:reconstruct:partitionPrecipitation:sizeMismatch');
end

function test_audit_segments_split_disjoint_masks(testCase)
    % Disjoint masks never produce a row spanning an unfilled interval.
    times = datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:7).';
    rows = icemodel.forcing.reconstruct.auditSegments(times, ...
       logical([1; 1; 0; 0; 1; 1; 1; 0]), "tair", "method", "detail");
    testCase.verifyEqual(numel(rows), 2);
    testCase.verifyEqual(rows{1}{2}, times(1));
    testCase.verifyEqual(rows{1}{3}, times(2));
    testCase.verifyEqual(rows{2}{2}, times(5));
    testCase.verifyEqual(rows{2}{3}, times(7));
    testCase.verifyEqual(rows{1}{7}, 'method');
    testCase.verifyEmpty(icemodel.forcing.reconstruct.auditSegments( ...
       times, false(8, 1), "tair", "method", "detail"));
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.auditSegments( ...
       times, false(7, 1), "tair", "method", "detail"), ...
       'icemodel:reconstruct:auditSegments:sizeMismatch');
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.auditSegments( ...
       times(1), true, "tair", "method", "detail"), ...
       'icemodel:reconstruct:auditSegments:missingCadence');
end

function test_empty_audit_timezone_accepts_final_tier(testCase)
    % A UTC empty composition audit can accept later proxy-segment rows.
    series = icemodel.test.fixtures.makeReconstructSeries();
    times = series.Properties.RowTimes;
    codes = icemodel.forcing.reconstruct.provenanceCodes();
    methods = struct('name', {}, 'code', {}, 'estimate', {}, ...
       'seasons', {}, 'buckets', {});
    composed = icemodel.forcing.reconstruct.reconstructSeries( ...
       series, struct('channel', "tair", 'methods', methods));
    testCase.verifyEqual(composed.audit.start_time.TimeZone, 'UTC');

    composed.series.tair(10) = NaN;
    composed.provenance.tair(10) = codes.missing;
     proxy = struct('series', timetable(times, 261 * ones(numel(times), 1), ...
        'VariableNames', {'tair'}), 'name', "mar", 'code_name', "mar");
     opts = icemodel.forcing.reconstruct.setopts(required_channels="tair");
     calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
        times, composed.series.tair, proxy.series.tair, "tair", ...
        fit_years=unique(year(times)).');
     calibration.source = "mar";
     channel_plan = struct('channel', "tair", ...
        'proxy_calibrations', struct('source', "mar", ...
        'parameters', calibration));
     plan = struct('channels', channel_plan);
     [returned, provenance, audit] = ...
        icemodel.forcing.reconstruct.lastResortProxies( ...
        composed.series, composed.provenance, composed.audit, ...
        proxy, codes, opts, plan=plan);
    testCase.verifyTrue(isfinite(returned.tair(10)));
    testCase.verifyEqual(provenance.tair(10), codes.mar);
    testCase.verifyEqual(height(audit), 1);
end

function test_reconstruct_per_channel_cap_override(testCase)
   % B3/D-21: a per-channel override replaces the global tier-1 cap for
   % that channel only, and the (0, ceiling] contract holds at setopts
   % and again at the direct engine boundary.
   series = icemodel.test.fixtures.makeReconstructSeries();
   gap = (4000:4002).';                          % 3 h interior run
   series.tair(gap) = NaN;
   series.rh(gap) = NaN;
   methods = struct('name', {}, 'code', {}, 'estimate', {}, ...
      'seasons', {}, 'buckets', {});
   channel_methods = [ ...
      struct('channel', "tair", 'methods', methods), ...
      struct('channel', "rh", 'methods', methods)];

   % Default cap (6 h) bridges both 3 h gaps in tier 1.
   composed = icemodel.forcing.reconstruct.reconstructSeries( ...
      series, channel_methods, interp_channels=["tair", "rh"]);
   testCase.verifyTrue(all(isfinite(composed.series.tair(gap))));
   testCase.verifyTrue(all(isfinite(composed.series.rh(gap))));

   % A 2 h tair override rejects the 3 h gap while rh keeps the default.
   composed = icemodel.forcing.reconstruct.reconstructSeries( ...
      series, channel_methods, interp_channels=["tair", "rh"], ...
      cap_hours_by_channel=struct('tair', 2));
   testCase.verifyTrue(all(isnan(composed.series.tair(gap))));
   testCase.verifyTrue(all(isfinite(composed.series.rh(gap))));

   % setopts refuses values outside a channel ceiling; the engine enforces
   % the same ceiling on its own.
   testCase.verifyError(@() icemodel.forcing.reconstruct.setopts( ...
      cap_hours_by_channel=struct('albedo', 31)), ...
      'icemodel:reconstruct:setopts:capHours');
   testCase.verifyError(@() icemodel.forcing.reconstruct.setopts( ...
      cap_hours_by_channel=struct('tair', 7)), ...
      'icemodel:reconstruct:setopts:capHours');
   testCase.verifyError(@() icemodel.forcing.reconstruct.reconstructSeries( ...
      series, channel_methods, interp_channels=["tair", "rh"], ...
      cap_hours_by_channel=struct('tair', 7)), ...
      'icemodel:reconstruct:fillShortGaps:capHours');
end

%% Fixtures
