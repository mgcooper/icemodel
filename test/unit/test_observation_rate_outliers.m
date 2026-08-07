function tests = test_observation_rate_outliers
   %TEST_OBSERVATION_RATE_OUTLIERS Cover the cross-year observation-rate flag.
   tests = functiontests(localfunctions);
end

function summary = buildSummary(rates_mwe_per_day, site_ids, statuses)
   % Build the minimal summary contract the diagnostic consumes. Every
   % site-year uses a 100-day window so the requested rate is exact.
   n = numel(rates_mwe_per_day);
   window_days = 100;
   summary = table( ...
      "case" + string(1:n)', ...
      string(site_ids(:)), ...
      (2010:2010 + n - 1)', ...
      string(statuses(:)), ...
      rates_mwe_per_day(:) * window_days, ...
      repmat(datetime(2010, 6, 1), n, 1), ...
      repmat(datetime(2010, 6, 1) + days(window_days), n, 1), ...
      'VariableNames', {'case_id', 'site_id', 'year', 'status', ...
      'observation_intact_mwe', 'window_start', 'window_end'});
end

function test_compressed_site_year_is_flagged_against_its_own_station(testCase)
   % One station with a compressed year: the low year must flag, and only it.
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   returned = icemodel.verification.observationRateOutliers( ...
      buildSummary([0.030, 0.032, 0.028, 0.009], repmat("KAN_M", 1, 4), ...
      repmat("completed", 1, 4)), policy);

   expected = [false; false; false; true];
   testCase.verifyEqual(returned.rate_outlier_flag, expected)
   testCase.verifyEqual(returned.observed_rate_mwe_per_day(4), 0.009, ...
      AbsTol=1e-12)
   testCase.verifyEqual(returned.station_median_rate_mwe_per_day(1), 0.029, ...
      AbsTol=1e-12)
   % The ratio is what makes the flag auditable rather than a bare boolean.
   testCase.verifyEqual(returned.rate_ratio_to_station_median(4), ...
      0.009 / 0.029, RelTol=1e-12)
end

function test_station_without_enough_years_is_never_flagged(testCase)
   % A two-year station has no family, so its median carries no information.
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   returned = icemodel.verification.observationRateOutliers( ...
      buildSummary([0.030, 0.001], repmat("KAN_U", 1, 2), ...
      repmat("completed", 1, 2)), policy);

   testCase.verifyFalse(any(returned.rate_outlier_flag))
   testCase.verifyEqual(returned.station_scored_year_count, [2; 2])
end

function test_low_melt_station_is_scored_against_itself(testCase)
   % A uniformly low-melt station must not be flagged merely for being low,
   % which is why each station is compared against its own median.
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   returned = icemodel.verification.observationRateOutliers( ...
      buildSummary([0.030, 0.031, 0.029, 0.005, 0.0052, 0.0048], ...
      [repmat("KAN_L", 1, 3), repmat("KAN_HIGH", 1, 3)], ...
      repmat("completed", 1, 6)), policy);

   testCase.verifyFalse(any(returned.rate_outlier_flag))
end

function test_incomplete_site_years_are_excluded(testCase)
   % Only completed site-years carry a comparable observed total.
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   returned = icemodel.verification.observationRateOutliers( ...
      buildSummary([0.030, 0.032, 0.028, 0.009], repmat("KAN_M", 1, 4), ...
      ["completed", "completed", "completed", "excluded"]), policy);

   testCase.verifyEqual(height(returned), 3)
   testCase.verifyFalse(any(returned.rate_outlier_flag))
end

function test_nonpositive_window_yields_no_rate(testCase)
   % A zero-length window cannot produce a rate, and must not produce a flag.
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   summary = buildSummary([0.030, 0.032, 0.028], repmat("KAN_M", 1, 3), ...
      repmat("completed", 1, 3));
   summary.window_end(3) = summary.window_start(3);
   returned = icemodel.verification.observationRateOutliers(summary, policy);

   testCase.verifyTrue(isnan(returned.observed_rate_mwe_per_day(3)))
   testCase.verifyFalse(any(returned.rate_outlier_flag))
end
