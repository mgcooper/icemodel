function tests = test_reconstruct_production
   %TEST_RECONSTRUCT_PRODUCTION Verify the planner, driver, and report builder.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % Fixture tree: staged-met layout, proxy cache, out/qa/report dirs.
   [~, ~, ~, ~, cleanup] = icemodel.test.helpers.bootstrapTestEnvironment();
   testCase.TestData.cleanup = cleanup;
   root = tempname;
   testCase.TestData.root = root;
   mkdir(fullfile(root, 'met', 'promice'));
   mkdir(fullfile(root, 'proxies'));
   mkdir(fullfile(root, 'qa'));
   mkdir(fullfile(root, 'figures'));
   mkdir(fullfile(root, 'report'));
end

function teardown(testCase)
   % Remove the fixture tree.
   if isfolder(testCase.TestData.root)
      rmdir(testCase.TestData.root, 's')
   end
   clear testCase.TestData.cleanup
end

%% stationMethodPlan

function test_plan_admits_near_donor_and_gates_far_donor(testCase)
   % A collocated, correlated donor is admitted; a donor beyond the 60 km /
   % 600 m geometry gate never enters the candidate pool.
   [target, donor] = syntheticPair(0.05, -20);        % ~6 km, -20 m
   heights = struct('present', true, 'records', ...
      struct('station', "AWS5", 'height_cm', 250));
   donor.series.Properties.UserData = struct( ...
      'sensor_heights', heights, 'aws_types', [0 1]);
   far = donor;
   far.location.lat_wgs84 = far.location.lat_wgs84 + 2;   % ~220 km away
   far.station = "far";

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      [donor, far], emptyProxies(), channels="tair", seed=11);

   names = string({plan.channels(1).methods.name});
   testCase.verifyTrue(any(contains(names, "donor:dsta")));
   testCase.verifyFalse(any(contains(names, "donor:far")));
   donor_methods = contains(names, "donor:dsta");
   testCase.verifyFalse(any(string( ...
      {plan.channels(1).methods(donor_methods).seasons}) == "all"));
   testCase.verifyTrue(all(contains(string( ...
      {plan.channels(1).methods.audit_context_id}), "tair:")));
   testCase.verifyTrue(all(arrayfun(@(m) isscalar(m.buckets), ...
      plan.channels(1).methods)));
   parameters = plan.channels(1).methods(find(donor_methods, 1)).parameters;
   testCase.verifyEqual(parameters.height_provenance.sensor_heights, heights);
   testCase.verifyEqual(parameters.height_provenance.aws_types, [0 1]);
   testCase.verifyEqual(string(plan.station), "tsta");
end

function test_plan_applies_donor_cap_after_skill_evaluation(testCase)
   % An unusable nearest donor must not hide a farther eligible donor with
   % validated skill when max_donors is one.
   [target, donor] = syntheticPair(0.05, -20);
   nearest = donor;
   nearest.station = "nearest";
   nearest.location.lat_wgs84 = target.location.lat_wgs84 + 0.001;
   nearest.series.tair(:) = NaN;

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      [nearest, donor], emptyProxies(), channels="tair", ...
      max_donors=1, seed=11);

   names = string({plan.channels(1).methods.name});
   testCase.verifyTrue(any(startsWith(names, "donor:dsta")));
   testCase.verifyFalse(any(startsWith(names, "donor:nearest")));
end

function test_direct_plan_excludes_target_from_donors(testCase)
   % Direct planner callers cannot use a GC-Net spelling of the target's
   % own finite observations as a perfect held-out donor — through EVERY
   % alias spelling of every cataloged station, including the GEUS
   % PROMICE continuation codes (POLICY A8).
   pairs = [ ...
      "dy2", "DYE_2"
      "cp1", "CP1"
      "nae", "NASA_E"
      "nau", "NASA_U"
      "nse", "NASA_SE"
      "sdl", "Saddle"
      "sdm", "SouthDome"
      "tun", "TUNU_N"];
   for k = 1:size(pairs, 1)
      [target, ~] = syntheticPair(0.05, -20);
      target.station = pairs(k, 1);
      self = target;
      self.station = pairs(k, 2);
      self.family = "gcnet";

      plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
         self, emptyProxies(), channels="tair", seed=11);

      names = string({plan.channels(1).methods.name});
      denials = string(plan.channels(1).denials.candidate);
      testCase.verifyFalse(any(contains([names(:); denials(:)], ...
         "donor:" + pairs(k, 2))), ...
         "self-donor escaped exclusion: " + pairs(k, 2));
   end
end

function test_direct_plan_keeps_ktransect_neighbor_distinct(testCase)
   % The approved AWS10/KAN_U crosswalk relation is a nearby donor hypothesis,
   % not a same-station alias to remove from candidate evaluation.
   [target, donor] = syntheticPair(0.001, 10);
   target.station = "kanu";
   donor.station = "AWS10";
   donor.family = "ktransect";

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, donor, ...
      emptyProxies(), channels="tair", seed=11);

   names = string({plan.channels(1).methods.name});
   denials = string(plan.channels(1).denials.candidate);
   testCase.verifyTrue(any(contains([names(:); denials(:)], ...
      "donor:AWS10")));
end

function test_plan_applies_donor_cap_per_validation_stratum(testCase)
   % Distinct seasonal donors receive independent budgets when max_donors is
   % one; skill in one season must not evict the only usable donor in another.
   [target, winter] = syntheticPair(0.05, -20);
   summer = winter;
   winter.station = "winter";
   summer.station = "summer";
   summer.location.lat_wgs84 = target.location.lat_wgs84 + 0.06;
   warm = ismember(month(target.series.Properties.RowTimes), 4:9);
   winter.series.tair(warm) = NaN;
   summer.series.tair(~warm) = NaN;
   jja = find(year(target.series.Properties.RowTimes) == 2019 ...
      & month(target.series.Properties.RowTimes) == 6, 12, 'first');
   target.series.tair(jja) = NaN;

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      [winter, summer], emptyProxies(), channels="tair", ...
      max_donors=1, knot_candidates=0, min_overlap_hours=1000, ...
      n_gaps=2, seed=11);

   names = string({plan.channels(1).methods.name});
   testCase.verifyTrue(any(startsWith(names, "donor:winter")));
   testCase.verifyTrue(any(startsWith(names, "donor:summer")));
end

function test_plan_denies_wind_donor_that_does_not_beat_mar(testCase)
   % Wind donors have a channel-specific MAR benchmark in addition to the
   % common persistence/climatology admission baseline.
   [target, donor] = syntheticPair(0.05, -20);
   truth = target.series.wspd;
   target.series.wspd(5000:5011) = NaN;
   target.series.wspd(30000:30143) = NaN;
   donor.series.wspd = truth + 0.5 * sin((1:numel(truth)).' / 20);
   proxy_series = target.series;
   proxy_series.wspd = truth;
   proxy = struct('series', proxy_series, 'name', "mar", ...
      'code_name', "mar");

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, donor, ...
      proxy, channels="wspd", n_gaps=2, seed=11);

   names = string({plan.channels(1).methods.name});
   testCase.verifyFalse(any(startsWith(names, "donor:dsta")));
   donor_denial = startsWith(plan.channels(1).denials.candidate, ...
      "donor:dsta");
   donor_reasons = plan.channels(1).denials.reasons(donor_denial);
   testCase.verifyTrue(any(contains( ...
      donor_reasons, "calibrated MAR wind")), strjoin(donor_reasons, " | "));
end

function test_plan_fits_without_selection_draw_leakage(testCase)
    % Changing only held-out selection observations must not change fitted
    % donor parameters or the climatology estimate.
    [target, donor] = syntheticPair(0.05, -20);
    seed = 11;
    first = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
       donor, emptyProxies(), channels="tair", seed=seed);
    census = icemodel.forcing.reconstruct.gapCensus(target.series, ...
       channels="tair", latitude=target.location.lat_wgs84, ...
       longitude=target.location.lon_wgs84);
    planned = census.runs.bucket >= 2 ...
       & census.runs.bucket < numel( ...
       icemodel.forcing.reconstruct.bucketEdges());
    stratum = census.runs(find(planned, 1), :);
    bucket = stratum.bucket;
    season = stratum.season;
    seasons = ["DJF", "MAM", "JJA", "SON"];
    season_index = find(seasons == season, 1);
    draw = icemodel.forcing.reconstruct.syntheticMissingness( ...
       target.series, "tair", census.runs, ...
       years=first.split.years_selection, ...
       seed=seed + 100 * bucket + 1000 * season_index, ...
       n_gaps=icemodel.forcing.reconstruct.setopts().plan_n_gaps, ...
       bucket=bucket, season=season);
    changed = target;
    changed.series.tair(draw.mask) = changed.series.tair(draw.mask) + 0.25;
    second = icemodel.forcing.reconstruct.stationMethodPlan(changed, ...
       donor, emptyProxies(), channels="tair", seed=seed);

    names1 = string({first.channels(1).methods.name});
    names2 = string({second.channels(1).methods.name});
    donor1 = find(startsWith(names1, "donor:dsta"), 1);
    donor2 = find(startsWith(names2, "donor:dsta"), 1);
    testCase.assertNotEmpty(donor1);
    testCase.assertNotEmpty(donor2);
    testCase.verifyEqual(first.channels(1).methods(donor1).parameters, ...
       second.channels(1).methods(donor2).parameters);
    climatology1 = find(names1 == "climatology", 1);
    climatology2 = find(names2 == "climatology", 1);
    testCase.assertNotEmpty(climatology1);
    testCase.assertNotEmpty(climatology2);
    testCase.verifyEqual(first.channels(1).methods(climatology1).estimate, ...
       second.channels(1).methods(climatology2).estimate);
end

function test_plan_supports_role_reversal(testCase)
   % The same planner fills the DONOR family from the PROMICE-style target
   % — the role-contract acceptance shape (DesignSpec Decision 9): swap the
   % structs and nothing about the call changes.
   [target, donor] = syntheticPair(0.05, -20);

   reversed = icemodel.forcing.reconstruct.stationMethodPlan(donor, ...
      target, emptyProxies(), channels="tair", seed=12);

   testCase.verifyEqual(string(reversed.station), "dsta");
   testCase.verifyGreaterThan(numel(reversed.channels(1).methods), 0);
end

function test_plan_refuses_methods_without_evaluation_year(testCase)
   % A one-year record has no disjoint whole-year evaluation holdout, so
   % selection draws may diagnose candidates but cannot authorize filling.
   [target, donor] = syntheticPair(0.05, -20);
   keep = year(target.series.Properties.RowTimes) == 2019;
   target.series = target.series(keep, :);
   donor.series = donor.series(keep, :);

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      donor, emptyProxies(), channels="tair", seed=12);

   testCase.verifyEmpty(plan.split.years_evaluation);
   testCase.verifyEmpty(plan.channels(1).methods);
   testCase.verifyTrue(any(contains(string( ...
      plan.channels(1).denials.reasons), ...
      "no disjoint evaluation year")));
end

function test_direct_plan_rejects_precipitation_channels(testCase)
   % Direct callers cannot bypass total-source precipitation adoption policy.
   [target, donor] = syntheticPair(0.05, -20);

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.stationMethodPlan(target, donor, ...
      emptyProxies(), channels=["tair", "ppt"], seed=12), ...
      'icemodel:reconstruct:stationMethodPlan:plannedPrecipitation');
end

function test_plan_split_uses_only_joint_core_record_years(testCase)
   % A peripheral channel cannot fabricate holdout years outside core support.
   [target, donor] = syntheticPair(0.05, -20);
   outside_core = year(target.series.Properties.RowTimes) > 2019;
   for channel = ["tair", "rh", "psfc"]
      target.series.(channel)(outside_core) = NaN;
   end

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      donor, emptyProxies(), channels=["tair", "lwd"], seed=12);

   testCase.verifyEqual(plan.split.years_selection, 2019);
   testCase.verifyEmpty(plan.split.years_evaluation);
   testCase.verifyTrue(all(arrayfun(@(entry) isempty(entry.methods), ...
      plan.channels)));
end

function test_plan_denies_candidate_without_common_baseline_support(testCase)
   % Candidate skill cannot be inferred when the baseline has no paired draw.
   [target, donor] = syntheticPair(0.05, -20);

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      donor, emptyProxies(), channels="tair", seed=12, ...
      climatology_min_support=1e9);

   candidates = string(plan.channels(1).denials.candidate);
   reasons = string(plan.channels(1).denials.reasons);
   donor_denial = startsWith(candidates, "donor:dsta");
   testCase.verifyTrue(any(contains(reasons(donor_denial), ...
      "common baseline support")), ...
      "planner denials: " + strjoin(candidates + "=" + reasons, " | "));
end

function test_plan_reports_donor_with_insufficient_overlap(testCase)
   % A geometrically eligible donor with no finite overlap remains visible
   % in candidate diagnostics instead of disappearing before admission.
   [target, donor] = syntheticPair(0.05, -20);
   donor.series.tair(:) = NaN;

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      donor, emptyProxies(), channels="tair", seed=19);

   denials = plan.channels(1).denials;
   testCase.verifyTrue(any(string(denials.candidate) == "donor:dsta"));
   testCase.verifyTrue(any(contains(string(denials.reasons), ...
      "overlap", 'IgnoreCase', true)));
end

function test_plan_rejects_zero_overlap_proxy_candidates(testCase)
   % Zero-overlap proxy and empirical corrections remain recorded for
   % provenance but cannot enter held-out competition as raw identities.
   [target, donor] = syntheticPair(0.05, -20);
   seed = 23;
   proxy_series = target.series;
   split = icemodel.forcing.reconstruct.validationSplit(2019:2021, ...
      station=target.station, seed=seed);
   in_selection = ismember(year(target.series.Properties.RowTimes), ...
      split.years_selection);
   target.series.lwd(in_selection) = NaN;
   proxy = struct('series', proxy_series, 'name', "mar", ...
      'code_name', "mar");

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      donor([]), proxy, channels="lwd", core_channels="tair", seed=seed);

   channel_plan = plan.channels(1);
   testCase.verifyEqual(numel(channel_plan.proxy_calibrations), 1);
   testCase.verifyEqual( ...
      channel_plan.proxy_calibrations.parameters.n_overlap, 0);
   candidate_names = string(channel_plan.denials.candidate);
   reasons = string(channel_plan.denials.reasons);
   testCase.verifyTrue(any(candidate_names == "proxy:mar" ...
      & contains(reasons, "no positive finite target overlap")));
   testCase.verifyTrue(any(candidate_names == "estimator:brutsaert" ...
      & contains(reasons, "no positive finite target overlap")));
   method_names = string({channel_plan.methods.name});
   testCase.verifyFalse(any(startsWith(method_names, "proxy:") ...
      | startsWith(method_names, "estimator:")));
end

function test_plan_fits_temperature_elevation_adjustment_on_overlap(testCase)
   % A donor beyond the elevation threshold fits its lapse correction from
   % selection-year overlap and persists that evidence with the method.
   [target, donor] = syntheticPair(0.05, -200);
   dz = target.location.elev_m - donor.location.elev_m;
   expected_lapse = -0.006;
   donor.series.tair = target.series.tair - expected_lapse * dz;

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      donor, emptyProxies(), channels="tair", seed=23);

   names = string({plan.channels(1).methods.name});
   method = plan.channels(1).methods(find(startsWith(names, ...
      "donor:dsta"), 1));
   testCase.verifyEqual(method.parameters.elevation_adjustment.mode, ...
      "overlap_fit");
   testCase.verifyEqual( ...
      method.parameters.elevation_adjustment.n_overlap_hours, ...
      method.parameters.elevation_adjustment.n_overlap * 0.25);
   testCase.verifyEqual( ...
      method.parameters.elevation_adjustment.lapse_rate, ...
      expected_lapse, 'AbsTol', 1e-12);
   testCase.verifyEqual(method.parameters.fit_years, ...
      plan.split.years_selection);
   testCase.verifyEqual(method.uncertainty, "not_provided");
end

function test_plan_covers_coarse_donor_support(testCase)
   % An hourly donor must cover the target's support-held 15-minute axis
   % (four samples per posting), not one exact-match sample per hour —
   % otherwise the coverage gate denies skilled donors for cadence reasons.
   [target, donor] = syntheticPair(0.05, -20);
   donor.series = retime(donor.series, 'hourly', 'mean');

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      donor, emptyProxies(), channels="tair", seed=14);

   names = string({plan.channels(1).methods.name});
   testCase.verifyTrue(any(contains(names, "donor:dsta")));
end

function test_plan_stops_donor_at_final_support_interval(testCase)
   % A donor posting supports one cadence interval, never the remaining
   % tail of a longer target record.
   [target, donor] = syntheticPair(0.05, -20);
   donor.series = donor.series(1:(end - 100), :);
   donor_end = donor.series.Properties.RowTimes(end);

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      donor, emptyProxies(), channels="tair", seed=14);

   names = string({plan.channels(1).methods.name});
   method = plan.channels(1).methods(find(contains(names, ...
      "donor:dsta"), 1));
   tail = target.series.Properties.RowTimes >= donor_end + minutes(15);
   testCase.verifyTrue(all(isnan(method.estimate(tail))));
end

function test_plan_uses_channel_specific_observed_mask(testCase)
   % A reconstructed RH sample must not suppress simultaneous observed
   % temperature from the same donor.
   [target, donor] = syntheticPair(0.05, -20);
   times = donor.series.Properties.RowTimes;
   donor.observed_mask = timetable(times, ...
      true(numel(times), 1), false(numel(times), 1), ...
      'VariableNames', {'tair', 'rh'});

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      donor, emptyProxies(), channels="tair", seed=14);

   names = string({plan.channels(1).methods.name});
   testCase.verifyTrue(any(contains(names, "donor:dsta")));
end

function test_plan_evaluates_elevation_bracketing_candidate(testCase)
   % Two donors on opposite sides of the target elevation produce the
   % explicit policy-required elevation-interpolation candidate.
   [target, low] = syntheticPair(0.03, -200);
   high = low;
   low.station = "low";
   high.station = "high";
   high.location.elev_m = target.location.elev_m + 200;
   low.series.tair = target.series.tair + 1.2;
   high.series.tair = target.series.tair - 1.2;

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, ...
      [low, high], emptyProxies(), channels="tair", seed=14);

   names = string({plan.channels(1).methods.name});
   bracket = find(startsWith(names, "bracket:"), 1);
   testCase.assertNotEmpty(bracket);
   parameters = plan.channels(1).methods(bracket).parameters;
   adjustment = parameters.elevation_adjustment;
   testCase.verifyEqual(adjustment.mode, "elevation_bracket");
end

function test_plan_uses_persistence_baseline_for_short_gaps(testCase)
   % Bucket-1 candidates are graded on blocked <=6 h draws against
   % persistence, while longer admissions retain the climatology baseline.
   [target, donor] = syntheticPair(0.03, -20);

   plan = icemodel.forcing.reconstruct.stationMethodPlan(target, donor, ...
      emptyProxies(), channels="tair", seed=14, n_gaps=3, ...
      min_overlap_hours=24);

   methods = plan.channels(1).methods;
   short = methods([methods.buckets] == 1);
   testCase.assertNotEmpty(short);
   baselines = arrayfun(@(method) ...
      string(method.selection.baseline), short);
   testCase.verifyEqual(baselines, repmat("persistence", size(baselines)));
   baseline_rmse = arrayfun(@(method) ...
      method.selection.baseline_rmse, short);
   testCase.verifyTrue(all(isfinite(baseline_rmse)));
   % Persistence and each candidate retain the same finite comparison support
   % when no other draw removes these deterministic pre-gap anchors.
   expected_common = round(arrayfun(@(method) ...
      method.selection.n * method.selection.coverage, short));
   comparison_n = arrayfun(@(method) ...
      method.selection.comparison_n, short);
   testCase.verifyEqual(comparison_n, expected_common);
end

function test_selection_metrics_ignore_evaluation_year_truth(testCase)
   % Evaluation-year variability cannot set the native seam scale used for
   % selection admission.
   [target, donor] = syntheticPair(0.03, -20);
   baseline = icemodel.forcing.reconstruct.stationMethodPlan(target, donor, ...
      emptyProxies(), channels="tair", seed=14, n_gaps=3, ...
      selection_fraction=0.34, min_overlap_hours=24);
   altered = target;
   is_evaluation = ismember(year(altered.series.Properties.RowTimes), ...
      baseline.split.years_evaluation) & isfinite(altered.series.tair);
   rows = find(is_evaluation);
   altered.series.tair(rows) = 260 + 45 * (-1) .^ (1:numel(rows)).';

   returned = icemodel.forcing.reconstruct.stationMethodPlan(altered, donor, ...
      emptyProxies(), channels="tair", seed=14, n_gaps=3, ...
      selection_fraction=0.34, min_overlap_hours=24);

   expected = baseline.channels(1).methods;
   actual = returned.channels(1).methods;
   testCase.assertNotEmpty(expected);
   expected_id = string({expected.name}) + "|" + string([expected.buckets]) ...
      + "|" + string({expected.seasons});
   actual_id = string({actual.name}) + "|" + string([actual.buckets]) ...
      + "|" + string({actual.seasons});
   testCase.verifyEqual(sort(actual_id), sort(expected_id));
   for k = 1:numel(expected)
      match = find(actual_id == expected_id(k), 1);
      testCase.verifyEqual(actual(match).selection, expected(k).selection);
   end
end

function test_blocked_plan_rejects_invalid_selection_fraction(testCase)
   % The no-core early return must enforce the same open-interval split contract.
   [target, donor] = syntheticPair(0.05, -20);
   target.series.tair(:) = NaN;

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.stationMethodPlan(target, donor([]), ...
      emptyProxies(), channels="tair", core_channels="tair", seed=11, ...
      selection_fraction=1), ...
      'MATLAB:validators:mustBeInRange');
end

%% fillPromiceStation

function test_selected_data_root_resolves_supported_layouts(testCase)
   % Canonical, compact, flat, and custom met paths share one resolver.
   root = string(testCase.TestData.root);
   [data_root, met_root] = ...
      icemodel.forcing.reconstruct.selectedDataRoot( ...
      fullfile(root, 'input', 'met', 'promice'));
   testCase.verifyEqual(string(data_root), root);
   testCase.verifyEqual(string(met_root), fullfile(root, 'input', 'met'));
   [flat_data, flat_met] = ...
      icemodel.forcing.reconstruct.selectedDataRoot( ...
      fullfile(root, 'input', 'met'));
   testCase.verifyEqual(string(flat_data), root);
   testCase.verifyEqual(string(flat_met), fullfile(root, 'input', 'met'));
   [compact_data, compact_met] = ...
      icemodel.forcing.reconstruct.selectedDataRoot( ...
      fullfile(root, 'met', 'promice'));
   testCase.verifyEqual(string(compact_data), root);
   testCase.verifyEqual(string(compact_met), fullfile(root, 'met'));
   [custom_data, custom_met] = ...
      icemodel.forcing.reconstruct.selectedDataRoot( ...
      fullfile(root, 'custom', 'promice'));
   testCase.verifyEqual(string(custom_data), fullfile(root, 'custom'));
   testCase.verifyEqual(string(custom_met), fullfile(root, 'custom'));
end

function test_explicit_donor_list_excludes_target_station(testCase)
   % A caller-supplied target token must not enter its own held-out donor pool.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, true);
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), donor_sites="tsta", ...
      use_ktransect=false, use_gcnet=false, write=false, opts=opts);

   methods = string({result.plan.channels(1).methods.name});
   denials = string(result.plan.channels(1).denials.candidate);
   testCase.verifyFalse(any(contains(lower([methods(:); denials(:)]), ...
      "donor:tsta")));
end

function test_fill_station_keeps_ktransect_in_selected_data_root(testCase)
   % Canonical and flat custom met roots must not fall back to repository
   % donors; each selected-root fixture intentionally lacks height proof.
   for layout = ["canonical", "flat", "canonical_flat"]
      data_root = fullfile(testCase.TestData.root, 'selected-' + layout);
      station_root = data_root;
      if ismember(layout, ["canonical", "canonical_flat"])
         station_root = fullfile(data_root, 'input');
      end
      met_dir = fullfile(station_root, 'met', 'promice');
      mkdir(met_dir);
      writeFixtureStation(station_root, "tsta", 0, 0, true);
      if layout == "canonical_flat"
         source_file = fullfile(met_dir, ...
            'met_tsta_promice_20200101_20211231_15m.mat');
         met_dir = fullfile(station_root, 'met');
         copyfile(source_file, met_dir);
      end

      kt_root = fullfile(data_root, 'eval', 'ktransect');
      case_root = fullfile(kt_root, 'aws9');
      mkdir(case_root);
      donor = syntheticMet(0, false);
      heights = struct('present', false, 'records', struct([]));
      expected_metadata = ktransectMetadata("AWS9", heights);
      donor.Properties.UserData = ktransectMetadata("AWS6", heights);
      targets = struct('data', donor, ...
         'metadata', donor.Properties.UserData);
      evaluation_file = fullfile(case_root, 'observations.mat');
      save(evaluation_file, 'targets');
      location = struct('lat_wgs84', 67.0, 'lon_wgs84', -48.8, ...
         'elev_m', 1200);
      evaluation_rel = fullfile('aws9', 'observations.mat');
      leg = ktransectLeg(evaluation_file, evaluation_rel, expected_metadata);
       entry = struct('site_id', 'AWS9', 'site_location', location, ...
          'colocation', struct('ktransect', leg));
       manifest = struct('cases', entry);
       writeJsonFixture(fullfile(kt_root, 'manifest.json'), manifest);

       testCase.verifyError(@() ...
         icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
         met_dir=met_dir, donor_sites=string.empty(1, 0), ...
         use_ktransect=true, use_gcnet=false, write=false), ...
         ['icemodel:reconstruct:fillPromiceStation:' ...
         'ktransectIdentityMismatch']);
       donor.Properties.UserData = expected_metadata;
       targets.data = donor;
       targets.metadata = expected_metadata;
       save(evaluation_file, 'targets');
       legacy_leg = ...
          ktransectLeg(evaluation_file, evaluation_rel, expected_metadata);
       legacy_leg = rmfield(legacy_leg, ...
          {'evaluation_size_bytes', 'evaluation_sha256'});
       entry.colocation.ktransect = legacy_leg;
       manifest.cases = entry;
      writeJsonFixture(fullfile(kt_root, 'manifest.json'), manifest);
      % POLICY A3: absent sensor-height provenance downgrades the donor
      % with a warning; it never aborts the fill. Reaching this warning
      % also proves the legacy leg passed the identity checks above.
      testCase.verifyWarning(@() ...
         icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
         met_dir=met_dir, donor_sites=string.empty(1, 0), ...
         use_ktransect=true, use_gcnet=false, write=false), ...
         ['icemodel:reconstruct:fillPromiceStation:' ...
         'missingKtransectHeightProvenance']);
   end
end

function test_fill_station_rejects_ktransect_path_escape(testCase)
   % A manifest-relative traversal must not load a donor outside the
   % caller-selected K-transect root.
   data_root = fullfile(testCase.TestData.root, 'ktransect-escape');
   met_dir = fullfile(data_root, 'met', 'promice');
   mkdir(met_dir);
   writeFixtureStation(data_root, "tsta", 0, 0, true);
   kt_root = fullfile(data_root, 'eval', 'ktransect');
   mkdir(kt_root);
   location = struct('lat_wgs84', 67.0, 'lon_wgs84', -48.8, ...
      'elev_m', 1200);
   leg = struct('evaluation_file', ...
      fullfile('..', '..', '..', 'outside.mat'));
   entry = struct('site_id', 'AWS9', 'site_location', location, ...
      'colocation', struct('ktransect', leg));
   writeJsonFixture(fullfile(kt_root, 'manifest.json'), ...
      struct('cases', entry));

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=met_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=true, use_gcnet=false, write=false), ...
      ['icemodel:reconstruct:fillPromiceStation:' ...
      'ktransectPathOutsideRoot']);
end

function test_fill_station_rejects_changed_ktransect_donor_bytes(testCase)
   % A manifest-pinned donor cannot be altered after staging, even when its
   % station, DOI, and height metadata remain unchanged.
   data_root = fullfile(testCase.TestData.root, 'ktransect-tamper');
   met_dir = fullfile(data_root, 'met', 'promice');
   mkdir(met_dir);
   writeFixtureStation(data_root, "tsta", 0, 0, true);
   kt_root = fullfile(data_root, 'eval', 'ktransect');
   case_root = fullfile(kt_root, 'aws9');
   mkdir(case_root);
   heights = struct('present', true, 'records', ...
      struct('station', "AWS9", 'height_cm', 250));
   metadata = ktransectMetadata("AWS9", heights);
   donor = syntheticMet(0, false);
   donor.Properties.UserData = metadata;
   targets = struct('data', donor, 'metadata', metadata);
   evaluation_file = fullfile(case_root, 'observations.mat');
   save(evaluation_file, 'targets');
   evaluation_rel = fullfile('aws9', 'observations.mat');
   leg = ktransectLeg(evaluation_file, evaluation_rel, metadata);
   location = struct('lat_wgs84', 67.0, 'lon_wgs84', -48.8, ...
      'elev_m', 1200);
   entry = struct('site_id', 'AWS9', 'site_location', location, ...
      'colocation', struct('ktransect', leg));
   writeJsonFixture(fullfile(kt_root, 'manifest.json'), ...
      struct('cases', entry));

   changed_after_staging = true;
   save(evaluation_file, 'changed_after_staging', '-append');

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=met_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=true, use_gcnet=false, write=false), ...
      ['icemodel:reconstruct:fillPromiceStation:' ...
      'ktransectArtifactIdentityMismatch']);
end

function test_fill_station_keeps_measured_ktransect_albedo_eligible(testCase)
   % A8 treats the screened ratio of coincident measured swu/swd as a
   % measurement product, so K-transect albedo reaches donor evaluation.
   data_root = fullfile(testCase.TestData.root, 'ktransect-albedo');
   met_dir = fullfile(data_root, 'met', 'promice');
   mkdir(met_dir);
   writeFixtureStation(data_root, "tsta", 0, 0, true);
   station_file = fullfile(met_dir, ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(station_file, 'met');
    met = S.met;
    met.albedo(30001:30144) = NaN;
    save(station_file, 'met');
    recordNativeMetIdentity(data_root, "tsta", station_file);

   kt_root = fullfile(data_root, 'eval', 'ktransect');
   case_root = fullfile(kt_root, 'aws9');
   mkdir(case_root);
   donor = syntheticMet(0, false);
   heights = struct('present', true, 'records', ...
      struct('station', "AWS9", 'height_cm', 250));
   metadata = ktransectMetadata("AWS9", heights);
   donor.Properties.UserData = metadata;
   targets = struct('data', donor, 'metadata', metadata);
   evaluation_file = fullfile(case_root, 'observations.mat');
   save(evaluation_file, 'targets');
   location = struct('lat_wgs84', 67.0, 'lon_wgs84', -48.8, ...
      'elev_m', 1200);
   evaluation_rel = fullfile('aws9', 'observations.mat');
   leg = ktransectLeg(evaluation_file, evaluation_rel, metadata);
    entry = struct('site_id', 'AWS9', 'site_location', location, ...
       'colocation', struct('ktransect', leg));
    manifest_file = fullfile(kt_root, 'manifest.json');
    mismatched = entry;
    mismatched.site_location.lat_wgs84 = ...
       mismatched.site_location.lat_wgs84 + 1;
    writeJsonFixture(manifest_file, struct('cases', mismatched));
    opts = icemodel.forcing.reconstruct.setopts( ...
       required_channels="albedo", core_channels="albedo", ...
       plan_channels="albedo", interp_channels=string.empty(1, 0), ...
       last_resort_proxies=false, plan_n_gaps=1, seed=9);

    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
       met_dir=met_dir, donor_sites=string.empty(1, 0), ...
       use_ktransect=true, use_gcnet=false, write=false, opts=opts), ...
       ['icemodel:reconstruct:fillPromiceStation:' ...
       'ktransectIdentityMismatch']);
    writeJsonFixture(manifest_file, struct('cases', entry));
    result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=met_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=true, use_gcnet=false, write=false, opts=opts);

   evaluated = lower([string({result.plan.channels.methods.name}), ...
      string(result.plan.channels.denials.candidate).']);
   testCase.verifyTrue(any(startsWith(evaluated, "donor:aws9")));
end

function test_fill_station_ignores_filled_artifact_in_native_directory(testCase)
   % An output colocated with native met cannot feed reconstructed values
   % back into planning as if they were PROMICE observations.
   data_root = fullfile(testCase.TestData.root, 'mixed-products');
   met_dir = fullfile(data_root, 'met', 'promice');
   mkdir(met_dir);
   writeFixtureStation(data_root, "tsta", 0, 0, true);
   native_file = fullfile(met_dir, ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(native_file, 'met');
   native_first = S.met.tair(1);
   met = S.met;
   met.tair(:) = 290;
   met.Properties.UserData.gapfill_product = 'promice_filled';
   save(fullfile(met_dir, ...
      'met_tsta_promice_filled_20190101_20221231_15m.mat'), 'met');
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=9);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=met_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, write=false, opts=opts);

   testCase.verifyEqual(result.filled.tair(1), native_first, ...
      'AbsTol', 1e-12);
end

function test_fill_station_requires_native_metadata_site(testCase)
   % A correctly named file without an authoritative metadata station cannot
   % cross the native-observation trust boundary.
   data_root = fullfile(testCase.TestData.root, 'native-identity');
   met_dir = fullfile(data_root, 'met', 'promice');
   mkdir(met_dir);
   writeFixtureStation(data_root, "tsta", 0, 0, false);
   filename = fullfile(met_dir, ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   met.Properties.UserData = rmfield(met.Properties.UserData, 'site');
   save(filename, 'met');

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=met_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, write=false), ...
       'icemodel:reconstruct:fillPromiceStation:nativeIdentityMismatch');
end

function test_fill_station_accepts_normalized_native_metadata_site(testCase)
   % Canonical compact tokens must match the separator-bearing station names
   % written by the PROMICE builder.
   data_root = fullfile(testCase.TestData.root, 'native-alias-identity');
   met_dir = fullfile(data_root, 'met', 'promice');
   mkdir(met_dir);
   writeFixtureStation(data_root, "tsta", 0, 0, false);
   filename = fullfile(met_dir, ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   met.Properties.UserData.site = "T_STA";
   save(filename, 'met');
   recordNativeMetIdentity(data_root, "tsta", filename);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=met_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, write=false);
   testCase.verifyEqual(result.site, "tsta");
end

function test_fill_station_rejects_unpinned_native_payload_change(testCase)
   % A staged MAT changed after import cannot become reconstruction truth.
   data_root = fullfile(testCase.TestData.root, 'native-payload-identity');
   met_dir = fullfile(data_root, 'met', 'promice');
   mkdir(met_dir);
   writeFixtureStation(data_root, "tsta", 0, 0, false);
   filename = fullfile(met_dir, ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   met.tair(1) = met.tair(1) + 1;
   save(filename, 'met');

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=met_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, write=false), ...
      ['icemodel:reconstruct:fillPromiceStation:' ...
      'nativeArtifactIdentityMismatch']);
end

function test_fill_station_validates_proxy_metadata_identity(testCase)
   % Proxy calibration requires both the target point and a source-family
   % metadata signature; a filename alone cannot establish either identity.
   data_root = fullfile(testCase.TestData.root, 'proxy-identity');
   met_dir = fullfile(data_root, 'met', 'promice');
   mkdir(met_dir);
   writeFixtureStation(data_root, "tsta", 0, 0, false);
   writeFixtureMar(data_root, "tsta");
   proxy_file = fullfile(data_root, 'met', 'mar3.11', ...
      'met_tsta_mar3.11_20200101_20211231_15m.mat');
   S = load(proxy_file, 'mar_met');
   mar_met = S.mar_met;
   mar_met.Properties.UserData.site = "dsta";
   save(proxy_file, 'mar_met');

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=met_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, write=false), ...
      'icemodel:reconstruct:acceptanceWindow:proxyIdentityMismatch');

   mar_met.Properties.UserData.site = "tsta";
   mar_met.Properties.UserData.source_files = ...
      "MARv3.12-ERA5-15km-2020.nc";
   save(proxy_file, 'mar_met');
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=met_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, write=false), ...
       'icemodel:reconstruct:acceptanceWindow:proxyIdentityMismatch');
end

function test_fill_station_loads_only_acceptance_selected_proxy(testCase)
   % A malformed-date proxy with a wider saved axis cannot bypass the exact
   % validated file set pinned by the acceptance window.
   data_root = fullfile(testCase.TestData.root, 'proxy-selection');
   met_dir = fullfile(data_root, 'met', 'promice');
   mkdir(met_dir);
   writeFixtureStation(data_root, "tsta", 0, 0, false);
   writeFixtureMar(data_root, "tsta");
   valid_file = fullfile(data_root, 'met', 'mar3.11', ...
      'met_tsta_mar3.11_20200101_20211231_15m.mat');
   S = load(valid_file, 'mar_met');
   mar_met = S.mar_met;
   mar_met.Properties.RowTimes = ...
      datetime(2010, 1, 1, 'TimeZone', 'UTC') ...
      + hours((0:height(mar_met) - 1).');
   save(fullfile(fileparts(valid_file), ...
      'met_tsta_mar3.11_malformed_15m.mat'), 'mar_met');
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=5);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=met_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, write=false, opts=opts);

   calibrations = result.plan.channels(1).proxy_calibrations;
   mar = find(string({calibrations.source}) == "mar", 1);
   testCase.assertNotEmpty(mar);
   testCase.verifyGreaterThan(calibrations(mar).parameters.n_overlap, 0);
end

function test_public_station_tokens_reject_path_components(testCase)
   % Public site and donor identifiers must fail before entering any glob,
   % manifest path, or output filename.
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("../escape", ...
      write=false), ...
      'icemodel:reconstruct:mustBeStationToken:invalidToken');
   data_root = fullfile(testCase.TestData.root, 'token-validation');
   met_dir = fullfile(data_root, 'met', 'promice');
   mkdir(met_dir);
   writeFixtureStation(data_root, "tsta", 0, 0, true);
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=met_dir, donor_sites="../escape", write=false), ...
      'icemodel:reconstruct:mustBeStationToken:invalidToken');
   testCase.verifyError(@() ...
      icemodel.verification.report.buildGapFillReport(sites="../escape", ...
      render=false, filled_dir=fullfile(data_root, 'met', ...
      'promice_filled'), qa_dir=fullfile(data_root, 'qa'), ...
      fig_dir=fullfile(data_root, 'figures'), ...
      report_dir=fullfile(data_root, 'report')), ...
      'icemodel:reconstruct:mustBeStationToken:invalidToken');
end

function test_fill_station_writes_canonical_artifacts(testCase)
   % The driver plans, fills, adopts MAR precip, and ships through the
   % canonical writer with provenance channels and readiness verdicts.
     root = testCase.TestData.root;
     writeFixtureStation(root, "tsta", 0, 0, true);
     writeFixtureStation(root, "dsta", 0.05, -20, false);
     native_file = fullfile(root, 'met', 'promice', ...
        'met_tsta_promice_20200101_20211231_15m.mat');
      native = load(native_file, 'met');
      met = native.met;
      met.Properties.UserData.site = "T_STA";
      precip_samples = [49; 53; 57];
      met.ppt(precip_samples) = [0.2; 0.3; 0.3];
      met.rainf(precip_samples(2)) = 0.1;
      met.snowf(precip_samples(3)) = 0.2;
      % An inconsistent finite native pair must preserve rain, discard
      % snow, then adopt the proxy total and reconstruct its exact
      % complement on the delivered quarter-hour support.
      met.rainf(205:208) = 0.04;
      met.snowf(205:208) = 0.04;
      % Gap the upward shortwave so the shipped product must exercise the
      % post-tier derivation and carry its provenance column (A7/B10).
      met.swu(5001:5012) = NaN;
      save(native_file, 'met');
     recordNativeMetIdentity(root, "tsta", native_file);
     writeFixtureMar(root, "tsta");
    completeFixtureMarPrecip(root, "tsta");
    writeFixtureMerraPrecipPatch(root, "tsta");

    opts = icemodel.forcing.reconstruct.setopts( ...
       seed=13, min_overlap_hours=1e9);
    result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
       met_dir=fullfile(root, 'met', 'promice'), ...
      out_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
       donor_sites="dsta", use_ktransect=false, use_gcnet=false, ...
       opts=opts);

    % Canonical writemet naming under the forcings token directory.
    testCase.verifyTrue(isfile(result.met_file));
    testCase.verifyTrue(contains(result.met_file, "promice_filled"));
    producer = jsondecode(fileread(fullfile(root, 'qa', 'plans', ...
       'tsta-report-inputs.json')));
    testCase.verifyEqual(string(producer.path_base), ...
       "selected_data_root");
    testCase.verifyTrue(all(arrayfun(@(record) ...
       ~java.io.File(char(record.path)).isAbsolute(), ...
       producer.artifacts)));
    testCase.verifyFalse(contains(string(jsonencode(producer)), root));
    S = load(result.met_file);
   met = S.met;
   testCase.verifyTrue(ismember("tair_provenance", ...
      string(met.Properties.VariableNames)));
    codes = icemodel.forcing.reconstruct.provenanceCodes();
     % Both-missing native phases stay missing at reconstruction (POLICY
     % A10/D-18); only single-missing complements derive by arithmetic.
     testCase.verifyTrue(isnan(met.rainf(precip_samples(1))));
     testCase.verifyTrue(isnan(met.snowf(precip_samples(1))));
     testCase.verifyEqual(met.rainf(precip_samples(2:3)), [0.1; 0.1], ...
        'AbsTol', 1e-12);
     testCase.verifyEqual(met.snowf(precip_samples(2:3)), [0.2; 0.2], ...
        'AbsTol', 1e-12);
    testCase.verifyTrue(any(string(result.audit.method) == ...
       "complement:total_minus_rain"));
    testCase.verifyTrue(any(string(result.audit.method) == ...
       "complement:total_minus_snow"));
    testCase.verifyTrue(any(met.tair_provenance == codes.observed));
    metadata = met.Properties.UserData;
    testCase.verifyEqual(string(metadata.site), "tsta");
    testCase.verifyEqual(string(metadata.gapfill_engine_version), ...
       string(icemodel.internal.version()));
    policy_file = fullfile(icemodel.internal.fullpath, 'icemodel', ...
       '+icemodel', '+forcing', '+reconstruct', 'POLICY.md');
     testCase.verifyEqual(string(metadata.gapfill_policy_sha256), ...
        icemodel.verification.setup.fileSha256(policy_file));
     testCase.verifyEqual(string(metadata.gapfill_donors), "dsta");
     testCase.verifyEqual(string(metadata.gapfill_channels), ...
        string({result.plan.channels.channel}));
   % The seeded 3 h tair gap fills in tier 1.
   testCase.verifyTrue(any(met.tair_provenance == codes.bounded_interp));
   % MAR precipitation adoption stamps code 4 on placeholder samples.
   testCase.verifyTrue(any(met.ppt_provenance == codes.mar));
    testCase.verifyTrue(ismember("boom_height_provenance", ...
       string(met.Properties.VariableNames)));
    testCase.verifyTrue(all(isnan(met.boom_height(5001:5012))));
    testCase.verifyTrue(all(met.boom_height_provenance(5001:5012) == ...
       codes.missing));
    % Union-driven provenance attachment (A7): the derived upward-shortwave
    % channel ships its provenance column, and the seeded gap carries the
    % algebraic derivation code with values equal to albedo * swd (B10).
    testCase.verifyTrue(ismember("swu_provenance", ...
       string(met.Properties.VariableNames)));
    testCase.verifyTrue(all(met.swu_provenance(5001:5012) == ...
       codes.derived_shortwave));
    derived = met.swu_provenance == codes.derived_shortwave;
    testCase.verifyEqual(met.swu(derived), ...
       met.albedo(derived) .* met.swd(derived), 'AbsTol', 1e-12);
    testCase.verifyTrue(all(isfinite(met.ppt(met.ppt_provenance == ...
       codes.mar))));
    adopted = met.ppt_provenance == codes.mar;
    testCase.verifyEqual(met.ppt(adopted), ...
       met.rainf(adopted) + met.snowf(adopted), 'AbsTol', 1e-12);
     testCase.verifyTrue(all(ismember(met.rainf_provenance(adopted), ...
        [codes.observed, codes.mar])));
     testCase.verifyTrue(all(met.snowf_provenance(adopted) == codes.mar));
     support = (205:208).';
     testCase.verifyEqual(met.rainf(support), 0.04 + zeros(4, 1), ...
        'AbsTol', 1e-12);
     testCase.verifyTrue(all( ...
        icemodel.forcing.helpers.precipitationConsistency( ...
        met.ppt(support), met.rainf(support), met.snowf(support))));
     testCase.verifyGreaterThan(numel(unique(met.ppt(support))), 1);
     testCase.verifyEqual(met.rainf_provenance(support), ...
        repmat(codes.observed, 4, 1));
     testCase.verifyEqual(met.ppt_provenance(support), ...
        repmat(codes.mar, 4, 1));
     testCase.verifyEqual(met.snowf_provenance(support), ...
        repmat(codes.mar, 4, 1));
     testCase.verifyFalse(any(met.ppt_provenance == codes.merra2));
    precip_audit = startsWith(string(result.audit.method), ...
       "proxy:mar:precip_adoption");
     testCase.verifyTrue(any(precip_audit));
    % Final audit rows must not retain the provisional tair outage after
    % the MAR last-resort tier has filled it.
    testCase.verifyFalse(any(string(result.audit.method) == "unfilled" ...
       & string(result.audit.channel) == "tair"));
    plan_ids = strings(0, 1);
    for c = 1:numel(result.plan.channels)
       method_ids = string( ...
          {result.plan.channels(c).methods.audit_context_id}).';
       plan_ids = [plan_ids; method_ids]; %#ok<AGROW>
    end
    used_plan_context = startsWith(string(result.audit.context_id), ...
       string(result.audit.channel) + ":candidate-");
    testCase.verifyTrue(all(ismember( ...
       string(result.audit.context_id(used_plan_context)), plan_ids)));
    fixed_ids = string({result.plan.fixed_methods.audit_context_id}).';
    fixed_names = string({result.plan.fixed_methods.name}).';
    used_fixed_context = ismember(string(result.audit.method), fixed_names);
    testCase.verifyTrue(any(used_fixed_context));
    testCase.verifyTrue(all(ismember( ...
       string(result.audit.context_id(used_fixed_context)), fixed_ids)));
    testCase.verifyEqual(height(result.plan.audit_contexts), ...
       numel(unique(string(result.audit.context_id))));
    testCase.verifyTrue(all(ismember(string(result.audit.context_id), ...
       string(result.plan.audit_contexts.context_id))));
    % The caller's extreme overlap requirement reaches the planner: no
    % donor transfer may survive by silently reverting to defaults.
    names = strings(0, 1);
    for c = 1:numel(result.plan.channels)
       names = [names; ...
          string({result.plan.channels(c).methods.name}).']; %#ok<AGROW>
    end
    testCase.verifyFalse(any(startsWith(names, "donor:")));
   % Native samples are byte-identical to the input where observed.
   native = load(fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat'));
   observed = met.tair_provenance == codes.observed;
   testCase.verifyEqual(met.tair(observed), native.met.tair(observed));
   % Readiness verdicts exist for both fixture years.
   testCase.verifyEqual(sort(result.readiness.year), [2020; 2021]);
   testCase.verifyTrue(all(ismember( ...
      string(result.readiness.verdict_icemodel), ...
      ["ready", "not_forcing_ready"])));
   testCase.verifyTrue(all(ismember( ...
      string(result.readiness.verdict_snowmodel), ...
      ["ready", "not_forcing_ready"])));
   % Plan summary JSON and readiness CSV sidecars exist.
   testCase.verifyTrue(isfile(fullfile(root, 'qa', 'plans', ...
      'tsta-plan-summary.json')));
   plan_sidecar = load(fullfile(root, 'qa', 'plans', 'tsta-plan.mat'), ...
      'flat_run_findings_record');
   testCase.verifyTrue(istable(plan_sidecar.flat_run_findings_record));
   testCase.verifyTrue(isfile(fullfile(root, 'qa', 'ledger', ...
      'tsta-readiness.csv')));
   testCase.verifyTrue(isfile(fullfile(root, 'qa', 'splits', ...
      'tsta-split.json')));
end

function test_fill_family_ships_family_filled_product(testCase)
   % A non-promice family (the IMAU staging shape: ud.station identity,
   % ud.site_location point, no raw-source pins, path-pin-only family
   % manifest) must fill through the same driver and ship the
   % <family>_filled product. Finite native albedo ships as observed —
   % the explicit non-promice albedo-provenance policy (bead
   % icemodel-g1n.49): only the PROMICE builder is known to have injected
   % legacy fills, so only promice fails closed on missing raw provenance.
   root = testCase.TestData.root;
   writeFixtureFamilyStation(root, "tsta", "imau", true);
   writeFixtureMar(root, "tsta");
   completeFixtureMarPrecip(root, "tsta");
   writeFixtureMerraPrecipPatch(root, "tsta");

   opts = icemodel.forcing.reconstruct.setopts( ...
      seed=13, min_overlap_hours=1e9);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      family="imau", met_dir=fullfile(root, 'met', 'imau'), ...
      out_dir=fullfile(root, 'met', 'imau_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, opts=opts);

   % The product ships under the family token with the family identity.
   testCase.verifyTrue(isfile(result.met_file));
   testCase.verifyTrue(contains(result.met_file, "imau_filled"));
   S = load(result.met_file);
   met = S.met;
   returned = string(met.Properties.UserData.gapfill_product);
   expected = "imau_filled";
   testCase.verifyEqual(returned, expected);
   % No raw replay and an all-false winter mask: every native albedo
   % sample is finite in this fixture and every one ships as observed.
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   testCase.verifyTrue(all(met.albedo_provenance == codes.observed));
   % Seeded gaps fill completely: 3 h gaps in tier 1 and 36 h gaps by the
   % MAR last resort, so the required state channels ship gap-free.
   for name = ["tair", "rh", "psfc", "lwd"]
      testCase.verifyTrue(all(isfinite(met.(name))), name);
   end
   testCase.verifyTrue(any(met.tair_provenance == codes.bounded_interp));
   testCase.verifyTrue(any(met.tair_provenance == codes.mar));
   % The readiness ledger sidecar ships under the same site key.
   testCase.verifyTrue(isfile(fullfile(root, 'qa', 'ledger', ...
      'tsta-readiness.csv')));
end

function test_fill_half_hour_family_restores_each_posting(testCase)
   % Guarded K-transect input represents each half-hour posting as two
   % held rows. Reconstruction must map each pair back to its own posting,
   % then disaggregate a filled posting across exactly those two rows.
   root = testCase.TestData.root;
   writeFixtureFamilyStation(root, "half", "ktransect", false, 1800);
   filename = fullfile(root, 'met', 'ktransect', ...
      'met_half_ktransect_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   source_row = 199;
   gap_time = met.Properties.RowTimes(source_row);
   met.tair(source_row:source_row + 1) = NaN;
   save(filename, 'met');

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=4);
   result = icemodel.forcing.reconstruct.fillPromiceStation("half", ...
      family="ktransect", met_dir=fullfile(root, 'met', 'ktransect'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   support = result.filled.Properties.RowTimes >= gap_time ...
      & result.filled.Properties.RowTimes < gap_time + minutes(30);
   before = result.filled.Properties.RowTimes ...
      >= gap_time - minutes(30) ...
      & result.filled.Properties.RowTimes < gap_time;
   after = result.filled.Properties.RowTimes ...
      >= gap_time + minutes(30) ...
      & result.filled.Properties.RowTimes < gap_time + minutes(60);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   testCase.verifyEqual(nnz(support), 2);
   testCase.verifyGreaterThan(numel(unique( ...
      result.filled.tair(support))), 1);
   expected_posting = (unique(result.filled.tair(before)) ...
      + unique(result.filled.tair(after))) / 2;
   testCase.verifyEqual(mean(result.filled.tair(support)), ...
      expected_posting, 'AbsTol', 1e-9);
   testCase.verifyTrue(all(result.provenance.tair(support) ...
      == codes.bounded_interp));
   observed = result.provenance.tair == codes.observed;
   testCase.verifyEqual(result.filled.tair(observed), met.tair(observed));
end

function test_fill_family_promice_default_output_unchanged(testCase)
   % family="promice" must reproduce the historical implicit default
   % exactly: the family kwarg is a pure generalization (bead
   % icemodel-g1n.49) and the promice production output is the reference
   % this guard protects.
   data_root = fullfile(testCase.TestData.root, 'family-regression');
   met_dir = fullfile(data_root, 'met', 'promice');
   mkdir(met_dir);
   writeFixtureStation(data_root, "tsta", 0, 0, true);
   writeFixtureMar(data_root, "tsta");
   opts = icemodel.forcing.reconstruct.setopts( ...
      seed=13, min_overlap_hours=1e9);

   expected = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=met_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, write=false, opts=opts);
   returned = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      family="promice", met_dir=met_dir, ...
      donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, write=false, opts=opts);

   % Identical filled values, provenance codes, audit rows, and readiness
   % verdicts prove the explicit family token changes nothing. isequaln
   % compares the NaN-bearing timetables with NaN treated as equal.
   testCase.verifyTrue(isequaln(returned.filled, expected.filled));
   testCase.verifyTrue(isequaln(returned.provenance, ...
      expected.provenance));
   testCase.verifyTrue(isequaln(returned.audit, expected.audit));
   testCase.verifyTrue(isequaln(returned.readiness, expected.readiness));
end

function test_fill_family_auto_donors_degrade_without_registry(testCase)
   % donor_sites="auto" for a family absent from the donor registry must
   % degrade to no station donors with a warning instead of crashing or
   % glob-loading another family's stations (bead icemodel-g1n.49).
   root = testCase.TestData.root;
   writeFixtureFamilyStation(root, "tsta", "imau", true);
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=9);

   result = testCase.verifyWarning(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      family="imau", met_dir=fullfile(root, 'met', 'imau'), ...
      donor_sites="auto", use_ktransect=false, use_gcnet=false, ...
      write=false, opts=opts), ...
      'icemodel:reconstruct:fillPromiceStation:noDonorInventory');
   testCase.verifyEqual(result.site, "tsta");
end

function test_fill_station_rolls_back_partial_publication(testCase)
   % A sidecar install failure must not leave a newly published met file.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   completeFixtureNativePrecip(root, "tsta");
   qa_dir = fullfile(root, 'qa-atomic');
   plans_dir = fullfile(qa_dir, 'plans');
   mkdir(plans_dir);
   mkdir(fullfile(qa_dir, 'ledger'));
   [status, message] = system(sprintf('chmod 500 "%s"', plans_dir));
   testCase.assertEqual(status, 0, message);
   permission_cleanup = onCleanup(@() system( ...
      sprintf('chmod 700 "%s"', plans_dir)));
   out_dir = fullfile(root, 'met', 'promice_filled-atomic');
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=5);

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), out_dir=out_dir, ...
      qa_dir=qa_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, opts=opts), ...
      'icemodel:reconstruct:fillPromiceStation:installFailed');

   testCase.verifyEmpty(dir(fullfile(out_dir, ...
      'met_tsta_promice_filled_*.mat')));
   testCase.verifyFalse(isfile(fullfile(qa_dir, 'splits', ...
      'tsta-split.json')));
   clear permission_cleanup
end

function test_winter_albedo_stamp_reenters_as_fillable(testCase)
   % Albedo samples carrying exactly the native builder's winter constant
   % in its stamp months are legacy fills: the driver masks them on
   % input, fills them through the tiers/last resort, and restores the
   % constant with the constant code only where nothing better fills —
   % none may ship labeled observed.
    root = testCase.TestData.root;
    writeFixtureStation(root, "tsta", 0, 0, true);
   % Stamp one winter month with the exact constant, as the native
   % builder does.
   f = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(f);
   met = S.met;
   stamp = icemodel.parameterLookup('promice_winter_albedo');
   in_dec = month(met.Properties.RowTimes) == 12 ...
      & year(met.Properties.RowTimes) == 2020;
   met.albedo(in_dec) = stamp;
   raw_times = met.Properties.RowTimes(1:4:end);
   raw_albedo = met.albedo(1:4:end);
   raw_albedo(month(raw_times) == 12 & year(raw_times) == 2020) = NaN;
   source_file = writeRawPromiceAlbedo( ...
      root, "tsta", raw_times, raw_albedo);
    met = pinRawSource(met, source_file);
    save(f, 'met');
    recordNativeMetIdentity(root, "tsta", f);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
       met_dir=fullfile(root, 'met', 'promice'), ...
       out_dir=fullfile(root, 'met', 'promice_filled'), ...
       qa_dir=fullfile(root, 'qa'), ...
       donor_sites=string.empty(1, 0), use_ktransect=false, ...
       use_gcnet=false, ...
       write=false, opts=icemodel.forcing.reconstruct.setopts( ...
       required_channels=["tair", "albedo"], core_channels="tair", ...
       plan_channels=["tair", "albedo"], interp_channels="tair", ...
       climatology_min_support=1e9, last_resort_proxies=false, ...
       plan_n_gaps=1, seed=13));

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   in_year = year(result.filled.Properties.RowTimes) == 2020;
   stamped = in_year ...
      & month(result.filled.Properties.RowTimes) == 12;
   % Every stamped sample is finite with NON-observed provenance.
   testCase.verifyTrue(all(isfinite(result.filled.albedo(stamped))));
    testCase.verifyTrue(all(result.provenance.albedo(stamped) ...
       == codes.constant));
    % The bridge audits under its own method; it meets the anchors
    % continuously so no seam note is expected here.
    bridge_audit = string(result.audit.method) == ...
       "winter_albedo_bridge";
    testCase.verifyTrue(any(contains(string( ...
       result.audit.detail(bridge_audit)), "seasonal bridge")));
    % The bridge interior never dips below the dry-snow floor (POLICY
    % B13/D-15a); the seam blend may taper the few edge samples toward
    % the adjoining observations.
    interior = stamped ...
       & day(result.filled.Properties.RowTimes) >= 5 ...
       & day(result.filled.Properties.RowTimes) <= 25;
    testCase.verifyTrue(all(result.filled.albedo(interior) ...
       >= stamp - 1e-12));
   % Genuinely observed albedo elsewhere keeps observed provenance.
    testCase.verifyTrue(any(result.provenance.albedo(~stamped) ...
       == codes.observed));
    testCase.verifyFalse(isfolder(fullfile(root, 'qa', 'splits')));
end

function test_window_bounds_product_to_staged_proxy_span(testCase)
   % POLICY A6/D-17: the staged proxy window clips the product. With MAR
   % staged for 2020 only, the 2020-2021 native record ships a 2020-only
   % product and ledger; the 2021 native samples never reconstruct.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false, false);
   writeFixtureMarWindow(root, "tsta", ...
      datetime(2020, 1, 1, 'TimeZone', 'UTC'), ...
      datetime(2020, 12, 31, 'TimeZone', 'UTC'));

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=7);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   times = result.filled.Properties.RowTimes;
   testCase.verifyTrue(all(year(times) == 2020));
   testCase.verifyEqual(result.readiness.year, 2020);

   % A record wholly outside the staged window has no product span.
   writeFixtureMarWindow(root, "tsta", ...
      datetime(2030, 1, 1, 'TimeZone', 'UTC'), ...
      datetime(2030, 12, 31, 'TimeZone', 'UTC'));
   delete(fullfile(root, 'met', 'mar3.11', ...
      'met_tsta_mar3.11_20200101_20201231_15m.mat'));
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts), ...
      'icemodel:reconstruct:fillPromiceStation:windowRecordDisjoint');

   % A production refusal transactionally retires every stale product
   % artifact while preserving unrelated station files.
   out_dir = fullfile(root, 'met', 'promice_filled');
   qa_dir = fullfile(root, 'qa-disjoint');
   mkdir(out_dir);
   mkdir(fullfile(qa_dir, 'plans'));
   mkdir(fullfile(qa_dir, 'ledger'));
   mkdir(fullfile(qa_dir, 'splits'));
   stale = [ ...
      string(fullfile(out_dir, ...
      'met_tsta_promice_filled_20200101_20211231_15m.mat')); ...
      string(fullfile(qa_dir, 'plans', 'tsta-plan.mat')); ...
      string(fullfile(qa_dir, 'plans', 'tsta-plan-summary.json')); ...
      string(fullfile(qa_dir, 'plans', 'tsta-report-inputs.json')); ...
      string(fullfile(qa_dir, 'ledger', 'tsta-readiness.csv')); ...
      string(fullfile(qa_dir, 'splits', 'tsta-split.json'))];
   for file = stale.'
      fid = fopen(file, 'w');
      fprintf(fid, 'stale');
      fclose(fid);
   end
   foreign = fullfile(out_dir, ...
      'met_keep_promice_filled_20200101_20211231_15m.mat');
   fid = fopen(foreign, 'w');
   fprintf(fid, 'keep');
   fclose(fid);
   invoke = @() icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), out_dir=out_dir, ...
      qa_dir=qa_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, opts=opts);
   testCase.verifyError(invoke, ...
      'icemodel:reconstruct:fillPromiceStation:windowRecordDisjoint');
   testCase.verifyFalse(any(arrayfun(@isfile, stale)));
   testCase.verifyTrue(isfile(foreign));
   % Repeating the refusal with no prior station files is still exact.
   testCase.verifyError(invoke, ...
      'icemodel:reconstruct:fillPromiceStation:windowRecordDisjoint');
end

function test_empty_proxy_window_refuses_product(testCase)
   % POLICY A6: no validated MAR/MERRA file means there is no fillable
   % product span. Donors and climatology must not extend that empty span.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false, false);
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      plan_n_gaps=1, seed=7);

   invoke = @() icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);
   testCase.verifyError(invoke, ...
      'icemodel:reconstruct:fillPromiceStation:noProxyWindow');
end

function test_invalid_proxy_window_retires_stale_product(testCase)
   % An A6 validation error occurs before product construction, but it must
   % still retire a prior station artifact so stale bytes cannot survive.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false, false);
   writeFixtureMar(root, "tsta");
   proxy_file = fullfile(root, 'met', 'mar3.11', ...
      'met_tsta_mar3.11_20200101_20211231_15m.mat');
   S = load(proxy_file, 'mar_met');
   mar_met = S.mar_met;
   mar_met.Properties.RowTimes = mar_met.Properties.RowTimes + minutes(7);
   save(proxy_file, 'mar_met');

   out_dir = fullfile(root, 'met', 'promice_filled');
   qa_dir = fullfile(root, 'qa-invalid-window');
   mkdir(out_dir);
   mkdir(fullfile(qa_dir, 'plans'));
   stale = [ ...
      string(fullfile(out_dir, ...
      'met_tsta_promice_filled_20200101_20211231_15m.mat')); ...
      string(fullfile(qa_dir, 'plans', 'tsta-plan.mat'))];
   for file = stale.'
      fid = fopen(file, 'w');
      fprintf(fid, 'stale');
      fclose(fid);
   end
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair");
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), out_dir=out_dir, ...
      qa_dir=qa_dir, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, opts=opts), ...
      'icemodel:reconstruct:acceptanceWindow:proxyWindowMismatch');
   testCase.verifyFalse(any(arrayfun(@isfile, stale)));
end

function test_proxy_window_inside_guarded_posting_keeps_boundary(testCase)
   % A proxy window starting at 00:15 still owns the guarded 00:00 hourly
   % posting. The shipped artifact starts at the exact proxy boundary.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false, false);
   t0 = datetime(2020, 1, 1, 0, 15, 0, 'TimeZone', 'UTC');
   t1 = datetime(2020, 12, 31, 'TimeZone', 'UTC');
   writeFixtureMarWindow(root, "tsta", t0, t1);
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=7);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);
   testCase.verifyEqual(result.filled.Properties.RowTimes(1), t0);
   testCase.verifyEqual(result.filled.Properties.RowTimes(end), ...
      t1 + hours(23) + minutes(45));
end

function test_window_union_includes_same_source_extender(testCase)
   % POLICY A6: a wider-duration file must not silently drop a staged
   % sibling that extends coverage (the production EGP case, where a
   % long 2012-2018-style file outlasted the record-start file and cost
   % the product its staged 2019). Here the 18-month anchor is widest
   % and a 12-month sibling extends past its end: both files pin, the
   % window is the union, and the driver fills from the extender span.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false, false);
   writeFixtureMarWindow(root, "tsta", ...
      datetime(2020, 1, 1, 'TimeZone', 'UTC'), ...
      datetime(2021, 6, 30, 'TimeZone', 'UTC'));
   % Two extenders that OVERLAP EACH OTHER beyond the anchor prove the
   % merge dedups against everything already combined, not just the
   % anchor span (overlapping partial restagings are the D-22 world).
   writeFixtureMarWindow(root, "tsta", ...
      datetime(2021, 4, 1, 'TimeZone', 'UTC'), ...
      datetime(2021, 10, 31, 'TimeZone', 'UTC'));
   writeFixtureMarWindow(root, "tsta", ...
      datetime(2021, 8, 1, 'TimeZone', 'UTC'), ...
      datetime(2021, 12, 31, 'TimeZone', 'UTC'));
   location = struct('lat_wgs84', 67.0, 'lon_wgs84', -48.8);
   [window, files] = icemodel.forcing.reconstruct.acceptanceWindow( ...
      "tsta", met_dir=fullfile(root, 'met', 'promice'), ...
      location=location);
   testCase.verifyEqual(window, ...
      [datetime(2020, 1, 1, 'TimeZone', 'UTC'), ...
      datetime(2021, 12, 31, 23, 45, 0, 'TimeZone', 'UTC')]);
   testCase.verifyEqual(numel(files), 3);

   % The driver's proxy union serves the extender-only span: a November
   % 2021 tair outage (outside the anchor file) still fills from MAR.
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   gap = met.Properties.RowTimes ...
      >= datetime(2021, 11, 1, 'TimeZone', 'UTC') ...
      & met.Properties.RowTimes < datetime(2021, 11, 3, 'TimeZone', 'UTC');
   met.tair(gap) = NaN;
   save(filename, 'met');
   recordNativeMetIdentity(root, "tsta", filename);
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      climatology_min_support=1e9, plan_n_gaps=1, seed=7);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   gap_out = result.filled.Properties.RowTimes ...
      >= datetime(2021, 11, 1, 'TimeZone', 'UTC') ...
      & result.filled.Properties.RowTimes ...
      < datetime(2021, 11, 3, 'TimeZone', 'UTC');
   testCase.verifyTrue(all(isfinite(result.filled.tair(gap_out))));
   testCase.verifyTrue(any(result.provenance.tair(gap_out) == codes.mar));
end

function test_modis_albedo_tier_outranks_rcm_last_resort(testCase)
   % POLICY A11/B12: a residual albedo gap adopts staged MODIS daily
   % albedo (code 7, its own audit method) ahead of the RCM last resort,
   % even when a MAR proxy could fill the same span; the daily value
   % attaches through the single daily->met-cadence helper.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   % A July week of missing albedo is beyond the tier-1 cap and outside
   % the winter mask, so only last-resort tiers can own it.
   gap = met.Properties.RowTimes >= datetime(2020, 7, 6, 'TimeZone', 'UTC') ...
      & met.Properties.RowTimes < datetime(2020, 7, 13, 'TimeZone', 'UTC');
   met.albedo(gap) = NaN;
   save(filename, 'met');
   recordNativeMetIdentity(root, "tsta", filename);
   writeFixtureMar(root, "tsta");
   writeFixtureModis(root, "tsta");

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels=["tair", "albedo"], core_channels="tair", ...
      plan_channels=["tair", "albedo"], interp_channels="tair", ...
      climatology_min_support=1e9, plan_n_gaps=1, seed=11);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      modis_dir=fullfile(root, 'userdata', 'modis'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   gap_out = result.filled.Properties.RowTimes ...
      >= datetime(2020, 7, 6, 'TimeZone', 'UTC') ...
      & result.filled.Properties.RowTimes ...
      < datetime(2020, 7, 13, 'TimeZone', 'UTC');
   % The whole gap fills from MODIS, never the MAR proxy.
   testCase.verifyTrue(all(isfinite(result.filled.albedo(gap_out))));
   testCase.verifyTrue(all(result.provenance.albedo(gap_out) ...
      == codes.modis));
   testCase.verifyTrue(any(string(result.audit.method) ...
      == "modis:daily_albedo"));
   % Away from any seam taper the constant daily value passes through
   % the cadence attachment exactly.
   mid = gap_out & day(result.filled.Properties.RowTimes) == 9;
   testCase.verifyEqual(result.filled.albedo(mid), ...
      0.6 * ones(nnz(mid), 1), 'AbsTol', 1e-9);
   % Native albedo elsewhere keeps observed provenance.
   testCase.verifyTrue(any(result.provenance.albedo(~gap_out) ...
      == codes.observed));
end

function test_boom_height_gaps_fail_closed_without_visit_registry(testCase)
   % Station handovers do not prove maintenance-bounded sensor geometry, so
   % a short gap remains missing in the PRODUCT; the runtime fallback chain
   % (measured -> interpolated -> nominal) owns geometry, so the gap never
   % grades a readiness verdict (POLICY A3).
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   met.boom_height(5001:5012) = NaN;
   metadata = met.Properties.UserData;
   metadata.station_transition_times = ...
       met.Properties.RowTimes(8005);
    met.Properties.UserData = metadata;
    save(filename, 'met');
    recordNativeMetIdentity(root, "tsta", filename);

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=5);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   testCase.verifyTrue(all(isnan(result.filled.boom_height(5001:5012))));
   testCase.verifyTrue(all(result.provenance.boom_height(5001:5012) == ...
      icemodel.forcing.reconstruct.provenanceCodes().missing));
   testCase.verifyFalse(any(string(result.audit.method) == ...
      "maintenance_bounded_linear"));
   testCase.verifyFalse(any(contains( ...
      string(result.readiness.reason_icemodel), "boom_height")));
   testCase.verifyFalse(any(contains( ...
      string(result.readiness.reason_snowmodel), "boom_height")));
end

function test_raw_albedo_mask_removes_every_legacy_fill(testCase)
   % A builder-marked albedo series uses its raw L3 finite mask, so both
   % winter constants and linear/end fills re-enter reconstruction as
   % missing while true hourly observations retain observed provenance.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
    raw_times = met.Properties.RowTimes(1:4:end);
    raw_albedo = nan(numel(raw_times), 1);
    stamp = icemodel.parameterLookup('promice_winter_albedo');
    raw_albedo(1000) = 0.7;
    raw_albedo(1100) = stamp;
    raw_albedo(5000) = 0.7;
    raw_albedo(5001) = 2;
    raw_albedo(5100) = 0;
    raw_albedo(5200) = 1;
   source_file = writeRawPromiceAlbedo(root, "tsta", ...
      raw_times, raw_albedo);

   % Simulate the legacy builder: every missing source value became finite,
   % the summer observation stays raw, and the winter observation is overwritten.
   met.albedo(:) = 0.6;
    winter_overwritten = met.Properties.RowTimes >= raw_times(1000) ...
       & met.Properties.RowTimes < raw_times(1000) + hours(1);
    native_constant_support = met.Properties.RowTimes >= raw_times(1100) ...
       & met.Properties.RowTimes < raw_times(1100) + hours(1);
    observed_support = met.Properties.RowTimes >= raw_times(5000) ...
      & met.Properties.RowTimes < raw_times(5000) + hours(1);
    invalid_support = met.Properties.RowTimes >= raw_times(5001) ...
       & met.Properties.RowTimes < raw_times(5001) + hours(1);
    zero_support = met.Properties.RowTimes >= raw_times(5100) ...
       & met.Properties.RowTimes < raw_times(5100) + hours(1);
    one_support = met.Properties.RowTimes >= raw_times(5200) ...
       & met.Properties.RowTimes < raw_times(5200) + hours(1);
    met.albedo(observed_support) = 0.7;
    met.albedo(zero_support) = 0;
    met.albedo(one_support) = 1;
    met.albedo(winter_overwritten | native_constant_support) = stamp;
   met = pinRawSource(met, source_file);
    met.Properties.UserData.albedo_policy = ...
       "albedo from PROMICE L3 albedo; fillPromiceAlbedo(fillwinter=1)";
    save(filename, 'met');
    recordNativeMetIdentity(root, "tsta", filename);

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="albedo", core_channels="tair", ...
      plan_channels="albedo", interp_channels=string.empty(1, 0), ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      qa_dir=fullfile(root, 'qa-raw-albedo'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   testCase.verifyTrue(all(result.provenance.albedo(observed_support) ...
      == codes.observed));
    testCase.verifyFalse(any(result.provenance.albedo(winter_overwritten) ...
       == codes.observed));
    testCase.verifyTrue(all(result.provenance.albedo(native_constant_support) ...
       == codes.observed));
    testCase.verifyEqual(result.filled.albedo(native_constant_support), ...
       repmat(stamp, nnz(native_constant_support), 1));
    testCase.verifyFalse(any(result.provenance.albedo(invalid_support) ...
       == codes.observed));
    testCase.verifyTrue(all(result.provenance.albedo( ...
       zero_support | one_support) == codes.observed));
   invented = find(~observed_support, 1);
   testCase.verifyNotEqual(result.provenance.albedo(invented), ...
      codes.observed);
end

function test_legacy_albedo_requires_raw_source(testCase)
   % Silently treating a legacy-filled albedo artifact as observed is less
   % safe than refusing a stale/missing raw-source provenance link.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   met.Properties.UserData.source_file = fullfile(root, 'missing.nc');
   met.Properties.UserData.albedo_policy = ...
      "albedo from PROMICE L3 albedo; fillPromiceAlbedo(fillwinter=1)";
   save(filename, 'met');

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false), ...
      'icemodel:reconstruct:fillPromiceStation:missingAlbedoSource');
end

function test_legacy_albedo_source_must_stay_in_selected_root(testCase)
   % A staged artifact copied into another data root cannot follow its stale
   % raw-source link back into the original root.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   original = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   selected = fullfile(root, 'selected');
   selected_met = fullfile(selected, 'met', 'promice');
   mkdir(selected_met);
   copyfile(original, selected_met);

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=selected_met, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, write=false), ...
      'icemodel:reconstruct:fillPromiceStation:albedoSourceOutsideRoot');
end

function test_legacy_albedo_source_rebinds_after_root_relocation(testCase)
   % A copied data tree resolves historical absolute raw-source provenance to
   % the unique same-named source under its new selected root.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   original = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(original, 'met');
   recorded_source = string(S.met.Properties.UserData.source_file);
   selected = fullfile(root, 'selected');
   selected_met = fullfile(selected, 'met', 'promice');
   selected_raw = fullfile(selected, 'raw');
   mkdir(selected_met);
   mkdir(selected_raw);
    copyfile(original, selected_met);
    copyfile(recorded_source, selected_raw);
    selected_file = fullfile(selected_met, ...
       'met_tsta_promice_20200101_20211231_15m.mat');
    recordNativeMetIdentity(selected, "tsta", selected_file);
    writeFixtureMar(selected, "tsta");

   returned = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=selected_met, donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, write=false);

   testCase.verifyEqual(string(returned.site), "tsta");
end

function test_raw_promice_replay_rejects_changed_source_bytes(testCase)
   % Aggregate support counts are insufficient identity: an unrelated raw-file
   % metadata edit must still invalidate replay of builder-derived masks.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   source_file = string(S.met.Properties.UserData.source_file);
   ncwriteatt(source_file, '/', 'post_stage_mutation', 'changed bytes');

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false), ...
       'icemodel:reconstruct:fillPromiceStation:rawSourceIdentityMismatch');
end

function test_promice_native_without_byte_identity_is_rejected(testCase)
   % A path-only legacy manifest cannot authorize PROMICE target bytes;
   % A1/D-22 requires refreshing the manifest's size+SHA-256 identity.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   met.Properties.UserData = rmfield(met.Properties.UserData, ...
      {'source_size_bytes', 'source_sha256'});
   save(filename, 'met');
   recordNativeMetIdentity(root, "tsta", filename);
   manifest_file = fullfile(root, 'eval', 'promice', 'manifest.json');
   manifest = jsondecode(fileread(manifest_file));
   manifest.cases.colocation.promice = struct('met_files', ...
      "promice/met_tsta_promice_20200101_20211231_15m.mat");
   writeJsonFixture(manifest_file, manifest);

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false), ...
      'icemodel:reconstruct:fillPromiceStation:missingNativeArtifactIdentity');
end

function test_finite_albedo_requires_builder_provenance(testCase)
   % Finite staged albedo without a declared builder policy cannot be
   % distinguished from legacy interpolation and must fail closed.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   metadata = met.Properties.UserData;
   metadata = rmfield(metadata, 'albedo_policy');
   met.Properties.UserData = metadata;
   save(filename, 'met');

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false), ...
      'icemodel:reconstruct:fillPromiceStation:missingAlbedoProvenance');
end

function test_estimated_native_lwd_reenters_as_fillable(testCase)
   % Builder-estimated longwave is a candidate, never an observation that
   % may train, validate, or pass provenance accounting as native.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   met.Properties.UserData.lwd_estimated = true;
    met.Properties.UserData.lwd_policy = ...
       "lwd estimated from tair and vapor pressure";
    save(filename, 'met');
    recordNativeMetIdentity(root, "tsta", filename);

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="lwd", core_channels="tair", ...
      plan_channels="lwd", interp_channels=string.empty(1, 0), ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      qa_dir=fullfile(root, 'qa-estimated-lwd'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   testCase.verifyFalse(any(result.provenance.lwd == codes.observed));
end

function test_estimated_lwd_policy_requires_consistent_flag(testCase)
   % An empirical-policy artifact cannot become native observation merely
   % because its lwd_estimated flag is stale or absent.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   met.Properties.UserData.lwd_policy = ...
      "empirical incoming longwave reconstruction";
   save(filename, 'met');

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false), ...
      'icemodel:reconstruct:fillPromiceStation:inconsistentLwdProvenance');
end

function test_builder_darkness_zero_reenters_reconstruction(testCase)
   % A physical zero inserted by the native builder is not a source
   % observation; replaying the raw source must route SWD through the formal
   % darkness method and SWU through the final algebraic derivation contract.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   raw_times = met.Properties.RowTimes(1:4:end);
    raw_swd = met.swd(1:4:end);
    raw_swu = met.swu(1:4:end);
    raw_swd(1) = NaN;
    raw_swd(2) = -1;
    raw_swu(1) = NaN;
    source_file = writeRawPromiceAlbedo(root, "tsta", raw_times, ...
       met.albedo(1:4:end), raw_swd, raw_swu);
   met.swd(1:4) = 0;
   met.swd(5:8) = 0;
   metadata = met.Properties.UserData;
   metadata.swd_raw_fallback_count = nnz(isfinite(raw_swd));
   metadata.swd_negative_clamped_count = 1;
   metadata.swd_darkness_zero_filled_count = 1;
    metadata.swu_raw_fallback_count = nnz(isfinite(raw_swu));
    metadata.swu_negative_clamped_count = 0;
    metadata.swu_darkness_zero_filled_count = 1;
    metadata.swd_source_file_observations_present = true;
    metadata.swu_source_file_observations_present = true;
    source_info = dir(source_file);
    metadata.source_file = source_file;
    metadata.source_size_bytes = source_info.bytes;
    metadata.source_sha256 = ...
       icemodel.verification.setup.fileSha256(source_file);
    met.Properties.UserData = metadata;
    save(filename, 'met');
    recordNativeMetIdentity(root, "tsta", filename);

   opts = icemodel.forcing.reconstruct.setopts( ...
       required_channels=["swd", "swu"], core_channels="swd", ...
       plan_channels=["swd", "swu"], interp_channels="swd", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      qa_dir=fullfile(root, 'qa-builder-darkness'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   testCase.verifyEqual(result.filled.swd(1:4), zeros(4, 1));
   testCase.verifyEqual(result.provenance.swd(1:4), ...
      repmat(codes.darkness, 4, 1));
   testCase.verifyEqual(result.provenance.swd(5:8), ...
      repmat(codes.clamped_shortwave, 4, 1));
    testCase.verifyEqual(result.provenance.swd(9:12), ...
       repmat(codes.raw_shortwave, 4, 1));
    testCase.verifyEqual(result.filled.swu(1:4), zeros(4, 1));
    testCase.verifyEqual(result.provenance.swu(1:4), ...
       repmat(codes.derived_shortwave, 4, 1));
    testCase.verifyEqual(result.provenance.swu(5:8), ...
       repmat(codes.raw_shortwave, 4, 1));
    testCase.verifyTrue(any(string(result.audit.channel) == "swd" ...
       & string(result.audit.method) == "darkness_zero"));
    testCase.verifyTrue(any(string(result.audit.channel) == "swu" ...
       & string(result.audit.method) == "derived_shortwave"));
end

function test_raw_fallback_shortwave_gated_by_validity(testCase)
   % A15's hard limits apply at the raw-pyranometer fallback tier too
   % (the fallback prefers raw over nothing, never impossible values —
   % A7): a raw swd far above the TOA ceiling (the GC-Net-legacy
   % evening-shifted class) and a raw swu above its paired swd must not
   % ship as code 13. The samples turn missing and re-enter
   % reconstruction while in-bounds raw fallback keeps its code.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   raw_times = met.Properties.RowTimes(1:4:end);
   raw_swd = met.swd(1:4:end);
   raw_swu = met.swu(1:4:end);
   % Poison one evening posting far above the ceiling and one daylight
   % posting whose upward flux exceeds the paired downwelling flux.
   swd_poison = find(raw_times == datetime(2020, 8, 15, 20, 0, 0, ...
      'TimeZone', 'UTC'));
   swu_poison = find(raw_times == datetime(2020, 8, 16, 14, 0, 0, ...
      'TimeZone', 'UTC'));
   raw_swd(swd_poison) = 950;
   raw_swu(swu_poison) = raw_swd(swu_poison) + 50;
   source_file = writeRawPromiceAlbedo(root, "tsta", raw_times, ...
      met.albedo(1:4:end), raw_swd, raw_swu);
   % Ship the poisoned raw selections over their staged quarter-hour
   % support and stamp the builder's source-selection counts so the
   % replay reproduces the fallback masks.
   swd_support = (4 * (swd_poison - 1) + 1):(4 * swd_poison);
   swu_support = (4 * (swu_poison - 1) + 1):(4 * swu_poison);
   met.swd(swd_support) = raw_swd(swd_poison);
   met.swu(swu_support) = raw_swu(swu_poison);
   metadata = met.Properties.UserData;
   metadata.swd_raw_fallback_count = nnz(isfinite(raw_swd));
   metadata.swd_negative_clamped_count = 0;
   metadata.swd_darkness_zero_filled_count = 0;
   metadata.swu_raw_fallback_count = nnz(isfinite(raw_swu));
   metadata.swu_negative_clamped_count = 0;
   metadata.swu_darkness_zero_filled_count = 0;
   metadata.swd_source_file_observations_present = true;
   metadata.swu_source_file_observations_present = true;
   source_info = dir(source_file);
   metadata.source_file = source_file;
   metadata.source_size_bytes = source_info.bytes;
   metadata.source_sha256 = ...
      icemodel.verification.setup.fileSha256(source_file);
   met.Properties.UserData = metadata;
   save(filename, 'met');
   recordNativeMetIdentity(root, "tsta", filename);

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels=["swd", "swu"], core_channels="swd", ...
      plan_channels=["swd", "swu"], interp_channels="swd", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      qa_dir=fullfile(root, 'qa-raw-gate'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   % Precondition: the poisoned value really exceeds the swd ceiling.
   testCase.verifyFalse(icemodel.forcing.reconstruct.physicalValidity( ...
      "swd", 950, raw_times(swd_poison), latitude=67, longitude=-48.8));
   % Nothing over the poisoned support ships as raw fallback, and
   % whatever ships stays inside the validity ceiling (or is missing).
   testCase.verifyTrue(all(result.provenance.swd(swd_support) ...
      ~= codes.raw_shortwave));
   swd_vals = result.filled.swd(swd_support);
   [~, ~, ceiling_limit] = ...
      icemodel.forcing.reconstruct.physicalValidity("swd", ...
      zeros(numel(swd_support), 1), ...
      result.filled.Properties.RowTimes(swd_support), ...
      latitude=67, longitude=-48.8);
   testCase.verifyTrue(all(isnan(swd_vals) ...
      | swd_vals <= ceiling_limit + 1e-9));
   testCase.verifyTrue(all(result.provenance.swu(swu_support) ...
      ~= codes.raw_shortwave));
   swu_vals = result.filled.swu(swu_support);
   testCase.verifyTrue(all(isnan(swu_vals) ...
      | swu_vals <= result.filled.swd(swu_support) + 1e-9));
   % The gate is surgical: in-bounds raw fallback keeps its code.
   testCase.verifyTrue(any(result.provenance.swd == codes.raw_shortwave));
   testCase.verifyTrue(any(result.provenance.swu == codes.raw_shortwave));
end

function test_missing_swu_derives_after_swd_and_albedo(testCase)
   % The production path preserves native upward shortwave and derives only
   % its missing samples after the operand channels are final.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   gap = (5001:5012).';
   expected = met.albedo(gap) .* met.swd(gap);
   met.swu(gap) = NaN;
    native_swu = met.swu(1);
    save(filename, 'met');
    recordNativeMetIdentity(root, "tsta", filename);

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="swu", core_channels="tair", ...
      plan_channels="swu", interp_channels=string.empty(1, 0), ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      qa_dir=fullfile(root, 'qa-derived-swu'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   testCase.verifyEqual(result.filled.swu(gap), expected, ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(result.provenance.swu(gap), ...
      repmat(codes.derived_shortwave, numel(gap), 1));
   testCase.verifyEqual(result.filled.swu(1), native_swu);
   testCase.verifyEqual(result.provenance.swu(1), codes.observed);
   swu_plan = result.plan.channels( ...
      string({result.plan.channels.channel}) == "swu");
   testCase.verifyEmpty(swu_plan.methods);
end

function test_staged_usr_fixture_ships_swu_product(testCase)
   % READ-BOUNDARY COVERAGE (POLICY A16/D-24): guarded staged artifacts
   % written before the rename still carry the pypromice source name usr.
   % The loader must rename on read so the product ships canonical swu and
   % the source name never survives past the boundary.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   % Recreate the legacy artifact: the staged channel under its source
   % name, with a gap so the boundary-renamed channel also exercises the
   % downstream derivation.
   met = renamevars(met, "swu", "usr");
   gap = (5001:5012).';
   expected = met.albedo(gap) .* met.swd(gap);
   met.usr(gap) = NaN;
   save(filename, 'met');
   recordNativeMetIdentity(root, "tsta", filename);

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="swu", core_channels="tair", ...
      plan_channels="swu", interp_channels=string.empty(1, 0), ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      qa_dir=fullfile(root, 'qa-legacy-usr'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   % Only the canonical name survives the read; usr dies at the boundary.
   names = string(result.filled.Properties.VariableNames);
   testCase.verifyTrue(ismember("swu", names));
   testCase.verifyFalse(ismember("usr", names));
   testCase.verifyFalse(ismember("usr", ...
      string(result.provenance.Properties.VariableNames)));
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   returned = result.filled.swu(gap);
   testCase.verifyEqual(returned, expected, 'AbsTol', 1e-12);
   testCase.verifyEqual(result.provenance.swu(gap), ...
      repmat(codes.derived_shortwave, numel(gap), 1));
   testCase.verifyEqual(result.provenance.swu(1), codes.observed);
end

function test_donor_legacy_albedo_requires_raw_source(testCase)
   % A corrupt donor source link is an integrity failure, not a reason to
   % silently omit that donor and continue with lower tiers.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, true);
   writeFixtureStation(root, "dsta", 0.05, -20, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_dsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   met.Properties.UserData.source_file = fullfile(root, 'missing.nc');
   met.Properties.UserData.albedo_policy = ...
      "albedo from PROMICE L3 albedo; fillPromiceAlbedo(fillwinter=1)";
   save(filename, 'met');

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), donor_sites="dsta", ...
      use_ktransect=false, use_gcnet=false, write=false), ...
      'icemodel:reconstruct:fillPromiceStation:missingAlbedoSource');
end

function test_readiness_counts_absent_core_channel_as_missing(testCase)
   % MATLAB's all(empty,2) is true; readiness must instead give a station
   % with an absent configured core channel zero native-core coverage.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, true);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
    met.psfc = [];
    save(filename, 'met');
    recordNativeMetIdentity(root, "tsta", filename);
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="psfc", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   testCase.verifyEqual(result.readiness.native_core_coverage, ...
      zeros(height(result.readiness), 1));
   testCase.verifyTrue(all(ismember( ...
      string(result.readiness.verdict_icemodel), ...
      ["ready", "not_forcing_ready"])));
   testCase.verifyTrue(all(ismember( ...
      string(result.readiness.verdict_snowmodel), ...
      ["ready", "not_forcing_ready"])));
   testCase.verifyFalse(any(result.readiness.worth_filling));
   testCase.verifyTrue(all(contains(result.readiness.triage_note, ...
      "native core coverage 0% below")));
   testCase.verifyEqual(result.readiness.n_admitted_methods, ...
      zeros(height(result.readiness), 1));
end

function test_readiness_splits_scalar_and_relational_violations(testCase)
   % Verdicts grade completeness plus SCALAR bounds only: the negative
   % precipitation sample flips the verdict, while the relational
   % exceedances (swd above the TOA ceiling, swu above swd) surface in
   % the diagnostic columns and never reach a verdict or its reason
   % (POLICY A15/D-28).
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   source_row = mod((1:height(met)).' - 1, 4) == 0;
   daylight = find(year(met.Properties.RowTimes) == 2020 ...
      & met.swd > 0 & source_row, 2, 'first');
   met.swd(daylight(1)) = 1e6;
   met.swu(daylight(2)) = met.swd(daylight(2)) + 1;
   met.ppt(:) = 0;
    met.ppt(daylight(2)) = -1;
    save(filename, 'met');
    recordNativeMetIdentity(root, "tsta", filename);
    proxy_file = fullfile(root, 'met', 'mar3.11', ...
       'met_tsta_mar3.11_20200101_20211231_15m.mat');
    S = load(proxy_file, 'mar_met');
    mar_met = S.mar_met;
    mar_met.ppt(daylight(2)) = NaN;
    mar_met.rainf(daylight(2)) = NaN;
    mar_met.snowf(daylight(2)) = NaN;
    save(proxy_file, 'mar_met');
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels=["swd", "swu", "ppt"], core_channels="tair", ...
      plan_channels=["swd", "swu"], interp_channels="swd", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   row = result.readiness(result.readiness.year == 2020, :);
   % The A5 IceModel set is fixed: caller planning options cannot add ppt
   % or swu to it, so the precipitation defect affects only snow readiness.
   testCase.verifyEqual(string(row.verdict_icemodel), "ready");
   testCase.verifyEqual(string(row.reason_icemodel), "");
   testCase.verifyFalse(contains(string(row.reason_icemodel), "swd"));
   testCase.verifyFalse(contains(string(row.reason_icemodel), "swu"));
   testCase.verifyEqual(row.worst_residual_invalid, 0);
   testCase.verifyEqual(string(row.verdict_snowmodel), ...
      "not_forcing_ready");
   testCase.verifyTrue(contains(string(row.reason_snowmodel), ...
      "snowfall input"));
   % Relational exceedances land in the diagnostic columns instead.
   testCase.verifyGreaterThan(row.worst_relational_invalid, 0);
   testCase.verifyTrue(contains(string(row.relational_diagnostic), ...
      "swd"));
   testCase.verifyTrue(contains(string(row.relational_diagnostic), ...
      "swu"));
   % The unpoisoned year still carries structural relational exceedances
   % (held hourly postings exceed the instantaneous TOA ceiling near
   % sunset — representation noise, not bad data), yet it grades READY:
   % that is exactly the D-28 point of demoting relational rules to
   % diagnostics.
   clean = result.readiness(result.readiness.year == 2021, :);
   testCase.verifyEqual(string(clean.verdict_icemodel), "ready");
   testCase.verifyGreaterThan(clean.worst_relational_invalid, 0);
end

function test_readiness_rejects_inconsistent_precipitation_phase(testCase)
   % A complete required-channel set is not runnable when its phase split
   % contradicts total precipitation.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   met.ppt(:) = 0;
   met.rainf(:) = 0;
   met.snowf(:) = 0;
   met.ppt(101) = 0.1;
    met.rainf(101) = 0.2;
    save(filename, 'met');
    recordNativeMetIdentity(root, "tsta", filename);
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   row = result.readiness(result.readiness.year == 2020, :);
   testCase.verifyEqual(string(row.verdict_icemodel), "ready");
   testCase.verifyEqual(string(row.verdict_snowmodel), ...
      "not_forcing_ready");
   testCase.verifyTrue(contains(string(row.reason_snowmodel), ...
      "snowfall input"));
   testCase.verifyEqual(result.filled.rainf(101), 0.2, 'AbsTol', 1e-12);
   testCase.verifyTrue(isnan(result.filled.ppt(101)));
   testCase.verifyTrue(isnan(result.filled.snowf(101)));
   testCase.verifyGreaterThan(row.worst_residual_snowmodel, 0);
end

function test_nonready_station_refuses_publication(testCase)
   % A station with staggered core support still returns a diagnostic
   % ledger under write=false, but a canonical product cannot be published.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   writeFixtureMar(root, "tsta");
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   midpoint = floor(height(met) / 2);
   met.tair((midpoint + 1):end) = NaN;
    met.rh(1:midpoint) = NaN;
    save(filename, 'met');
    recordNativeMetIdentity(root, "tsta", filename);
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels=["tair", "rh"], ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);

   out_dir = fullfile(root, 'met', 'promice_filled');
   qa_dir = fullfile(root, 'qa-no-core');
   mkdir(out_dir);
   mkdir(fullfile(qa_dir, 'plans'));
   stale = [ ...
      string(fullfile(out_dir, ...
      'met_tsta_promice_filled_20200101_20211231_15m.mat')); ...
      string(fullfile(qa_dir, 'plans', 'tsta-plan.mat'))];
   for file = stale.'
      fid = fopen(file, 'w');
      fprintf(fid, 'stale');
      fclose(fid);
   end
   invoke = @() icemodel.forcing.reconstruct.fillPromiceStation( ...
      "tsta", met_dir=fullfile(root, 'met', 'promice'), ...
      out_dir=out_dir, qa_dir=qa_dir, ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, opts=opts);
   testCase.verifyError(invoke, ...
      'icemodel:reconstruct:fillPromiceStation:notForcingReady');
   testCase.verifyFalse(any(arrayfun(@isfile, stale)));

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);
   testCase.verifyEqual(result.met_file, "");
   testCase.verifyTrue(all( ...
      result.readiness.verdict_icemodel == "not_forcing_ready"));
   testCase.verifyTrue(all( ...
      result.readiness.verdict_snowmodel == "not_forcing_ready"));
   testCase.verifyFalse(any(result.readiness.worth_filling));
   testCase.verifyTrue(all(contains(result.readiness.triage_note, ...
      "native core coverage 0% below")));
   testCase.verifyEqual(result.readiness.n_admitted_methods, ...
      zeros(height(result.readiness), 1));
   testCase.verifyEqual(result.plan.split.selection_fraction, ...
      opts.selection_fraction);
end

function test_flat_run_screen_excludes_working_native_samples(testCase)
   % A corroborated buried-sensor run is reported and removed from the
   % reconstruction working copy without changing the staged artifact.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   saved = load(filename, 'met');
   met = saved.met;
   times = met.Properties.RowTimes;
   buried = times >= datetime(2020, 3, 25, 'TimeZone', 'UTC') ...
      & times < datetime(2020, 3, 29, 'TimeZone', 'UTC');
   met.tair(buried) = 250 + 0.05 * sin( ...
      2 * pi * hour(times(buried)) / 24);
   met.rh(buried) = 85;
   met.swd(buried) = 0.5;
   sigma = icemodel.physicalConstant('SB');
   met.lwd(buried) = sigma * met.tair(buried) .^ 4 + 1;
   staged_tair = met.tair;
   save(filename, 'met');
   recordNativeMetIdentity(root, "tsta", filename);
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   testCase.verifyEqual(height(result.flat_run_findings), 1);
   testCase.verifyTrue(contains(result.flat_run_findings.channels, ...
      "tair"));
   in_result = result.filled.Properties.RowTimes >= ...
      result.flat_run_findings.start_time ...
      & result.filled.Properties.RowTimes <= ...
      result.flat_run_findings.end_time;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   testCase.verifyFalse(any( ...
      result.provenance.tair(in_result) == codes.observed));
   unchanged = load(filename, 'met');
   testCase.verifyEqual(unchanged.met.tair, staged_tair);
end

function test_gcnet_donor_requires_channel_origin_flags(testCase)
   % Origin masks remain per channel, and an unflagged LRin channel cannot
   % silently treat reconstructed values as native observations.
   filename = writeGcnetDonorFixture(testCase.TestData.root);

   donor = icemodel.forcing.helpers.readGcnetDonor(string(filename));

   testCase.verifyEqual(donor.station, "TEST");
   testCase.verifyTrue(all(donor.observed_mask.tair));
   testCase.verifyFalse(donor.observed_mask.rh(1));
   testCase.verifyFalse(any(donor.observed_mask.lwd));
   bad_file = fullfile(testCase.TestData.root, 'bad_surface.nc');
   fid = fopen(bad_file, 'w');
   fprintf(fid, 'not netcdf');
   fclose(fid);
   testCase.verifyEmpty( ...
      icemodel.forcing.helpers.readGcnetDonor(string(bad_file)));
end

function test_gcnet_donor_real_file_contract(testCase)
   % The REAL Vandecrux surface files carry a drifting fractional-day
   % time coordinate and no location attributes; the reader must return
   % the exact hourly row-index axis (shared gcnetHourlyAxis convention)
   % and catalog coordinates for the station token.
   root = testCase.TestData.root;
   filename = fullfile(root, 'CP1_surface.nc');
   nccreate(filename, 'time', 'Dimensions', {'time', 4});
   % Fractional days with drift: the third step is deliberately short.
   ncwrite(filename, 'time', [36000; 36000.0416667; 36000.083; ...
      36000.1249999]);
   ncwriteatt(filename, 'time', 'units', 'days since 1900-1-1 0:0:0');
   for name = ["Ta_2m", "RH_2m", "WS_10m", "SRin", "LRin"]
      nccreate(filename, name, 'Dimensions', {'time', 4});
      ncwrite(filename, name, [1; 2; 3; 4]);
      nccreate(filename, name + "_origin", 'Dimensions', {'time', 4});
      ncwrite(filename, name + "_origin", zeros(4, 1));
   end

   donor = icemodel.forcing.helpers.readGcnetDonor(string(filename));

   testCase.verifyEqual(donor.station, "CP1");
   t = donor.series.Properties.RowTimes;
   % Row-index hourly axis anchored at the first timestamp, drift gone.
   expected = t(1) + hours(0:3).';
   testCase.verifyEqual(t, expected);
   % Coordinates come from the station catalog (the dataset's Dataverse
   % metadata), not from absent file attributes.
   testCase.verifyEqual(donor.location.lat_wgs84, 69.88, 'AbsTol', 1e-9);
   testCase.verifyEqual(donor.location.lon_wgs84, -46.99, 'AbsTol', 1e-9);
end

function test_boom_channel_materializes_when_native_lacks_it(testCase)
   % POLICY A3: a native product with no boom channel still ships an
   % all-missing boom_height column so the runtime loader's provenance
   % gate passes and every sample lands on the nominal rung.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = removevars(S.met, 'boom_height');
   save(filename, 'met');
   recordNativeMetIdentity(root, "tsta", filename);
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=3);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   testCase.verifyTrue(ismember("boom_height", ...
      string(result.filled.Properties.VariableNames)));
   testCase.verifyTrue(all(isnan(result.filled.boom_height)));
   testCase.verifyTrue(all(result.provenance.boom_height == ...
      codes.missing));
end

function test_last_resort_adopts_bounded_proxy_values(testCase)
   % The last-resort tier adopts aligned, in-bounds proxy values for
   % residual missing required-channel samples in catalog order, stamps
   % the proxy code, appends one audit row per channel and source, and
   % never touches finite samples or out-of-bounds proxy values.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:99)).';
   n = numel(times);
   tair = 260 + zeros(n, 1);
   tair(41:60) = NaN;                            % residual outage
   filled = timetable(times, tair, 'VariableNames', {'tair'});
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   code = repmat(codes.observed, n, 1);
   code(41:60) = codes.missing;
   provenance = timetable(times, code, 'VariableNames', {'tair'});
   audit = table('Size', [0 7], 'VariableTypes', {'cellstr', ...
      'datetime', 'datetime', 'double', 'cellstr', 'cellstr', ...
      'cellstr'}, ...
      'VariableNames', {'channel', 'start_time', 'end_time', ...
      'duration_hours', 'method', 'detail', 'context_id'});
   % The engine's audit datetimes are zoned; the empty schema must match
   % or appended rows cannot concatenate.
   audit.start_time.TimeZone = 'UTC';
   audit.end_time.TimeZone = 'UTC';

   % The preferred proxy covers most of the outage, with two poisoned
   % samples. A second proxy covers only those two; it must not be mixed
   % into the same outage.
   p_tair = 262 + zeros(n, 1);
   p_tair(45) = 500;
   p_tair(46) = NaN;
   proxy_series = timetable(times, p_tair, 'VariableNames', {'tair'});
   proxies = struct('series', proxy_series, 'name', "mar", ...
      'code_name', "mar");
   second = nan(n, 1);
   second([45 46]) = 261;
    proxies(2) = struct('series', timetable(times, second, ...
       'VariableNames', {'tair'}), 'name', "merra2", ...
       'code_name', "merra2");
    opts = icemodel.forcing.reconstruct.setopts();
    native = timetable(times, 260 + zeros(n, 1), ...
       'VariableNames', {'tair'});
    calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
       times, native.tair, p_tair, "tair", fit_years=2020);
    calibration.source = "mar";
    channel_plan = struct('channel', "tair", ...
       'proxy_calibrations', struct('source', "mar", ...
       'parameters', calibration));
    plan = struct('channels', channel_plan);

    [returned, prov, audit_out, denials] = ...
       icemodel.forcing.reconstruct.lastResortProxies(filled, ...
       provenance, audit, proxies, codes, opts, native=native, plan=plan);

   % Adopted samples carry the proxy value and code; poisoned samples
   % stay missing; finite samples are untouched.
   expected = [41:44, 47:60].';
    testCase.verifyTrue(all(isfinite(returned.tair(expected))));
    testCase.verifyEqual(returned.tair([41; 60]), [260; 260], ...
       'AbsTol', 1e-12);
     testCase.verifyEqual(returned.tair(47:54), ...
        260 + zeros(8, 1), 'AbsTol', 1e-12);
   testCase.verifyEqual(prov.tair(expected), ...
      repmat(codes.mar, numel(expected), 1));
   testCase.verifyTrue(all(isnan(returned.tair([45 46]))));
   testCase.verifyEqual(prov.tair([45; 46]), ...
      repmat(codes.missing, 2, 1));
   testCase.verifyEqual(returned.tair(1:40), tair(1:40));
   % Two rows record the two actually contiguous adopted segments.
   row = audit_out(strcmp(audit_out.method, ...
      'proxy:mar:last_resort'), :);
     testCase.verifyEqual(height(row), 2);
     testCase.verifyEqual(row.duration_hours, [4; 14]);
     testCase.verifyTrue(all(contains(string(row.detail), ...
        "whole-outage source")));
   testCase.verifyFalse(any(strcmp(audit_out.method, ...
      'proxy:merra2:last_resort')));
   % Final-tier denial notes: the poisoned samples the chosen source
   % cannot supply carry the whole-outage refusal so the driver's audit
   % reconciliation ships the actual last-tier cause.
   testCase.verifyEqual(denials.tair([45; 46]), repmat( ...
      "last-resort source lacks valid values (whole-outage policy)", ...
      2, 1));
   testCase.verifyEqual(unique(denials.tair(expected)), "");
end

function test_last_resort_rejects_post_blend_outside_validity(testCase)
   % A twilight proxy zero is valid before tapering; blending it toward
   % the prior daylight anchor pushes it above the local ceiling. The
   % post-blend candidate must remain missing because A15/B6 applies the
   % candidate bound after every estimate adjustment.
   times = (datetime(2020, 3, 20, 'TimeZone', 'UTC') + hours(0:47)).';
   toa = icemodel.forcing.reconstruct.toaIrradiance(times, 0, 0);
   gap = find(toa == 0 & [false; toa(1:end - 1) > 0], 1);
   swd = zeros(numel(times), 1);
   swd(gap - 1) = 0.9 * toa(gap - 1);
   swd(gap) = NaN;
   filled = timetable(times, swd, 'VariableNames', {'swd'});
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   code = repmat(codes.observed, numel(times), 1);
   code(gap) = codes.missing;
   provenance = timetable(times, code, 'VariableNames', {'swd'});
   audit = table('Size', [0 7], 'VariableTypes', {'cellstr', ...
      'datetime', 'datetime', 'double', 'cellstr', 'cellstr', 'cellstr'}, ...
      'VariableNames', {'channel', 'start_time', 'end_time', ...
      'duration_hours', 'method', 'detail', 'context_id'});
   audit.start_time.TimeZone = 'UTC';
   audit.end_time.TimeZone = 'UTC';
   proxy_swd = 0.5 * toa;
   proxy_swd(gap) = 0;
   proxy = struct('series', timetable(times, proxy_swd, ...
       'VariableNames', {'swd'}), 'name', "mar", 'code_name', "mar");
    calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
       times, swd, proxy.series.swd, "swd", fit_years=2020, ...
       target_toa=toa);
    calibration.source = "mar";
    channel_plan = struct('channel', "swd", ...
       'proxy_calibrations', struct('source', "mar", ...
       'parameters', calibration));
    plan = struct('channels', channel_plan);

    [returned, returned_provenance, returned_audit, denials] = ...
       icemodel.forcing.reconstruct.lastResortProxies( ...
       filled, provenance, audit, proxy, codes, ...
       icemodel.forcing.reconstruct.setopts(required_channels="swd"), ...
       latitude=0, longitude=0, plan=plan);

   testCase.verifyTrue(isnan(returned.swd(gap)));
   testCase.verifyEqual(returned_provenance.swd(gap), codes.missing);
   testCase.verifyEqual(height(returned_audit), 0);
   testCase.verifyEqual(denials.swd(gap), "post-blend seam rejection");
end

function test_last_resort_uses_swd_elevation_bins(testCase)
   % D-28: last resort must apply the planner's elevation-binned SWD
   % correction before candidate validation. The legacy seasonal scalar
   % leaves this civil-twilight run above its ceiling; the binned ratio
   % makes it admissible.
   times = (datetime(2020, 3, 20, 'TimeZone', 'UTC') ...
      + minutes(0:15:(24 * 60 - 15))).';
   raw_proxy = 100 + zeros(numel(times), 1);
   corrected_proxy = 20 + zeros(numel(times), 1);
   interval = median(diff(times));
   elevation = ...
      icemodel.forcing.helpers.intervalMaximumSolarElevation( ...
      times, 0, 0, interval);
   toa = icemodel.forcing.reconstruct.toaIrradiance(times, 0, 0);
   raw_valid = icemodel.forcing.reconstruct.physicalValidity( ...
      "swd", raw_proxy, times, latitude=0, longitude=0, ...
      toa=toa, elevation=elevation);
   corrected_valid = icemodel.forcing.reconstruct.physicalValidity( ...
      "swd", corrected_proxy, times, latitude=0, longitude=0, ...
      toa=toa, elevation=elevation);
   candidate = ~raw_valid & corrected_valid;
   pair = find(candidate(1:end - 1) & candidate(2:end), 1);
   testCase.assertNotEmpty(pair);
   gap = [pair; pair + 1];

   swd = corrected_proxy;
   swd(gap) = NaN;
   filled = timetable(times, swd, 'VariableNames', {'swd'});
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   code = repmat(codes.observed, numel(times), 1);
   code(gap) = codes.missing;
   provenance = timetable(times, code, 'VariableNames', {'swd'});
   audit = table('Size', [0 7], 'VariableTypes', {'cellstr', ...
      'datetime', 'datetime', 'double', 'cellstr', 'cellstr', 'cellstr'}, ...
      'VariableNames', {'channel', 'start_time', 'end_time', ...
      'duration_hours', 'method', 'detail', 'context_id'});
   audit.start_time.TimeZone = 'UTC';
   audit.end_time.TimeZone = 'UTC';
   proxy = struct('series', timetable(times, raw_proxy, ...
      'VariableNames', {'swd'}), 'name', "mar", 'code_name', "mar");

   bands = icemodel.forcing.reconstruct.solarElevationBands();
   ratios = 0.2 * ones(1, numel(bands.calibration_bin_edges_deg) - 1);
   calibration = struct('channel', "swd", 'mode', "multiplicative", ...
      'corrections', struct('DJF', 2, 'MAM', 2, 'JJA', 2, 'SON', 2), ...
      'n_overlap', 100, ...
      'bin_edges_deg', bands.calibration_bin_edges_deg, ...
      'binned_corrections', struct('DJF', ratios, 'MAM', ratios, ...
      'JJA', ratios, 'SON', ratios));
   plan = struct('channels', struct('channel', "swd", ...
      'proxy_calibrations', struct('source', "mar", ...
      'parameters', calibration)));

   [returned, returned_provenance] = ...
      icemodel.forcing.reconstruct.lastResortProxies( ...
      filled, provenance, audit, proxy, codes, ...
      icemodel.forcing.reconstruct.setopts(required_channels="swd"), ...
      latitude=0, longitude=0, plan=plan);

   testCase.verifyEqual(returned.swd(gap), corrected_proxy(gap), ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(returned_provenance.swd(gap), ...
      repmat(codes.mar, numel(gap), 1));
end

function test_last_resort_admits_zero_overlap_identity_low_confidence(testCase)
   % A source with no target overlap carries no calibration evidence, but
   % refusing it left megasample in-bounds gaps unfilled; POLICY A11/D-25
   % admits the identity values as last resort with a low-confidence
   % audit stamp instead.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:9)).';
   tair = 260 + zeros(numel(times), 1);
   filled = timetable(times, tair, 'VariableNames', {'tair'});
   filled.tair(5) = NaN;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   provenance = timetable(times, repmat(codes.observed, numel(times), 1), ...
      'VariableNames', {'tair'});
   provenance.tair(5) = codes.missing;
   audit = table('Size', [0 7], 'VariableTypes', {'cellstr', ...
      'datetime', 'datetime', 'double', 'cellstr', 'cellstr', 'cellstr'}, ...
      'VariableNames', {'channel', 'start_time', 'end_time', ...
      'duration_hours', 'method', 'detail', 'context_id'});
   audit.start_time.TimeZone = 'UTC';
   audit.end_time.TimeZone = 'UTC';
   proxy = struct('series', timetable(times, 265 + zeros(numel(times), 1), ...
      'VariableNames', {'tair'}), 'name', "mar", 'code_name', "mar");
   calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
      times, nan(numel(times), 1), proxy.series.tair, "tair", fit_years=2020);
   calibration.source = "mar";
   plan = struct('channels', struct('channel', "tair", ...
      'proxy_calibrations', struct('source', "mar", ...
      'parameters', calibration)));

   [returned, returned_provenance, returned_audit] = ...
      icemodel.forcing.reconstruct.lastResortProxies( ...
      filled, provenance, audit, proxy, codes, ...
      icemodel.forcing.reconstruct.setopts(required_channels="tair"), ...
      plan=plan);

   testCase.verifyTrue(isfinite(returned.tair(5)));
   testCase.verifyEqual(returned_provenance.tair(5), codes.mar);
   testCase.verifyEqual(height(returned_audit), 1);
   testCase.verifyTrue(contains(string(returned_audit.detail(1)), ...
      "low confidence"));
end

function test_last_resort_rejects_taper_beyond_lwd_bound(testCase)
   % The cen Apr 2016 class: a genuine native observation ABOVE the lwd
   % ceiling anchors the seam, and the taper toward it lifts otherwise
   % valid proxy fills out of bounds. A15/B6 requires those post-blend
   % samples to remain missing.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:99)).';
   n = numel(times);
   lwd_bounds = icemodel.forcing.reconstruct.physicalBounds("lwd");
   lwd = 300 + zeros(n, 1);
   lwd(40) = lwd_bounds(2) + 60;   % out-of-bounds observed seam anchor
   lwd(41:52) = NaN;               % residual outage
   filled = timetable(times, lwd, 'VariableNames', {'lwd'});
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   code = repmat(codes.observed, n, 1);
   code(41:52) = codes.missing;
   provenance = timetable(times, code, 'VariableNames', {'lwd'});
   audit = table('Size', [0 7], 'VariableTypes', {'cellstr', ...
      'datetime', 'datetime', 'double', 'cellstr', 'cellstr', 'cellstr'}, ...
      'VariableNames', {'channel', 'start_time', 'end_time', ...
      'duration_hours', 'method', 'detail', 'context_id'});
   audit.start_time.TimeZone = 'UTC';
   audit.end_time.TimeZone = 'UTC';
   proxy = struct('series', timetable(times, 300 + zeros(n, 1), ...
      'VariableNames', {'lwd'}), 'name', "mar", 'code_name', "mar");
   native = timetable(times, 300 + zeros(n, 1), 'VariableNames', {'lwd'});

   [returned, returned_provenance, returned_audit, denials] = ...
      icemodel.forcing.reconstruct.lastResortProxies( ...
      filled, provenance, audit, proxy, codes, ...
      icemodel.forcing.reconstruct.setopts(required_channels="lwd"), ...
      native=native);

   % Out-of-bounds taper samples remain missing; admissible later samples
   % may still fill and must stay inside the ceiling.
   gap = (41:52).';
   rejected = isnan(returned.lwd(gap));
   testCase.verifyTrue(any(rejected));
   testCase.verifyTrue(any(~rejected));
   testCase.verifyTrue(all(returned.lwd(gap(~rejected)) ...
      <= lwd_bounds(2)));
   testCase.verifyEqual(returned_provenance.lwd(gap(rejected)), ...
      repmat(codes.missing, nnz(rejected), 1));
   testCase.verifyEqual(denials.lwd(gap(rejected)), ...
      repmat("post-blend seam rejection", nnz(rejected), 1));
   testCase.verifyFalse(any(contains(string(returned_audit.detail), ...
      "seam capped at validity bound")));
end

function test_last_resort_clamps_calibrated_rh_into_bounds(testCase)
   % D-27: an overlap correction can push humidity above 100 percent;
   % that excess is calibration arithmetic, not physics, so the corrected
   % candidate clamps into the registry bounds and the audit notes the
   % clamp instead of refusing megasample gaps.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:199)).';
   n = numel(times);
   rh = 95 + zeros(n, 1);
   rh(150:174) = NaN;              % residual outage
   filled = timetable(times, rh, 'VariableNames', {'rh'});
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   code = repmat(codes.observed, n, 1);
   code(150:174) = codes.missing;
   provenance = timetable(times, code, 'VariableNames', {'rh'});
   audit = table('Size', [0 7], 'VariableTypes', {'cellstr', ...
      'datetime', 'datetime', 'double', 'cellstr', 'cellstr', 'cellstr'}, ...
      'VariableNames', {'channel', 'start_time', 'end_time', ...
      'duration_hours', 'method', 'detail', 'context_id'});
   audit.start_time.TimeZone = 'UTC';
   audit.end_time.TimeZone = 'UTC';
   % Overlap bias is +10 (observed 95 vs proxy 85), so the outage proxy
   % values of 95 correct to 105 and must clamp to the 100 ceiling.
   proxy_rh = 85 + zeros(n, 1);
   proxy_rh(150:174) = 95;
   proxy = struct('series', timetable(times, proxy_rh, ...
      'VariableNames', {'rh'}), 'name', "mar", 'code_name', "mar");
   native = timetable(times, 95 + zeros(n, 1), 'VariableNames', {'rh'});
   calibration = icemodel.forcing.reconstruct.fitProxyCalibration( ...
      times, filled.rh, proxy_rh, "rh", fit_years=2020);
   calibration.source = "mar";
   plan = struct('channels', struct('channel', "rh", ...
      'proxy_calibrations', struct('source', "mar", ...
      'parameters', calibration)));

   [returned, returned_provenance, returned_audit] = ...
      icemodel.forcing.reconstruct.lastResortProxies( ...
      filled, provenance, audit, proxy, codes, ...
      icemodel.forcing.reconstruct.setopts(required_channels="rh"), ...
      native=native, plan=plan);

   % Mid-outage samples beyond the seam-taper reach hold the clamped
   % ceiling exactly; nothing anywhere exceeds the registry bounds.
   rh_bounds = icemodel.forcing.reconstruct.physicalBounds("rh");
   gap = (150:174).';
   testCase.verifyTrue(all(isfinite(returned.rh(gap))));
   testCase.verifyTrue(all(returned.rh(gap) <= rh_bounds(2) + 1e-9));
   testCase.verifyEqual(returned.rh(160:164), 100 + zeros(5, 1), ...
      'AbsTol', 1e-9);
   testCase.verifyEqual(returned_provenance.rh(gap), ...
      repmat(codes.mar, numel(gap), 1));
   testCase.verifyTrue(contains(string(returned_audit.detail(1)), ...
      "calibration clamped to bounds"));
end

function test_last_resort_denies_span_without_usable_source(testCase)
   % A span no catalog source covers records the no-usable-source denial
   % note for audit reconciliation and leaves the samples missing.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:9)).';
   n = numel(times);
   tair = 260 + zeros(n, 1);
   tair(5) = NaN;
   filled = timetable(times, tair, 'VariableNames', {'tair'});
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   code = repmat(codes.observed, n, 1);
   code(5) = codes.missing;
   provenance = timetable(times, code, 'VariableNames', {'tair'});
   audit = table('Size', [0 7], 'VariableTypes', {'cellstr', ...
      'datetime', 'datetime', 'double', 'cellstr', 'cellstr', 'cellstr'}, ...
      'VariableNames', {'channel', 'start_time', 'end_time', ...
      'duration_hours', 'method', 'detail', 'context_id'});
   audit.start_time.TimeZone = 'UTC';
   audit.end_time.TimeZone = 'UTC';
   % The catalog source carries a different channel entirely, so the
   % outage has no usable values from any source.
   proxy = struct('series', timetable(times, 80 + zeros(n, 1), ...
      'VariableNames', {'rh'}), 'name', "mar", 'code_name', "mar");

   [returned, returned_provenance, returned_audit, denials] = ...
      icemodel.forcing.reconstruct.lastResortProxies( ...
      filled, provenance, audit, proxy, codes, ...
      icemodel.forcing.reconstruct.setopts(required_channels="tair"));

   testCase.verifyTrue(isnan(returned.tair(5)));
   testCase.verifyEqual(returned_provenance.tair(5), codes.missing);
   testCase.verifyEmpty(returned_audit);
   expected = "no usable last-resort source";
   testCase.verifyEqual(denials.tair(5), expected);
   testCase.verifyEqual(unique(denials.tair([1:4, 6:10].')), "");
end

function test_last_resort_rejects_lone_island_adoption(testCase)
   % One admitted posting flanked by unfillable samples renders as an
   % isolated spike pair after the 15-minute expansion (the KAN_M
   % 2011-03-22 class); it is a rendering artifact, not information, so
   % the single-posting run is refused with a denial note while a
   % two-posting run in the same outage survives.
   times = (datetime(2020, 1, 1, 'TimeZone', 'UTC') + hours(0:99)).';
   n = numel(times);
   tair = 260 + zeros(n, 1);
   tair(41:60) = NaN;                            % residual outage
   filled = timetable(times, tair, 'VariableNames', {'tair'});
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   code = repmat(codes.observed, n, 1);
   code(41:60) = codes.missing;
   provenance = timetable(times, code, 'VariableNames', {'tair'});
   audit = table('Size', [0 7], 'VariableTypes', {'cellstr', ...
      'datetime', 'datetime', 'double', 'cellstr', 'cellstr', 'cellstr'}, ...
      'VariableNames', {'channel', 'start_time', 'end_time', ...
      'duration_hours', 'method', 'detail', 'context_id'});
   audit.start_time.TimeZone = 'UTC';
   audit.end_time.TimeZone = 'UTC';
   % The source covers one isolated mid-outage posting plus a separate
   % two-posting run; the rest of the outage stays uncoverable.
   p_tair = nan(n, 1);
   p_tair(50) = 262;
   p_tair([55; 56]) = 262;
   proxy = struct('series', timetable(times, p_tair, ...
      'VariableNames', {'tair'}), 'name', "mar", 'code_name', "mar");
   native = timetable(times, 260 + zeros(n, 1), ...
      'VariableNames', {'tair'});

   [returned, returned_provenance, ~, denials] = ...
      icemodel.forcing.reconstruct.lastResortProxies( ...
      filled, provenance, audit, proxy, codes, ...
      icemodel.forcing.reconstruct.setopts(required_channels="tair"), ...
      native=native);

   % The lone island is refused with its own denial; the two-posting run
   % adopts normally.
   testCase.verifyTrue(isnan(returned.tair(50)));
   testCase.verifyEqual(returned_provenance.tair(50), codes.missing);
   testCase.verifyEqual(denials.tair(50), ...
      "lone single-posting adoption between missing neighbors");
   testCase.verifyTrue(all(isfinite(returned.tair([55; 56]))));
   testCase.verifyEqual(returned_provenance.tair([55; 56]), ...
      repmat(codes.mar, 2, 1));
end

function test_native_loader_uses_saved_time_coverage(testCase)
   % The selected native file is the widest timetable, not the largest MAT.
   root = testCase.TestData.root;
   met_dir = fullfile(root, 'met', 'promice');
   times = (datetime(2019, 1, 1, 'TimeZone', 'UTC'): ...
      hours(1):datetime(2020, 12, 31, 23, 0, 0, ...
      'TimeZone', 'UTC')).';
   met = timetable(times, 260 + zeros(numel(times), 1), ...
      'VariableNames', {'tair'});
   met.Properties.UserData = struct('site', "wide", ...
      'lat', 67, 'lon', -48.8, 'elev', 1200);
   wide_file = fullfile(met_dir, ...
      'met_wide_promice_20190101_20201231_15m.mat');
   save(wide_file, 'met');
   recordNativeMetIdentity(root, "wide", wide_file);

   % Incompressible padding makes this shorter record the larger file.
   short_met = met(1:10, :);
   stream = RandStream('mt19937ar', 'Seed', 9);
   padding = rand(stream, 2e5, 1);
   short_file = fullfile(met_dir, ...
      'met_wide_promice_20200101_20200110_15m.mat');
   save(short_file, 'short_met', 'padding');
   recordNativeMetIdentity(root, "wide", short_file);
   writeFixtureMarWindow(root, "wide", ...
      datetime(2019, 1, 1, 'TimeZone', 'UTC'), ...
      datetime(2020, 12, 31, 'TimeZone', 'UTC'));

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, seed=2);
   result = icemodel.forcing.reconstruct.fillPromiceStation("wide", ...
      met_dir=met_dir, out_dir=fullfile(root, 'out'), ...
      qa_dir=fullfile(root, 'qa-wide'), donor_sites=string.empty(1, 0), ...
      use_ktransect=false, use_gcnet=false, write=false, opts=opts);
   testCase.verifyEqual(result.filled.Properties.RowTimes, times);
   testCase.verifyEqual(result.readiness.year, [2019; 2020]);
end

function test_fill_station_rejects_unverified_quarter_hour_axis(testCase)
   % An unstamped quarter-hour artifact cannot redefine native PROMICE cadence.
   root = testCase.TestData.root;
   writeFixtureStation(root, "nostamp", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_nostamp_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   metadata = met.Properties.UserData;
   metadata = rmfield(metadata, {'met_resample_policy', ...
      'met_resample_source_cadence_seconds'});
   met.Properties.UserData = metadata;
   save(filename, 'met');
   recordNativeMetIdentity(root, "nostamp", filename);

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("nostamp", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false), ...
      'icemodel:reconstruct:fillPromiceStation:unverifiedNativeCadence');
end

function test_fill_station_rejects_malformed_cadence_stamp(testCase)
   % Array-valued provenance must fail with the cadence contract error.
   root = testCase.TestData.root;
   sites = ["badpolicy", "badcadence"];

   % Exercise each scalar field guard without accepting ambiguous metadata.
   for site = sites
      writeFixtureStation(root, site, 0, 0, false);
      filename = fullfile(root, 'met', 'promice', sprintf( ...
         'met_%s_promice_20200101_20211231_15m.mat', site));
      S = load(filename, 'met');
      met = S.met;
      metadata = met.Properties.UserData;
      if site == "badpolicy"
         metadata.met_resample_policy = [ ...
            "interval_start_zero_order_hold", "ambiguous"];
      else
         metadata.met_resample_source_cadence_seconds = [3600, 3600];
      end
      met.Properties.UserData = metadata;
      save(filename, 'met');
      recordNativeMetIdentity(root, site, filename);

      testCase.verifyError(@() ...
         icemodel.forcing.reconstruct.fillPromiceStation(site, ...
         met_dir=fullfile(root, 'met', 'promice'), ...
         donor_sites=string.empty(1, 0), use_ktransect=false, ...
         use_gcnet=false, write=false), ...
         'icemodel:reconstruct:fillPromiceStation:unverifiedNativeCadence');
   end
end

function test_fill_station_reconstructs_hourly_then_restores_support(testCase)
   % Guarded PROMICE input represents hourly postings as four held rows.
   % The engine fills one hourly value, then restores the quarter-hour
   % support: the FILLED posting disaggregates mean-preservingly over its
   % four quarter-hours while OBSERVED postings remain exact held copies
   % and provenance codes repeat (POLICY D-30/A7).
   root = testCase.TestData.root;
   writeFixtureStation(root, "hour", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_hour_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
    hourly = S.met(1:4:end, :);
    hourly.Properties.DimensionNames{1} = 'Time';
    hourly.Properties.UserData = rmfield(hourly.Properties.UserData, ...
       {'met_resample_policy', 'met_resample_source_cadence_seconds'});
    gap_time = hourly.Properties.RowTimes(100);
   hourly.tair(100) = NaN;
    met = icemodel.forcing.helpers.resampleMetTimestep(hourly, "15m");
    met.ppt(:) = 0;
    met.rainf(:) = 0;
    met.snowf(:) = 0;
    save(filename, 'met');
    recordNativeMetIdentity(root, "hour", filename);
    writeFixtureMar(root, "hour");

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=4);
   result = icemodel.forcing.reconstruct.fillPromiceStation("hour", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      out_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa-hour'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, opts=opts);

   support = result.filled.Properties.RowTimes >= gap_time ...
      & result.filled.Properties.RowTimes < gap_time + hours(1);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   testCase.verifyEqual(height(result.filled), height(met));
   testCase.verifyTrue(all(isfinite(result.filled.tair(support))));
   % D-30: the filled posting's quarter-hours are a smooth curve (not
   % four held copies) whose mean equals the hourly reconstruction — the
   % linear midpoint of the two observed neighbor postings — exactly.
   quarter = result.filled.tair(support);
   testCase.verifyGreaterThan(numel(unique(quarter)), 1);
   before_support = result.filled.Properties.RowTimes ...
      >= gap_time - hours(1) ...
      & result.filled.Properties.RowTimes < gap_time;
   after_support = result.filled.Properties.RowTimes ...
      >= gap_time + hours(1) ...
      & result.filled.Properties.RowTimes < gap_time + hours(2);
   % Observed postings stay exact held copies (A7): four equal samples.
   testCase.verifyEqual(numel(unique( ...
      result.filled.tair(before_support))), 1);
   testCase.verifyEqual(numel(unique( ...
      result.filled.tair(after_support))), 1);
   expected_posting = (unique(result.filled.tair(before_support)) ...
      + unique(result.filled.tair(after_support))) / 2;
   testCase.verifyEqual(mean(quarter), expected_posting, 'AbsTol', 1e-9);
   % Observed samples are untouched by the disaggregation, and no sample
   % anywhere leaves the scalar registry bounds.
   observed = result.provenance.tair == codes.observed;
   testCase.verifyEqual(result.filled.tair(observed), met.tair(observed));
   tair_bounds = icemodel.forcing.reconstruct.physicalBounds("tair");
   testCase.verifyTrue(all(result.filled.tair >= tair_bounds(1) ...
      & result.filled.tair <= tair_bounds(2)));
   testCase.verifyTrue(all(result.provenance.tair(support) ...
      == codes.bounded_interp));
   row = result.audit(string(result.audit.method) == "bounded_interp" ...
      & string(result.audit.channel) == "tair", :);
   testCase.verifyEqual(row.duration_hours, 1);

   % The report must retain the hourly audit's one-hour end-exclusive span
   % after the delivered product has returned to quarter-hour support.
   report = icemodel.verification.report.buildGapFillReport( ...
      sites="hour", render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa-hour'), ...
      fig_dir=fullfile(root, 'figures-hour'), ...
      report_dir=fullfile(root, 'report-hour'));
   report_row = report.figure_ledger( ...
      report.figure_ledger.channel == "tair" ...
      & report.figure_ledger.method == "bounded_interp", :);
   testCase.verifyEqual(hours(report_row.gap_end - report_row.gap_start), ...
      report_row.duration_hours, 'AbsTol', 1e-12);
end

function test_fill_station_retries_swd_on_delivered_axis(testCase)
   % An hourly low-sun posting can fail the coarse source-axis seam test
   % even though its four delivered quarters admit a smooth D-32 bridge.
   % The final delivered-axis retry must close that residual without
   % changing either observed anchor or leaving a stale unfilled audit row.
   root = testCase.TestData.root;
   writeFixtureStation(root, "hour", 0, 0, false);
   filename = fullfile(root, 'met', 'promice', ...
      'met_hour_promice_20200101_20211231_15m.mat');
   S = load(filename, 'met');
   met = S.met;
   gap_time = datetime(2020, 10, 12, 20, 0, 0, 'TimeZone', 'UTC');
   gap = met.Properties.RowTimes >= gap_time ...
      & met.Properties.RowTimes < gap_time + hours(1);
   before = met.Properties.RowTimes >= gap_time - hours(1) ...
      & met.Properties.RowTimes < gap_time;
   after = met.Properties.RowTimes >= gap_time + hours(1) ...
      & met.Properties.RowTimes < gap_time + hours(2);
   met.swd(before) = 70;
   met.swd(gap) = NaN;
   met.swd(after) = 0;
   met.swu(gap) = NaN;
   save(filename, 'met');
   recordNativeMetIdentity(root, "hour", filename);
   writeFixtureMar(root, "hour");

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="swd", core_channels="swd", ...
      plan_channels="swd", interp_channels="swd", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=4);
   result = icemodel.forcing.reconstruct.fillPromiceStation("hour", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   testCase.verifyTrue(all(isfinite(result.filled.swd(gap))));
   testCase.verifyEqual(result.provenance.swd(gap), ...
      repmat(codes.bounded_interp, nnz(gap), 1));
   testCase.verifyEqual(result.filled.swu(gap), ...
      result.filled.albedo(gap) .* result.filled.swd(gap), ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(result.provenance.swu(gap), ...
      repmat(codes.derived_shortwave, nnz(gap), 1));
   testCase.verifyEqual(result.filled.swd(before), met.swd(before));
   testCase.verifyEqual(result.filled.swd(after), met.swd(after));
   delivered = string(result.audit.channel) == "swd" ...
      & contains(string(result.audit.detail), "delivered-axis");
   testCase.verifyTrue(any(delivered));
   testCase.verifyFalse(any(string(result.audit.channel) == "swd" ...
      & string(result.audit.method) == "unfilled"));
   testCase.verifyFalse(any(string(result.audit.channel) == "swu" ...
      & string(result.audit.method) == "unfilled"));
end

function test_unfilled_audit_reports_final_tier_denial(testCase)
   % Residual unfilled audit rows must state the LAST tier's actual
   % refusal, not the stale provisional reason recorded before the
   % last-resort tier ran.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, true);
   writeFixtureMar(root, "tsta");
   writeFixtureMerraPrecipPatch(root, "tsta");
   % Strip MAR's humidity so the seeded rh outages meet a source with no
   % usable values: the tier must deny, and the denial must ship.
   mar_file = fullfile(root, 'met', 'mar3.11', ...
      'met_tsta_mar3.11_20200101_20211231_15m.mat');
   S = load(mar_file, 'mar_met');
   mar_met = S.mar_met;
   mar_met.rh(:) = NaN;
   save(mar_file, 'mar_met');
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels=["tair", "rh"], core_channels="tair", ...
      plan_channels=["tair", "rh"], interp_channels="tair", ...
      min_overlap_hours=1e9, climatology_min_support=1e9, ...
      plan_n_gaps=1, seed=5);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   % The rh outages remain, and their audit rows now carry the tier's
   % actual denial cause appended to the provisional reason. The rh-only
   % short outage meets no usable source at all, while the joint long
   % outage has a chosen source that cannot supply rh — both causes ship.
   unfilled = result.audit(string(result.audit.method) == "unfilled" ...
      & string(result.audit.channel) == "rh", :);
   testCase.verifyGreaterThan(height(unfilled), 0);
   testCase.verifyTrue(all(contains(string(unfilled.detail), ...
      "final tier:")));
   testCase.verifyTrue(any(contains(string(unfilled.detail), ...
      "last-resort source lacks valid values")));
   testCase.verifyTrue(any(contains(string(unfilled.detail), ...
      "no usable last-resort source")));
   % MAR fills the tair outages, so reconciliation drops their rows.
   testCase.verifyFalse(any(string(result.audit.method) == "unfilled" ...
      & string(result.audit.channel) == "tair"));
end

function test_post_final_interpolation_closes_proxy_sliver(testCase)
   % A long native outage is too large for tier 1. Last resort fills all
   % but two hourly postings, exposing a short interior sliver that the
   % D-32 second pass can now bridge without changing source values.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, true);
   writeFixtureMar(root, "tsta");
   writeFixtureMerraPrecipPatch(root, "tsta");

   native_file = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(native_file, 'met');
   source_times = S.met.Properties.RowTimes(1:4:end);
   source_missing = isnan(S.met.tair(1:4:end));
   edges = diff([false; source_missing; false]);
   starts = find(edges == 1);
   stops = find(edges == -1) - 1;
   [~, longest] = max(stops - starts + 1);
   first = starts(longest) + 2;
   sliver_times = source_times(first:first + 1);

   % Remove exactly the two proxy postings from the otherwise complete
   % source. The rest of the original outage remains a MAR last-resort
   % fill and becomes the second pass's finite anchor context.
   mar_file = fullfile(root, 'met', 'mar3.11', ...
      'met_tsta_mar3.11_20200101_20211231_15m.mat');
   S = load(mar_file, 'mar_met');
   mar_met = S.mar_met;
   poison = mar_met.Properties.RowTimes >= sliver_times(1) ...
      & mar_met.Properties.RowTimes < sliver_times(end) + hours(1);
   mar_met.tair(poison) = NaN;
   save(mar_file, 'mar_met');

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      min_overlap_hours=1e9, climatology_min_support=1e9, ...
      plan_n_gaps=1, seed=5);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   delivered = result.filled.Properties.RowTimes >= sliver_times(1) ...
      & result.filled.Properties.RowTimes < sliver_times(end) + hours(1);
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   testCase.verifyTrue(all(isfinite(result.filled.tair(delivered))));
   testCase.verifyEqual(result.provenance.tair(delivered), ...
      repmat(codes.bounded_interp, nnz(delivered), 1));
   post_final = string(result.audit.method) == "bounded_interp" ...
      & contains(string(result.audit.detail), "post-final residual");
   testCase.verifyTrue(any(post_final));
   testCase.verifyFalse(any(string(result.audit.method) == "unfilled" ...
      & string(result.audit.channel) == "tair"));
end

function test_merra_precip_adoption_derives_rainf_complement(testCase)
   % Staged MERRA-2 ships ppt and snowf but no rainf channel; adoption
   % must derive rain as the exact nonnegative complement ppt - snowf so
   % the shipped product carries a full source split (POLICY D-31).
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   % A few finite native precip observations anchor the seam step scale
   % like the canonical-artifacts fixture does.
   native_file = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   native = load(native_file, 'met');
   met = native.met;
   met.ppt([49; 53; 57]) = [0.2; 0.3; 0.3];
   save(native_file, 'met');
   recordNativeMetIdentity(root, "tsta", native_file);
   writeFixtureMar(root, "tsta");
   % Remove MAR's precipitation so the policy order falls through to the
   % staged MERRA-2 source for the placeholder outage.
   mar_file = fullfile(root, 'met', 'mar3.11', ...
      'met_tsta_mar3.11_20200101_20211231_15m.mat');
   S = load(mar_file, 'mar_met');
   mar_met = S.mar_met;
   mar_met.ppt(:) = NaN;
   mar_met.rainf(:) = NaN;
   mar_met.snowf(:) = NaN;
   save(mar_file, 'mar_met');
   writeFixtureMerraPrecipSplit(root, "tsta");
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      last_resort_proxies=false, plan_n_gaps=1, seed=6);

   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   codes = icemodel.forcing.reconstruct.provenanceCodes();
   adopted = result.provenance.ppt == codes.merra2;
   testCase.verifyTrue(any(adopted));
   % The adoption ships a FULL split: both phases finite and nonnegative.
   testCase.verifyTrue(all(isfinite(result.filled.rainf(adopted))));
   testCase.verifyTrue(all(isfinite(result.filled.snowf(adopted))));
   testCase.verifyTrue(all(result.filled.rainf(adopted) >= 0));
   testCase.verifyTrue(all(result.filled.snowf(adopted) >= 0));
   % Per-sample complement identity holds across every adopted sample,
   % including the guarded source hour.
   guard_support = (501:504).';
   testCase.verifyEqual(result.filled.rainf(adopted) ...
      + result.filled.snowf(adopted), result.filled.ppt(adopted), ...
      'AbsTol', 1e-9);
   testCase.verifyTrue(all(result.provenance.rainf(adopted) ...
      == codes.merra2));
   % The inconsistent staged sample (snowf above ppt) exercises the
   % nonnegative guard: zero rain over the whole posting support, never
   % negative rain, and every delivered sample keeps the split identity.
   testCase.verifyTrue(all(adopted(guard_support)));
   testCase.verifyEqual(result.filled.rainf(guard_support), ...
      zeros(4, 1), 'AbsTol', 1e-8);
   testCase.verifyEqual(result.filled.rainf(guard_support) ...
      + result.filled.snowf(guard_support), ...
      result.filled.ppt(guard_support), 'AbsTol', 1e-9);
   testCase.verifyTrue(any(startsWith(string(result.audit.method), ...
      "proxy:merra2:precip_adoption")));
end

function test_precip_seam_blend_cannot_lower_total_below_native_phase(testCase)
   % A proxy total may pass the native-phase veto before seam blending but
   % fall below a preserved phase afterward. The observation wins and that
   % total sample remains missing rather than publishing rainf > ppt.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   native_file = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(native_file, 'met');
   met = S.met;
   met.ppt([49; 57]) = 0;
   met.rainf(53) = 0.08;
   conflict_time = met.Time(53);
   save(native_file, 'met');
   recordNativeMetIdentity(root, "tsta", native_file);
   writeFixtureMar(root, "tsta");

   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      plan_n_gaps=1, seed=6);
   result = icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts);

   support = result.filled.Time >= conflict_time ...
      & result.filled.Time < conflict_time + hours(1);
   testCase.verifyTrue(all(isnan(result.filled.ppt(support))));
   testCase.verifyEqual(result.filled.rainf(support), ...
      0.08 * ones(nnz(support), 1), 'AbsTol', 1e-12);
   finite_total = isfinite(result.filled.ppt);
   testCase.verifyFalse(any(isfinite(result.filled.rainf) ...
      & finite_total & result.filled.rainf > result.filled.ppt));
   testCase.verifyFalse(any(isfinite(result.filled.snowf) ...
      & finite_total & result.filled.snowf > result.filled.ppt));
end

function test_negative_native_rain_refuses_instead_of_proxy_replacement(testCase)
   % A finite corrected tipping-bucket observation is immutable. If it is
   % invalid, publication must refuse rather than erase and replace it.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   native_file = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(native_file, 'met');
   met = S.met;
   met.rainf(53) = -0.01;
   save(native_file, 'met');
   recordNativeMetIdentity(root, "tsta", native_file);
   writeFixtureMar(root, "tsta");
   opts = icemodel.forcing.reconstruct.setopts( ...
      required_channels="tair", core_channels="tair", ...
      plan_channels="tair", interp_channels="tair", ...
      plan_n_gaps=1, seed=6);

   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      donor_sites=string.empty(1, 0), use_ktransect=false, ...
      use_gcnet=false, write=false, opts=opts), ...
      'icemodel:reconstruct:fillPromiceStation:invalidPrecipitation');
   reloaded = load(native_file, 'met');
   testCase.verifyEqual(reloaded.met.rainf(53), -0.01);
end

function test_fill_station_errors_without_native_met(testCase)
   % A missing staged native met fails loudly.
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.fillPromiceStation("nope", ...
      met_dir=fullfile(testCase.TestData.root, 'met', 'promice'), ...
      use_ktransect=false, use_gcnet=false, write=false), ...
      'icemodel:reconstruct:fillPromiceStation:missingNativeMet');
end

%% buildGapFillReport

function test_report_builds_ledgered_figures_and_qmd(testCase)
   % The generator renders windowed figures, reconciles the ledger, and
   % writes the fixed-structure QMD from saved artifacts only.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, true);
   writeFixtureStation(root, "dsta", 0.05, -20, false);
   writeFixtureStation(root, "egp", 0.1, -30, false);
   tsta_file = fullfile(root, 'met', 'promice', ...
      'met_tsta_promice_20200101_20211231_15m.mat');
   S = load(tsta_file, 'met');
   met = S.met;
   for first = [10001, 11001, 12001, 13001]
      met.tair(first:first + 3) = NaN;
   end
   % A four-day corroborated sensor-failure signature exercises the
   % production exclusion and the report's systematic native-QC category.
   buried = met.Properties.RowTimes >= ...
      datetime(2020, 3, 25, 'TimeZone', 'UTC') ...
      & met.Properties.RowTimes < ...
      datetime(2020, 3, 29, 'TimeZone', 'UTC');
   met.tair(buried) = 250 + 0.05 * sin( ...
      2 * pi * hour(met.Properties.RowTimes(buried)) / 24);
   met.rh(buried) = 85;
   met.swd(buried) = 0.5;
   sigma = icemodel.physicalConstant('SB');
   met.lwd(buried) = sigma * met.tair(buried) .^ 4 + 1;
   save(tsta_file, 'met');
   recordNativeMetIdentity(root, "tsta", tsta_file);
   egp_file = fullfile(root, 'met', 'promice', ...
      'met_egp_promice_20200101_20211231_15m.mat');
   S = load(egp_file, 'met');
   met = S.met;
   partial_day = dateshift(met.Properties.RowTimes, 'start', 'day') == ...
      datetime(2020, 6, 1, 'TimeZone', 'UTC');
   complete_day = dateshift(met.Properties.RowTimes, 'start', 'day') == ...
      datetime(2020, 6, 2, 'TimeZone', 'UTC');
   met.swd(partial_day) = 900;
   met.swd(find(partial_day, 1)) = NaN;
    met.swd(complete_day) = 800;
    save(egp_file, 'met');
    recordNativeMetIdentity(root, "egp", egp_file);
   writeFixtureMar(root, "tsta");
   completeFixtureMarPrecip(root, "tsta");
   writeFixtureMerraPrecipPatch(root, "tsta");
   completeFixtureNativePrecip(root, "dsta");
   completeFixtureNativePrecip(root, "egp");
   report_opts = icemodel.forcing.reconstruct.setopts(seed=13, ...
      rmse_improvement=0.17, min_variability_ratio=0.63, ...
      max_variability_ratio=1.37, min_coverage=0.82);
   icemodel.forcing.reconstruct.fillPromiceStation("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      out_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      donor_sites="dsta", use_ktransect=false, use_gcnet=false, ...
      opts=report_opts);
   icemodel.forcing.reconstruct.fillPromiceStation("dsta", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      out_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      donor_sites="dsta", use_ktransect=false, ...
      use_gcnet=false, opts=report_opts);
   egp_result = icemodel.forcing.reconstruct.fillPromiceStation("egp", ...
      met_dir=fullfile(root, 'met', 'promice'), ...
      out_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      donor_sites="dsta", use_ktransect=false, ...
      use_gcnet=false, opts=report_opts);
   % Canonical metadata aliases identify the same public station token.
   tsta_filled = dir(fullfile(root, 'met', 'promice_filled', ...
      'met_tsta_promice_filled_*_15m.mat'));
   S = load(fullfile(tsta_filled(1).folder, tsta_filled(1).name), 'met');
   met = S.met;
   metadata = met.Properties.UserData;
   metadata.site = "T_STA";
   met.Properties.UserData = metadata;
   save(fullfile(tsta_filled(1).folder, tsta_filled(1).name), 'met');
   refreshReportInputArtifact(fullfile(root, 'qa', 'plans', ...
      'tsta-report-inputs.json'), "filled", root);
   S = load(egp_result.met_file, 'met');
   met = S.met;
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   met.swd_provenance(complete_day) = codes.clamped_shortwave;
   save(egp_result.met_file, 'met');
   refreshReportInputArtifact(fullfile(root, 'qa', 'plans', ...
      'egp-report-inputs.json'), "filled", root);
   % A proxy window staged after production is intentionally unpinned and
   % must not change the report's producer-recorded acceptance window.
   extra_proxy = 1;
   save(fullfile(root, 'met', 'mar3.11', ...
      'met_tsta_mar3.11_20300101_20301231_15m.mat'), 'extra_proxy');
   ghost = readtable(fullfile(root, 'qa', 'ledger', ...
      'tsta-readiness.csv'), 'TextType', 'string');
   ghost.site(:) = "ghost";
   writetable(ghost, fullfile(root, 'qa', 'ledger', ...
      'ghost-readiness.csv'));
   % Unrelated figures from a previous subset render are outside this
   % requested cohort and must not poison ledger reconciliation.
   fid = fopen(fullfile(root, 'figures', 'ghost_tair_stale.png'), 'w');
   fclose(fid);
   fid = fopen(fullfile(root, 'figures', 'tsta_tair_stale.png'), 'w');
   fclose(fid);
   fid = fopen(fullfile(root, 'report', 'gapfill-report.html'), 'w');
   fprintf(fid, 'stale html');
   fclose(fid);
   mkdir(fullfile(root, 'report', 'gapfill-report_files'));
   fid = fopen(fullfile(root, 'report', 'gapfill-report_files', 'stale.js'), 'w');
   fprintf(fid, 'stale asset');
   fclose(fid);

   report = icemodel.verification.report.buildGapFillReport( ...
      render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report'));

   testCase.verifyTrue(isfile(report.qmd_file));
   testCase.verifyTrue(isfile(report.detail_qmd_file));
   qmd_text = string(fileread(report.qmd_file));
   testCase.verifyThat(qmd_text, ...
      matlab.unittest.constraints.ContainsSubstring( ...
      "format:" + newline + "  pdf:" + newline));
   testCase.verifyThat(qmd_text, ...
      matlab.unittest.constraints.ContainsSubstring( ...
      "  html:" + newline + "    toc: true" + newline));
   testCase.verifyFalse(isfile(fullfile(root, 'report', ...
      'gapfill-report.html')));
   testCase.verifyFalse(isfolder(fullfile(root, 'report', ...
      'gapfill-report_files')));
   testCase.verifyTrue(isfile(fullfile(root, 'qa', 'plans', ...
      'tsta-report-inputs.json')));
   testCase.verifyGreaterThan(height(report.figure_ledger), 0);
   tsta_tair = report.summary.site == "tsta" ...
      & report.summary.channel == "tair";
   tsta_tair_figures = report.figure_ledger.site == "tsta" ...
      & report.figure_ledger.channel == "tair" ...
      & report.figure_ledger.method ~= "overview" ...
      & ~startsWith(report.figure_ledger.method, "interpretation:");
   testCase.verifyEqual(report.summary.figures(tsta_tair), ...
      nnz(tsta_tair_figures));
   detail_figures = report.figure_ledger.method ~= "overview" ...
      & ~startsWith(report.figure_ledger.method, "interpretation:");
   testCase.verifyTrue(any(report.figure_ledger.segments_shown( ...
      detail_figures) > 1));
   testCase.verifyTrue(all(report.figure_ledger.segments_shown( ...
      tsta_tair_figures) <= 6));
   testCase.verifyEqual(height(unique(report.figure_ledger(:, ...
      {'site', 'channel', 'method'}), 'rows')), ...
      height(report.figure_ledger));
   testCase.verifyEqual(numel(unique(report.figure_ledger.file)), ...
      height(report.figure_ledger));
   % POLICY A14/D-19: every station leads with one full-period overview
   % figure, present on disk and first in that station's ledger rows.
   % D-31: overview figures live in the dedicated overview/ subfolder.
   for overview_site = ["tsta", "dsta", "egp"]
      overview_rows = report.figure_ledger.site == overview_site ...
         & report.figure_ledger.method == "overview";
      testCase.verifyEqual(nnz(overview_rows), 1);
      testCase.verifyEqual( ...
         string(report.figure_ledger.channel(overview_rows)), "all");
      testCase.verifyTrue(isfile(fullfile(root, 'figures', ...
         'overview', overview_site + "_overview.png")));
      site_rows = find(report.figure_ledger.site == overview_site);
      testCase.verifyEqual( ...
         string(report.figure_ledger.method(site_rows(1))), "overview");
   end
   testCase.verifyFalse(isfile(fullfile(root, 'figures', ...
      'tsta_tair_stale.png')));
   testCase.verifyTrue(isfile(fullfile(root, 'figures', ...
      'ghost_tair_stale.png')));
   text = string(fileread(report.qmd_file));
   for section = ["# Executive Summary", "# Background", "# Methods", ...
         "# Results", "# Summary", "# Appendices"]
      testCase.verifyTrue(contains(text, section));
   end
   detail_text = string(fileread(report.detail_qmd_file));
   % D-31 layout: the overview and detail appendix sections each carry
   % station navigation and anchors in their owning reports.
   testCase.verifyTrue(contains(text, "- [TSTA](#overview-tsta)"));
   testCase.verifyTrue(contains(text, "- [DSTA](#overview-dsta)"));
   testCase.verifyTrue(contains(text, "- [EGP](#overview-egp)"));
   testCase.verifyFalse(contains(text, "- [TSTA](#detail-tsta)"));
   testCase.verifyTrue(contains(text, "### TSTA {#overview-tsta}"));
   testCase.verifyTrue(contains(detail_text, "- [TSTA](#detail-tsta)"));
   testCase.verifyTrue(contains(detail_text, "## TSTA {#detail-tsta}"));
   testCase.verifyTrue(contains(detail_text, "## DSTA {#detail-dsta}"));
   testCase.verifyTrue(contains(detail_text, "## EGP {#detail-egp}"));
   results_text = extractBefore(extractAfter(text, "# Results"), ...
      "## Summary tables");
   testCase.verifyTrue(contains(results_text, "[!["));
    testCase.verifyFalse(contains(results_text, newline + "!["));
    testCase.verifyTrue(contains(text, '"blend_hours"'));
    testCase.verifyTrue(contains(text, '"donor_sites"'));
    testCase.verifyTrue(contains(text, '"use_ktransect"'));
    testCase.verifyTrue(contains(text, '"use_gcnet"'));
   testCase.verifyFalse(contains(text, "../figures/figures/"));
   testCase.verifyTrue(contains(text, "../figures/"));
   % Every preview links from the report that owns its display class.
   for k = 1:height(report.figure_ledger)
      file = string(report.figure_ledger.file(k));
      figure_path = "../figures/" + file;
      is_detail = string(report.figure_ledger.method(k)) ~= "overview" ...
         && ~startsWith(string(report.figure_ledger.method(k)), ...
         "interpretation:");
      if is_detail
         testCase.verifyTrue(contains(detail_text, ...
            "](" + figure_path + ")](" + figure_path + ")"));
         testCase.verifyFalse(contains(text, figure_path));
      else
         testCase.verifyTrue(contains(text, ...
            "](" + figure_path + ")](" + figure_path + ")"));
         testCase.verifyFalse(contains(detail_text, figure_path));
      end
   end
   testCase.verifyFalse(contains(text, "Candidate decisions:"));
   testCase.verifyFalse(contains(text, "Channel summary:"));
   testCase.verifyTrue(contains(detail_text, ...
      "Admitted candidate decisions (compact view):"));
   testCase.verifyTrue(contains(detail_text, "Channel summary:"));
   testCase.verifyTrue(contains(text, "For the selected cohort"));
   testCase.verifyTrue(contains(text, "a 17% RMSE improvement"));
   testCase.verifyTrue(contains(text, ...
      "over persistence for gaps up to 6 h and station day-of-year climatology"));
   testCase.verifyTrue(contains(text, "within 0.63--1.37"));
    testCase.verifyTrue(contains(text, "and 82% coverage"));
    testCase.verifyTrue(contains(text, ...
       "No staged-proxy acceptance window is available for"));
    for csv = ["gapfill_summary.csv", "gapfill_method_diagnostics.csv", ...
          "gapfill_interpretation_catalog.csv", ...
          "gapfill_residual_gaps.csv", ...
          "gapfill_readiness_blockers.csv", ...
          "gapfill_absent_products.csv", ...
          "promice_filled_readiness.csv", "gapfill_figure_ledger.csv"]
       testCase.verifyTrue(contains(text, "../qa/" + csv));
    end
   testCase.verifyFalse(contains(text, ...
      "covers the previously blocking KANM/KANL"));
   testCase.verifyTrue(all(ismember(["segments_shown", ...
      "has_before_context", "has_after_context"], ...
      string(report.figure_ledger.Properties.VariableNames))));
   testCase.verifyEqual( ...
      hours(report.figure_ledger.gap_end ...
      - report.figure_ledger.gap_start), ...
      report.figure_ledger.duration_hours, 'AbsTol', 1e-12);
   testCase.verifyFalse(any(report.figure_ledger.method == "native_context"));
   testCase.verifyEqual(height(report.interpretation_catalog), 7);
   sensor_case = report.interpretation_catalog.category ...
      == "native_sensor_flat_run";
   testCase.verifyTrue(report.interpretation_catalog.present(sensor_case));
   testCase.verifyEqual(report.interpretation_catalog.site(sensor_case), ...
      "tsta");
   testCase.verifyTrue(contains(text, ...
      "## Scientific interpretation and decision-sensitive outcomes"));
   egp_swd = report.summary.site == "egp" ...
      & report.summary.channel == "swd";
   expected_figures = nnz(report.figure_ledger.site == "egp" ...
      & report.figure_ledger.channel == "swd" ...
      & report.figure_ledger.method ~= "overview" ...
      & ~startsWith(report.figure_ledger.method, "interpretation:"));
   testCase.verifyEqual(report.summary.figures(egp_swd), expected_figures);
   testCase.verifyTrue(isfile(fullfile(root, 'qa', ...
      'gapfill_method_diagnostics.csv')));
   testCase.verifyGreaterThan(height(report.method_diagnostics), 0);
   testCase.verifyTrue(all(ismember(["decision", "season", ...
      "selection_baseline", "selection_baseline_rmse", ...
      "uncertainty_status", "selection_sigma1_coverage", ...
      "selection_sigma2_coverage", "evaluation_rmse", ...
      "evaluation_variability_ratio", "evaluation_sigma1_coverage", ...
      "evaluation_sigma2_coverage"], ...
      string(report.method_diagnostics.Properties.VariableNames))));
   short_admitted = report.method_diagnostics.decision == "admitted" ...
      & report.method_diagnostics.bucket == 1;
   testCase.verifyTrue(any(short_admitted));
   testCase.verifyTrue(all(report.method_diagnostics.selection_baseline( ...
      short_admitted) == "persistence"));
   testCase.verifyTrue(any(report.summary.site == "dsta" ...
      & report.summary.channel == "tair"));
   dsta_ppt = report.summary(report.summary.site == "dsta" ...
      & report.summary.channel == "ppt", :);
   testCase.verifyEqual(dsta_ppt.residual_missing_pct, 0);
   testCase.verifyEqual(dsta_ppt.unfilled_segments, 0);
   dsta_residual = report.residual_gaps( ...
      report.residual_gaps.site == "dsta" ...
      & report.residual_gaps.channel == "ppt", :);
   testCase.verifyEmpty(dsta_residual);
   testCase.verifyTrue(contains(text, ...
      "## Why anything remains unfilled"));
   testCase.verifyTrue(contains(text, "| dsta |"));
   testCase.verifyTrue(any(report.summary.site == "tsta" ...
      & report.summary.channel == "albedo"));
   % Every A6 product uses its pinned 2020--2021 proxy window; the
   % unpinned 2030 file above cannot change TSTA.
   combined = readtable(fullfile(root, 'qa', ...
      'promice_filled_readiness.csv'), 'TextType', 'string');
   testCase.verifyTrue(all(ismember(["window_start", "window_end", ...
      "policy_verdict"], string(combined.Properties.VariableNames))));
   testCase.verifyEqual(combined.policy_verdict, ...
      combined.verdict_icemodel);
   testCase.verifyEqual(unique(combined.window_start), ...
      datetime(2020, 1, 1));
   testCase.verifyEqual(unique(combined.window_end), ...
      datetime(2021, 12, 31, 23, 45, 0));
   testCase.verifyFalse(any(combined.site == "ghost"));

   % A byte-valid path substitution must still fail semantic station/product
   % identity checks without replacing the published report.
   ledger_file = fullfile(root, 'qa', 'gapfill_figure_ledger.csv');
   % An explicit subset owns only its requested station. An unrelated native
   % station awaiting production must not be invented as this report's
   % absence or block the subset render.
   writeFixtureStation(root, "xsta", 0.15, -15, false);
   subset = icemodel.verification.report.buildGapFillReport( ...
      sites="tsta", render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report'));
   testCase.verifyEmpty(subset.absent_products);
   qmd_hash = icemodel.verification.setup.fileSha256(subset.qmd_file);
   ledger_hash = icemodel.verification.setup.fileSha256(ledger_file);
   tsta_manifest_file = fullfile(root, 'qa', 'plans', ...
      'tsta-report-inputs.json');
   outside_dir = fullfile(root, 'outside-report-inputs');
   mkdir(outside_dir);
   % Every manifest role is byte-identical after copying, but its selected
   % role root still controls whether the report may trust it.
   for escaped_role = ["native", "filled", "plan", "readiness", ...
         "proxy_window"]
      escaped_manifest = jsondecode(fileread(tsta_manifest_file));
      escaped_index = find(string({escaped_manifest.artifacts.role}) ...
         == escaped_role, 1);
       escaped_record = escaped_manifest.artifacts(escaped_index);
       [~, ~, extension] = fileparts(escaped_record.path);
       escaped_file = fullfile(outside_dir, escaped_role + extension);
       copyfile(fullfile(root, escaped_record.path), escaped_file);
       escaped_manifest.artifacts(escaped_index).path = fullfile( ...
          "outside-report-inputs", escaped_role + extension);
      writeJsonFixture(tsta_manifest_file, escaped_manifest);
      testCase.verifyError(@() ...
         icemodel.verification.report.buildGapFillReport(sites="tsta", ...
         render=false, ...
         filled_dir=fullfile(root, 'met', 'promice_filled'), ...
         qa_dir=fullfile(root, 'qa'), ...
         fig_dir=fullfile(root, 'figures'), ...
         report_dir=fullfile(root, 'report')), ...
         'icemodel:report:buildGapFillReport:inputPathOutsideRoot');
      escaped_manifest.artifacts(escaped_index) = escaped_record;
      writeJsonFixture(tsta_manifest_file, escaped_manifest);
   end

   identity_manifest = jsondecode(fileread(tsta_manifest_file));
   native_index = find(string({identity_manifest.artifacts.role}) ...
      == "native");
   original_native = identity_manifest.artifacts(native_index);
    identity_manifest.artifacts(native_index).path = fullfile( ...
       'met', 'promice', ...
       'met_dsta_promice_20200101_20211231_15m.mat');
    writeJsonFixture(tsta_manifest_file, identity_manifest);
    refreshReportInputArtifact(tsta_manifest_file, "native", root);
   testCase.verifyError(@() ...
      icemodel.verification.report.buildGapFillReport(sites="tsta", ...
      render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report')), ...
      'icemodel:report:buildGapFillReport:artifactIdentityMismatch');
   identity_manifest.artifacts(native_index) = original_native;
   writeJsonFixture(tsta_manifest_file, identity_manifest);

   % A hash-valid proxy from another station remains the wrong scientific
   % input even though it is inside the selected met root.
   identity_manifest = jsondecode(fileread(tsta_manifest_file));
   proxy_index = find(string({identity_manifest.artifacts.role}) ...
      == "proxy_window", 1);
   original_proxy = identity_manifest.artifacts(proxy_index);
    identity_manifest.artifacts(proxy_index).path = fullfile( ...
       'met', 'mar3.11', ...
       'met_dsta_mar3.11_20200101_20211231_15m.mat');
    writeJsonFixture(tsta_manifest_file, identity_manifest);
    refreshReportInputArtifact(tsta_manifest_file, "proxy_window", root);
   testCase.verifyError(@() ...
      icemodel.verification.report.buildGapFillReport(sites="tsta", ...
      render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report')), ...
      'icemodel:report:buildGapFillReport:inputIdentityMismatch');
   identity_manifest.artifacts(proxy_index) = original_proxy;
   writeJsonFixture(tsta_manifest_file, identity_manifest);
   testCase.verifyEqual( ...
      icemodel.verification.setup.fileSha256(report.qmd_file), qmd_hash);
   testCase.verifyEqual( ...
      icemodel.verification.setup.fileSha256(ledger_file), ledger_hash);

   % A late Quarto failure occurs after every staged CSV/QMD/figure exists;
   % none of those staged bytes may replace the published report.
   old_path = getenv('PATH');
   path_cleanup = onCleanup(@() setenv('PATH', old_path));
   setenv('PATH', '');
   testCase.verifyError(@() ...
      icemodel.verification.report.buildGapFillReport(sites="dsta", ...
      render=true, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report')), ...
      'icemodel:report:buildGapFillReport:renderFailed');
   setenv('PATH', old_path);
   clear path_cleanup
   testCase.verifyEqual( ...
      icemodel.verification.setup.fileSha256(report.qmd_file), qmd_hash);
   testCase.verifyEqual( ...
      icemodel.verification.setup.fileSha256(ledger_file), ledger_hash);

   % Producer-pinned donor-pool controls must agree across a report cohort.
   plan_file = fullfile(root, 'qa', 'plans', 'dsta-plan.mat');
   saved_plan = load(plan_file, 'plan_record', 'audit_record');
   original_plan = saved_plan.plan_record;
   plan_record = original_plan;
   plan_record.reconstruction_options.donor_sites = "other";
   audit_record = saved_plan.audit_record;
   save(plan_file, 'plan_record', 'audit_record');
   manifest_file = fullfile(root, 'qa', 'plans', ...
      'dsta-report-inputs.json');
    refreshReportInputArtifact(manifest_file, "plan", root);
   testCase.verifyError(@() ...
      icemodel.verification.report.buildGapFillReport( ...
      sites=["tsta", "dsta"], render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report')), ...
      'icemodel:report:buildGapFillReport:inconsistentPolicy');
   plan_record = original_plan;
   save(plan_file, 'plan_record', 'audit_record');
    refreshReportInputArtifact(manifest_file, "plan", root);
   % An option outside the headline prose is still rendered in the complete
   % Methods record and therefore must be producer-identical.
   plan_record = original_plan;
   plan_record.reconstruction_options.blend_hours = 3;
   save(plan_file, 'plan_record', 'audit_record');
    refreshReportInputArtifact(manifest_file, "plan", root);
   testCase.verifyError(@() ...
      icemodel.verification.report.buildGapFillReport( ...
      sites=["tsta", "dsta"], render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report')), ...
      'icemodel:report:buildGapFillReport:inconsistentPolicy');
   plan_record = original_plan;
   save(plan_file, 'plan_record', 'audit_record');
    refreshReportInputArtifact(manifest_file, "plan", root);

   % A report rerun must fail before changing outputs if any consumed
   % artifact no longer matches the producer-recorded byte identity.
   readiness_file = fullfile(root, 'qa', 'ledger', ...
      'tsta-readiness.csv');
   fid = fopen(readiness_file, 'a');
   fprintf(fid, '\n');
   fclose(fid);
   testCase.verifyError(@() ...
      icemodel.verification.report.buildGapFillReport(render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report')), ...
      'icemodel:report:buildGapFillReport:inputIdentityMismatch');
   testCase.verifyEqual( ...
      icemodel.verification.setup.fileSha256(report.qmd_file), qmd_hash);
end

function test_acceptance_window_derives_from_staged_proxies(testCase)
   % The policy window is the union span of identity-validated staged proxy
   % files; even a non-widest stale artifact must be rejected.
   root = testCase.TestData.root;
   writeFixtureStation(root, "tsta", 0, 0, false);
   writeFixtureMar(root, "tsta");
   location = struct('lat_wgs84', 67.0, 'lon_wgs84', -48.8, ...
      'elev_m', 1200);

   [returned, files] = ...
      icemodel.forcing.reconstruct.acceptanceWindow("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), location=location);
   expected = [datetime(2020, 1, 1, 'TimeZone', 'UTC'), ...
      datetime(2021, 12, 31, 23, 45, 0, 'TimeZone', 'UTC')];
   testCase.verifyEqual(returned, expected);
   testCase.verifyEqual(numel(files), 1);
    testCase.verifyTrue(endsWith(files, ...
       'met_tsta_mar3.11_20200101_20211231_15m.mat'));

    % The same proxy remains discoverable from a flat selected met root.
    proxy_file = files(1);
    flat_file = fullfile(root, 'met', ...
       'met_tsta_mar3.11_20200101_20211231_15m.mat');
    copyfile(proxy_file, flat_file);
    delete(proxy_file);

    % A malformed remnant in the preferred source directory must not hide
    % the valid flat-layout fallback artifact.
    malformed_file = fullfile(fileparts(proxy_file), ...
       'met_tsta_mar3.11_malformed_15m.mat');
    S = load(flat_file, 'mar_met');
    mar_met = S.mar_met;
    save(malformed_file, 'mar_met');
    [flat_window, flat_files] = ...
       icemodel.forcing.reconstruct.acceptanceWindow("tsta", ...
       met_dir=fullfile(root, 'met'), location=location);
    testCase.verifyEqual(flat_window, expected);
    testCase.verifyEqual(string(flat_files), string(flat_file));
    delete(malformed_file);

    % Filename bounds cannot claim days absent from the saved timetable.
    valid_file = flat_files(1);
     S = load(valid_file, 'mar_met');
     mar_met = S.mar_met;
     complete_mar = mar_met;

     % Filename tokens identify boundary dates; the returned window pins the
     % exact timetable endpoints rather than inventing the missing hours.
     mar_met = complete_mar(3:end - 2, :);
     save(valid_file, 'mar_met');
     partial_window = ...
        icemodel.forcing.reconstruct.acceptanceWindow("tsta", ...
        met_dir=fullfile(root, 'met'), location=location);
     testCase.verifyEqual(partial_window, ...
        mar_met.Properties.RowTimes([1, end]).');

     mar_met = mar_met(day(mar_met.Properties.RowTimes) ~= 1 ...
        | month(mar_met.Properties.RowTimes) ~= 1 ...
        | year(mar_met.Properties.RowTimes) ~= 2020, :);
    save(valid_file, 'mar_met');
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.acceptanceWindow("tsta", ...
       met_dir=fullfile(root, 'met'), location=location), ...
       'icemodel:reconstruct:acceptanceWindow:proxyWindowMismatch');
    mar_met = complete_mar;
    save(valid_file, 'mar_met');
    mar_met(100, :) = [];
    save(valid_file, 'mar_met');
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.acceptanceWindow("tsta", ...
       met_dir=fullfile(root, 'met'), location=location), ...
       'icemodel:reconstruct:acceptanceWindow:proxyWindowMismatch');
     mar_met = complete_mar;
     save(valid_file, 'mar_met');
     mar_met.Properties.RowTimes = ...
        mar_met.Properties.RowTimes + minutes(7);
     save(valid_file, 'mar_met');
     testCase.verifyError(@() ...
        icemodel.forcing.reconstruct.acceptanceWindow("tsta", ...
        met_dir=fullfile(root, 'met'), location=location), ...
        'icemodel:reconstruct:acceptanceWindow:proxyWindowMismatch');
     mar_met = complete_mar;
     save(valid_file, 'mar_met');

     % A one-row same-day artifact has matching endpoint dates but cannot
     % prove continuous 15-minute support.
     delete(valid_file);
     one_row_file = fullfile(fileparts(valid_file), ...
        'met_tsta_mar3.11_20200101_20200101_15m.mat');
     mar_met = complete_mar(1, :);
     save(one_row_file, 'mar_met');
     testCase.verifyError(@() ...
        icemodel.forcing.reconstruct.acceptanceWindow("tsta", ...
        met_dir=fullfile(root, 'met'), location=location), ...
        'icemodel:reconstruct:acceptanceWindow:proxyWindowMismatch');
     delete(one_row_file);
     mar_met = complete_mar;
     save(valid_file, 'mar_met');

    % Adjacent same-source files are all pinned: the widest anchors and
    % the sibling extending before it joins, so the window covers every
    % staged span reconstruction loads (POLICY A6).
    delete(valid_file);
    first_file = fullfile(fileparts(valid_file), ...
       'met_tsta_mar3.11_20200101_20200630_15m.mat');
    mar_met = complete_mar(complete_mar.Properties.RowTimes < ...
       datetime(2020, 7, 1, 'TimeZone', 'UTC'), :);
    save(first_file, 'mar_met');
    second_file = fullfile(fileparts(valid_file), ...
       'met_tsta_mar3.11_20200701_20211231_15m.mat');
    mar_met = complete_mar(complete_mar.Properties.RowTimes >= ...
       datetime(2020, 7, 1, 'TimeZone', 'UTC'), :);
    save(second_file, 'mar_met');
    [selected_window, selected_files] = ...
       icemodel.forcing.reconstruct.acceptanceWindow("tsta", ...
       met_dir=fullfile(root, 'met'), location=location);
    testCase.verifyEqual(selected_window, [ ...
       datetime(2020, 1, 1, 'TimeZone', 'UTC'), ...
       datetime(2021, 12, 31, 23, 45, 0, 'TimeZone', 'UTC')]);
    testCase.verifyEqual(sort(string(selected_files)), ...
       sort([string(first_file); string(second_file)]));

    % Adjacent filename dates cannot conceal a sub-day gap between the
    % exact timetable supports.
    S = load(first_file, 'mar_met');
    first_mar = S.mar_met;
    mar_met = first_mar(1:end - 24, :);
    save(first_file, 'mar_met');
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.acceptanceWindow("tsta", ...
       met_dir=fullfile(root, 'met'), location=location), ...
       'icemodel:reconstruct:acceptanceWindow:proxyCoverageGap');
    mar_met = first_mar;
    save(first_file, 'mar_met');

    % Individually complete files cannot define one scalar window across an
    % unstaged day in their union.
    delete(second_file);
    second_file = fullfile(fileparts(valid_file), ...
       'met_tsta_mar3.11_20200702_20211231_15m.mat');
    mar_met = complete_mar(complete_mar.Properties.RowTimes >= ...
       datetime(2020, 7, 2, 'TimeZone', 'UTC'), :);
    save(second_file, 'mar_met');
    testCase.verifyError(@() ...
       icemodel.forcing.reconstruct.acceptanceWindow("tsta", ...
       met_dir=fullfile(root, 'met'), location=location), ...
       'icemodel:reconstruct:acceptanceWindow:proxyCoverageGap');
    delete(first_file);
    delete(second_file);
    mar_met = complete_mar;
    save(valid_file, 'mar_met');

    none = icemodel.forcing.reconstruct.acceptanceWindow("ghost", ...
      met_dir=fullfile(root, 'met', 'promice'), location=location);
   testCase.verifyTrue(all(isnat(none)));

   S = load(valid_file, 'mar_met');
   mar_met = S.mar_met;
   mar_met.Properties.UserData.site = "other";
   stale_file = fullfile(fileparts(valid_file), ...
      'met_tsta_mar3.11_20200101_20201231_15m.mat');
   save(stale_file, 'mar_met');
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.acceptanceWindow("tsta", ...
      met_dir=fullfile(root, 'met', 'promice'), location=location), ...
      'icemodel:reconstruct:acceptanceWindow:proxyIdentityMismatch');
   testCase.verifyError(@() ...
      icemodel.forcing.reconstruct.acceptanceWindow("../escape", ...
      met_dir=fullfile(root, 'met', 'promice'), location=location), ...
      'icemodel:reconstruct:mustBeStationToken:invalidToken');
end

function test_report_errors_without_products(testCase)
   % No filled products means no report, loudly.
   root = testCase.TestData.root;
   testCase.verifyError(@() ...
      icemodel.verification.report.buildGapFillReport(render=false, ...
      filled_dir=fullfile(root, 'met', 'promice_filled'), ...
      qa_dir=fullfile(root, 'qa'), ...
      fig_dir=fullfile(root, 'figures'), ...
      report_dir=fullfile(root, 'report')), ...
       'icemodel:report:buildGapFillReport:noFilledProducts');
end

function test_shell_quote_preserves_literal_metacharacters(testCase)
   % Quarto paths containing quotes or shell substitutions must reach the
   % process as one literal argument without expansion.
   value = "path with ' quote;$(printf injected);`printf bad`";
   command = "printf %s " + ...
      icemodel.verification.helpers.shellQuote(value);

   [status, returned] = system(command);

   testCase.verifyEqual(status, 0);
   testCase.verifyEqual(string(returned), value);
end

%% Fixtures

function [target, donor] = syntheticPair(dlat, delev)
   %SYNTHETICPAIR In-memory target/donor structs with correlated channels.
   % Both carry seeded gaps so either can play the target role (a gapless
   % target correctly admits nothing - there is nothing to validate), and
   % both share a seeded weather anomaly: without irreducible day-to-day
   % variability, climatology reproduces a deterministic fixture exactly
   % and no donor could ever clear the improvement gate.
   % Three fixture years keep the one-year fitting-overlap requirement off
   % the threshold edge for every split seed.
   anomaly = weatherAnomaly();
   target = struct('series', syntheticMet(0, true, anomaly, 2019, 3), ...
      'station', "tsta", ...
      'location', struct('lat_wgs84', 67.0, 'lon_wgs84', -48.8, ...
      'elev_m', 1200));
   donor = struct('series', ...
      syntheticMet(0.4, true, 0.9 * anomaly, 2019, 3), ...
      'station', "dsta", ...
      'family', "promice", ...
      'location', struct('lat_wgs84', 67.0 + dlat, 'lon_wgs84', -48.8, ...
      'elev_m', 1200 + delev), 'observed_mask', []);
   % The reversed plan needs the target back in donor-struct shape.
   target.family = "promice";
   target.observed_mask = [];
end

function anomaly = weatherAnomaly()
   %WEATHERANOMALY Seeded synoptic-scale temperature anomaly (K).
   % A bounded random walk at daily scale, support-held to 15 minutes, so
   % consecutive samples stay smooth while day-to-day values are genuinely
   % unpredictable by day-of-year climatology.
   stream = RandStream('mt19937ar', 'Seed', 7);
   n_days = 1097;
   daily = cumsum(randn(stream, n_days, 1));
   daily = 3 * (daily - mean(daily)) / std(daily);
   anomaly = repelem(daily, 96);
end

function proxies = emptyProxies()
   %EMPTYPROXIES Typed empty proxy pool.
   proxies = struct('series', {}, 'name', {}, 'code_name', {});
end

function filename = writeGcnetDonorFixture(root)
   %WRITEGCNETDONORFIXTURE Save a minimal origin-flagged GC-Net file.
   filename = fullfile(root, 'TEST_surface.nc');
   nccreate(filename, 'time', 'Dimensions', {'time', 3});
   ncwrite(filename, 'time', [0; 1; 2]);
   ncwriteatt(filename, 'time', 'units', ...
      'hours since 2020-1-1 0:0:0');
   names = ["Ta_2m", "RH_2m", "WS_10m", "SRin", "LRin"];
   for k = 1:numel(names)
      nccreate(filename, names(k), 'Dimensions', {'time', 3});
      ncwrite(filename, names(k), [1; 2; 3]);
      if names(k) ~= "LRin"
         nccreate(filename, names(k) + "_origin", ...
            'Dimensions', {'time', 3});
         origin = [0; 0; 0];
         if names(k) == "RH_2m"
            origin(1) = 1;
         end
         ncwrite(filename, names(k) + "_origin", origin);
      end
   end
   ncwriteatt(filename, '/', 'latitude', 67);
   ncwriteatt(filename, '/', 'longitude', -48);
   ncwriteatt(filename, '/', 'elevation', 1000);
end

function met = syntheticMet(offset, with_gaps, anomaly, start_year, n_years)
   %SYNTHETICMET N years of 15-minute synthetic met with optional gaps.
   if nargin < 4
      start_year = 2020;
   end
   if nargin < 5
      n_years = 2;
   end
   times = (datetime(start_year, 1, 1, 'TimeZone', 'UTC'): ...
      minutes(15):datetime(start_year + n_years - 1, 12, 31, 23, 45, 0, ...
      'TimeZone', 'UTC')).';
   n = numel(times);
   if nargin < 3
      anomaly = zeros(n, 1);
   end
   % Callers hand a whole-day anomaly at least as long as the axis; slice
   % to the axis so leap-year day counts cannot desynchronize the fixture.
   anomaly = anomaly(1:n);
   doy = day(times, 'dayofyear');
   hod = hour(times) + minute(times) / 60;
   tair = 255 + offset + anomaly + ...
      15 * sin(2 * pi * (doy - 30) / 366) + ...
      2 * sin(2 * pi * hod / 24);
   rh = min(100, max(5, 80 + 10 * sin(2 * pi * doy / 366)));
   wspd = max(0.5, 6 + 2 * sin(2 * pi * doy / 366));
   psfc = 88000 + 500 * sin(2 * pi * doy / 366);
   elevation = icemodel.forcing.helpers.solarElevation(times, 67.0, -48.8);
   swd = max(0, 600 * sind(max(elevation(:), 0)));
   albedo = min(0.95, max(0.4, 0.8 - 0.2 * sin(2 * pi * doy / 366)));
   swu = albedo .* swd;
   lwd = 230 + 60 * sin(2 * pi * (doy - 30) / 366);
   ppt = nan(n, 1);
   rainf = nan(n, 1);
   snowf = nan(n, 1);
   boom_height = 2.5 + 0.15 * sin(2 * pi * doy / 366);
   met = timetable(times, tair, rh, wspd, psfc, swd, swu, lwd, albedo, ...
      ppt, rainf, snowf, boom_height, ...
      'VariableNames', {'tair', 'rh', 'wspd', 'psfc', 'swd', 'swu', ...
      'lwd', 'albedo', 'ppt', 'rainf', 'snowf', 'boom_height'});
   if with_gaps
      % Seeded gaps: one tier-1-sized (3 h) and one donor-sized (36 h) per
      % core channel, placed away from the record edges.
      for name = ["tair", "rh", "psfc", "lwd"]
         met.(name)(5001:5012) = NaN;            % 3 hourly postings
         met.(name)(30001:30144) = NaN;          % 36 hourly postings
      end
      met.boom_height(5001:5012) = NaN;          % maintenance-safe short gap
      met.boom_height(30001:30144) = NaN;        % remains readiness-blocking
   end
end

function writeFixtureStation(root, site, dlat, delev, with_gaps, with_proxy)
   %WRITEFIXTURESTATION Save one production-viable staged-met fixture.
   % Boundary side (POLICY A16/D-24): this fixture emulates a CANONICAL
   % staged artifact, so its upward shortwave ships as swu. The legacy
   % pre-rename artifact side (a staged usr channel) is covered once, in
   % test_staged_usr_fixture_ships_swu_product, by renaming this fixture's
   % swu back to the source name before the fill.
   if nargin < 6
      with_proxy = true;
   end
   source = syntheticMet(0.4 * (site == "dsta"), with_gaps);
   source.Properties.DimensionNames{1} = 'Time';
   met = icemodel.forcing.helpers.resampleMetTimestep( ...
      source(1:4:end, :), "15m");
   raw_times = met.Properties.RowTimes(1:4:end);
   source_file = writeRawPromiceAlbedo(root, site, raw_times, ...
      met.albedo(1:4:end));
    source_info = dir(source_file);
    ud = struct('site', site, 'lat', 67.0 + dlat, 'lon', -48.8, ...
       'elev', 1200 + delev, 'source_file', source_file, ...
       'met_resample_policy', ...
       met.Properties.UserData.met_resample_policy, ...
       'met_resample_source_cadence_seconds', ...
       met.Properties.UserData.met_resample_source_cadence_seconds, ...
       'source_size_bytes', source_info.bytes, ...
      'source_sha256', ...
      icemodel.verification.setup.fileSha256(source_file), ...
      'albedo_policy', ...
      "albedo from PROMICE L3 albedo; fillPromiceAlbedo(fillwinter=1)");
    met.Properties.UserData = ud;
    filename = fullfile(root, 'met', 'promice', sprintf( ...
       'met_%s_promice_20200101_20211231_15m.mat', site));
    save(filename, 'met');
    recordNativeMetIdentity(root, site, filename);
    if with_proxy
       writeFixtureMar(root, site, dlat);
    end
 end

function writeFixtureFamilyStation(root, site, family, with_gaps, ...
      cadence_seconds)
   %WRITEFIXTUREFAMILYSTATION Save one viable non-promice staged fixture.
   % Emulates the IMAU staging shape (bead icemodel-g1n.49): identity rides
   % ud.station, top-level metadata carries lat_wgs84/lon_wgs84 but no
   % elev_m so the complete point lives only in ud.site_location, and no
   % raw-source or albedo-policy builder pins exist. The family manifest
   % pins the staged path only, mirroring the production data/eval/imau
   % manifest form.
   if nargin < 5
      cadence_seconds = 3600;
   end
   source = syntheticMet(0, with_gaps);
   source.Properties.DimensionNames{1} = 'Time';
   source_stride = cadence_seconds / 900;
   met = icemodel.forcing.helpers.resampleMetTimestep( ...
      source(1:source_stride:end, :), "15m");
   ud = struct('station', upper(site), ...
      'site_location', struct('lat_wgs84', 67.0, 'lon_wgs84', -48.8, ...
      'elev_m', 1200), ...
      'lat_wgs84', 67.0, 'lon_wgs84', -48.8, ...
      'met_resample_policy', ...
      met.Properties.UserData.met_resample_policy, ...
      'met_resample_source_cadence_seconds', ...
      met.Properties.UserData.met_resample_source_cadence_seconds);
   met.Properties.UserData = ud;
   family_dir = fullfile(root, 'met', char(family));
   if ~isfolder(family_dir)
      mkdir(family_dir);
   end
   filename = fullfile(family_dir, sprintf( ...
      'met_%s_%s_20200101_20211231_15m.mat', site, family));
   save(filename, 'met');
   recordFamilyMetPathIdentity(site, family, filename);
   writeFixtureMar(root, site);
end

function recordFamilyMetPathIdentity(site, family, filename)
   %RECORDFAMILYMETPATHIDENTITY Pin one staged path in a family manifest.
   % Mirrors the production non-promice manifests (e.g. data/eval/imau):
   % the colocation.<family> leg pins the met_files path without per-file
   % hashes, exercising the legacy path-identity branch of
   % verifyNativeMetIdentity (bead icemodel-g1n.49).
   [met_dir, stem, ext] = fileparts(filename);
   data_root = ...
      icemodel.forcing.reconstruct.selectedDataRoot(string(met_dir));
   manifest_dir = fullfile(data_root, 'eval', char(family));
   if ~isfolder(manifest_dir)
      mkdir(manifest_dir);
   end
   [~, family_folder] = fileparts(met_dir);
   leg = struct('met_files', string(fullfile(string(family_folder), ...
      string(stem) + string(ext))));
   entry = struct('case_id', site, ...
      'colocation', struct(char(family), leg));
   manifest = struct('dataset_family', string(family), 'cases', entry);
   writeJsonFixture(fullfile(manifest_dir, 'manifest.json'), manifest);
end

function met = pinRawSource(met, source_file)
   %PINRAWSOURCE Refresh fixture provenance after replacing raw bytes.
   source_info = dir(source_file);
   met.Properties.UserData.source_file = source_file;
   met.Properties.UserData.source_size_bytes = source_info.bytes;
   met.Properties.UserData.source_sha256 = ...
      icemodel.verification.setup.fileSha256(source_file);
end

function writeFixtureMar(root, site, dlat)
   %WRITEFIXTUREMAR Save one synthetic staged MAR met with finite precip.
   % The driver reads proxies from the canonical sibling directory
   % met/mar3.11 using the staged filename convention.
   if nargin < 3
      dlat = 0;
   end
   mar_met = syntheticMet(1.0, false);
   mar_met.ppt = 0.1 + 0.05 * sin((1:height(mar_met)).' / 500);
   % Real staged MAR ships rainf/snowf beside the total; adoption rescales
   % this source split onto the tapered total (POLICY A10/D-18).
   mar_met.rainf = 0.25 * mar_met.ppt;
   mar_met.snowf = 0.75 * mar_met.ppt;
   mar_met.ppt(100:101) = NaN;
   mar_met.rainf(100:101) = NaN;
   mar_met.snowf(100:101) = NaN;
   mar_met.Properties.UserData = struct('site', site, ...
      'lat_wgs84', 67.0 + dlat, 'lon_wgs84', -48.8, ...
      'mar_qc_status', "applied", ...
      'source_files', "MARv3.11-ERA5-15km-2020.nc");
   mar_dir = fullfile(root, 'met', 'mar3.11');
   if ~isfolder(mar_dir)
      mkdir(mar_dir);
   end
   save(fullfile(mar_dir, sprintf( ...
      'met_%s_mar3.11_20200101_20211231_15m.mat', site)), 'mar_met');
end

function writeFixtureMarWindow(root, site, t0, t1)
   %WRITEFIXTUREMARWINDOW Save a staged MAR met spanning one custom window.
   % Satisfies the acceptanceWindow contract: 15-minute cadence, endpoint
   % dates matching the filename window, and pinned identity metadata.
   times = (t0:minutes(15):t1 + hours(23) + minutes(45)).';
   mar_met = timetable(times, ...
      260 + zeros(numel(times), 1), 0.1 + zeros(numel(times), 1), ...
      'VariableNames', {'tair', 'ppt'});
   mar_met.Properties.UserData = struct('site', site, ...
      'lat_wgs84', 67.0, 'lon_wgs84', -48.8, ...
      'mar_qc_status', "applied", ...
      'source_files', "MARv3.11-ERA5-15km-window.nc");
   mar_dir = fullfile(root, 'met', 'mar3.11');
   if ~isfolder(mar_dir)
      mkdir(mar_dir);
   end
   save(fullfile(mar_dir, sprintf('met_%s_mar3.11_%s_%s_15m.mat', ...
      site, string(t0, 'yyyyMMdd'), string(t1, 'yyyyMMdd'))), 'mar_met');
end

function writeFixtureModis(root, site)
   %WRITEFIXTUREMODIS Save one synthetic staged MODIS daily-albedo artifact.
   % Matches the canonical writeuserdata schema the staging layer produces:
   % a daily 'Data' timetable with an 'albedo' column plus top-level
   % artifact metadata carrying the native cadence stamp.
   days_axis = (datetime(2020, 1, 1):caldays(1):datetime(2021, 12, 31)).';
   Data = timetable(days_axis, 0.6 * ones(numel(days_axis), 1), ...
      'VariableNames', {'albedo'});
   artifact_metadata = struct('site', site, ...
      'artifact_cadence_seconds', 86400);
   modis_dir = fullfile(root, 'userdata', 'modis');
   if ~isfolder(modis_dir)
      mkdir(modis_dir);
   end
   save(fullfile(modis_dir, sprintf( ...
      '%s_modis_20200101_20211231_86400s.mat', site)), ...
      'Data', 'artifact_metadata');
end

function writeFixtureMerraPrecipSplit(root, site)
   %WRITEFIXTUREMERRAPRECIPSPLIT Save staged MERRA-2 precip without rainf.
   % Mirrors the staged MERRA-2 channel inventory (POLICY D-31): the
   % corrected total and snowfall ship, but no rain channel does, so
   % adoption must derive rain as the exact complement ppt - snowf.
   axis_met = syntheticMet(0, false);
   times = axis_met.Properties.RowTimes;
   n = numel(times);
   ppt = 0.2 + 0.1 * sin((1:n).' / 700);
   snowf = 0.75 * ppt;
   % One inconsistent sample ON THE HOURLY SOURCE GRID (rows 1, 5, 9, ...)
   % exercises the nonnegative guard: staged snowfall exceeding the total
   % must derive zero rain, never negative.
   snowf(501) = ppt(501) + 0.05;
   merra_met = timetable(times, ppt, snowf, ...
      'VariableNames', {'ppt', 'snowf'});
   merra_met.Properties.UserData = struct('site', site, ...
      'lat_wgs84', 67.0, 'lon_wgs84', -48.8, ...
      'merra_source_time_coordinate', "native_at_reader", ...
      'merra_time_relabel_policy', ...
      "time_averaged_center_to_interval_start", ...
      'merra_time_upsample_policy', ...
      "zero_order_hold_over_declared_support", ...
      'merra_collection_support_hours', ...
      struct('slv', 1, 'rad', 1, 'flx', 1, 'glc', 3));
   merra_dir = fullfile(root, 'met', 'merra2');
   if ~isfolder(merra_dir)
      mkdir(merra_dir);
   end
   save(fullfile(merra_dir, sprintf( ...
      'met_%s_merra2_20200101_20211231_15m.mat', site)), 'merra_met');
end

function writeFixtureMerraPrecipPatch(root, site)
   %WRITEFIXTUREMERRAPRECIPPATCH Save a partial fallback that must not mix.
   mar_met = syntheticMet(0, false);
   merra_met = timetable(mar_met.Properties.RowTimes, ...
      nan(height(mar_met), 1), 'VariableNames', {'ppt'});
   merra_met.ppt(100:101) = 0.2;
   merra_met.Properties.UserData = struct('site', site, ...
      'lat_wgs84', 67.0, 'lon_wgs84', -48.8, ...
      'merra_source_time_coordinate', "native_at_reader", ...
      'merra_time_relabel_policy', ...
      "time_averaged_center_to_interval_start", ...
      'merra_time_upsample_policy', ...
      "zero_order_hold_over_declared_support", ...
      'merra_collection_support_hours', ...
      struct('slv', 1, 'rad', 1, 'flx', 1, 'glc', 3));
   merra_dir = fullfile(root, 'met', 'merra2');
   if ~isfolder(merra_dir)
      mkdir(merra_dir);
   end
   save(fullfile(merra_dir, sprintf( ...
      'met_%s_merra2_20200101_20211231_15m.mat', site)), 'merra_met');
end

function completeFixtureMarPrecip(root, site)
   %COMPLETEFIXTUREMARPRECIP Make publication fixtures snow-model ready.
   filename = fullfile(root, 'met', 'mar3.11', sprintf( ...
      'met_%s_mar3.11_20200101_20211231_15m.mat', site));
   S = load(filename, 'mar_met');
   mar_met = S.mar_met;
   missing = ~isfinite(mar_met.ppt);
   mar_met.ppt(missing) = 0.2;
   mar_met.rainf(missing) = 0.05;
   mar_met.snowf(missing) = 0.15;
   save(filename, 'mar_met');
end

function completeFixtureNativePrecip(root, site)
   %COMPLETEFIXTURENATIVEPRECIP Make a no-proxy fixture snow-model ready.
   filename = fullfile(root, 'met', 'promice', sprintf( ...
      'met_%s_promice_20200101_20211231_15m.mat', site));
   S = load(filename, 'met');
   met = S.met;
   met.ppt(:) = 0;
   met.rainf(:) = 0;
   met.snowf(:) = 0;
   save(filename, 'met');
   recordNativeMetIdentity(root, site, filename);
end

function metadata = ktransectMetadata(station, heights)
   %KTRANSECTMETADATA Minimal pinned source identity for donor fixtures.
   child = struct('year', 2010, 'doi', "10.1594/PANGAEA.950093");
   location = struct('lat_wgs84', 67.0, 'lon_wgs84', -48.8, ...
      'elev_m', 1200);
   metadata = struct('station', station, ...
      'doi', "10.1594/PANGAEA.950093", ...
      'bundle_doi', "10.1594/PANGAEA.947483", ...
      'license', "CC-BY-4.0", 'children', child, ...
      'sensor_heights', heights, 'site_location', location);
end

function leg = ktransectLeg(filename, relative_file, metadata)
   %KTRANSECTLEG Pin one donor fixture's bytes and source metadata.
   info = dir(filename);
   leg = struct('evaluation_file', relative_file, ...
      'evaluation_size_bytes', info.bytes, ...
      'evaluation_sha256', ...
      icemodel.verification.setup.fileSha256(filename), ...
      'doi', metadata.doi, 'bundle_doi', metadata.bundle_doi, ...
      'license', metadata.license, 'children', metadata.children);
end

function refreshReportInputArtifact(manifest_file, role, data_root)
   %REFRESHREPORTINPUTARTIFACT Re-pin one deliberately mutated test artifact.
   manifest = jsondecode(fileread(manifest_file));
   match = find(string({manifest.artifacts.role}) == role, 1);
    pathname = string(fullfile(data_root, ...
       manifest.artifacts(match).path));
   info = dir(pathname);
   manifest.artifacts(match).bytes = info.bytes;
   manifest.artifacts(match).sha256 = ...
      icemodel.verification.setup.fileSha256(pathname);
   writeJsonFixture(manifest_file, manifest);
end

function recordNativeMetIdentity(~, site, filename)
   %RECORDNATIVEMETIDENTITY Pin one fixture in its producer manifest.
   [met_dir, stem, ext] = fileparts(filename);
   data_root = icemodel.forcing.reconstruct.selectedDataRoot(string(met_dir));
   manifest_file = fullfile(data_root, 'eval', 'promice', 'manifest.json');
   manifest_dir = fileparts(manifest_file);
   if ~isfolder(manifest_dir)
      mkdir(manifest_dir);
   end
   info = dir(filename);
   [~, family] = fileparts(met_dir);
   identity = struct('file', string(fullfile(string(family), ...
      string(stem) + string(ext))), ...
      'size_bytes', info.bytes, 'sha256', ...
      icemodel.verification.setup.fileSha256(filename));
   leg = struct('met_file_identities', identity);
   entry = struct('case_id', site, ...
      'colocation', struct('promice', leg));
   if isfile(manifest_file)
      manifest = jsondecode(fileread(manifest_file));
      cases = manifest.cases;
      case_ids = string({cases.case_id});
      match = case_ids == site;
      if any(match)
         current = cases(match);
         identities = current.colocation.promice.met_file_identities;
         declared = string({identities.file});
         same_file = declared == string(identity.file);
         if any(same_file)
            identities(same_file) = identity;
         else
            identities(end + 1) = identity;
         end
         current.colocation.promice.met_file_identities = identities;
         cases(match) = current;
      else
         cases(end + 1) = entry;
      end
      manifest.cases = cases;
   else
      manifest = struct('dataset_family', "promice", 'cases', entry);
   end
   writeJsonFixture(manifest_file, manifest);
end

function writeJsonFixture(filename, value)
   %WRITEJSONFIXTURE Write one compact JSON fixture with guaranteed closure.
   fid = fopen(filename, 'w');
   cleaner = onCleanup(@() fclose(fid));
   fprintf(fid, '%s', jsonencode(value));
   clear cleaner
end

function filename = writeRawPromiceAlbedo(root, site, times, albedo, swd, swu)
   %WRITERAWPROMICEALBEDO Save the minimal raw-L3 albedo source contract.
   source_dir = fullfile(root, 'raw', 'hour');
   if ~isfolder(source_dir)
      mkdir(source_dir);
   end
   filename = fullfile(source_dir, upper(site) + "_hour.nc");
   if isfile(filename)
      delete(filename);
   end
   n = numel(times);
   nccreate(filename, 'time', 'Dimensions', {'time', n});
   ncwrite(filename, 'time', hours(times - times(1)));
   ncwriteatt(filename, 'time', 'units', ...
      "hours since " + string(times(1), 'yyyy-MM-dd HH:mm:ss'));
   nccreate(filename, 'albedo', 'Dimensions', {'time', n});
   ncwrite(filename, 'albedo', albedo);
    if nargin >= 5
       nccreate(filename, 'dsr', 'Dimensions', {'time', n});
       ncwrite(filename, 'dsr', swd);
    end
    if nargin >= 6
       nccreate(filename, 'usr', 'Dimensions', {'time', n});
       ncwrite(filename, 'usr', swu);
    end
   ncwriteatt(filename, '/', 'site_id', upper(site));
   ncwriteatt(filename, '/', 'latitude', 67);
   ncwriteatt(filename, '/', 'longitude', -48.8);
   ncwriteatt(filename, '/', 'altitude', 1200);
end
