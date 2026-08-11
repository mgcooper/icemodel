function tests = test_ablation_performance_metrics
   %TEST_ABLATION_PERFORMANCE_METRICS Verify the model-versus-observation scoring.
   tests = functiontests(localfunctions);
end

function test_perfect_diagnostic_scores_exactly(testCase)
   % A diagnostic equal to the observations must score zero error and unit
   % efficiency, which pins the sign convention and the metric definitions.

   results = syntheticResults();
   [per_case, aggregate] = ...
      icemodel.verification.ablationPerformanceMetrics(results);

   rows = per_case.case_id == "perfect" ...
      & per_case.diagnostic == "model_melt_mwe" ...
      & per_case.density_kg_m3 == scoringDensity();
   testCase.verifyTrue(per_case.scored(rows));
   testCase.verifyEqual(per_case.endpoint_error_mwe(rows), 0, 'AbsTol', 1e-12);
   testCase.verifyEqual(per_case.bias_mwe(rows), 0, 'AbsTol', 1e-12);
   testCase.verifyEqual(per_case.mae_mwe(rows), 0, 'AbsTol', 1e-12);
   testCase.verifyEqual(per_case.rmse_mwe(rows), 0, 'AbsTol', 1e-12);
   testCase.verifyEqual(per_case.nse(rows), 1, 'AbsTol', 1e-12);

   % The aggregate carries one row per diagnostic and density, so its height
   % follows the policy band rather than a hand-maintained count.
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   n_densities = numel(policy.effective_density_kg_m3);
   n_diagnostics = numel(unique(aggregate.diagnostic));
   testCase.verifyEqual(height(aggregate), n_diagnostics * n_densities);
   testCase.verifyEqual(numel(unique(aggregate.density_kg_m3)), n_densities);
end

function test_over_prediction_has_positive_signed_error(testCase)
   % Positive endpoint error and bias must mean the model exceeded the
   % observation, because the report states the direction in prose.

   results = syntheticResults();
   per_case = icemodel.verification.ablationPerformanceMetrics(results);

   % A constant offset shifts the level but not the rate, so it shows up in
   % the endpoint error and cancels out of the increments; the two are
   % reported separately.
   rows = per_case.case_id == "perfect" ...
      & per_case.diagnostic == "model_runoff_mwe" ...
      & per_case.density_kg_m3 == scoringDensity();
   testCase.verifyEqual(per_case.endpoint_error_mwe(rows), 0.05, ...
      'AbsTol', 1e-12);
   testCase.verifyEqual(per_case.bias_mwe(rows), 0, 'AbsTol', 1e-12);
   testCase.verifyEqual(per_case.rmse_mwe(rows), 0, 'AbsTol', 1e-12);

   % A diagnostic that ablates 20 percent too fast does bias the increments.
   rows = per_case.case_id == "perfect" ...
      & per_case.diagnostic == "model_ablation_proxy_mwe" ...
      & per_case.density_kg_m3 == scoringDensity();
   testCase.verifyGreaterThan(per_case.endpoint_error_mwe(rows), 0);
   testCase.verifyGreaterThan(per_case.bias_mwe(rows), 0);
   testCase.verifyGreaterThanOrEqual(per_case.mae_mwe(rows), 0);
end

function test_unscorable_cases_are_retained_with_a_reason(testCase)
   % An incomplete case must appear with an explicit reason. The metrics must
   % not drop it, because a dropped case looks the same as an agreeing case.

   results = syntheticResults();
   [per_case, aggregate] = ...
      icemodel.verification.ablationPerformanceMetrics(results);

   policy = icemodel.verification.namelists.promiceAblationPolicy();
   n_densities = numel(policy.effective_density_kg_m3);
   n_diagnostics = numel(unique(aggregate.diagnostic));

   % The unscorable case still occupies one row per diagnostic and density.
   rows = per_case.case_id == "unavailable";
   testCase.verifyEqual(nnz(rows), n_diagnostics * n_densities);
   testCase.verifyFalse(any(per_case.scored(rows)));
   testCase.verifyTrue(all(strlength(per_case.reason(rows)) > 0));
   testCase.verifyTrue(all(aggregate.n_excluded == 1));
   testCase.verifyTrue(all(aggregate.n_scored == 2));
end

function test_pooled_rmse_weights_by_sample_count(testCase)
   % Pooling must weight a long site-year more than a short one, otherwise a
   % two-sample case would carry the same headline weight as a full season.

   results = syntheticResults();
   [per_case, aggregate] = ...
      icemodel.verification.ablationPerformanceMetrics(results);

   rows = per_case.scored & per_case.diagnostic == "model_melt_mwe" ...
      & per_case.density_kg_m3 == scoringDensity();
   weights = per_case.n_increments(rows);
   expected = sqrt(sum(weights .* per_case.rmse_mwe(rows) .^ 2) ...
      / sum(weights));
   pooled = aggregate.pooled_rmse_mwe( ...
      aggregate.diagnostic == "model_melt_mwe" ...
      & aggregate.density_kg_m3 == scoringDensity());
   testCase.verifyEqual(pooled, expected, 'AbsTol', 1e-12);

   % Guard against the assertion being vacuous: pooling must span more than
   % one case, the weights must differ, and the error must be nonzero.
   testCase.verifyGreaterThan(nnz(rows), 1);
   testCase.verifyGreaterThan(numel(unique(weights)), 1);
   testCase.verifyGreaterThan(max(per_case.rmse_mwe(rows)), 0);

   % An unweighted mean must NOT reproduce the pooled value, otherwise the
   % weighting is untested.
   unweighted = sqrt(mean(per_case.rmse_mwe(rows) .^ 2));
   testCase.verifyNotEqual(round(pooled, 12), round(unweighted, 12));
end

function density = scoringDensity()
   %SCORINGDENSITY Dense endpoint of the policy band, used by the fixtures.
   %
   % This function reads the value from the policy instead of stating a
   % literal. A change to the band then fails the tests, and it does not
   % leave every fixture filter matching zero rows.

   policy = icemodel.verification.namelists.promiceAblationPolicy();
   density = policy.effective_density_kg_m3(end);
end

function results = syntheticResults()
   %SYNTHETICRESULTS Build two scorable cases and one unscorable saved case.
   %
   % The two scorable cases differ in length and in error so that pooling by
   % sample count is observable; with a single case, or with a case whose
   % error is zero, any weighting rule produces the same pooled value.

   % Build from geometric lowering, since that is what the metrics convert.
   % Use the policy's dense endpoint so the expected values are exact.
   time = datetime(2019, 7, 1, 'TimeZone', 'UTC') + hours(0:4)';
   observation_lowering_m = (0:4)' * 0.1;
   ro_liq = icemodel.physicalConstant('ro_liq');
   observed = observation_lowering_m * scoringDensity() / ro_liq;

   % One diagnostic matches exactly, one is offset by a constant, one ablates
   % too slowly, and one ablates 20 percent too fast.
   model_melt_mwe = observed;
   model_runoff_mwe = observed + 0.05;
   model_solid_loss_mwe = observed * 0.5;
   model_ablation_proxy_mwe = observed * 1.2;
   observation_intact_mwe = observed;
   aligned = timetable(observation_intact_mwe, observation_lowering_m, ...
      model_melt_mwe, model_runoff_mwe, model_solid_loss_mwe, ...
      model_ablation_proxy_mwe, 'RowTimes', time);

   % A shorter, worse case. Its melt diagnostic ablates too fast, which is a
   % RATE error: the metrics score per-step increments, so a constant offset
   % would cancel and leave the pooled RMSE untestable.
   short_time = datetime(2020, 7, 1, 'TimeZone', 'UTC') + hours(0:2)';
   short_lowering_m = (0:2)' * 0.1;
   short_observed = short_lowering_m * scoringDensity() / ro_liq;
   observation_intact_mwe = short_observed;
   observation_lowering_m = short_lowering_m;
   model_melt_mwe = short_observed * 1.3;
   model_runoff_mwe = short_observed + 0.05;
   model_solid_loss_mwe = short_observed * 0.5;
   model_ablation_proxy_mwe = short_observed * 1.2;
   short_aligned = timetable(observation_intact_mwe, ...
      observation_lowering_m, model_melt_mwe, model_runoff_mwe, ...
      model_solid_loss_mwe, model_ablation_proxy_mwe, ...
      'RowTimes', short_time);

   results = struct( ...
      'case_id', {'perfect', 'short', 'unavailable'}, ...
      'site_id', {'KAN_M', 'KAN_U', 'KAN_L'}, ...
      'year', {2019, 2020, 2019}, ...
      'status', {'completed', 'completed', 'unavailable'}, ...
      'aligned', {aligned, short_aligned, timetable.empty(0, 0)});
end
