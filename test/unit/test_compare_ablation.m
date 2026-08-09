function tests = test_compare_ablation
   %TEST_COMPARE_ABLATION Verify the PROMICE ablation comparison.
   tests = functiontests(localfunctions);
end

function test_policy_is_fixed_and_reuses_budget_registry(testCase)
   % Policy values must remain predeclared and share production field names.
   policy = icemodel.verification.namelists.promiceAblationPolicy();

   testCase.verifyEqual(policy.numerical.atol_mwe, 1e-10);
   testCase.verifyEqual(policy.numerical.rtol, 1e-8);
   testCase.verifyEqual(policy.numerical.q_floor_mwe, 1e-12);
   testCase.verifyEqual(policy.scientific.signal_floor_mwe, 0.025);
   testCase.verifyEqual(policy.scientific.materiality_limit, 0.05);
   % The band is a documented modelling decision, not a magic constant, so
   % assert its structure rather than restating the endpoints.
   testCase.verifyEqual(policy.effective_density_kg_m3, [600, 870, 900]);
   testCase.verifyEqual(policy.effective_density_reference_kg_m3, 870);
   testCase.verifyTrue(ismember(policy.effective_density_reference_kg_m3, ...
      policy.effective_density_kg_m3));
   testCase.verifyLessThan(max(policy.effective_density_kg_m3), ...
      icemodel.physicalConstant('ro_ice'));
   testCase.verifyEqual(policy.model_output_cadence_seconds, 3600);
   readiness = ...
      icemodel.verification.namelists.promiceAblationReadiness();
   testCase.verifyEqual(policy.observation_cadence_seconds, ...
      readiness.observation_cadence_seconds);
   testCase.verifyEqual(policy.observation_field, readiness.target_field);
   testCase.verifyEqual(policy.support_flag_fields, ...
      readiness.support_flag_fields);
   testCase.verifyEqual(policy.direct_zero_flag_fields, ...
      readiness.direct_zero_flag_fields);
   testCase.verifyEqual(policy.datum_break_flag_fields, ...
      readiness.datum_break_flag_fields);
   testCase.verifyEqual(policy.ordinary_gap_flag_fields, ...
      readiness.ordinary_gap_flag_fields);
   testCase.verifyEqual(policy.metadata_only_flag_fields, ...
      readiness.metadata_only_flag_fields);
   testCase.verifyEqual(policy.required_observation_fields, ...
      [readiness.target_field, readiness.snow_variable, ...
      readiness.support_flag_fields]);
   testCase.verifyEqual(policy.required_model_fields, ...
      string([icemodel.namelists.budgetoutputs('all'), ...
      icemodel.namelists.cumulativeoutputs()]));
   testCase.verifyEqual(policy.snow_continuity_threshold_m, 0.05);
   testCase.verifyEqual(policy.ice_exposure_threshold_m, 0.01);
   testCase.verifySubstring(policy.model_interval_convention, "[Time,Time+dt)");
   testCase.verifySubstring(policy.rationale.density, "600--900");
   testCase.verifySubstring(policy.effective_density_role, ...
      "not an intact glacier-ice density");
end

function test_comparison_rebases_and_closes_budgets(testCase)
   % Bundled inputs must integrate every direct posting while excluding the
   % future interval that begins at the final observation timestamp.
   [observations, model, increment] = makeInputs();

   [summary, aligned, diagnostics, policy] = ...
      icemodel.verification.compareAblation(observations, model);

   testCase.verifyEqual(summary.classification, "within_materiality");
   testCase.verifyTrue(summary.physical_comparable);
   testCase.verifyEqual(summary.eligible_sample_count, 5);
   testCase.verifyEqual(height(aligned), 5);
   testCase.verifyEqual(aligned.observation_lowering_m, (0:0.1:0.4)', ...
      AbsTol=1e-12);
   testCase.verifyEqual(aligned.model_solid_loss_mwe, (0:4)' * increment, ...
      AbsTol=1e-12, RelTol=1e-12);
   testCase.verifyEqual(aligned.model_melt_mwe, (0:4)' * increment, ...
      AbsTol=1e-12, RelTol=1e-12);
   testCase.verifyEqual(aligned.model_runoff_mwe, (0:4)' * increment / 2, ...
      AbsTol=1e-12, RelTol=1e-12);
   testCase.verifyEqual(diagnostics.excluded.total_unique_excluded, 0);
   testCase.verifyTrue(all(diagnostics.identities.passed));
   testCase.verifyTrue(any(diagnostics.materiality.material));
   testCase.verifyEqual(diagnostics.top_deletion.count, 2);
   testCase.verifyEqual(diagnostics.top_deletion.height_m, 0.2);

   % Exported mass is reported separately from the quantized geometry, and the
   % aligned surface-loss series accumulates only the exported mass.
   testCase.verifyEqual(diagnostics.top_deletion.export_solid_mwe, 0.06, ...
      AbsTol=1e-12);
   testCase.verifyEqual(diagnostics.top_deletion.export_liquid_mwe, 0.02, ...
      AbsTol=1e-12);
   % The prefix is the state at each eligible time, so the two export rows in
   % the ledger appear at the first and third aligned steps.
   testCase.verifyEqual(aligned.model_surface_mass_loss_mwe, ...
      [0; 0.04; 0.04; 0.08; 0.08], AbsTol=1e-12);
   testCase.verifyFalse(any(diagnostics.effective_density.rigorous_bound));
   testCase.verifyEqual(summary.policy_version, policy.version);
end

function test_signed_density_band_orders_bounds_and_preserves_endpoints(testCase)
   % Negative cumulative lowering reverses which density endpoint is the
   % numeric lower value without changing the density-specific sensitivity.
   [observations, model] = makeInputs();
   observations.data.ablation = 7 - (0:0.1:0.4)';

   [~, aligned, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   [~, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   policy = icemodel.verification.namelists.promiceAblationPolicy();
   endpoint_values = aligned.observation_lowering_m ...
      .* policy.effective_density_kg_m3 ./ ro_liq;
   testCase.verifyEqual(aligned.observation_lower_mwe, ...
      min(endpoint_values, [], 2), AbsTol=1e-14)
   testCase.verifyEqual(aligned.observation_upper_mwe, ...
      max(endpoint_values, [], 2), AbsTol=1e-14)
   testCase.verifyTrue(all(aligned.observation_lower_mwe ...
      <= aligned.observation_upper_mwe))
   testCase.verifyEqual( ...
      diagnostics.effective_density.observation_sensitivity_mwe, ...
      -0.4 * icemodel.verification.namelists.promiceAblationPolicy() ...
      .effective_density_kg_m3(:) ./ ro_liq, AbsTol=1e-14)
end

function test_trace_snow_is_censored_without_splitting_window(testCase)
   % Trace snow above the exposed-ice threshold removes only that posting.
   [observations, model] = makeInputs();
   observations.data.snow_depth(3) = 0.02;

   [summary, aligned, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   testCase.verifyEqual(summary.window_start, observations.data.Time(1));
   testCase.verifyEqual(summary.window_end, observations.data.Time(end));
   testCase.verifyEqual(summary.eligible_sample_count, 4);
   testCase.verifyFalse(any(aligned.Time == observations.data.Time(3)));
   testCase.verifyEqual(diagnostics.excluded.snow_censored, 1);
end

function test_negative_snow_is_unknown_not_censored(testCase)
   % Negative snow is invalid support rather than an exposed-ice observation.
   [observations, model] = makeInputs();
   observations.data.snow_depth(3) = -0.01;

   [summary, aligned, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   testCase.verifyEqual(summary.eligible_sample_count, 4);
   testCase.verifyFalse(any(aligned.Time == observations.data.Time(3)));
   testCase.verifyEqual(diagnostics.excluded.unknown_snow_depth, 1);
   testCase.verifyEqual(diagnostics.excluded.snow_censored, 0);
end

function test_shared_snow_classifier_partitions_boundary_values(testCase)
   % Negative/nonfinite, exposed, and censored values are mutually exclusive.
   snow = [-0.01; NaN; Inf; 0; 0.01; 0.02];
   [ice_exposed, snow_censored, unknown_snow] = ...
      icemodel.verification.helpers.classifySnowDepth(snow, 0.01);

   testCase.verifyEqual(ice_exposed, ...
      [false; false; false; true; true; false]);
   testCase.verifyEqual(snow_censored, ...
      [false; false; false; false; false; true]);
   testCase.verifyEqual(unknown_snow, ...
      [true; true; true; false; false; false]);
   testCase.verifyTrue(all( ...
      double(ice_exposed) + double(snow_censored) ...
      + double(unknown_snow) == 1));
end

function test_leading_flags_are_excluded_without_bridging(testCase)
   % Flagged postings may trim the leading edge when the remaining direct
   % support is contiguous; no flagged posting lies between its endpoints.
   [observations, model, increment] = makeInputs();
   observations.data.surface_height_flag(1) = 1;
   observations.data.station_transition_flag(2) = 1;
   observations.data.step_detected_flag(3) = 1;
   observations.data.step_correctable_flag(3) = 1;

   [summary, aligned, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   testCase.verifyEqual(summary.eligible_sample_count, 2);
   testCase.verifyEqual(height(aligned), 2);
   testCase.verifyEqual(aligned.model_solid_loss_mwe, [0; increment], ...
      AbsTol=1e-12);
   testCase.verifyEqual(diagnostics.excluded.gap_bridged, 1);
   testCase.verifyEqual(diagnostics.excluded.station_transition, 1);
   testCase.verifyEqual(diagnostics.excluded.unresolved_step, 1);
   testCase.verifyEqual( ...
      diagnostics.excluded.step_correctable_but_unresolved, 1);
end

function test_a_negative_correctable_flag_still_counts_as_flagged(testCase)
   % step_correctable_but_unresolved used to test the raw flag column with
   % > 0, so a malformed negative posting read as unflagged here while the
   % readiness writer read the same posting as flagged and the two ledgers
   % disagreed for the same site-year. classifyObservationSupport now governs
   % both with one finite-and-nonzero rule, and a negative posting counts.

   [observations, model] = makeInputs();
   observations.data.surface_height_flag(1) = 1;
   observations.data.station_transition_flag(2) = 1;
   observations.data.step_detected_flag(3) = 1;
   observations.data.step_correctable_flag(3) = -1;

   [~, ~, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   testCase.verifyEqual(diagnostics.excluded.unresolved_step, 1);
   testCase.verifyEqual( ...
      diagnostics.excluded.step_correctable_but_unresolved, 1);
end

function test_a_nonfinite_correctable_flag_does_not_count_as_flagged(testCase)
   % The other half of the same rule: NaN is missing data, not a raised flag,
   % so the unresolved step is counted but not as correctable.

   [observations, model] = makeInputs();
   observations.data.surface_height_flag(1) = 1;
   observations.data.station_transition_flag(2) = 1;
   observations.data.step_detected_flag(3) = 1;
   observations.data.step_correctable_flag(3) = NaN;

   [~, ~, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   testCase.verifyEqual(diagnostics.excluded.unresolved_step, 1);
   testCase.verifyEqual( ...
      diagnostics.excluded.step_correctable_but_unresolved, 0);
end

function test_separate_domain_terms_prevent_cancellation(testCase)
   % Equal cloned-bottom addition and merge export close aggregate remeshing,
   % but each material solid-domain term must retain its own scientific scenario.
   [observations, model] = makeInputs();
   model.data.mass_budget_cloned_bottom_solid_mwe(2) = 0.1;
   model.data.mass_budget_merge_export_solid_mwe(2) = 0.1;
   model.data.mass_budget_cloned_bottom_solid_gross_mwe(2) = 0.1;
   model.data.mass_budget_merge_export_solid_gross_mwe(2) = 0.1;

   [summary, ~, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   testCase.verifyTrue(summary.physical_comparable);
   testCase.verifyEqual(summary.classification, "non_identifiable");
   names = diagnostics.scenarios.scenario;
   testCase.verifyTrue(any(names == "merge_delete_solid"));
   testCase.verifyTrue(any(names == "cloned_bottom_solid"));
   separate = ismember(names, ["merge_delete_solid", "cloned_bottom_solid"]);
   testCase.verifyTrue(all(diagnostics.scenarios.material(separate)));
   testCase.verifyTrue(all(diagnostics.scenarios.changes_classification( ...
      separate)));
end

function test_native_gross_survives_hourly_cancellation(testCase)
   % Native event gross must remain visible when every hourly signed
   % remesh term has cancelled to zero.
   [observations, model] = makeInputs();
   model.data.mass_budget_remesh_solid_gross_mwe(2) = 0.2;
   model.data.mass_budget_cloned_bottom_solid_gross_mwe(2) = 0.1;
   model.data.mass_budget_merge_export_solid_gross_mwe(2) = 0.1;

   [summary, ~, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   testCase.verifyTrue(summary.physical_comparable);
   material = diagnostics.materiality( ...
      diagnostics.materiality.channel == "remesh_solid", :);
   testCase.verifyEqual(material.signed_net_mwe, 0);
   testCase.verifyEqual(material.gross_mwe, 0.2);
   testCase.verifyTrue(material.material);
   identity = diagnostics.identities( ...
      diagnostics.identities.identity == "remesh_solid", :);
   testCase.verifyEqual(identity.normalization, 0.2, AbsTol=1e-12);
end

function test_liquid_storage_remains_diagnostic_only(testCase)
   % Retained meltwater is material state evidence, but without a demonstrated
   % observation operator it cannot arithmetically correct solid loss.
   [observations, model] = makeInputs();
   phase_liquid = model.data.mass_budget_phase_liquid_mwe;
   model.data.mass_budget_remesh_liquid_mwe(:) = 0;
   model.data.mass_budget_merge_export_liquid_mwe(:) = 0;
   model.data.mass_budget_remesh_liquid_gross_mwe(:) = 0;
   model.data.mass_budget_merge_export_liquid_gross_mwe(:) = 0;
   model.data.mass_budget_liquid_start_mwe = ...
      1 + [0; cumsum(phase_liquid(1:end - 1))];
   model.data.mass_budget_liquid_end_mwe = ...
      model.data.mass_budget_liquid_start_mwe + phase_liquid;
   model.data.mass_budget_liquid_storage_gross_mwe = ...
      abs(phase_liquid);

   [summary, ~, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   liquid = diagnostics.materiality( ...
      diagnostics.materiality.channel == "endpoint_liquid_storage", :);
   testCase.verifyTrue(liquid.material);
   testCase.verifyEqual(summary.classification, "within_materiality");
   testCase.verifyFalse(any(contains( ...
      diagnostics.scenarios.scenario, "liquid")));
end

function test_condensation_overflow_remains_diagnostic_only(testCase)
   % Overflow is liquid outflow and vapor-energy accounting, not solid loss.
   [observations, model] = makeInputs();
   [Lv, ro_liq] = icemodel.physicalConstant('Lv', 'ro_liq');
   overflow = 0.1;
   model.data.mass_budget_condensation_overflow_mwe(2) = overflow;
   model.data.mass_budget_condensation_overflow_gross_mwe(2) = overflow;
   model.data.mass_budget_vapor_potential_j_m2(2) = ...
      ro_liq * Lv * overflow;
   model.data.mass_budget_vapor_potential_gross_j_m2(2) = ...
      ro_liq * Lv * overflow;

   [summary, ~, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   overflow_row = diagnostics.materiality( ...
      diagnostics.materiality.channel == "condensation_overflow", :);
   testCase.verifyTrue(overflow_row.material);
   testCase.verifyEqual(summary.classification, "within_materiality");
   testCase.verifyFalse(any( ...
      diagnostics.scenarios.scenario == "condensation_overflow"));
end

function test_physical_dry_vapor_closes_fixed_mwe_identity(testCase)
   % Dry-vapor postings must close on the physical intrinsic-density basis used
   % by both production solver kernels and fixed-reference MWE diagnostics.

   [observations, model] = makeInputs();
   [Ls, Lv, ro_ice, ro_liq] = ...
      icemodel.physicalConstant('Ls', 'Lv', 'ro_ice', 'ro_liq');
   d_pevp = -0.01;
   d_psbl = d_pevp * (Lv * ro_liq) / (Ls * ro_ice);
   vapor_solid_mwe = ro_ice / ro_liq * d_psbl;
   event_row = 2;

   % Post one dry-sublimation event and carry its storage change through every
   % later checkpoint so both vapor energy and solid storage remain closed.
   model.data.mass_budget_vapor_solid_mwe(event_row) = vapor_solid_mwe;
   model.data.mass_budget_vapor_solid_gross_mwe(event_row) = ...
      abs(vapor_solid_mwe);
   model.data.mass_budget_vapor_potential_j_m2(event_row) = ...
      ro_liq * Lv * d_pevp;
   model.data.mass_budget_vapor_potential_gross_j_m2(event_row) = ...
      abs(ro_liq * Lv * d_pevp);
   model.data.mass_budget_solid_end_mwe(event_row:end) = ...
      model.data.mass_budget_solid_end_mwe(event_row:end) + vapor_solid_mwe;
   model.data.mass_budget_solid_start_mwe(event_row + 1:end) = ...
      model.data.mass_budget_solid_start_mwe(event_row + 1:end) ...
      + vapor_solid_mwe;
   model.data.mass_budget_solid_storage_gross_mwe(event_row) = ...
      model.data.mass_budget_solid_storage_gross_mwe(event_row) ...
      + abs(vapor_solid_mwe);

   [summary, ~, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);
   vapor = diagnostics.identities( ...
      diagnostics.identities.identity == "vapor_energy", :);
   testCase.verifyTrue(vapor.passed);
   testCase.verifyEqual(vapor.residual, 0, 'AbsTol', 1e-8);
   testCase.verifyTrue(summary.physical_comparable);
end

function test_exact_window_normalization_ignores_storage_path_length(testCase)
   % The closure scale is |Delta q| plus physical-flux gross. Storage
   % path length is a separate materiality diagnostic and must not inflate it.
   [observations, model, increment] = makeInputs();
   model.data.mass_budget_solid_storage_gross_mwe(:) = 1e6;

   [~, ~, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   identity = diagnostics.identities( ...
      diagnostics.identities.identity == "solid_storage", :);
   testCase.verifyEqual(identity.normalization, 8 * increment, ...
      AbsTol=1e-12, RelTol=1e-12);
end

function test_forcing_step_residuals_cannot_cancel(testCase)
   % Equal and opposite row defects retain a zero window residual but must fail
   % the step gate instead of producing a physically comparable window.
   [observations, model] = makeInputs();
   model.data.mass_budget_solid_end_mwe(2) = ...
      model.data.mass_budget_solid_end_mwe(2) + 1e-4;
   model.data.mass_budget_solid_end_mwe(3) = ...
      model.data.mass_budget_solid_end_mwe(3) - 1e-4;

   [summary, ~, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   identity = diagnostics.identities( ...
      diagnostics.identities.identity == "solid_storage", :);
   testCase.verifyEqual(identity.residual, 0, AbsTol=1e-12);
   testCase.verifyTrue(identity.window_passed);
   testCase.verifyFalse(identity.step_passed);
   testCase.verifyEqual(identity.failed_step_count, 2);
   testCase.verifyFalse(summary.physical_comparable);
end

function test_phase_mass_uses_combined_forcing_step_increment(testCase)
   % Opposing phase residuals across rows must not cancel. Compensating remesh
   % terms keep the other storage and domain identities closed for isolation.
   [observations, model] = makeInputs();
   phase_error = [0; 1e-4; -1e-4; 0; 0; 0];
   model.data.mass_budget_phase_liquid_mwe = ...
      model.data.mass_budget_phase_liquid_mwe + phase_error;
   model.data.mass_budget_remesh_liquid_mwe = ...
      model.data.mass_budget_remesh_liquid_mwe - phase_error;
   model.data.mass_budget_merge_export_liquid_mwe = ...
      -model.data.mass_budget_remesh_liquid_mwe;
   model.data.mass_budget_phase_liquid_gross_mwe = ...
      abs(model.data.mass_budget_phase_liquid_mwe);
   model.data.mass_budget_remesh_liquid_gross_mwe = ...
      abs(model.data.mass_budget_remesh_liquid_mwe);
   model.data.mass_budget_merge_export_liquid_gross_mwe = ...
      abs(model.data.mass_budget_merge_export_liquid_mwe);

   [~, ~, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   identity = diagnostics.identities( ...
      diagnostics.identities.identity == "phase_mass", :);
   testCase.verifyEqual(identity.residual, 0, AbsTol=1e-12);
   testCase.verifyTrue(identity.window_passed);
   testCase.verifyFalse(identity.step_passed);
   testCase.verifyEqual(identity.failed_step_count, 2);
end

function test_unverified_endpoint_scenarios_never_govern(testCase)
   % Adversarial caller values remain visible but cannot supersede the central
   % classification without validated observation-artifact provenance.
   [observations, model] = makeInputs();

   [summary, ~, diagnostics] = icemodel.verification.compareAblation( ...
      observations, model, endpoint_deficit_kg_m2=[-100, 100]);

   testCase.verifyEqual(summary.base_classification, "within_materiality");
   testCase.verifyEqual(summary.classification, "within_materiality");
   testCase.verifyEqual(summary.endpoint_deficit_scenario_count, 2);
   testCase.verifyEmpty(summary.non_identifiable_reasons);
   endpoint_rows = diagnostics.scenarios.role == "endpoint";
   testCase.verifyFalse(any(diagnostics.scenarios.credible(endpoint_rows)));
   testCase.verifyTrue(any(diagnostics.scenarios.changes_classification( ...
      endpoint_rows)));
end

function test_failed_identity_blocks_physical_comparison(testCase)
   % A storage residual above gross tolerance must take precedence over a
   % visually matching cumulative series.
   [observations, model] = makeInputs();
   model.data.mass_budget_solid_end_mwe(end - 1) = ...
      model.data.mass_budget_solid_end_mwe(end - 1) + 1e-4;

   [summary, ~, diagnostics] = ...
      icemodel.verification.compareAblation(observations, model);

   testCase.verifyEqual(summary.classification, ...
      "not_physically_comparable");
   testCase.verifyFalse(summary.physical_comparable);
   testCase.verifyFalse(diagnostics.identities.passed(1));
end

function test_unmatched_observation_has_directional_bias(testCase)
   % A closed comparison outside the fixed five-percent limit reports bias when
   % no credible scenario changes that conclusion.
   [observations, model] = makeInputs();
   observations.data.ablation(end) = observations.data.ablation(end) + 0.2;

   summary = icemodel.verification.compareAblation(observations, model);

   testCase.verifyEqual(summary.base_classification, "model_low");
   testCase.verifyEqual(summary.classification, "model_low");
end

function test_raw_timetables_honor_explicit_window(testCase)
   % Raw timetables and both explicit endpoint arguments share the bundle path.
   [observations, model, increment] = makeInputs();
   t = observations.data.Time;

   [summary, aligned] = icemodel.verification.compareAblation( ...
      observations.data, model.data, window_start=t(2), window_end=t(5));

   testCase.verifyEqual(summary.window_start, t(2));
   testCase.verifyEqual(summary.window_end, t(5));
   testCase.verifyEqual(aligned.model_solid_loss_mwe(end), 3 * increment, ...
      AbsTol=1e-12);
end

function test_required_observation_field_has_stable_error(testCase)
   % Incomplete staged observations must fail before temporal selection.
   [observations, model] = makeInputs();
   observations.data = removevars(observations.data, 'surface_height_flag');

   testCase.verifyError(@() icemodel.verification.compareAblation( ...
      observations, model), ...
      'icemodel:verification:compareAblation:missingObservationField');
end

function test_required_model_field_has_stable_error(testCase)
   % Omitting one closure term must prevent physical comparability claims.
   [observations, model] = makeInputs();
   model.data = removevars(model.data, 'mass_budget_vapor_solid_mwe');

   testCase.verifyError(@() icemodel.verification.compareAblation( ...
      observations, model), ...
      'icemodel:verification:compareAblation:missingModelField');
end

function test_missing_internal_model_row_has_stable_error(testCase)
   % A missing interval-start row must not compress elapsed time into the next
   % surviving increment or appear to retain complete physical support.
   [observations, model] = makeInputs();
   model.data(3, :) = [];

   testCase.verifyError(@() icemodel.verification.compareAblation( ...
      observations, model), ...
      'icemodel:verification:compareAblation:noncontiguousModelSupport');
end

function test_shifted_model_grid_has_stable_error(testCase)
   % Uniform spacing is insufficient when interval starts are shifted off the
   % native UTC output boundaries.
   [observations, model] = makeInputs();
   observations.data.Time = observations.data.Time + minutes(30);
   model.data.Time = model.data.Time + minutes(30);

   testCase.verifyError(@() icemodel.verification.compareAblation( ...
      observations, model), ...
      'icemodel:verification:compareAblation:offGridModelTime');
end

function test_missing_internal_observation_preserves_cumulative_total(testCase)
   % Exact cumulative endpoints remain comparable across a missing plotting
   % marker; no increment or rate is inferred across the omitted posting.
   [observations, model, increment] = makeInputs();
   observations.data(3, :) = [];

   [summary, aligned] = ...
      icemodel.verification.compareAblation(observations, model);

   testCase.verifyEqual(height(aligned), 4);
   testCase.verifyEqual(aligned.model_solid_loss_mwe, ...
      [0; increment; 3 * increment; 4 * increment], AbsTol=1e-12);
   testCase.verifyEqual(summary.model_solid_loss_mwe, 4 * increment, ...
      AbsTol=1e-12);
   testCase.verifyEqual(summary.observation_lowering_m, 0.4, AbsTol=1e-12);
end

function test_gap_flagged_internal_observation_is_marker_only(testCase)
   % An ordinary surface-height gap removes one marker but does not invalidate
   % later cumulative levels on the same datum.
   [observations, model, increment] = makeInputs();
   observations.data.surface_height_flag(3) = 1;

   [summary, aligned] = ...
      icemodel.verification.compareAblation(observations, model);

   testCase.verifyEqual(height(aligned), 4);
   testCase.verifyEqual(summary.model_solid_loss_mwe, 4 * increment, ...
      AbsTol=1e-12);
   testCase.verifyEqual(summary.observation_lowering_m, 0.4, AbsTol=1e-12);
end

function test_unknown_internal_gap_flags_are_marker_only(testCase)
   % Unknown surface-gap or correction metadata cannot be an endpoint marker,
   % but neither changes the station datum by itself.
   fields = ["surface_height_flag", "step_correctable_flag"];
   values = [NaN, Inf];
   for k = 1:numel(fields)
      [observations, model, increment] = makeInputs();
      flag = observations.data.(fields(k));
      flag(3) = values(k);
      observations.data.(fields(k)) = flag;
      [summary, aligned] = ...
         icemodel.verification.compareAblation(observations, model);
      testCase.verifyEqual(height(aligned), 4);
      testCase.verifyEqual(summary.model_solid_loss_mwe, 4 * increment, ...
         AbsTol=1e-12);
   end
end

function test_internal_datum_breaks_have_stable_error(testCase)
   % Finite events and nonfinite unknown values in either datum flag invalidate
   % every cumulative window that crosses their posting.
   fields = ["station_transition_flag", "step_detected_flag", ...
      "station_transition_flag", "step_detected_flag"];
   values = [1, 1, NaN, Inf];
   for k = 1:numel(fields)
      [observations, model] = makeInputs();
      flag = observations.data.(fields(k));
      flag(3) = values(k);
      observations.data.(fields(k)) = flag;
      testCase.verifyError(@() icemodel.verification.compareAblation( ...
         observations, model), ...
         'icemodel:verification:compareAblation:observationDatumBreak');
   end
end

function test_bad_input_and_nonnumeric_fields_have_stable_errors(testCase)
   % The public boundary rejects unsupported containers and nonnumeric payloads.
   [observations, model] = makeInputs();
   testCase.verifyError(@() icemodel.verification.compareAblation( ...
      42, model), 'icemodel:verification:compareAblation:badInput');

   observations.data.ablation = cellstr(string(observations.data.ablation));
   testCase.verifyError(@() icemodel.verification.compareAblation( ...
      observations, model), ...
      'icemodel:verification:compareAblation:nonNumericField');
end

function test_time_window_and_support_errors_are_stable(testCase)
   % Duplicate times, disjoint bounds, and one eligible point each fail with a
   % distinct error id so readiness code can classify the reason.
   [observations, model] = makeInputs();
   duplicate = observations;
   duplicate.data.Time(2) = duplicate.data.Time(1);
   testCase.verifyError(@() icemodel.verification.compareAblation( ...
      duplicate, model), ...
      'icemodel:verification:compareAblation:duplicateTime');

   testCase.verifyError(@() icemodel.verification.compareAblation( ...
      observations, model, ...
      window_start=datetime(2030, 1, 1, 'TimeZone', 'UTC'), ...
      window_end=datetime(2030, 1, 2, 'TimeZone', 'UTC')), ...
      'icemodel:verification:compareAblation:invalidWindow');

   observations.data.surface_height_flag(2:end) = 1;
   testCase.verifyError(@() icemodel.verification.compareAblation( ...
      observations, model), ...
      'icemodel:verification:compareAblation:insufficientSupport');
end

function test_nonfinite_ledger_and_endpoint_errors_are_stable(testCase)
   % A hidden nonfinite model increment and a nonfinite endpoint scenario must
   % be rejected rather than dropped from the scientific accounting.
   [observations, model] = makeInputs();
   model.data.mass_budget_unapplied_vapor_j_m2(3) = NaN;
   testCase.verifyError(@() icemodel.verification.compareAblation( ...
      observations, model), ...
      'icemodel:verification:compareAblation:nonfiniteModelWindow');

   [observations, model] = makeInputs();
   testCase.verifyError(@() icemodel.verification.compareAblation( ...
      observations, model, endpoint_deficit_kg_m2=NaN), ...
      'icemodel:verification:compareAblation:badEndpointDeficit');
end

function test_negative_gross_has_stable_error(testCase)
   % Gross is an absolute accepted-event total and cannot be negative.
   [observations, model] = makeInputs();
   model.data.mass_budget_phase_solid_gross_mwe(2) = -1;

   testCase.verifyError(@() icemodel.verification.compareAblation( ...
      observations, model), ...
      'icemodel:verification:compareAblation:negativeGross');
end

function [observations, model, increment] = makeInputs()
   %MAKEINPUTS Build a closed synthetic five-hour ablation comparison.
   time = datetime(2026, 7, 1, 'TimeZone', 'UTC') + hours((0:4)');
   n = numel(time);
   [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   increment = 0.1 * ro_ice / ro_liq;

   % The observation carries an arbitrary installation offset to prove that
   % independent rebasing, rather than absolute height, controls the result.
   ablation = 7 + (0:0.1:0.4)';
   observations.data = timetable(ablation, zeros(n, 1), zeros(n, 1), ...
      zeros(n, 1), zeros(n, 1), zeros(n, 1), 'RowTimes', time, ...
      'VariableNames', {'ablation', 'snow_depth', 'surface_height_flag', ...
      'station_transition_flag', 'step_detected_flag', ...
      'step_correctable_flag'});
   observations.format = 'timeseries';

   % Initialize every required diagnostic from the canonical registry, then
   % populate one physically closed interval-start melt/export ledger. The final
   % row has a large future-interval increment that must not enter [t0,t4).
   model_time = [time(1) - hours(1); time];
   n_model = numel(model_time);
   model.data = timetable('RowTimes', model_time);
   fields = icemodel.namelists.budgetoutputs('all');
   for k = 1:numel(fields)
      model.data.(fields{k}) = zeros(n_model, 1);
   end
   solid_loss = [0; increment * ones(n - 1, 1); 100 * increment];
   phase_solid = -solid_loss;
   phase_liquid = -phase_solid;
   remesh_liquid = -phase_liquid;
   solid_start = 10 - [0; cumsum(solid_loss(1:end - 1))];
   solid_end = solid_start - solid_loss;
   model.data.mass_budget_solid_start_mwe = solid_start;
   model.data.mass_budget_solid_end_mwe = solid_end;
   model.data.mass_budget_liquid_start_mwe(:) = 1;
   model.data.mass_budget_liquid_end_mwe(:) = 1;
   model.data.mass_budget_phase_solid_mwe = phase_solid;
   model.data.mass_budget_phase_liquid_mwe = phase_liquid;
   model.data.mass_budget_remesh_liquid_mwe = remesh_liquid;
   model.data.mass_budget_merge_export_liquid_mwe = -remesh_liquid;
   model.data.mass_budget_solid_storage_gross_mwe = solid_loss;
   model.data.mass_budget_phase_solid_gross_mwe = abs(phase_solid);
   model.data.mass_budget_phase_liquid_gross_mwe = abs(phase_liquid);
   model.data.mass_budget_remesh_liquid_gross_mwe = abs(remesh_liquid);
   model.data.mass_budget_merge_export_liquid_gross_mwe = ...
      abs(remesh_liquid);
   model.data.mass_budget_top_deletion_count([2, 4]) = 1;
   model.data.mass_budget_top_deletion_height_m([2, 4]) = 0.1;
   model.data.mass_budget_top_deletion_count(end) = 99;
   model.data.mass_budget_top_deletion_height_m(end) = 9.9;
   % Exported mass is deliberately unequal to the quantized cell height so a
   % comparator that confused the two would fail.
   model.data.mass_budget_top_export_solid_mwe([2, 4]) = 0.03;
   model.data.mass_budget_top_export_liquid_mwe([2, 4]) = 0.01;
   model.data.mass_budget_top_export_solid_mwe(end) = 5.5;
   model.data.mass_budget_interior_merge_count(3) = 1;
   model.data.melt = (0:n_model - 1)' * increment;
   model.data.runoff = (0:n_model - 1)' * increment / 2;
   model.data.freeze = (0:n_model - 1)' * increment / 4;
   model.data.dlayer = (0:n_model - 1)' * 0.02;
   model.format = 'timeseries';
end
