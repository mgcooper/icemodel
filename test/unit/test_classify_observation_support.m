function tests = test_classify_observation_support
   %TEST_CLASSIFY_OBSERVATION_SUPPORT Flag-support rules for observation rows.
   %
   % icemodel.verification.helpers.classifyObservationSupport judges whether a
   % PROMICE observation row is supported; the comparator, the readiness
   % writer, the evaluation runner, and the report builder all call it. These
   % tests pin the rules at flag-list lengths above one, since a per-row
   % collapse and a per-matrix count only agree when every list has exactly
   % one member.

   tests = functiontests(localfunctions);
end

function policy = fixturePolicy()
   %FIXTUREPOLICY Build a policy with multi-member flag lists.
   %
   % The real policy lists mostly have one member today. Testing at one member
   % cannot distinguish a per-row collapse from a per-matrix count, so the
   % fixture uses two.

   policy = struct( ...
      'support_flag_fields', ["s1", "s2", "g1", "g2", "m1", "m2"], ...
      'direct_zero_flag_fields', ["s1", "s2"], ...
      'datum_break_flag_fields', ["d1", "d2"], ...
      'ordinary_gap_flag_fields', ["g1", "g2"], ...
      'metadata_only_flag_fields', ["m1", "m2"], ...
      'station_transition_flag_field', "d1", ...
      'unresolved_step_flag_field', "d2");
end

function fields = fixtureFields()
   %FIXTUREFIELDS Column names matching fixturePolicy, target column first.

   fields = ["ablation", "s1", "s2", "g1", "g2", "m1", "m2", "d1", "d2"];
end

function support = classify(values)
   %CLASSIFY Run the helper over the fixture policy and column order.

   support = icemodel.verification.helpers.classifyObservationSupport( ...
      values, fixtureFields(), "ablation", fixturePolicy());
end

function test_clean_row_is_directly_supported(testCase)
   % A finite target with every flag finite and zero is direct support.

   returned = classify([1.5, 0, 0, 0, 0, 0, 0, 0, 0]);
   verifyTrue(testCase, returned.target_finite)
   verifyTrue(testCase, returned.quality_finite)
   verifyTrue(testCase, returned.direct_flags_zero)
   verifyTrue(testCase, returned.flag_clean)
   verifyTrue(testCase, returned.datum_intact)
   verifyFalse(testCase, returned.gap_flagged)
   verifyFalse(testCase, returned.metadata_flagged)
   verifyFalse(testCase, returned.station_transition)
   verifyFalse(testCase, returned.unresolved_step)
end

function test_a_nonfinite_target_is_not_flag_clean(testCase)
   % flag_clean must carry the target's finiteness, not just the flags, so a
   % row with a missing target is never reported as flag_clean.

   returned = classify([NaN, 0, 0, 0, 0, 0, 0, 0, 0]);
   verifyFalse(testCase, returned.target_finite)
   verifyFalse(testCase, returned.flag_clean)

   % The flags themselves are still clean; only the target is missing.
   verifyTrue(testCase, returned.quality_finite)
   verifyTrue(testCase, returned.direct_flags_zero)
end

function test_a_nonfinite_direct_flag_is_not_a_zero_flag(testCase)
   % NaN == 0 is already false, so all(values == 0) alone rejects this row.
   % The isfinite term matters for the set-flag masks below, where without it
   % a NaN would read as a raised flag. Testing it here keeps
   % direct_flags_zero correct on its own rather than via quality_finite.

   returned = classify([1.5, NaN, 0, 0, 0, 0, 0, 0, 0]);
   verifyFalse(testCase, returned.quality_finite)
   verifyFalse(testCase, returned.direct_flags_zero)
   verifyFalse(testCase, returned.flag_clean)
end

function test_a_set_flag_in_the_second_list_member_is_seen(testCase)
   % A flag set only in the second member of a multi-member list must still
   % register as raised.

   gap_second = classify([1.5, 0, 0, 0, 1, 0, 0, 0, 0]);
   verifyTrue(testCase, gap_second.gap_flagged)

   metadata_second = classify([1.5, 0, 0, 0, 0, 0, 1, 0, 0]);
   verifyTrue(testCase, metadata_second.metadata_flagged)
end

function test_nonfinite_is_not_a_set_flag_but_negative_is(testCase)
   % A NaN or an Inf is missing data, not a raised flag. Counting it as one
   % inflates every exclusion count these masks feed.

   nonfinite = classify([1.5, 0, 0, NaN, 0, NaN, 0, NaN, NaN]);
   verifyFalse(testCase, nonfinite.gap_flagged)
   verifyFalse(testCase, nonfinite.metadata_flagged)
   verifyFalse(testCase, nonfinite.station_transition)
   verifyFalse(testCase, nonfinite.unresolved_step)

   % A finite negative posting is malformed rather than missing, and counts
   % as set; excluding the row is the conservative direction.
   negative = classify([1.5, 0, 0, -1, 0, 0, 0, 0, 0]);
   verifyTrue(testCase, negative.gap_flagged)
end

function test_datum_intact_requires_every_datum_flag_finite_and_zero(testCase)
   % datum_intact is the only mask that is true by absence of a break, so a
   % nonfinite datum flag has to break it rather than be ignored.

   set_flag = classify([1.5, 0, 0, 0, 0, 0, 0, 1, 0]);
   verifyFalse(testCase, set_flag.datum_intact)
   verifyTrue(testCase, set_flag.station_transition)

   nonfinite = classify([1.5, 0, 0, 0, 0, 0, 0, NaN, 0]);
   verifyFalse(testCase, nonfinite.datum_intact)

   % Nonfinite is not a raised transition flag, only a broken datum.
   verifyFalse(testCase, nonfinite.station_transition)
end

function test_the_two_datum_roles_are_selected_by_name(testCase)
   % station_transition and unresolved_step come from differently ordered
   % policy lists. Positional indexing would swap them; this row sets only the
   % step flag and must not report a station transition.

   returned = classify([1.5, 0, 0, 0, 0, 0, 0, 0, 1]);
   verifyTrue(testCase, returned.unresolved_step)
   verifyFalse(testCase, returned.station_transition)
end

function test_masks_are_row_shaped_for_a_multi_row_matrix(testCase)
   % Every consumer indexes these masks against timetable rows, so the shape
   % has to follow the input rather than collapse to a scalar.

   values = [ ...
      1.5, 0, 0, 0, 0, 0, 0, 0, 0; ...
      NaN, 0, 0, 0, 0, 0, 0, 0, 0; ...
      2.5, 0, 0, 1, 0, 0, 0, 0, 0];
   returned = classify(values);
   expected = [true; false; true];
   verifyEqual(testCase, returned.target_finite, expected)
   verifyEqual(testCase, returned.gap_flagged, [false; false; true])
end

function test_a_missing_column_is_rejected(testCase)
   % A missing column must raise, rather than let an all() or any() collapse
   % compute over fewer columns than the caller asked for.

   verifyError(testCase, @() ...
      icemodel.verification.helpers.classifyObservationSupport( ...
      zeros(1, 3), ["ablation", "s1", "s2"], "ablation", fixturePolicy()), ...
      'icemodel:verification:classifyObservationSupport:missingField')
end

function test_the_real_policy_flag_lists_resolve(testCase)
   % The fixture uses invented field names, so pin that the shipped policy's
   % own lists still name columns the helper can find.

   policy = icemodel.verification.namelists.promiceAblationReadiness();
   fields = icemodel.verification.helpers.observationSupportFields( ...
      policy.target_field, policy);
   returned = icemodel.verification.helpers.classifyObservationSupport( ...
      zeros(2, numel(fields)), fields, policy.target_field, policy);
   verifyTrue(testCase, all(returned.direct_flags_zero))
   verifyTrue(testCase, all(returned.datum_intact))
   verifyFalse(testCase, any(returned.gap_flagged))
end
