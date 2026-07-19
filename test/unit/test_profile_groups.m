function tests = test_profile_groups
   %TEST_PROFILE_GROUPS Verify profile identity/date grouping contracts.
   tests = functiontests(localfunctions);
end

function test_empty_table_has_no_groups(testCase)
   % Empty profile tables must not invent an ungrouped profile.
   groups = icemodel.verification.helpers.profileGroups(table());
   testCase.verifyEmpty(groups)
end

function test_mar_profile_id_and_date_define_groups(testCase)
   % MAR rows with repeated depths remain separated by explicit id/date.
   t1 = datetime(2019, 5, 10, 12, 0, 0, 'TimeZone', 'UTC');
   t2 = datetime(2019, 6, 10, 12, 0, 0, 'TimeZone', 'UTC');
   profiles = table(["mar_a"; "mar_a"; "mar_b"; "mar_b"], ...
      [t1; t1; t2; t2], [0; 1; 0; 1], [350; 500; 360; 510], ...
      'VariableNames', {'profile_id', 'datetime', 'depth', 'density'});

   groups = icemodel.verification.helpers.profileGroups(profiles);

   testCase.verifyEqual(numel(groups), 2)
   testCase.verifyEqual([groups.id]', ["mar_a"; "mar_b"])
   testCase.verifyEqual(height(groups(1).data), 2)
   testCase.verifyEqual(groups(2).datetime, t2)
end

function test_sumup_stable_provenance_defines_groups(testCase)
   % Depth, measurement, reference, and method ids do not split one physical
   % SUMup profile; distinct name keys and dates remain separate profiles.
   t = datetime(2010, 5, 8, 'TimeZone', 'UTC');
   t2 = t + days(1);
   profiles = table(repmat("core", 6, 1), [1; 1; 2; 2; 1; 1], ...
      [10; 11; 10; 10; 10; 11], [1; 2; 1; 1; 1; 2], ...
      repmat(70, 6, 1), repmat(-40, 6, 1), ...
      [t; t; t; t; t2; t2], [0; 1; 0; 1; 0; 1], ...
      (1:6)', [300; 400; 350; 450; 325; 425], ...
      'VariableNames', {'name', 'name_key', 'reference_key', 'method_key', ...
      'latitude', 'longitude', 'datetime', 'midpoint', 'measurement_id', ...
      'density'});

   groups = icemodel.verification.helpers.profileGroups(profiles);

   testCase.verifyEqual(numel(groups), 3)
   testCase.verifyEqual([groups.id]', ...
      ["name_key|1"; "name_key|2"; "name_key|1"])
   testCase.verifyEqual(arrayfun(@(g) height(g.data), groups), [2; 2; 2])
   testCase.verifyEqual([groups.datetime]', [t; t; t2])
end

function test_sumup_legacy_name_and_location_fallbacks(testCase)
   % Older tables without name_key group by resolved name; unnamed tables use
   % location so adjacent source profiles on one date remain distinct.
   t = datetime(2010, 5, 8, 'TimeZone', 'UTC');
   named = table(["a"; "a"; "b"], nan(3, 1), repmat(t, 3, 1), ...
      [0; 1; 0], 'VariableNames', ...
      {'name', 'name_key', 'datetime', 'depth'});
   groups = icemodel.verification.helpers.profileGroups(named);
   testCase.verifyEqual([groups.id]', ["name|a"; "name|b"])
   testCase.verifyEqual(arrayfun(@(g) height(g.data), groups), [2; 1])

   located = table(strings(3, 1), nan(3, 1), [70; 70; 71], ...
      [-40; -40; -41], repmat(t, 3, 1), [0; 1; 0], ...
      'VariableNames', {'name', 'name_key', 'latitude', 'longitude', ...
      'datetime', 'depth'});
   groups = icemodel.verification.helpers.profileGroups(located);
   testCase.verifyEqual(numel(groups), 2)
   testCase.verifyEqual(arrayfun(@(g) height(g.data), groups), [2; 1])
end

function test_undated_legacy_table_stays_one_group(testCase)
   % A legacy single-profile table remains usable without fabricated dates.
   profile = table([0; 1], [300; 500], ...
      'VariableNames', {'depth', 'density'});

   groups = icemodel.verification.helpers.profileGroups(profile);

   testCase.verifyEqual(numel(groups), 1)
   testCase.verifyEqual(groups.id, "profile")
   testCase.verifyTrue(isnat(groups.datetime))
   testCase.verifyEqual(groups.data, profile)
end

function test_timetable_row_times_define_profile_dates(testCase)
   % Canonical SUMup temperature timetables use Time rather than an explicit
   % datetime variable, and must still remain separated by UTC date.
   t1 = datetime(2012, 5, 1, 'TimeZone', 'UTC');
   t2 = datetime(2013, 5, 1, 'TimeZone', 'UTC');
   Time = [t1; t1; t2; t2];
   profile_id = repmat("core", 4, 1);
   depth = [0; 1; 0; 1];
   temperature = [-15; -14; -13; -12];
   profiles = timetable(Time, profile_id, depth, temperature);

   groups = icemodel.verification.helpers.profileGroups(profiles);

   testCase.verifyEqual(numel(groups), 2)
   testCase.verifyEqual([groups.datetime]', [t1; t2])
   testCase.verifyEqual(arrayfun(@(g) height(g.data), groups), [2; 2])
end

function test_non_tabular_profiles_are_rejected(testCase)
   % Profile grouping requires canonical table/timetable metadata.
   testCase.verifyError( ...
      @() icemodel.verification.helpers.profileGroups([1; 2]), ...
      'icemodel:verification:profileGroups:badInput')
end

function test_timestamp_resolution_preserves_subdaily_model_profiles(testCase)
   % Model histories can retain separate profiles within one UTC date while
   % observation grouping continues to use the default calendar-date support.
   t0 = datetime(2012, 5, 1, 'TimeZone', 'UTC');
   times = [t0; t0; t0 + hours(12); t0 + hours(12)];
   profiles = table(repmat("model", 4, 1), times, [0; 1; 0; 1], ...
      [300; 500; 310; 510], 'VariableNames', ...
      {'profile_id', 'datetime', 'depth', 'density'});

   daily = icemodel.verification.helpers.profileGroups(profiles);
   timestamped = icemodel.verification.helpers.profileGroups( ...
      profiles, time_resolution="timestamp");

   testCase.verifyEqual(numel(daily), 1)
   testCase.verifyEqual(numel(timestamped), 2)
   testCase.verifyEqual(arrayfun(@(g) height(g.data), timestamped), [2; 2])
end
