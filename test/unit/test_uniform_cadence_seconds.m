function tests = test_uniform_cadence_seconds
   %TEST_UNIFORM_CADENCE_SECONDS Cover the single cadence-derivation helper.
   %
   % artifactCadenceMatches and artifactMetadata delegate the entire
   % is-this-regular decision here, so every branch needs a direct test: a
   % regression in the short-payload branch would let a writer stamp a cadence
   % on a payload that has none.
   tests = functiontests(localfunctions);
end

function value = payload(time)
   %PAYLOAD Build a minimal timetable carrying only a time coordinate.
   value = timetable(zeros(numel(time), 1), 'RowTimes', time, ...
      'VariableNames', {'tair'});
end

function test_regular_cadence_is_returned_in_seconds(testCase)
   % The ordinary case: a regularly sampled payload reports its spacing.
   time = datetime(2019, 7, 1) + minutes(0:15:60)';
   returned = icemodel.forcing.helpers.uniformCadenceSeconds(payload(time));

   expected = 900;
   testCase.verifyEqual(returned, expected)
end

function test_fewer_than_two_rows_has_no_cadence(testCase)
   % A single sample defines no spacing, and an empty payload defines none
   % either. Both must report NaN rather than inventing a cadence.
   one_row = payload(datetime(2019, 7, 1));
   testCase.verifyTrue(isnan( ...
      icemodel.forcing.helpers.uniformCadenceSeconds(one_row)))

   empty_rows = payload(datetime.empty(0, 1));
   testCase.verifyTrue(isnan( ...
      icemodel.forcing.helpers.uniformCadenceSeconds(empty_rows)))
end

function test_irregular_spacing_has_no_cadence(testCase)
   % One gap is enough to disqualify the payload, because a stamped cadence
   % would claim regularity the samples do not have.
   time = datetime(2019, 7, 1) ...
      + minutes([0, 15, 30, 75, 90])';
   testCase.verifyTrue(isnan( ...
      icemodel.forcing.helpers.uniformCadenceSeconds(payload(time))))
end

function test_median_defines_the_cadence_not_the_first_step(testCase)
   % A corrupt leading sample must not define the accepted cadence. The first
   % step here is wrong, so the payload is irregular and must be rejected
   % rather than accepted at the leading step's spacing.
   time = datetime(2019, 7, 1) + minutes([0, 5, 20, 35, 50])';
   testCase.verifyTrue(isnan( ...
      icemodel.forcing.helpers.uniformCadenceSeconds(payload(time))))
end

function test_round_off_within_tolerance_is_still_regular(testCase)
   % Datetime arithmetic leaves sub-microsecond noise, which must not be read
   % as a real gap.
   time = datetime(2019, 7, 1) + seconds([0, 900, 1800, 2700])' ...
      + seconds([0, 1e-9, 0, -1e-9])';
   returned = icemodel.forcing.helpers.uniformCadenceSeconds(payload(time));

   testCase.verifyEqual(returned, 900, AbsTol=1e-6)
end
