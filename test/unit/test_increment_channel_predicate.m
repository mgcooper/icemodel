function tests = test_increment_channel_predicate
   %TEST_INCREMENT_CHANNEL_PREDICATE Cover the shared increment predicate.
   %
   % retimeHourlyFixedStep and postprocess both decide aggregation with this
   % predicate, and they pass different container types, so the return shape
   % matters as much as the rule.
   tests = functiontests(localfunctions);
end

function test_increment_channels_are_recognized(testCase)
   % Any df_ channel and the recovery count are per-step totals that must sum.
   returned = icemodel.isIncrementChannel('df_liq');
   testCase.verifyTrue(returned)
   testCase.verifyTrue(icemodel.isIncrementChannel('df_rof'))
   testCase.verifyTrue(icemodel.isIncrementChannel("df_evp"))
   testCase.verifyTrue(icemodel.isIncrementChannel('cpl_recovery_count'))
end

function test_non_increment_channels_are_rejected(testCase)
   % State and diagnostic channels average, so they must not match.
   testCase.verifyFalse(icemodel.isIncrementChannel('Tsfc'))
   testCase.verifyFalse(icemodel.isIncrementChannel('melt'))

   % A name that merely contains df_ is not a df_ channel.
   testCase.verifyFalse(icemodel.isIncrementChannel('x_df_liq'))
end

function test_cell_input_returns_a_logical_array(testCase)
   % postprocess passes a cell column and ORs the result with an errH test,
   % so the return must be a logical of the same shape.
   names = {'Tsfc'; 'df_liq'; 'errH'; 'df_rof'};
   returned = icemodel.isIncrementChannel(names);

   expected = [false; true; false; true];
   testCase.verifyEqual(returned, expected)
   testCase.verifySize(returned, size(names))
end

function test_string_array_input_preserves_shape(testCase)
   % retimeHourlyFixedStep indexes a variable list with the mask, so a row
   % input must produce a row mask.
   names = ["df_liq", "Tsfc", "df_evp"];
   returned = icemodel.isIncrementChannel(names);

   expected = [true, false, true];
   testCase.verifyEqual(returned, expected)
end

function test_postprocess_sum_rule_keeps_the_errH_exception(testCase)
   % errH is a residual that sums even though it is not a df_ channel. The
   % predicate must not absorb that exception, or the one caller that needs
   % it would lose it.
   testCase.verifyFalse(icemodel.isIncrementChannel('errH'))
end
