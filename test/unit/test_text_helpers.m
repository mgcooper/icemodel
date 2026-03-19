function tests = test_text_helpers
   %TEST_TEXT_HELPERS Test text and blank-input helper behavior.
   tests = functiontests(localfunctions);
end

function test_isblanktext_matches_text_only_semantics(testCase)
   verifyTrue(testCase, isblanktext(''));
   verifyTrue(testCase, isblanktext(""));
   verifyTrue(testCase, isblanktext(string.empty()));

   verifyFalse(testCase, isblanktext([]));
   verifyFalse(testCase, isblanktext(missing));
   verifyFalse(testCase, isblanktext('abc'));
end

function test_isblankinput_preserves_empty_placeholder_behavior(testCase)
   verifyTrue(testCase, icemodel.isblankinput([]));
   verifyTrue(testCase, icemodel.isblankinput(''));
   verifyTrue(testCase, icemodel.isblankinput(""));
   verifyFalse(testCase, icemodel.isblankinput(0));
end

function test_isblankpart_handles_path_part_placeholders(testCase)
   verifyTrue(testCase, icemodel.isblankpart(''));
   verifyTrue(testCase, icemodel.isblankpart(""));
   verifyTrue(testCase, icemodel.isblankpart({}));

   verifyFalse(testCase, icemodel.isblankpart([]));
   verifyFalse(testCase, icemodel.isblankpart({'run001'}));
end

function test_istextlike_alias_preserves_ischarlike_behavior(testCase)
   cases = { ...
      'a', ...
      "", ...
      {'a', "b"}, ...
      {'a', 1}};

   returned = cellfun(@(value) ischarlike(value) == istextlike(value), cases);
   verifyTrue(testCase, all(returned));
   verifyFalse(testCase, istextlike("", 'nontrivial'));
   verifyFalse(testCase, ischarlike("", 'nontrivial'));
end

function test_setopts_treats_empty_optional_inputs_as_omitted(testCase)
   opts = icemodel.setopts('icemodel', 'behar', 2016, 'kanm', [], [], []);

   verifyEqual(testCase, opts.userdata, 'kanm');
   verifyEqual(testCase, opts.uservars, 'albedo');
   verifyEqual(testCase, opts.testname, '');
end
