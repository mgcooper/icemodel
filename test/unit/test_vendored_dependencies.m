function tests = test_vendored_dependencies
   %TEST_VENDORED_DEPENDENCIES Test icemodel's vendored matfunclib copies.
   %
   % icemodel/dependencies holds copies of matfunclib functions so icemodel
   % runs without matfunclib on the path. The 2026-09-07 refresh
   % (matfunclib-juq.23) took each copy from its matfunclib source. Two
   % results of that refresh are icemodel-specific and are not covered by
   % matfunclib's own tests, so they are covered here.
   %
   % Path note: the suite bootstrap calls icemodel.dependencies(), which
   % runs addpath(genpath(...)) on the sibling matfunclib checkout when it
   % is present. addpath prepends, so matfunclib's originals would shadow
   % the vendored copies and these tests would exercise the wrong files.
   % setupOnce puts the dependencies folder back in front, and each test
   % asserts the function under test resolves there.
   tests = functiontests(localfunctions);
end

function setupOnce(testCase)
   %SETUPONCE Put the vendored dependencies folder ahead of matfunclib.
   here = fileparts(mfilename('fullpath'));
   repo = fileparts(fileparts(here));
   depdir = fullfile(repo, 'icemodel', 'dependencies');
   assertTrue(testCase, isfolder(depdir), ...
      sprintf('Expected the vendored dependencies folder at %s', depdir));

   original = path();
   testCase.addTeardown(@() path(original));
   addpath(depdir, '-begin');

   testCase.TestData.depdir = depdir;
end

function verifyResolvesToVendoredCopy(testCase, funcname)
   %VERIFYRESOLVESTOVENDOREDCOPY Guard against testing matfunclib's copy.
   returned = which(funcname);
   expected = testCase.TestData.depdir;
   verifySubstring(testCase, returned, expected, sprintf( ...
      '%s resolved to %s, not the vendored copy.', funcname, returned));
end

function test_mustBeStruct_accepts_a_struct(testCase)
   %TEST_MUSTBESTRUCT_ACCEPTS_A_STRUCT One-argument validator call passes.
   %
   % Every icemodel call site is an arguments-block validator under
   % +icemodel/+netcdf, which calls mustBeStruct with one input, so the
   % one-argument branch is the live one.
   verifyResolvesToVendoredCopy(testCase, 'mustBeStruct');

   verifyWarningFree(testCase, @() mustBeStruct(struct('a', 1)));
   verifyWarningFree(testCase, @() mustBeStruct(struct('a', {1, 2})));
   verifyWarningFree(testCase, @() mustBeStruct(struct([])));
end

function test_mustBeStruct_rejects_a_nonstruct(testCase)
   %TEST_MUSTBESTRUCT_REJECTS_A_NONSTRUCT One-argument branch raises.
   %
   % The vendored copy resolves the caller from dbstack rather than calling
   % matfunclib's mcallername, which icemodel does not vendor. The error
   % must still carry matfunclib's identifier so callers that catch it are
   % unaffected by that deviation.
   verifyResolvesToVendoredCopy(testCase, 'mustBeStruct');

   verifyError(testCase, @() mustBeStruct(42), ...
      'custom:validators:expectedStructInput');
   verifyError(testCase, @() mustBeStruct('text'), ...
      'custom:validators:expectedStructInput');
   verifyError(testCase, @() mustBeStruct({1, 2}), ...
      'custom:validators:expectedStructInput');
end

function test_mustBeStruct_names_the_calling_function(testCase)
   %TEST_MUSTBESTRUCT_NAMES_THE_CALLING_FUNCTION dbstack lookup resolves.
   %
   % The message interpolates the caller name in upper case. Calling
   % through a local helper puts a known name on the stack, so the lookup
   % is observable rather than incidental.
   verifyResolvesToVendoredCopy(testCase, 'mustBeStruct');

   try
      callMustBeStructFromAHelper();
      verifyFail(testCase, 'mustBeStruct accepted a non-struct.');
   catch err
      verifyEqual(testCase, err.identifier, ...
         'custom:validators:expectedStructInput');
      verifySubstring(testCase, err.message, ...
         upper('callMustBeStructFromAHelper'));
   end
end

function test_mustBeStruct_two_argument_form_names_its_caller(testCase)
   %TEST_MUSTBESTRUCT_TWO_ARGUMENT_FORM_NAMES_ITS_CALLER Explicit name.
   %
   % The refresh added the one- and two-argument forms. With two inputs the
   % caller supplies the name, so dbstack must not be consulted and the
   % supplied name must reach the message.
   verifyResolvesToVendoredCopy(testCase, 'mustBeStruct');

   verifyWarningFree(testCase, ...
      @() mustBeStruct(struct('a', 1), 'suppliedName'));

   try
      mustBeStruct(42, 'suppliedName');
      verifyFail(testCase, 'mustBeStruct accepted a non-struct.');
   catch err
      verifyEqual(testCase, err.identifier, ...
         'custom:validators:expectedStructInput');
      verifySubstring(testCase, err.message, upper('suppliedName'));
      % The supplied name must win: the helper frame must not appear.
      verifyEmpty(testCase, ...
         regexp(err.message, upper('test_vendored'), 'once'));
   end
end

function test_mustBeStruct_four_argument_form_still_works(testCase)
   %TEST_MUSTBESTRUCT_FOUR_ARGUMENT_FORM_STILL_WORKS Arity is preserved.
   %
   % The four-argument form matfunclib documents keeps its own identifier,
   % which embeds the caller name, and its own message.
   verifyResolvesToVendoredCopy(testCase, 'mustBeStruct');

   verifyError(testCase, @() mustBeStruct(42, 'myfunc', 'S', 1), ...
      'custom:myfunc:expectedStructInput');
   verifyWarningFree(testCase, ...
      @() mustBeStruct(struct('a', 1), 'myfunc', 'S', 1));
end

function test_islogicalscalar_matches_its_name(testCase)
   %TEST_ISLOGICALSCALAR_MATCHES_ITS_NAME Newly vendored helper behaves.
   %
   % parseoptarg calls islogicalscalar on its documented logical-default
   % path. icemodel did not vendor it before the refresh, so that path
   % raised an undefined-function error.
   verifyResolvesToVendoredCopy(testCase, 'islogicalscalar');

   verifyTrue(testCase, islogicalscalar(true));
   verifyTrue(testCase, islogicalscalar(false));

   verifyFalse(testCase, islogicalscalar([true, false]));
   verifyFalse(testCase, islogicalscalar(1));
   verifyFalse(testCase, islogicalscalar('true'));
   verifyFalse(testCase, islogicalscalar(logical.empty()));
end

function test_parseoptarg_accepts_a_logical_default(testCase)
   %TEST_PARSEOPTARG_ACCEPTS_A_LOGICAL_DEFAULT The closure is complete.
   %
   % A logical DEFAULTOPT is the call path that reaches islogicalscalar.
   % Before that helper was vendored, both calls below raised an
   % undefined-function error. This is a regression test for the vendored
   % set's completeness, not for parseoptarg's own logic, which matfunclib
   % owns.
   verifyResolvesToVendoredCopy(testCase, 'parseoptarg');

   % No element of ARGS is in VALIDOPTS, so OPT falls back to DEFAULTOPT.
   [returned, args, nargs] = parseoptarg({1, 2}, {'flag'}, false);
   verifyFalse(testCase, returned);
   verifyEqual(testCase, args, {1, 2});
   verifyEqual(testCase, nargs, 2);

   % A match sets OPT true and removes the flag from ARGS.
   [returned, args, nargs] = parseoptarg({'flag', 1}, {'flag'}, false);
   verifyTrue(testCase, returned);
   verifyEqual(testCase, args, {1});
   verifyEqual(testCase, nargs, 1);
end

function test_withcd_license_text_is_not_corrupted(testCase)
   %TEST_WITHCD_LICENSE_TEXT_IS_NOT_CORRUPTED Guard the repaired file.
   %
   % A find-and-replace once rewrote Y to msg and X to cmd inside this
   % file's BSD license text, producing COPmsgRIGHT and EcmdPRESS. The
   % 2026-09-07 refresh repaired it. This test fails if a similar
   % substitution runs again.
   here = fileparts(mfilename('fullpath'));
   repo = fileparts(fileparts(here));
   target = fullfile(repo, 'icemodel', '+icemodel', '+internal', ...
      'private', 'withcd.m');
   assertTrue(testCase, isfile(target), ...
      sprintf('Expected withcd.m at %s', target));

   contents = fileread(target);
   verifyEmpty(testCase, regexp(contents, 'COPmsgRIGHT', 'once'));
   verifyEmpty(testCase, regexp(contents, 'EcmdPRESS', 'once'));
   verifySubstring(testCase, contents, 'THE COPYRIGHT HOLDERS');
end

%% local functions

function callMustBeStructFromAHelper()
   %CALLMUSTBESTRUCTFROMAHELPER Put a known name on the stack.
   %
   % mustBeStruct reads the caller from dbstack, so the assertion needs a
   % named frame between the test and the validator.
   mustBeStruct(42);
end
