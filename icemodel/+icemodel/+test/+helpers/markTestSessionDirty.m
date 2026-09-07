function markTestSessionDirty(label)
   %MARKTESTSESSIONDIRTY Record that this MATLAB session ran a test suite.
   %
   %  icemodel.test.helpers.markTestSessionDirty("run_unit_suite")
   %
   % Formal performance timing is sensitive to what ran earlier in the same
   % session: JIT state, persistents, and heap layout all carry over. Each
   % suite runner calls this at start so a later in-session formal perf run
   % can refuse a contaminated session (see
   % icemodel.test.helpers.assertCleanPerfSession). The record is an
   % environment variable so `matlab -batch` one-shot sessions start clean
   % and subprocesses inherit their parent's history.
   %
   % See also: icemodel.test.helpers.testSessionActivity,
   %  icemodel.test.helpers.assertCleanPerfSession

   name = 'ICEMODEL_TEST_SESSION_ACTIVITY';
   prior = getenv(name);

   % Append rather than overwrite so the record lists every suite run.
   if isempty(prior)
      setenv(name, char(label));
   else
      setenv(name, [prior ',' char(label)]);
   end
end
