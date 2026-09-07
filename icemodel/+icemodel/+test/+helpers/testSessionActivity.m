function activity = testSessionActivity()
   %TESTSESSIONACTIVITY Return the suites this MATLAB session has run.
   %
   %  activity = icemodel.test.helpers.testSessionActivity()
   %
   % Returns a string row vector of runner labels recorded by
   % icemodel.test.helpers.markTestSessionDirty, oldest first. An empty
   % string array means the session is clean for formal timing.
   %
   % See also: icemodel.test.helpers.markTestSessionDirty,
   %  icemodel.test.helpers.assertCleanPerfSession

   raw = getenv('ICEMODEL_TEST_SESSION_ACTIVITY');
   if isempty(raw)
      activity = strings(1, 0);
   else
      activity = string(strsplit(raw, ','));
   end
end
