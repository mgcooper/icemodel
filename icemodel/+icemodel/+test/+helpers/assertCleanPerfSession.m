function assertCleanPerfSession(isolation)
   %ASSERTCLEANPERFSESSION Refuse an in-session formal run in a dirty session.
   %
   %  icemodel.test.helpers.assertCleanPerfSession(isolation)
   %
   % A formal in-session perf run in a session that already ran other
   % suites produces timings contaminated by that history (the V5
   % helper-extraction revert came from such a reading). Process isolation
   % is immune because every case runs in a fresh subprocess. This check
   % errors, rather than warns, because a contaminated formal verdict is
   % worse than no verdict.
   %
   % See also: icemodel.test.helpers.markTestSessionDirty,
   %  icemodel.test.helpers.testSessionActivity

   % Subprocess cases start clean, so a dirty parent does not matter.
   if isolation == "process"
      return
   end

   activity = icemodel.test.helpers.testSessionActivity();
   if isempty(activity)
      return
   end

   error('icemodel:test:perf:contaminatedSession', ...
      ['this session already ran: %s. Formal in-session timings would ' ...
      'be contaminated. Use run_perf_suite(isolation="process"), or ' ...
      'run from a fresh MATLAB session.'], strjoin(activity, ', '))
end
