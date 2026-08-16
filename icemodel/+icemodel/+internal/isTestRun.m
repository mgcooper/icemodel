function tf = isTestRun()
   %ISTESTRUN True when the MATLAB test framework is on the call stack.
   %
   %  tf = icemodel.internal.isTestRun()
   %
   % Use this to suppress console output that helps an interactive user but
   % buries the pass and fail marks a test runner prints. A loader that
   % announces each file, or a plot helper that warns about an absent
   % variable, is useful at the prompt and is noise across a suite.
   %
   % Do not use this to change a computed result. Output that a caller reads,
   % and any error, must not depend on whether a test is running.
   %
   % Output
   %   tf - True inside a matlab.unittest run [logical].
   %
   % See also: icemodel.test.helpers.printFilePath,
   %  icemodel.verification.plotscatter

   % dbstack reports the framework's own files by full path. Matching the
   % installed testframework folder detects any runner, so runtests, a
   % TestRunner built by run_unit_suite, and a direct run(suite) all count.
   stack = dbstack('-completenames');
   if isempty(stack)
      tf = false;
      return
   end
   framework = fullfile(matlabroot, 'toolbox', 'matlab', 'testframework');
   tf = any(startsWith({stack.file}, framework));
end
