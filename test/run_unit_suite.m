function results = run_unit_suite(kwargs)
   %RUN_UNIT_SUITE Run icemodel unit tests with folder-based discovery.
   %
   %  results = run_unit_suite()
   %  results = run_unit_suite(selector="test_met_contracts.m")
   %  results = run_unit_suite(selector="contracts")
   %  results = run_unit_suite(debug=true)
   %
   % SELECTOR may be:
   %  - empty: run the whole unit suite
   %  - a file name under test/unit/
   %  - a subfolder under test/unit/
   %  - an absolute/relative existing file or folder path

   arguments
      kwargs.selector string = string.empty()
      kwargs.debug (1, 1) logical = false
      kwargs.stop_on_failure (1, 1) logical = false
      kwargs.verbosity (1, :) string {mustBeMember(kwargs.verbosity, ...
         ["terse", "concise", "detailed"])} = "detailed"
   end

   import matlab.unittest.TestRunner
   import matlab.unittest.TestSuite
   import matlab.unittest.Verbosity
   import matlab.unittest.plugins.StopOnFailuresPlugin

   thisdir = fileparts(mfilename('fullpath'));
   addpath(fullfile(fileparts(thisdir), 'icemodel'))
   addpath(thisdir)
   unitdir = fullfile(thisdir, 'unit');

   suite = buildUnitSuite(unitdir, kwargs.selector);
   runner = TestRunner.withTextOutput('Verbosity', ...
      mapVerbosity(kwargs.verbosity));

   if kwargs.debug || kwargs.stop_on_failure
      runner.addPlugin(StopOnFailuresPlugin)
   end

   results = runner.run(suite);
end

function suite = buildUnitSuite(unitdir, selector)

   import matlab.unittest.TestSuite

   if isempty(selector) || all(strlength(selector) == 0)
      suite = TestSuite.fromFolder(unitdir, 'IncludingSubfolders', true);
      return
   end

   target = resolveSelector(unitdir, selector);
   if exist(target, 'dir') == 7
      suite = TestSuite.fromFolder(target, 'IncludingSubfolders', true);
   elseif exist(target, 'file') == 2
      suite = TestSuite.fromFile(target);
   else
      error('unit test selector does not resolve to a file/folder: %s', target)
   end
end

function target = resolveSelector(unitdir, selector)

   selector = char(selector);
   if exist(selector, 'file') == 2 || exist(selector, 'dir') == 7
      target = selector;
      return
   end

   target = fullfile(unitdir, selector);
   if exist(target, 'file') == 2 || exist(target, 'dir') == 7
      return
   end

   if exist([target '.m'], 'file') == 2
      target = [target '.m'];
      return
   end
end

function v = mapVerbosity(name)

   import matlab.unittest.Verbosity

   switch char(name)
      case 'terse'
         v = Verbosity.Terse;
      case 'concise'
         v = Verbosity.Concise;
      case 'detailed'
         v = Verbosity.Detailed;
      otherwise
         error('unsupported verbosity: %s', name)
   end
end
