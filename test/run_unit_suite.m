function results = run_unit_suite(options)
   %RUN_UNIT_SUITE Run icemodel unit tests with folder-based discovery.
   %
   %  results = run_unit_suite()
   %  results = run_unit_suite(selector="test_met_contracts")
   %  results = run_unit_suite(selector="contracts")
   %  results = run_unit_suite(debug=true)
   %
   % SELECTOR may be:
   %  - empty: run the whole unit suite
   %  - a file name under test/unit/
   %  - a subfolder under test/unit/
   %  - an absolute/relative existing file or folder path

   arguments
      options.selector (1, 1) string = ""
      options.debug (1, 1) logical = false
      options.stop_on_failure (1, 1) logical = false
      options.verbosity (1, :) string ...
         {icemodel.validators.mustBeTestVerbosityName(options.verbosity)} ...
         = "concise" % "terse" "concise" "detailed"
   end

   import matlab.unittest.TestRunner
   import matlab.unittest.TestSuite
   import matlab.unittest.Verbosity
   import matlab.unittest.plugins.StopOnFailuresPlugin

   % Resolve the unit test directory.
   thisdir = fileparts(mfilename('fullpath'));
   unitdir = fullfile(thisdir, 'unit');

   % Keep the cleanup handle in scope so the caller's config is restored
   % when this entrypoint returns.
   [~, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(); %#ok<ASGLU>

   % Build the requested suite and configure the text runner once.
   suite = buildUnitSuite(unitdir, options.selector);
   runner = TestRunner.withTextOutput('Verbosity', ...
      mapVerbosity(options.verbosity));

   % Add a StopOnFailuresPlugin so failures can be inspected interactively.
   if options.debug || options.stop_on_failure
      runner.addPlugin(StopOnFailuresPlugin)
   end

   % Run the suite.
   results = runner.run(suite);

   % Print the results to the screen
   if options.verbosity == "detailed"
      for n = 1:numel(results)
         if results(n).Passed
            fprintf('Passed Test %s:%s\n', int2str(n), results(n).Name);
         else
            fprintf('Failed Test %s:%s\n', int2str(n), results(n).Name);
         end
      end
   end
end

function suite = buildUnitSuite(unitdir, selector)
   %BUILDUNITSUITE Resolve the requested unit-test selector into a suite.

   import matlab.unittest.TestSuite

   % Build the full unit suite when no selector was provided.
   if isblanktext(selector)
      suite = TestSuite.fromFolder(unitdir, 'IncludingSubfolders', true);
      return
   end

   % Resolve the selector to one file or folder under test/unit/.
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
   %RESOLVESELECTOR Map a selector string to a unit-test file or folder.

   selector = char(selector);
   if exist(selector, 'dir') == 7
      target = selector;
      return
   end

   if exist(selector, 'file') == 2
      pathhit = which(selector);
      if ~isempty(pathhit)
         target = pathhit;
      else
         target = selector;
      end
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
   %MAPVERBOSITY Convert a friendly name into MATLAB's verbosity enum.

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
