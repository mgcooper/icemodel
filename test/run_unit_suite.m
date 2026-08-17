function results = run_unit_suite(options)
   %RUN_UNIT_SUITE Run icemodel unit tests with folder-based discovery.
   %
   %  results = run_unit_suite()
   %  results = run_unit_suite(selector="test_met_contracts")
   %  results = run_unit_suite(selector="contracts")
   %  results = run_unit_suite(debug=true)
   %  results = run_unit_suite(progress_log="/tmp/unit_progress.log")
   %
   % SELECTOR may be:
   %  - empty: run the whole unit suite
   %  - a file name under test/unit/
   %  - a subfolder under test/unit/
   %  - an absolute/relative existing file or folder path
   %
   % The suite runs one test file at a time. A progress line prints
   % before and after each file (to stdout in a desktop session, to
   % stderr otherwise), and a per-file wall-clock table
   % prints at the end. Set PROGRESS_LOG to also append each progress line
   % to that file with a per-line open/write/close, so a run that hangs
   % still leaves a log whose last "..." line names the file that never
   % finished.

   arguments (Input)

      options.selector (1, 1) string ...
         = ""

      options.debug (1, 1) logical ...
         = false

      options.stop_on_failure (1, 1) logical ...
         = false

      options.verbosity (1, :) string ...
         {icemodel.validators.mustBeTestVerbosityName(options.verbosity)} ...
         = "concise" % "terse" "concise" "detailed"

      options.progress_log (1, 1) string ...
         = ""
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
      icemodel.test.helpers.bootstrapTestEnvironment();

   % Fail fast when the progress log cannot be opened for append: a
   % silently missing log would defeat its hang-diagnosis purpose.
   if strlength(options.progress_log) > 0
      fid = fopen(options.progress_log, 'a');
      if fid == -1
         error('icemodel:test:runUnitSuite:progressLogOpenFailed', ...
            'cannot open progress_log for append: %s', options.progress_log)
      end
      fclose(fid);
   end

   % Build the requested suite and configure the text runner once.
   suite = buildUnitSuite(unitdir, options.selector);
   runner = TestRunner.withTextOutput('Verbosity', ...
      mapVerbosity(options.verbosity));

   % Add a StopOnFailuresPlugin so failures can be inspected interactively.
   if options.debug || options.stop_on_failure
      runner.addPlugin(StopOnFailuresPlugin)
   end

   % Run the suite one test file at a time. File boundaries are fixture
   % boundaries, so per-file runs produce the same results as one
   % whole-suite run, and the loop makes a long run observable: stdout
   % buffers under matlab -batch when redirected, so without unbuffered
   % progress lines a slow run looks identical to a hung one (see
   % icemodel-6qo and icemodel-e9j).
   groups = fileKeys(suite);
   names = unique(groups, 'stable');
   ngroups = numel(names);
   chunks = cell(ngroups, 1);
   durations = zeros(ngroups, 1);
   for g = 1:ngroups
      members = suite(groups == names(g));
      progress(sprintf('[%d/%d] %s (%d tests) ...', ...
         g, ngroups, names(g), numel(members)), options.progress_log);
      tstart = tic;
      chunks{g} = runner.run(members);
      durations(g) = toc(tstart);
      progress(sprintf('[%d/%d] %s: %d passed, %d failed, %d incomplete in %.1f s', ...
         g, ngroups, names(g), sum([chunks{g}.Passed]), ...
         sum([chunks{g}.Failed]), sum([chunks{g}.Incomplete]), ...
         durations(g)), options.progress_log);
   end

   % Concatenate per-file results so the return value keeps the shape and
   % order a single whole-suite run produces. The empty seed keeps the
   % TestResult type when the suite resolves to zero files.
   results = [matlab.unittest.TestResult.empty(1, 0), chunks{:}];

   % Print per-file wall-clock, slowest first, so the files that dominate
   % the suite runtime are identifiable from any run.
   [sorted_durations, order] = sort(durations, 'descend');
   fprintf('\nPer-file wall-clock (slowest first):\n');
   for g = 1:ngroups
      fprintf('%9.1f s  %s\n', sorted_durations(g), names(order(g)));
   end

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

   % Restore the caller's config once the suite has finished. An early error
   % also restores it, because MATLAB deletes the cleanup object when the
   % scope ends.
   delete(suite_cleanup)
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

function keys = fileKeys(suite)
   %FILEKEYS Map each suite element to its containing test-file name.

   % A test name is '<file>/<test>' with optional parameterization after
   % the test token. Appending '/' before extractBefore guarantees a hit
   % even if a name ever arrives without a separator.
   keys = extractBefore(string({suite.Name}) + "/", "/");
end

function progress(line, logfile)
   %PROGRESS Emit one runner progress line to the console and optional log.

   % stderr is unbuffered under matlab -batch, so redirected runs stay
   % observable in real time. The desktop styles stderr as an error (with
   % an "Explain Error" button), so interactive sessions print to
   % stdout instead, which the desktop flushes immediately anyway. The
   % per-line open/write/close on the log file makes each line durable,
   % so a hard hang still leaves the name of the file that started but
   % never finished (the icemodel-e9j recipe).
   if usejava('desktop')
      stream = 1;
   else
      stream = 2;
   end
   fprintf(stream, '%s %s\n', ...
      char(datetime('now', 'Format', 'HH:mm:ss')), line);
   if strlength(logfile) > 0
      fid = fopen(logfile, 'a');
      if fid == -1
         % The runner verified this path opens before the suite started,
         % so a failure here is transient; surface it on the console
         % stream instead of dropping the line without notice.
         fprintf(stream, 'progress log append failed: %s\n', logfile);
      else
         fprintf(fid, '%s %s\n', ...
            char(datetime('now', 'Format', 'HH:mm:ss')), line);
         fclose(fid);
      end
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
