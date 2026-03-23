function profdir = baselineProfilerDir(pathname)
   %BASELINEPROFILERDIR Return the profiler-artifact folder for a baseline file.
   %
   %  profdir = icemodel.test.helpers.baselineProfilerDir(pathname)

   arguments
      pathname {mustBeTextScalar}
   end

   % Strip the extension so the profiler directory matches the baseline stem.
   [~, stem] = fileparts(char(pathname));

   % Keep profiler bundles in the repo-local test tree keyed by baseline stem.
   profdir = fullfile(icemodel.getpath('test'), 'baselines', 'profiler', stem);
end
