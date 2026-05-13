function archiveManagedBaseline(pathname, kind)
   %ARCHIVEMANAGEDBASELINE Archive an existing rolling baseline before overwrite.
   %
   %  icemodel.test.helpers.archiveManagedBaseline(pathname, "perf")
   %
   % The rolling baseline files are mutable acceptance targets. Before
   % replacing one, archive the prior managed MAT file and any saved profiler
   % artifacts under test/baselines/archive/ so older accepted states remain
   % available for later inspection.

   arguments
      pathname {mustBeTextScalar}
      kind (1, :) string {mustBeMember(kind, ["perf", "regression"])}
   end

   pathname = char(pathname);
   if ~isfile(pathname)
      return
   end

   % Create a timestamped archive folder for this overwrite event.
   timestamp = string(datetime('now', 'TimeZone', 'UTC', ...
      'Format', 'yyyyMMdd-HHmmss'));
   archivedir = fullfile(icemodel.getpath('test'), 'baselines', 'archive', ...
      char(kind), char(timestamp));
   if ~isfolder(archivedir)
      mkdir(archivedir);
   end

   [~, filename, ext] = fileparts(pathname);
   copyfile(pathname, fullfile(archivedir, [filename ext]));

   % Carry the profiler sidecar folder along with the baseline MAT file.
   src_profdir = icemodel.test.helpers.baselineProfilerDir(pathname);
   if isfolder(src_profdir)
      copyfile(src_profdir, fullfile(archivedir, [filename '_profiler']));
   end
end
