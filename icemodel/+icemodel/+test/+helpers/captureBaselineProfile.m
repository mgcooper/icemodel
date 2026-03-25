function [profile_summary, profile_meta, profile_artifacts] = ...
      captureBaselineProfile(kind, cases, output_file, kwargs)
   %CAPTUREBASELINEPROFILE Save a build-time profiler report alongside a baseline.
   %
   %  [profile_summary, profile_meta, profile_artifacts] = ...
   %     icemodel.test.helpers.captureBaselineProfile("perf", cases, output_file)
   %
   % The accepted baseline build should stay deterministic and timing-focused.
   % This helper therefore reruns the accepted workflow *after* the build and
   % captures profiler diagnostics as separate artifacts that can be archived
   % with the managed baseline.

   arguments
      kind (1, :) string {mustBeMember(kind, ["perf", "regression"])}
      cases table
      output_file {mustBeTextScalar}
      kwargs.history_size (1, 1) double {mustBeInteger, mustBePositive} = ...
         25000000
   end

   % Replace any previous managed profile bundle for this baseline target.
   profdir = icemodel.test.helpers.baselineProfilerDir(output_file);
   if exist(profdir, 'dir') == 7
      rmdir(profdir, 's');
   end
   mkdir(profdir);

   % Profile the accepted workflow in a second pass so the saved baseline
   % values are not affected by profiler overhead or altered execution order.
   profile clear
   profile('-historysize', kwargs.history_size);
   profile on

   % Re-run the accepted workflow for the selected suite type.
   switch kind
      case "perf"
         profilePerfBuild(cases);
      case "regression"
         profileRegressionBuild(cases);
   end

   profile off
   info = profile('info');
   profsave(info, profdir);
   html_entry = resolveProfileHtmlEntry(profdir);
   pruneProfileHtmlBundle(profdir, html_entry);

   % Keep both the raw profile info and a compact summary for later tools.
   info_file = fullfile(profdir, 'profile_info.mat');
   save(info_file, 'info');

   profile_summary = summarizeProfileInfo(info);
   profile_meta = struct();
   profile_meta.kind = kind;
   profile_meta.history_size = kwargs.history_size;
   profile_meta.case_ids = string(cases.case_id);
   profile_meta.timestamp_utc = datetime('now', 'TimeZone', 'UTC');
   profile_meta.matlab_version = string(version);
   profile_meta.host = string(computer);

   profile_artifacts = struct();
   profile_artifacts.dir = string(profdir);
   profile_artifacts.index_file = string(html_entry);
   profile_artifacts.info_file = string(info_file);
end

function profilePerfBuild(cases)
   %PROFILEPERFBUILD Re-run the managed perf cases under MATLAB profiling.

   for icase = 1:height(cases)
      fprintf('Profiling perf case %d/%d: %s\n', ...
         icase, height(cases), cases.case_id(icase))
      opts_run = icemodel.test.helpers.setModelOptsForCase(cases(icase, :));
      icemodel.test.helpers.runSmbModel(opts_run);
   end
end

function profileRegressionBuild(cases)
   %PROFILEREGRESSIONBUILD Re-run regression cases under MATLAB profiling.

   runoff_ref = icemodel.test.helpers.loadReference("runoff");

   for icase = 1:height(cases)
      c = cases(icase, :);
      fprintf('Profiling regression case %d/%d: %s\n', ...
         icase, height(cases), c.case_id)
      opts_run = icemodel.test.helpers.setModelOptsForCase(c);
      [ice1, ice2] = icemodel.test.helpers.runSmbModel(opts_run);
      [ice1, ~] = icemodel.postprocess( ...
         ice1, ice2, opts_run, opts_run.output_years);
      ridx = icemodel.test.helpers.findRunoffReferenceRow(runoff_ref, c);
      met = icemodel.test.helpers.loadProcessedMetForOutputYears(opts_run);
      if isempty(ridx)
         refrow = [];
      else
         refrow = runoff_ref(ridx, :);
      end
      icemodel.test.helpers.summarizeIce1Metrics(ice1, met, refrow);
   end
end

function html_entry = resolveProfileHtmlEntry(profdir)
   %RESOLVEPROFILEHTMLENTRY Locate the main profiler HTML entrypoint.

   html_entry = fullfile(profdir, 'file0.html');
   if exist(html_entry, 'file') == 2
      return
   end

   listing = dir(fullfile(profdir, '*.html'));
   if isempty(listing)
      html_entry = "";
   else
      html_entry = fullfile(listing(1).folder, listing(1).name);
   end
end

function pruneProfileHtmlBundle(profdir, html_entry)
   %PRUNEPROFILEHTMLBUNDLE Keep only the top-level profiler HTML page.

   html_entry = char(html_entry);
   if isempty(html_entry)
      return
   end
   listing = dir(fullfile(profdir, '*.html'));
   for i = 1:numel(listing)
      pathname = fullfile(listing(i).folder, listing(i).name);
      if strcmp(pathname, html_entry)
         continue
      end
      delete(pathname);
   end
end

function summary = summarizeProfileInfo(info)
   %SUMMARIZEPROFILEINFO Convert `profile('info')` output into a compact table.

   if ~isfield(info, 'FunctionTable') || isempty(info.FunctionTable)
      summary = table();
      return
   end

   % Read the profiler function table once before flattening its rows into a
   % compact, sortable struct array.
   ft = info.FunctionTable;
   total_time = sum([ft.TotalTime], 'omitnan');
   rows = repmat(struct( ...
      'FunctionName', "", ...
      'CompleteName', "", ...
      'FileName', "", ...
      'Type', "", ...
      'NumCalls', 0, ...
      'TotalTime', 0, ...
      'TotalRecursiveTime', 0, ...
      'FractionOfTotal', 0), numel(ft), 1);

   for i = 1:numel(ft)
      rows(i).FunctionName = string(ft(i).FunctionName);
      rows(i).CompleteName = string(ft(i).CompleteName);
      rows(i).FileName = string(ft(i).FileName);
      rows(i).Type = string(ft(i).Type);
      rows(i).NumCalls = ft(i).NumCalls;
      rows(i).TotalTime = ft(i).TotalTime;
      rows(i).TotalRecursiveTime = ft(i).TotalRecursiveTime;
      if total_time > 0
         rows(i).FractionOfTotal = ft(i).TotalTime / total_time;
      end
   end

   % Sort by total time so the dominant hotspots appear first.
   summary = struct2table(rows);
   summary = sortrows(summary, 'TotalTime', 'descend');
end
