function opts = caseSetopts(case_manifest, kwargs)
   %CASESETOPTS Build standard icemodel OPTS for one verification case.
   %
   %  opts = icemodel.verification.helpers.caseSetopts(case_manifest)
   %  opts = icemodel.verification.helpers.caseSetopts(case_manifest, ...
   %     startdate=..., enddate=..., dt=900)
   %
   %  Builds an opts struct via the standard icemodel.setopts contract using
   %  the manifest's case_id as both sitename and forcings, and the years
   %  touched by manifest.comparison_window as simyears. The standard chain
   %  (configureRun + createMetFileNames + loadmet) then resolves the staged
   %  met files at <ICEMODEL_INPUT_PATH>/met/met_<case>_<case>_<year>_<dt>.mat
   %  with no special-case bypass.
   %
   %  Inputs
   %    case_manifest : struct  Resolved case manifest from listcases.
   %
   %  Name-value
   %    startdate / enddate : datetime, default NaT — explicit window
   %        narrowing. When NaT, the manifest comparison_window is used.
   %    dt : double seconds, default derived from manifest.native_timestep
   %        ("hourly" -> 3600). Override only when re-sampling at runtime.
   %    smbmodel : string, default "icemodel"
   %    testname : string, default "verification"
   %
   %  Returns
   %    opts : struct  Standard model options. opts.startdate / opts.enddate
   %                   carry the comparison window so loadmet subsets the
   %                   staged calendar-year met file.
   %
   % See also: icemodel.setopts, icemodel.createMetFileNames,
   %  icemodel.verification.runIcemodelSnowCandidate

   arguments
      case_manifest (1, 1) struct
      kwargs.smbmodel  (1, 1) string = "icemodel"
      kwargs.testname  (1, 1) string = "verification"
      kwargs.dt        (1, 1) double = nativeTimestepToSeconds(case_manifest)
      kwargs.startdate = NaT('TimeZone', 'UTC')
      kwargs.enddate   = NaT('TimeZone', 'UTC')
   end

   case_id = string(case_manifest.case_id);
   [window_start, window_end] = resolveWindow(case_manifest, ...
      kwargs.startdate, kwargs.enddate);

   % Years touched by the comparison window become simyears. Standard
   % loadmet (per icemodel.loadmet) loads each calendar-year met file
   % under simyears and concatenates; opts.startdate / opts.enddate then
   % subset that concatenated met to the actual comparison window.
   simyears = unique(year([window_start; window_end]))';

   opts = icemodel.setopts(char(kwargs.smbmodel), char(case_id), simyears, ...
      char(case_id), [], [], char(kwargs.testname), false, false, ...
      'dt', kwargs.dt, ...
      'startdate', window_start, ...
      'enddate', window_end);
end

function dt_seconds = nativeTimestepToSeconds(case_manifest)
   %NATIVETIMESTEPTOSECONDS Map the manifest's native_timestep label to dt seconds.
   if isfield(case_manifest, 'native_timestep')
      label = lower(string(case_manifest.native_timestep));
   else
      label = "";
   end
   switch label
      case "hourly"
         dt_seconds = 3600;
      case {"15-min", "15 min", "15m", "15-minute", "15 minute"}
         dt_seconds = 900;
      otherwise
         % Default to hourly. The verification suite stages ESM-SnowMIP
         % met at native hourly resolution; runtime resampling to 15-min
         % is a future architectural change tracked separately.
         dt_seconds = 3600;
   end
end

function [window_start, window_end] = resolveWindow(manifest, ...
      startdate, enddate)
   %RESOLVEWINDOW Resolve the canonical comparison window for the case.

   if ~isnat(startdate) && ~isnat(enddate)
      window_start = icemodel.verification.setup.ensureUtc(startdate);
      window_end   = icemodel.verification.setup.ensureUtc(enddate);
      return
   end

   % Manifest's comparison_window is the staged-data window. Use it as
   % the canonical default.
   window_start = icemodel.verification.setup.ensureUtc( ...
      manifest.comparison_window.start);
   window_end   = icemodel.verification.setup.ensureUtc( ...
      manifest.comparison_window.end);
end
