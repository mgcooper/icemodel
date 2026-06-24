function metfname = createMetFileNames(opts)
   %CREATEMETFILENAMES Create icemodel met file names for model OPTS.
   %
   %  metfname = icemodel.createMetFileNames(opts)
   %
   % Naming forms (the window form is preferred when both are available):
   %   Per-year (simyears only):
   %     met_<met>_<forcings>_<YYYY>_<dt>.mat   (one file per simyear)
   %
   %   Window-stamped (opts.startdate / opts.enddate set):
   %     met_<met>_<forcings>_<YYYYMMDD>_<YYYYMMDD>_<dt>.mat   (single file)
   %
   % When the run window is set but no exact-window file is staged, the
   % function looks for a multi-year staged file whose encoded period
   % brackets the requested window and returns that filename instead.
   % This lets a single full-period met file serve many shorter sub-runs
   % (e.g. an ESM-SnowMIP site's full 1994-2014 file covers any
   % single-year smoke window).
   %
   % See also: icemodel.setopts icemodel.configureRun

   sitename = opts.sitename;
   forcings = opts.forcings;
   simyears = opts.simyears;

   % Deal with the case where met-station forcing data (as opposed to gridded
   % climate model forcing data) is requested for a nearby catchment by
   % replacing the catchment name in the metfile with the met station name.
   % For example, if sitename=="behar" and forcingdata=="kanm", this sets the
   % metfile name to met_kanm_kanm_YYYY rather than met_behar_kanm_YYYY, to
   % negate the need to create a second (identical) met_behar_kanm_YYYY file.
   if strcmpi(forcings, 'kanl') ...
         && ismember(sitename, {'ak4', 'upperbasin'})
      metname = 'kanl';
   elseif strcmpi(forcings, 'kanm') ...
         && ismember(sitename, {'slv1', 'slv2', 'behar'})
      metname = 'kanm';
   else
      metname = sitename;
   end

   switch opts.dt
      case 900
         dtstr = '15m.mat';
      case 3600
         dtstr = '1hr.mat';
      otherwise
         error('unsupported dt for met file naming: %g', opts.dt)
   end

   % Prefer the window form when both startdate and enddate are set
   % (regardless of whether simyears is also set). The window form is
   % the authoritative encoding for opts.startdate / opts.enddate runs,
   % including verification-suite multi-year met files.
   has_window = isfield(opts, 'startdate') && isfield(opts, 'enddate') ...
      && ~isempty(opts.startdate) && ~isempty(opts.enddate) ...
      && ~any(isnat(opts.startdate)) && ~any(isnat(opts.enddate));

   if has_window
      start_stamp = char(opts.startdate, 'yyyyMMdd');
      end_stamp = char(opts.enddate, 'yyyyMMdd');
      exact_name = ['met_' metname '_' forcings '_' ...
         start_stamp '_' end_stamp '_' dtstr];

      enclosing = findEnclosingMetFile(opts, metname, forcings, ...
         opts.startdate, opts.enddate, dtstr, exact_name);
      if ~isempty(enclosing)
         metfname = {enclosing};
      else
         metfname = {exact_name};
      end
      return
   end

   metfname = cell(1, numel(simyears));
   for n = 1:numel(simyears)
      simyear = num2str(simyears(n));
      metfname{n} = ['met_' metname '_' forcings '_' simyear '_' dtstr];
   end
end

function name = findEnclosingMetFile(opts, metname, forcings, ...
      startdate, enddate, dtstr, exact_name)
   % Locate a staged multi-year met file whose encoded YYYYMMDD-YYYYMMDD
   % period contains the requested run window. Returns '' when no match.
   % Delegates the glob/parse/bracket to the shared primitive
   % icemodel.forcing.helpers.findEnclosingWindowFile (also used by
   % icemodel.loadmet for userdata files).

   name = '';
   if ~isfield(opts, 'pathinput') || isempty(opts.pathinput)
      return
   end
   met_base = fullfile(opts.pathinput, 'met');
   prefix = ['met_' metname '_' forcings];
   % Look in the per-source subfolder (met/<forcings>/) first, then the flat
   % met dir, so a staged enclosing window file is found in either layout.
   for d = {fullfile(met_base, char(forcings)), met_base}
      met_dir = d{1};
      if ~isfolder(met_dir)
         continue
      end
      if isfile(fullfile(met_dir, exact_name))
         return   % exact file exists; caller uses exact_name (path resolved later)
      end
      enclosing = icemodel.forcing.helpers.findEnclosingWindowFile(met_dir, ...
         prefix, ['_' dtstr], startdate, enddate);
      if strlength(enclosing) > 0
         name = char(enclosing);
         return
      end
   end
end
