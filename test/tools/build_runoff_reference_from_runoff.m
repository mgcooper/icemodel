function RunoffReference = build_runoff_reference_from_runoff(kwargs)
   %BUILD_RUNOFF_REFERENCE_FROM_RUNOFF Build runoff reference baseline from
   %runoff repo.
   %
   %  RunoffReference = build_runoff_reference_from_runoff()
   %  RunoffReference = build_runoff_reference_from_runoff(simyear=2016)
   %
   % With no SIMYEAR input, build references for all valid observation years at
   % each site. With SIMYEAR specified, build only rows for sites where that
   % year is valid.
   %
   % This utility is intended for one-time/manual baseline generation. It
   % requires access to the sibling `runoff` repository and its data layout.

   arguments (Input)
      kwargs.simyear double = double.empty()
      kwargs.runoff_root string = string.empty()
      kwargs.output_file string = string.empty()
   end

   % Deal out arguments.
   [simyear, runoff_root, output_file] = deal(kwargs.simyear, ...
      kwargs.runoff_root, kwargs.output_file);

   % Validate the optional year selector before building any paths.
   if ~isempty(simyear)
      mustBeInteger(simyear)
      mustBePositive(simyear)
      if ~isscalar(simyear)
         error('simyear must be scalar when provided')
      end
   end

   % Add the canonical repo/test source tree before resolving repo-local
   % paths for the runoff reference build. Keep the cleanup handle in scope
   % so path/config state is restored when this entrypoint returns.
   [rootdir, ~, ~, ~, suite_cleanup] = ...
      icemodel.test.helpers.bootstrapTestEnvironment(configure_paths=false);
   assert(isa(suite_cleanup, 'onCleanup'));
   testdir = icemodel.getpath('test');

   % Resolve the runoff repo root and default output file.
   if isblanktext(runoff_root)
      runoff_root = fullfile(fileparts(rootdir), 'runoff');
   end
   if isblanktext(output_file)
      if isempty(simyear)
         output_file = fullfile(testdir, 'references', ...
            'runoff_reference.mat');
      else
         output_file = fullfile(testdir, 'references', ...
            sprintf('runoff_reference_%d.mat', simyear));
      end
   end

   % Ensure the sibling runoff repo exists before touching its data or paths.
   if exist(char(runoff_root), 'dir') ~= 7
      error('runoff root not found: %s', char(runoff_root))
   end

   % Run the sibling loaders under the runoff repo's own data config, then
   % restore the caller's path/config state on exit.
   runoff_cleanup = bootstrapRunoffEnvironment(runoff_root); %#ok<NASGU>

   % Enumerate all supported runoff reference sites and forcing families.
   sites = {'behar', 'ak4', 'slv1', 'slv2', 'upperbasin'};
   families = {'local', 'mar'};

   % Build one reference row per site/year/forcing-family combination.
   rows = struct([]);
   k = 0;
   for is = 1:numel(sites)
      sitename = sites{is};
      siteopts = setBasinOpts('sitename', sitename, 'smbmodel', 'icemodel');
      simyears = getReferenceYears(siteopts, simyear);

      for iy = 1:numel(simyears)
         simyear_i = simyears(iy);
         for ifam = 1:numel(families)
            family = families{ifam};
            if strcmp(family, 'local')
               forcings = localForcingForSite(sitename);
            else
               forcings = 'mar';
            end

            % Load the requested comparison window and external references.
            [t1_req, t2_req] = getPlotWindow(siteopts, simyear_i);

            Discharge = loadDischarge( ...
               sitename, t1_req, t2_req, units="m3");

            Catchment = loadCatchment( ...
               sitename, int2str(simyear_i));

            [t1_ref, t2_ref] = getReferenceWindow( ...
               Discharge, sitename, t1_req, t2_req);

            [Mar, Merra, Racmo] = loadRunoff( ...
               sitename, {'mar', 'merra', 'racmo'}, Discharge.Time, ...
               filetype="Data", racmo_version=siteopts.racmo_version);

            obs_final = getObservedFinal( ...
               Discharge, t1_ref, t2_ref);

            mar_final = getRunoffFinal( ...
               Mar, Catchment.med.ease.area, t1_ref, t2_ref);

            merra_final = getRunoffFinal( ...
               Merra, Catchment.med.ease.area, t1_ref, t2_ref);

            racmo_final = getRunoffFinal( ...
               Racmo, Catchment.med.ease.area, t1_ref, t2_ref);

            % Save the derived scalar comparison context for this row.
            k = k + 1;
            rows(k).sitename = string(sitename);
            rows(k).forcings = string(forcings);
            rows(k).simyear = simyear_i;
            rows(k).t1 = t1_ref;
            rows(k).t2 = t2_ref;
            rows(k).area_med_m2 = getArea(Catchment, 'med');
            rows(k).area_min_m2 = getArea(Catchment, 'min');
            rows(k).area_max_m2 = getArea(Catchment, 'max');
            rows(k).obs_final_m3 = obs_final;
            rows(k).mar_final_m3 = mar_final;
            rows(k).merra_final_m3 = merra_final;
            rows(k).racmo_final_m3 = racmo_final;
            rows(k).notes = getNotes(simyear_i);
         end
      end
   end

   if k == 0
      error('No runoff reference rows generated for the requested settings')
   end

   % Convert the assembled rows and write the static reference file.
   RunoffReference = struct2table(rows);
   outdir = fileparts(char(output_file));
   if exist(outdir, 'dir') ~= 7
      mkdir(outdir);
   end
   save(char(output_file), 'RunoffReference');
end

function f = localForcingForSite(site)
   %LOCALFORCINGFORSITE Map runoff sites to their local forcing families.
   switch lower(site)
      case {'behar', 'slv1', 'slv2', 'kanm'}
         f = 'kanm';
      case {'ak4', 'upperbasin', 'kanl'}
         f = 'kanl';
      otherwise
         error('no local forcing mapping for site: %s', site)
   end
end

function simyears = getReferenceYears(opts, simyear)
   %GETREFERENCEYEARS Resolve the runoff-reference years to build.
   if isempty(simyear)
      simyears = opts.simyears(:).';
   elseif any(opts.simyears == simyear)
      simyears = simyear;
   else
      simyears = [];
   end
end

function [t1, t2] = getPlotWindow(opts, simyear)
   %GETPLOTWINDOW Resolve the requested reference comparison window.
   try
      t1 = datetime(simyear, opts.m1, opts.d1, opts.h1, 0, 0, 'TimeZone', 'UTC');
      t2 = datetime(simyear, opts.m2, opts.d2, opts.h2, 0, 0, 'TimeZone', 'UTC');
   catch
      t1 = datetime(simyear, 6, 1, 0, 0, 0, 'TimeZone', 'UTC');
      t2 = datetime(simyear, 9, 1, 0, 0, 0, 'TimeZone', 'UTC');
   end
end

function [t1, t2] = getReferenceWindow(Discharge, sitename, t1_req, t2_req)
   %GETREFERENCEWINDOW Clip the requested window to valid observations.
   iobs = find(~isnan(Discharge.QM));
   if isempty(iobs)
      error('No valid discharge observations found for site: %s', sitename)
   end

   if ismember(sitename, {'slv1', 'slv2'})
      t1 = Discharge.Time(iobs(1));
      t2 = Discharge.Time(iobs(end));
   else
      t1 = max(t1_req, Discharge.Time(iobs(1)));
      t2 = min(t2_req, Discharge.Time(iobs(end)));
   end
end

function x = getObservedFinal(Discharge, t1, t2)
   %GETOBSERVEDFINAL Compute the retained observed runoff change.
   ikeep = isbetween(Discharge.Time, t1, t2);
   q = Discharge.QM(ikeep);
   if isempty(q)
      x = nan;
      return
   end
   q0 = q(find(~isnan(q), 1, 'first'));
   q = q - q0;
   x = q(end);
end

function x = getRunoffFinal(Data, area_m2, t1, t2)
   %GETRUNOFFFINAL Compute the retained modeled runoff change in cubic meters.
   x = nan;

   if ~istimetable(Data) || ~ismember('runoff', Data.Properties.VariableNames)
      return
   end

   ikeep = isbetween(Data.Time, t1, t2);
   r = Data.runoff(ikeep);
   if isempty(r)
      return
   end

   r = area_m2 * cumsum(r);
   r = r - r(1);
   x = r(end);
end

function note = getNotes(simyear)
   %GETNOTES Build the provenance note for one runoff reference row.
   note = "Derived via loadDischarge + " + ...
      "loadRunoff(Data) using prep_runoff window logic";
   note = note + "; reference year " + simyear;
end

function a = getArea(Catchment, size_name)
   %GETAREA Read one catchment area scalar from the runoff struct.
   a = nan;
   if isfield(Catchment, size_name) ...
         && isfield(Catchment.(size_name), 'ease') ...
         && isfield(Catchment.(size_name).ease, 'area')
      a = Catchment.(size_name).ease.area;
   end
end

function cleanup = bootstrapRunoffEnvironment(runoff_root)
   %BOOTSTRAPRUNOFFENVIRONMENT Install runoff function paths and data config.

   % Snapshot the caller's current icemodel environment and MATLAB path.
   [cfg_prev, extra_names, extra_values, old_path] = snapshotEnvironment();

   % Expose the runoff function tree before calling its loaders.
   addpath(genpath(fullfile(char(runoff_root), 'functions')))

   % Point the shared ICEMODEL_* paths at the runoff repo's managed data.
   data_root = fullfile(char(runoff_root), 'data', 'icemodel');
   icemodel.config( ...
      'ICEMODEL_DATA_PATH', data_root, ...
      'ICEMODEL_INPUT_PATH', fullfile(data_root, 'input'), ...
      'ICEMODEL_OUTPUT_PATH', fullfile(data_root, 'output'), ...
      'ICEMODEL_EVAL_PATH', fullfile(data_root, 'eval'), ...
      'ICEMODEL_USERDATA_PATH', fullfile(data_root, 'input', 'userdata'));

   cleanup = onCleanup(@() restoreEnvironment( ...
      cfg_prev, extra_names, extra_values, old_path));
end

function [cfg_prev, extra_names, extra_values, old_path] = snapshotEnvironment()
   %SNAPSHOTENVIRONMENT Capture the caller's path/config state for restore.

   % Snapshot the current ICEMODEL_* config through the canonical getter.
   cfg_prev = icemodel.config('getenv', true);

   % Preserve the extra env vars that runoff startup code may use or mutate.
   extra_names = ["GEUSRUNOFFPATH", "ICEMODEL_PROJECT_PATH"];
   extra_values = cell(size(extra_names));
   for n = 1:numel(extra_names)
      extra_values{n} = getenv(extra_names(n));
   end

   % Restore the full MATLAB path so transient runoff paths do not leak.
   old_path = path();
end

function restoreEnvironment(cfg_prev, extra_names, extra_values, old_path)
   %RESTOREENVIRONMENT Restore the caller's original path/config state.

   % Restore the original MATLAB path before returning to the caller.
   path(old_path)

   % Restore the caller's previous ICEMODEL_* config exactly as it was.
   names = fieldnames(cfg_prev);
   for n = 1:numel(names)
      setenv(names{n}, cfg_prev.(names{n}));
   end

   % Restore the non-ICEMODEL extras that the runoff workflow may use.
   for n = 1:numel(extra_names)
      setenv(extra_names(n), extra_values{n});
   end
end
