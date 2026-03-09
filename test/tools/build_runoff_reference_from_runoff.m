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

   % Resolve the repo root, runoff repo root, and default output file.
   rootdir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
   if strlength(runoff_root) == 0
      runoff_root = fullfile(fileparts(rootdir), 'runoff');
   end
   if strlength(output_file) == 0
      if isempty(simyear)
         output_file = fullfile(rootdir, 'test', 'references', ...
            'runoff_reference.mat');
      else
         output_file = fullfile(rootdir, 'test', 'references', ...
            sprintf('runoff_reference_%d.mat', simyear));
      end
   end

   % Ensure the sibling runoff repo exists, then expose its functions.
   if exist(char(runoff_root), 'dir') ~= 7
      error('runoff root not found: %s', runoff_root)
   end

   addpath(genpath(fullfile(char(runoff_root), 'functions')));

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
            rows(k).sitename = string(sitename); %#ok<AGROW>
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
   RunoffReference = struct2table(rows); %#ok<NASGU>
   outdir = fileparts(char(output_file));
   if exist(outdir, 'dir') ~= 7
      mkdir(outdir);
   end
   save(char(output_file), 'RunoffReference');
end

function f = localForcingForSite(site)
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
   if isempty(simyear)
      simyears = opts.simyears(:).';
   elseif any(opts.simyears == simyear)
      simyears = simyear;
   else
      simyears = [];
   end
end

function [t1, t2] = getPlotWindow(opts, simyear)
   try
      t1 = datetime(simyear, opts.m1, opts.d1, opts.h1, 0, 0, 'TimeZone', 'UTC');
      t2 = datetime(simyear, opts.m2, opts.d2, opts.h2, 0, 0, 'TimeZone', 'UTC');
   catch
      t1 = datetime(simyear, 6, 1, 0, 0, 0, 'TimeZone', 'UTC');
      t2 = datetime(simyear, 9, 1, 0, 0, 0, 'TimeZone', 'UTC');
   end
end

function [t1, t2] = getReferenceWindow(Discharge, sitename, t1_req, t2_req)
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
   note = "Derived via loadDischarge + " + ...
      "loadRunoff(Data) using prep_runoff window logic";
   note = note + "; reference year " + simyear;
end

function a = getArea(Catchment, size_name)
   a = nan;
   if isfield(Catchment, size_name) ...
         && isfield(Catchment.(size_name), 'ease') ...
         && isfield(Catchment.(size_name).ease, 'area')
      a = Catchment.(size_name).ease.area;
   end
end
