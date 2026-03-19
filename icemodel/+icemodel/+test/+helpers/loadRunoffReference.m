function ref = loadRunoffReference(simyear, pathname)
   %LOADRUNOFFREFERENCE Load runoff reference baseline table.
   %
   %  ref = icemodel.test.helpers.loadRunoffReference()
   %  ref = icemodel.test.helpers.loadRunoffReference(2016)
   %  ref = icemodel.test.helpers.loadRunoffReference(2016, pathname)
   %
   % Returns:
   %  ref - table with columns:
   %    sitename, forcings, simyear, t1, t2, area_med_m2, area_min_m2,
   %    area_max_m2, obs_final_m3, mar_final_m3, merra_final_m3, racmo_final_m3,
   %    notes

   if nargin < 1
      simyear = [];
   end
   if nargin < 2 || isempty(pathname)
      pathname = resolveDefaultPaths(icemodel.getpath('test'), simyear);
   end

   if isempty(pathname)
      error('runoff reference file not found in test/references')
   end

   if iscell(pathname)
      refs = cellfun(@(p) normalizeReference(loadReference(p)), pathname, ...
         'UniformOutput', false);
      ref = vertcat(refs{:});
   else
      if exist(pathname, 'file') ~= 2
         error('runoff reference file not found: %s', pathname)
      end
      ref = normalizeReference(loadReference(pathname));
   end

   if ~isempty(simyear)
      ref = ref(ref.simyear == simyear, :);
   end
end

function pathname = resolveDefaultPaths(testdir, simyear)
   %RESOLVEDEFAULTPATHS Resolve the default runoff-reference file(s) to load.
   refdir = fullfile(testdir, 'references');
   pathname = fullfile(refdir, 'runoff_reference.mat');
   if exist(pathname, 'file') == 2
      return
   end

   if ~isempty(simyear)
      pathname = fullfile(refdir, ...
         sprintf('runoff_reference_%d.mat', simyear));
      if exist(pathname, 'file') == 2
         return
      end
   end

   files = dir(fullfile(refdir, 'runoff_reference_*.mat'));
   if isempty(files)
      pathname = '';
   else
      pathname = fullfile({files.folder}, {files.name});
   end
end

function raw = loadReference(pathname)
   %LOADREFERENCE Load one runoff reference MAT file into a table.
   raw = icemodel.test.helpers.loadSavedTable(pathname, ["RunoffReference", "ref"]);
end

function ref = normalizeReference(raw)
   %NORMALIZEREFERENCE Normalize saved runoff reference columns and types.
   ref = raw;

   required = ["sitename", "forcings", "simyear", "t1", "t2", ...
      "area_med_m2", "area_min_m2", "area_max_m2", "obs_final_m3", ...
      "mar_final_m3", "merra_final_m3", "racmo_final_m3", "notes"];

   for i = 1:numel(required)
      v = required(i);
      if ~ismember(v, string(ref.Properties.VariableNames))
         switch v
            case {"sitename", "forcings", "notes"}
               ref.(v) = repmat("", height(ref), 1);
            case {"t1", "t2"}
               ref.(v) = repmat(NaT, height(ref), 1);
            otherwise
               ref.(v) = nan(height(ref), 1);
         end
      end
   end

   ref.sitename = string(ref.sitename);
   ref.forcings = string(ref.forcings);
   ref.notes = string(ref.notes);

   if ~isdatetime(ref.t1)
      ref.t1 = datetime(ref.t1, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
   end
   if ~isdatetime(ref.t2)
      ref.t2 = datetime(ref.t2, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
   end

   ref = movevars(ref, required, 'Before', 1);
end
