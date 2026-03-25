function ref = loadReference(name, kwargs)
   %LOADREFERENCE Load a test reference table.
   %
   %  ref = icemodel.test.helpers.loadReference("runoff")
   %  ref = icemodel.test.helpers.loadReference("runoff", simyear=2016)
   %  ref = icemodel.test.helpers.loadReference("runoff", ...
   %     filename="/path/to/file.mat")
   %
   % The loaded file path is printed to the console.

   arguments
      name (1, :) string {mustBeMember(name, ["runoff"])} ...
         = "runoff"

      kwargs.simyear double ...
         = []

      kwargs.filename string {mustBeTextScalarOrEmpty} ...
         = string.empty()
   end

   simyear = kwargs.simyear;

   if ~isblanktext(kwargs.filename)
      pathname = char(kwargs.filename);
   else
      pathname = char(icemodel.test.helpers.referenceFilePath( ...
         name, simyear=simyear));
   end

   switch name
      case "runoff"
         ref = loadRunoff(pathname, simyear);
   end

   % Print the loaded filename to the console.
   icemodel.test.helpers.printFilePath(pathname, "load");
end

function ref = loadRunoff(pathname, simyear)
   %LOADRUNOFF Load and normalize a runoff reference file.

   if iscell(pathname)
      refs = cellfun(@(p) normalizeRunoff(loadOneRunoff(p)), pathname, ...
         'UniformOutput', false);
      ref = vertcat(refs{:});
   else
      if exist(pathname, 'file') ~= 2
         error('icemodel:test:referenceNotFound', ...
            'Reference file not found: %s', pathname)
      end
      ref = normalizeRunoff(loadOneRunoff(pathname));
   end

   if ~isempty(simyear)
      ref = ref(ref.simyear == simyear, :);
   end
end

function raw = loadOneRunoff(pathname)
   %LOADONERUNOFF Load one runoff reference MAT file.
   raw = icemodel.test.helpers.loadSavedTable( ...
      pathname, ["RunoffReference", "ref"]);
end

function ref = normalizeRunoff(raw)
   %NORMALIZERUNOFF Normalize saved runoff reference columns and types.
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
      ref.t1 = datetime(ref.t1, 'ConvertFrom', 'datenum', ...
         'TimeZone', 'UTC');
   end
   if ~isdatetime(ref.t2)
      ref.t2 = datetime(ref.t2, 'ConvertFrom', 'datenum', ...
         'TimeZone', 'UTC');
   end

   ref = movevars(ref, required, 'Before', 1);
end
