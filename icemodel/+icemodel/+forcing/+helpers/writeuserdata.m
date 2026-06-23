function filenames = writeuserdata(Data, site, source, kwargs)
   %WRITEUSERDATA Save a Data timetable as met-swap userdata files.
   %
   %  filenames = icemodel.forcing.helpers.writeuserdata(Data, site, source)
   %  filenames = ... writeuserdata(_, outdir=..., naming="window")
   %
   % Saves DATA under the icemodel met-swap ("userdata") naming convention.
   % naming="yearly" (default): one file per calendar year,
   %
   %    <site>_<source>_<YYYY>.mat
   %
   % naming="window": one file spanning the full time axis,
   %
   %    <site>_<source>_<YYYYMMDD>_<YYYYMMDD>.mat
   %
   % icemodel.loadmet resolves either form (a window file bracketing the run
   % year is preferred; otherwise the per-year file).
   %
   % consumed by icemodel.loadmet when opts.userdata / opts.uservars
   % request that variables of the met file be swapped with the
   % corresponding columns of the userdata file. Each .mat file holds
   % one variable named Data, a timetable carrying location metadata as
   % table CustomProperties (X, Y, Lat, Lon, Elev, Slope, ScalarUnits).
   %
   % OUTDIR defaults to icemodel.getpath('userdata') (demo/data/input/
   % userdata when the demo or test config is active) and is created
   % when it does not exist.
   %
   % Inputs
   %  Data   - timetable of evaluation/forcing variables with location
   %           CustomProperties attached
   %  site   - site name encoded in the filename (e.g. "kanm")
   %  source - data source encoded in the filename (e.g. "merra")
   %
   % Outputs
   %  filenames - string column of the full paths written
   %
   % See also: icemodel.forcing.helpers.writemet, icemodel.loadmet,
   %  icemodel.setopts

   arguments
      Data timetable
      site (1, 1) string
      source (1, 1) string
      kwargs.outdir (1, 1) string = ""
      kwargs.requiremetadata (1, 1) logical = true
      kwargs.naming (1, 1) string ...
         {mustBeMember(kwargs.naming, ["window", "yearly"])} = "yearly"
   end

   if kwargs.requiremetadata
      requireLocationMetadata(Data)
   end

   outdir = kwargs.outdir;
   if outdir == ""
      outdir = string(icemodel.getpath('userdata'));
   end
   if ~isfolder(outdir)
      mkdir(outdir)
   end

   switch kwargs.naming
      case "window"
         % One file spanning the full time axis. The encoded YYYYMMDD-YYYYMMDD
         % period lets icemodel.loadmet pick the file bracketing a run year.
         t1 = char(min(Data.Time), 'yyyyMMdd');
         t2 = char(max(Data.Time), 'yyyyMMdd');
         filenames = fullfile(outdir, ...
            sprintf('%s_%s_%s_%s.mat', site, source, t1, t2));
         savedata(filenames, Data)
      case "yearly"
         years_present = unique(year(Data.Time));
         filenames = strings(numel(years_present), 1);
         for n = 1:numel(years_present)
            yyyy = years_present(n);
            filenames(n) = fullfile(outdir, ...
               sprintf('%s_%s_%d.mat', site, source, yyyy));
            savedata(filenames(n), Data(year(Data.Time) == yyyy, :))
         end
   end
end

%% Local functions
function requireLocationMetadata(Data)
   %REQUIRELOCATIONMETADATA Assert the Data-file CustomProperties contract.
   needed = ["X", "Y", "Lat", "Lon", "Elev", "Slope", "ScalarUnits"];
   have = string(fieldnames(Data.Properties.CustomProperties));
   missing = setdiff(needed, have);
   if ~isempty(missing)
      error('icemodel:forcing:writeuserdata:missingMetadata', ...
         'Data is missing CustomProperties metadata: %s', ...
         strjoin(missing, ', '));
   end
end

function savedata(filename, Data)
   %SAVEDATA Save DATA to FILENAME as a variable named Data.
   S.Data = Data;
   save(filename, '-struct', 'S')
end
