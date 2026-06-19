function filenames = writemet(met, site, forcings, kwargs)
   %WRITEMET Validate and save an icemodel met file.
   %
   %  filenames = icemodel.forcing.helpers.writemet(met, site, forcings)
   %  filenames = ... writemet(_, outdir=..., naming="yearly", validate=false)
   %
   % Saves MET (a timetable satisfying the icemodel met contract) under
   % the standard naming convention so icemodel.createMetFileNames and
   % icemodel.loadmet resolve it without special cases:
   %
   %  naming="window" (default): one file spanning the full time axis,
   %     met_<site>_<forcings>_<YYYYMMDD>_<YYYYMMDD>_<dt>.mat
   %
   %  naming="yearly": one file per calendar year (legacy form),
   %     met_<site>_<forcings>_<YYYY>_<dt>.mat
   %
   % OUTDIR defaults to fullfile(icemodel.getpath('input'), 'met'), the
   % met directory of the active icemodel workspace (demo/data when the
   % demo or test config is active). The directory is created when it
   % does not exist. The saved .mat file holds one variable named met.
   %
   % Inputs
   %  met      - timetable of forcing variables (see
   %             icemodel.forcing.helpers.metvariables)
   %  site     - site name encoded in the filename (e.g. "kanm")
   %  forcings - forcing source encoded in the filename (e.g. "mar")
   %
   % Outputs
   %  filenames - string column of the full paths written
   %
   % See also: icemodel.forcing.helpers.metfilename,
   %  icemodel.forcing.helpers.validatemet, icemodel.loadmet

   arguments
      met timetable
      site (1, 1) string
      forcings (1, 1) string
      kwargs.outdir (1, 1) string = ""
      kwargs.naming (1, 1) string {mustBeMember(kwargs.naming, ...
         ["window", "yearly"])} = "window"
      kwargs.validate (1, 1) logical = true
   end

   if kwargs.validate
      icemodel.forcing.helpers.validatemet(met)
   end

   outdir = kwargs.outdir;
   if outdir == ""
      outdir = string(fullfile(icemodel.getpath('input'), 'met'));
   end
   if ~isfolder(outdir)
      mkdir(outdir)
   end

   % Infer the timestep (seconds) from the (validated, uniform) time axis.
   dt = seconds(met.Time(2) - met.Time(1));

   switch kwargs.naming
      case "window"
         name = icemodel.forcing.helpers.metfilename(site, forcings, ...
            met.Time(1), met.Time(end), dt);
         filenames = fullfile(outdir, name);
         savemet(filenames, met)

      case "yearly"
         years_present = unique(year(met.Time));
         filenames = strings(numel(years_present), 1);
         for n = 1:numel(years_present)
            yyyy = years_present(n);
            name = icemodel.forcing.helpers.metfilename(site, forcings, ...
               yyyy, [], dt);
            filenames(n) = fullfile(outdir, name);
            savemet(filenames(n), met(year(met.Time) == yyyy, :))
         end
   end
end

%% Local functions
function savemet(filename, met)
   %SAVEMET Save MET to FILENAME as a variable named met.
   S.met = met;
   save(filename, '-struct', 'S')
end
