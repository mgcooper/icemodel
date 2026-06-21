function coverage = promiceSourceCoverage(models, dirs)
   %PROMICESOURCECOVERAGE Probe on-disk year coverage of each forcing source.
   %
   %  coverage = icemodel.verification.setup.promiceSourceCoverage( ...
   %     models, dirs)
   %
   % Probes the available calendar-year span of each requested forcing source
   % directly from the files on disk, so the staging driver can decouple each
   % leg's comparison window from a single hardcoded period. The probe never
   % fabricates coverage: a source that is absent or unreadable returns an
   % empty year set with a reason.
   %
   % Inputs
   %  models  string vector subset of ["promice","mar","merra","racmo"].
   %  dirs    struct with the (possibly empty) source directory for each
   %          model: dirs.promice, dirs.mar, dirs.merra, dirs.racmo. Empty
   %          fields fall back to each builder's reference layout (the same
   %          defaults the builders use).
   %
   % Outputs
   %  coverage  struct keyed by model. Each entry has:
   %     years     double row vector of available calendar years (may be empty)
   %     year_min  first available year (NaN when none)
   %     year_max  last available year (NaN when none)
   %     reason    "" when coverage was found, else why it is empty
   %
   % RACMO carries TWO physical archives on disk (the FGRN11 surface
   % no_subsurf run and the subsurface run) whose spans differ; this probe
   % reports the UNION span of the per-variable RACMO files found, which is
   % the span buildRacmoData can actually deliver from that directory.
   %
   % See also: icemodel.verification.setup.importPromiceSites

   arguments
      models (1, :) string
      dirs (1, 1) struct
   end

   coverage = struct();

   if ismember("mar", models)
      coverage.mar = probeMar(getfielddefault(dirs, "mar"));
   end
   if ismember("merra", models)
      coverage.merra = probeMerra(getfielddefault(dirs, "merra"));
   end
   if ismember("racmo", models)
      coverage.racmo = probeRacmo(getfielddefault(dirs, "racmo"));
   end
end

%% Local functions
function value = getfielddefault(s, name)
   if isfield(s, name)
      value = string(s.(name));
   else
      value = "";
   end
end

function entry = emptyCoverage(reason)
   entry = struct('years', [], 'year_min', NaN, 'year_max', NaN, ...
      'reason', string(reason));
end

function entry = yearsCoverage(years)
   years = unique(years(:)');
   entry = struct('years', years, 'year_min', min(years), ...
      'year_max', max(years), 'reason', "");
end

function entry = probeMar(source_dir)
   %PROBEMAR On-disk MAR yearly-file coverage (one *-YYYY.nc per year).
   if source_dir == ""
      source_dir = string(fullfile(icemodel.getpath('data'), 'forcing', 'mar'));
      if ~isfolder(source_dir)
         source_dir = "/Volumes/S03/DATA/greenland/mar3p11/RUH2";
      end
   end
   if ~isfolder(source_dir)
      entry = emptyCoverage(sprintf('MAR directory not found: %s', source_dir));
      return
   end
   files = dir(fullfile(source_dir, 'MARv3.11*.nc'));
   if isempty(files)
      files = dir(fullfile(source_dir, '*.nc'));
   end
   tokens = regexp({files.name}, '-(\d{4})\.nc$', 'tokens', 'once');
   have = ~cellfun(@isempty, tokens);
   if ~any(have)
      entry = emptyCoverage(sprintf('no MAR *-YYYY.nc files in %s', source_dir));
      return
   end
   years = cellfun(@(t) str2double(t{1}), tokens(have));
   entry = yearsCoverage(years);
end

function entry = probeMerra(source_dir)
   %PROBEMERRA On-disk MERRA-2 coverage (intersection across collections).
   if source_dir == ""
      source_dir = string(fullfile(icemodel.getpath('data'), 'forcing', 'merra2'));
      if ~isfolder(source_dir)
         source_dir = "/Volumes/S03/DATA/merra2/1hrly/ncfiles";
      end
   end
   if ~isfolder(source_dir)
      entry = emptyCoverage(sprintf( ...
         'MERRA-2 directory not found: %s', source_dir));
      return
   end
   % buildMerraData needs every collection (flx/rad/slv/glc); the deliverable
   % span is the INTERSECTION of the years present in each collection.
   collections = ["flx", "rad", "slv", "glc"];
   present = collections(arrayfun( ...
      @(c) isfolder(fullfile(source_dir, c)), collections));
   if isempty(present)
      entry = emptyCoverage(sprintf( ...
         'no MERRA-2 collection subdirectories in %s', source_dir));
      return
   end
   years = [];
   for c = present
      files = dir(fullfile(source_dir, c, '*_Nx.*.nc4*'));
      tokens = regexp({files.name}, '_Nx\.(\d{4})\d{4}\.', 'tokens', 'once');
      have = ~cellfun(@isempty, tokens);
      cyears = unique(cellfun(@(t) str2double(t{1}), tokens(have)));
      if c == present(1)
         years = cyears;
      else
         years = intersect(years, cyears);
      end
   end
   if isempty(years)
      entry = emptyCoverage(sprintf( ...
         'no common MERRA-2 years across collections in %s', source_dir));
      return
   end
   entry = yearsCoverage(years);
end

function entry = probeRacmo(source_dir)
   %PROBERACMO On-disk RACMO coverage (union of FGRN11 archive spans).
   if source_dir == ""
      source_dir = string(fullfile(icemodel.getpath('data'), 'forcing', 'racmo'));
      if ~isfolder(source_dir)
         source_dir = "/Volumes/S03/DATA/greenland/racmo2p3/surface";
      end
   end
   if ~isfolder(source_dir)
      entry = emptyCoverage(sprintf('RACMO directory not found: %s', source_dir));
      return
   end
   files = dir(fullfile(source_dir, '*.RACMO23p3_*.nc'));
   if isempty(files)
      entry = emptyCoverage(sprintf( ...
         'no RACMO *.RACMO23p3_*.nc files in %s', source_dir));
      return
   end
   % Filenames encode the archive span as FGRN11_<y0>_<y1>; the deliverable
   % span is the union of the per-file spans present.
   tokens = regexp({files.name}, 'FGRN11_(\d{4})_(\d{4})', 'tokens', 'once');
   have = ~cellfun(@isempty, tokens);
   if ~any(have)
      entry = emptyCoverage(sprintf( ...
         'no FGRN11_<y0>_<y1> span in RACMO files under %s', source_dir));
      return
   end
   years = [];
   for tok = tokens(have)
      y0 = str2double(tok{1}{1});
      y1 = str2double(tok{1}{2});
      years = [years, y0:y1]; %#ok<AGROW>
   end
   entry = yearsCoverage(years);
end
