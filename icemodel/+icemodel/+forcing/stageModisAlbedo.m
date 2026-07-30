function info = stageModisAlbedo(sites, kwargs)
   %STAGEMODISALBEDO Stage per-site GEUS MODIS daily albedo userdata artifacts.
   %
   %  info = icemodel.forcing.stageModisAlbedo(sites)
   %  info = ... stageModisAlbedo(sites, modis_dir=..., met_dir=..., ...
   %     met_source="promice", outdir=..., years=..., overwrite=false)
   %
   % Extracts the GEUS MODIS C6 daily albedo (Greenland_Reflectivity_
   % <YYYY>_5km_C6.nc) at each requested station and writes one
   % window-stamped userdata artifact per site through the canonical
   % icemodel.forcing.helpers.writeuserdata path (per-source subfolder
   % userdata/modis/, timetable named Data, native daily cadence, top-level
   % artifact_cadence_seconds). This is the staging half of the gap-fill
   % reconstruction MODIS tier (reconstruct/POLICY.md B12, D-15); runtime
   % attachment onto a met axis is the job of the single conversion helper
   % icemodel.forcing.modisToMetCadence, never of the model runtime.
   %
   % Extraction method: single NEAREST 5 km cell per station (readGeusModis
   % method "nearest"). The 5 km cell already integrates a footprint far
   % larger than the AWS-local surface, and a neighbourhood mean would blend
   % across the sharp ice-margin albedo gradients (ocean/tundra/bare-ice
   % mixing) that many PROMICE lowland stations sit on, biasing the series.
   % Nearest-cell also matches the default point sampling of the gridded
   % builders' MODIS channel, so staged and built MODIS resolve identically.
   % The method, selected cell indices, cell centre, and query-to-cell
   % distance are recorded in each artifact's metadata.
   %
   % Source loop: each yearly NetCDF is opened ONCE and every station is
   % extracted from that single covering read (readGeusModis point rows),
   % and each file is byte-pinned (size + SHA-256) into every artifact's
   % metadata so a changed source is detectable (POLICY A1 spirit).
   %
   % Physical validity: finite extracted values outside the reconstruction
   % albedo bounds (icemodel.forcing.reconstruct.physicalBounds("albedo"))
   % are masked to NaN before writing, so artifacts never carry values the
   % admission gates would reject.
   %
   % Station identity: each site's location is copied from its staged met
   % artifact met_<site>_<met_source>_*.mat so the MODIS artifact is
   % colocated with the met product it will reconstruct. The source-light
   % top-level artifact_metadata record (lat/lon/elev, the gap-filled
   % PROMICE convention) is preferred; met files carrying the full location
   % CustomProperties contract (X, Y, Lat, Lon, Elev, Slope) resolve from
   % the met timetable instead.
   %
   % Inputs
   %  sites - string row of met-site tokens (e.g. ["kanm" "kanl"])
   %
   % Name-value
   %  modis_dir  - directory of yearly GEUS reflectivity NetCDFs. Default ""
   %               resolves the ICEMODEL_MODIS_DIR environment variable;
   %               unresolved is an error (no silent partial staging).
   %  met_dir    - directory of staged met artifacts supplying station
   %               identity. Default: <input>/met/<met_source>.
   %  met_source - met filename source token (default "promice")
   %  outdir     - userdata root passed to writeuserdata (default "" lets
   %               the writer resolve icemodel.getpath('userdata'))
   %  years      - calendar years to stage. Default [] stages every year
   %               discovered in modis_dir; explicitly requested years must
   %               all be present or the call errors.
   %  overwrite  - replace existing artifacts (default false: canonical
   %               writer reuse/conflict rules apply)
   %
   % Outputs
   %  info - per-site struct row: site, filename, n_days, n_finite,
   %         first_finite, last_finite, coverage_years, cell_start,
   %         distance_m
   %
   % See also: icemodel.forcing.readGeusModis,
   %  icemodel.forcing.modisToMetCadence,
   %  icemodel.forcing.helpers.writeuserdata,
   %  icemodel.forcing.reconstruct.physicalBounds

   arguments
      sites (1, :) string {mustBeNonempty}
      kwargs.modis_dir (1, 1) string = ""
      kwargs.met_dir (1, 1) string = ""
      kwargs.met_source (1, 1) string = "promice"
      kwargs.outdir (1, 1) string = ""
      kwargs.years (1, :) double = []
      kwargs.overwrite (1, 1) logical = false
   end

   % Station tokens enter globs and output paths below, so reject malformed
   % identifiers before resolving or scanning any filesystem location.
   icemodel.forcing.reconstruct.mustBeStationToken(sites)

   % Resolve the source and met directories up front so every failure is a
   % clear precondition error rather than a mid-loop surprise.
   modis_dir = resolveModisDir(kwargs.modis_dir);
   met_dir = kwargs.met_dir;
   if met_dir == ""
      met_dir = string(fullfile(icemodel.getpath('input'), 'met', ...
         char(kwargs.met_source)));
   end

   % Inventory the yearly source files; the on-disk inventory (not a
   % hard-coded window) defines the default staged coverage.
   [files, file_years] = sourceInventory(modis_dir, kwargs.years);

   % Station identity comes from the staged met artifacts so the MODIS
   % artifact carries the exact colocated CustomProperties contract.
   locations = metLocations(sites, met_dir, kwargs.met_source);
   points = [[locations.lat_wgs84]', [locations.lon_wgs84]'];

   % One strictly regular daily axis spanning the staged years: dates with no
   % source coverage stay NaN, which keeps the artifact cadence uniform (so
   % artifact_cadence_seconds stamps) and the window stamp honest.
   t0 = datetime(min(file_years), 1, 1, 'TimeZone', 'UTC');
   t1 = datetime(max(file_years), 12, 31, 'TimeZone', 'UTC');
   Time = (t0:days(1):t1)';
   nsites = numel(sites);
   albedo = nan(numel(Time), nsites);

   % Read each yearly file ONCE, extracting all stations from one covering
   % hyperslab, and byte-pin the file (size + sha256) for the artifacts.
   source_files = repmat(struct('file', '', 'year', 0, ...
      'size_bytes', 0, 'sha256', ''), numel(files), 1);
   selection = [];
   for n = 1:numel(files)
      [year_albedo, year_time, year_selection] = ...
         icemodel.forcing.readGeusModis(files(n), points, "nearest");

      % The 5 km grid must be identical across the yearly files; a mismatch
      % means the directory mixes products and the cell pinning is invalid.
      if isempty(selection)
         selection = year_selection;
      elseif ~isequal(year_selection(1).grid_size, ...
            selection(1).grid_size) ...
            || string(year_selection(1).grid_sha256) ...
            ~= string(selection(1).grid_sha256)
         error('icemodel:forcing:stageModisAlbedo:gridMismatch', ...
            'source grid of %s differs from earlier yearly files', files(n))
      end

      % Place this file's days onto the full axis by date. Both axes are
      % daily UTC datetimes, so the offsets are exact integers. The C6 files
      % hold a fixed 366-day record, so a non-leap file carries one trailing
      % day belonging to the next calendar year: interior overhangs are
      % overwritten by the next file's native day (files are processed in
      % chronological order) and days beyond the staged axis are dropped.
      rows = round(days(year_time - Time(1))) + 1;
      keep = rows >= 1 & rows <= numel(Time);
      albedo(rows(keep), :) = year_albedo(keep, :);

      source_files(n) = struct( ...
         'file', char(files(n)), ...
         'year', file_years(n), ...
         'size_bytes', getfield(dir(files(n)), 'bytes'), ...
         'sha256', char(icemodel.verification.setup.fileSha256(files(n))));
   end

   % Mask finite values outside the reconstruction albedo bounds (SSOT; do
   % not restate the limits here) so staged artifacts are gate-clean.
   bounds = icemodel.forcing.reconstruct.physicalBounds("albedo");
   outside = isfinite(albedo) & (albedo < bounds(1) | albedo > bounds(2));
   albedo(outside) = NaN;

   % Write one window-stamped native-daily artifact per site through the
   % canonical writer, carrying the extraction/pinning provenance.
   axis_years = unique(year(Time))';
   info = repmat(struct('site', "", 'filename', "", 'n_days', 0, ...
      'n_finite', 0, 'first_finite', NaT('TimeZone', 'UTC'), ...
      'last_finite', NaT('TimeZone', 'UTC'), 'coverage_years', [], ...
      'cell_start', [0 0], 'distance_m', NaN), 1, nsites);
   for s = 1:nsites
      values = albedo(:, s);
      Data = timetable(Time, values, 'VariableNames', {'albedo'});

      % Reuse the canonical location attachment so the CustomProperties
      % schema exactly matches every other Data builder.
      Data = icemodel.forcing.helpers.attachLocationMetadata( ...
         Data, locations(s).location);

      % Coverage years are the axis years holding at least one finite value
      % for THIS site, classified by the canonical coverage contract.
      finite_here = isfinite(values);
      coverage_years = axis_years(arrayfun( ...
         @(y) any(finite_here & year(Time) == y), axis_years));
      metadata = icemodel.forcing.helpers.geusModisCoverageMetadata( ...
         axis_years, coverage_years);

      % Extraction + pinning provenance: how each value was selected, from
      % which bytes, under which bounds, and how to attach it downstream.
      metadata.extraction_method = selection(s).method;
      metadata.extraction_cell_start = selection(s).start;
      metadata.extraction_cell_count = selection(s).count;
      metadata.extraction_cell_lat = selection(s).cell_lat;
      metadata.extraction_cell_lon = selection(s).cell_lon;
      metadata.extraction_cell_distance_m = selection(s).distance_m;
      metadata.source_grid_size = selection(s).grid_size;
      metadata.source_grid_sha256 = selection(s).grid_sha256;
      metadata.albedo_physical_bounds = bounds;
      metadata.albedo_bounds_masked_count = sum(outside(:, s));
      metadata.met_identity_file = char(locations(s).met_file);
      metadata.source_files = source_files;
      metadata.attachment_helper = 'icemodel.forcing.modisToMetCadence';
      Data.Properties.UserData = metadata;

      % Canonical writer: window naming, explicit native (daily) cadence.
      filename = icemodel.forcing.helpers.writeuserdata(Data, sites(s), ...
         "modis", outdir=kwargs.outdir, naming="window", dt_out="", ...
         overwrite=kwargs.overwrite);

      % Compact per-site summary for the caller's coverage report.
      info(s).site = sites(s);
      info(s).filename = string(filename);
      info(s).n_days = numel(values);
      info(s).n_finite = sum(finite_here);
      if any(finite_here)
         info(s).first_finite = Time(find(finite_here, 1, 'first'));
         info(s).last_finite = Time(find(finite_here, 1, 'last'));
      end
      info(s).coverage_years = coverage_years;
      info(s).cell_start = selection(s).start;
      info(s).distance_m = selection(s).distance_m;
   end
end

function modis_dir = resolveModisDir(modis_dir)
   %RESOLVEMODISDIR Require an explicit or environment-configured source dir.

   % The caller's explicit path wins; otherwise the documented environment
   % override. There is no hard-coded volume fallback here: the source
   % location is site/machine configuration, not code.
   if modis_dir == ""
      modis_dir = string(getenv("ICEMODEL_MODIS_DIR"));
   end
   if modis_dir == "" || ~isfolder(modis_dir)
      error('icemodel:forcing:stageModisAlbedo:modisDirNotFound', ...
         ['MODIS source directory not found. Pass modis_dir=... or set ' ...
         'ICEMODEL_MODIS_DIR (got "%s").'], modis_dir)
   end
end

function [files, file_years] = sourceInventory(modis_dir, years)
   %SOURCEINVENTORY Discover and pin the yearly reflectivity files to stage.

   listing = dir(fullfile(modis_dir, 'Greenland_Reflectivity_*_C6.nc'));
   if isempty(listing)
      error('icemodel:forcing:stageModisAlbedo:noSourceFiles', ...
         'no Greenland_Reflectivity_*_C6.nc files under %s', modis_dir)
   end

   % Parse each file year from its basename, the same token readGeusModis
   % trusts for the time axis. Loop element-wise: regexp over a string array
   % nests its token output differently for one versus many names.
   file_years = zeros(1, numel(listing));
   for n = 1:numel(listing)
      tok = regexp(listing(n).name, '_(\d{4})_', 'tokens', 'once');
      if isempty(tok)
         error('icemodel:forcing:stageModisAlbedo:unparsableSourceName', ...
            'cannot parse a file year from %s under %s', ...
            listing(n).name, modis_dir)
      end
      file_years(n) = str2double(tok{1});
   end
   if numel(unique(file_years)) ~= numel(file_years)
      error('icemodel:forcing:stageModisAlbedo:ambiguousSourceYears', ...
         'multiple source files share a year under %s', modis_dir)
   end

   % An explicit year request is a contract: every requested year must have
   % a source file, so a mount hiccup cannot silently stage a subset.
   if ~isempty(years)
      missing = setdiff(years, file_years);
      if ~isempty(missing)
         error('icemodel:forcing:stageModisAlbedo:sourceYearMissing', ...
            'no source file for requested year(s) %s under %s', ...
            strjoin(string(missing), ', '), modis_dir)
      end
      keep = ismember(file_years, years);
      listing = listing(keep);
      file_years = file_years(keep);
   end

   % Stage in chronological order so the axis fill and pin list read plainly.
   [file_years, order] = sort(file_years);
   listing = listing(order);
   files = string(fullfile({listing.folder}, {listing.name}));
end

function locations = metLocations(sites, met_dir, met_source)
   %METLOCATIONS Copy each site's location contract from its met artifact.

   locations = repmat(struct('met_file', "", 'lat_wgs84', NaN, ...
      'lon_wgs84', NaN, 'location', struct()), 1, numel(sites));
   for s = 1:numel(sites)
      % Exactly one met artifact per site keeps the identity unambiguous.
      pattern = sprintf('met_%s_%s_*.mat', sites(s), met_source);
      match = dir(fullfile(met_dir, pattern));
      if isempty(match)
         error('icemodel:forcing:stageModisAlbedo:metArtifactNotFound', ...
            'no met artifact %s under %s', pattern, met_dir)
      end
      if numel(match) > 1
         error('icemodel:forcing:stageModisAlbedo:metArtifactAmbiguous', ...
            '%d met artifacts match %s under %s', ...
            numel(match), pattern, met_dir)
      end
      met_file = string(fullfile(match.folder, match.name));

      % The met artifact's saved location is the colocation contract; the
      % MODIS artifact copies it verbatim so colocation identity checks
      % (artifactIdentityMatches) compare equal points.
      location = metArtifactLocation(met_file);
      locations(s).met_file = met_file;
      locations(s).lat_wgs84 = location.lat_wgs84;
      locations(s).lon_wgs84 = location.lon_wgs84;
      locations(s).location = location;
   end
end

function location = metArtifactLocation(met_file)
   %METARTIFACTLOCATION Resolve the location identity saved in a met artifact.

   % Prefer the source-light top-level artifact_metadata record (the
   % gap-filled PROMICE builders save lat/lon/elev there) so the multi-year
   % met timetable never has to load just for a point identity; fall back to
   % the full location CustomProperties contract of the met timetable.
   names = string(who('-file', met_file));
   if ismember("artifact_metadata", names)
      loaded = load(met_file, 'artifact_metadata');
      location = locationFromMetadata(loaded.artifact_metadata);
      if ~isempty(location)
         return
      end
   end
   if ismember("met", names)
      loaded = load(met_file, 'met');
      if istimetable(loaded.met)
         location = locationFromCustomProperties( ...
            loaded.met.Properties.CustomProperties);
         if ~isempty(location)
            return
         end
      end
   end
   error('icemodel:forcing:stageModisAlbedo:badMetLocation', ...
      'met artifact %s lacks a finite saved location identity', met_file)
end

function location = locationFromMetadata(metadata)
   %LOCATIONFROMMETADATA Location from a saved artifact metadata record.

   % Accept both canonical spellings: the builder convention (lat_wgs84/
   % lon_wgs84/elev_m) and the PROMICE met convention (lat/lon/elev). The
   % result follows the attachLocationMetadata/projectLocation schema;
   % without saved x/y, projectLocation later derives EPSG:3413 from lat/lon.
   location = struct.empty;
   if ~isstruct(metadata)
      return
   end
   lat = firstFiniteField(metadata, ["lat_wgs84", "lat"]);
   lon = firstFiniteField(metadata, ["lon_wgs84", "lon"]);
   if isnan(lat) || isnan(lon)
      return
   end
   location = struct('lat_wgs84', lat, 'lon_wgs84', lon, ...
      'elev_m', firstFiniteField(metadata, ["elev_m", "elev"]));
end

function location = locationFromCustomProperties(custom)
   %LOCATIONFROMCUSTOMPROPERTIES Location from the CustomProperties contract.

   % Field names follow the canonical attachLocationMetadata/projectLocation
   % location schema; an incomplete or nonfinite contract yields empty so the
   % caller can fail with one clear identity error.
   location = struct.empty;
   needed = ["X", "Y", "Lat", "Lon", "Elev", "Slope"];
   have = string(fieldnames(custom));
   if ~all(ismember(needed, have)) ...
         || ~(isscalar(custom.Lat) && isfinite(custom.Lat)) ...
         || ~(isscalar(custom.Lon) && isfinite(custom.Lon))
      return
   end
   location = struct( ...
      'x_epsg3413', custom.X, ...
      'y_epsg3413', custom.Y, ...
      'lat_wgs84', custom.Lat, ...
      'lon_wgs84', custom.Lon, ...
      'elev_m', custom.Elev, ...
      'slope', custom.Slope);
end

function value = firstFiniteField(metadata, fields)
   %FIRSTFINITEFIELD First finite scalar among aliased metadata fields.

   % NaN signals "not recorded"; the caller decides whether that is fatal
   % (lat/lon) or acceptable (elevation may be legitimately unknown).
   value = NaN;
   for field = fields
      if isfield(metadata, field)
         candidate = metadata.(char(field));
         if isnumeric(candidate) && isscalar(candidate) ...
               && isfinite(candidate)
            value = double(candidate);
            return
         end
      end
   end
end
