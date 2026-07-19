function [profiles, status, dynamic_qa] = readMarDensitySnapshots( ...
      requested_dates, grid_index, kwargs)
   %READMARDENSITYSNAPSHOTS Read requested MAR RO1 density profiles.
   %
   %  [profiles, status, dynamic_qa] = ...
   %     icemodel.forcing.helpers.readMarDensitySnapshots( ...
   %     requested_dates, [grid_i grid_j], source_dir=...)
   %
   % Reads only the requested UTC calendar dates from MAR yearly RUH2 files.
   % The public output is the source-provided RO1 fixed-depth density product
   % on exact OUTLAY levels and bounds. It is grouped by profile_id and source
   % datetime, uses positive-down metre depths and kg m-3 density, and carries
   % explicit model-output, grid-sampling, version, and file provenance.
   %
   % This dedicated reader never passes profile arrays through readMar3p11's
   % generic level-by-time flattening path. It never substitutes a nearby date
   % or extrapolates beyond source coverage. Reduced sources, non-ice cells,
   % and unavailable dates return explicit status rows without a profile.
   % Native TIME metadata is checked against exact YYYY/MM/DD/HH/MIN channels;
   % the lossy float32 packed DATE channel is deliberately ignored.
   %
   % DZSN1/ROSN1 are read only for internal QA through marDynamicProfileQa.
   % They are not exposed as a public profile and never modify RO1.
   %
   % Inputs
   %   requested_dates : datetimes; time-of-day is mapped to the same UTC date.
   %   grid_index      : native MAR [row column] index matching LON/LAT arrays.
   %
   % Required name-value
   %   source_dir : directory holding unique MAR *-YYYY.nc files.
   %
   % Optional name-value
   %   site_id            : stable site token used in profile_id.
   %   sample_method      : grid-sampling provenance; currently "nearest".
   %   requested_location : requested [lat lon] WGS84 coordinates.

   arguments
      requested_dates datetime
      grid_index (1, 2) double {mustBeFinite, mustBeInteger, mustBePositive}
      kwargs.source_dir (1, 1) string
      kwargs.site_id (1, 1) string = ""
      kwargs.sample_method (1, 1) string = "nearest"
      kwargs.requested_location (1, 2) double = [NaN NaN]
   end

   % Missing request timestamps cannot identify an exact source snapshot.
   if any(isnat(requested_dates))
      error('icemodel:forcing:readMarDensitySnapshots:missingRequestedDate', ...
         'requested_dates must not contain NaT values.')
   end

   % Normalize requests to sorted unique UTC calendar dates because each MAR
   % daily snapshot has one public profile identity.
   requested_dates = ensureUtc(requested_dates(:));
   requested_days = unique(dateshift(requested_dates, 'start', 'day'));
   profiles = emptyProfiles();
   status = emptyStatus(requested_days);
   dynamic_qa = emptyDynamicQa();
   if isempty(requested_days)
      profiles = stampProfileMetadata(profiles);
      return
   end

   % A blank provenance method would make the sampled grid identity ambiguous.
   if strlength(strtrim(kwargs.sample_method)) == 0
      error('icemodel:forcing:readMarDensitySnapshots:blankSampleMethod', ...
         'sample_method must be nonblank.')
   end
   if kwargs.sample_method ~= "nearest"
      error('icemodel:forcing:readMarDensitySnapshots:unsupportedSampleMethod', ...
         ['A native [row column] index proves only nearest-cell sampling; ' ...
         'natural/polygon profiles require the forcing colocation collapse.'])
   end

   % Missing source roots are an optional-profile status, not a forcing error.
   source_dir = kwargs.source_dir;
   if ~isfolder(source_dir)
      status.status(:) = "source_missing";
      status.detail(:) = "MAR source directory is unavailable";
      profiles = stampProfileMetadata(profiles);
      return
   end

   % Accumulate compact per-date tables so repeated depths remain separated by
   % profile_id/date rather than being collapsed across the requested archive.
   profile_parts = cell(numel(requested_days), 1);
   qa_parts = cell(numel(requested_days), 1);
   profile_count = 0;
   qa_count = 0;
   request_years = year(requested_days);
   years = unique(request_years);
   for y = reshape(years, 1, [])
      % One yearly file is required for an unambiguous exact-date selection.
      request_rows = find(request_years == y);
      [filename, file_status, file_detail] = locateYearFile(source_dir, y);
      if file_status ~= "ok"
         status.status(request_rows) = file_status;
         status.detail(request_rows) = file_detail;
         continue
      end
      status.source_file(request_rows) = filename;

      % Missing profile variables are accepted for reduced MAR sources; other
      % malformed schema details fail closed as invalid_source statuses.
      try
         schema = inspectSource(filename);
      catch exception
         status.status(request_rows) = "invalid_source";
         status.detail(request_rows) = string(exception.identifier) + ": " ...
            + string(exception.message);
         continue
      end
      if ~schema.profile_available
         status.status(request_rows) = "profile_not_available";
         status.detail(request_rows) = schema.detail;
         continue
      end

      % A grid mismatch is explicit for every requested date in this file; the
      % reader never clips indices to the nearest valid cell.
      if any(grid_index > schema.grid_shape)
         status.status(request_rows) = "grid_out_of_range";
         status.detail(request_rows) = compose( ...
            "requested [%d %d], source grid [%d %d]", grid_index(1), ...
            grid_index(2), schema.grid_shape(1), schema.grid_shape(2));
         continue
      end

      % Resolve the sampled cell once because static surface type and location
      % apply to every daily profile in the same yearly file.
      try
         surface_type = readPointVariable(filename, schema.srf_name, ...
            schema, grid_index, 1, "");
         sampled_lat = readPointVariable(filename, schema.lat_name, ...
            schema, grid_index, 1, "");
         sampled_lon = readPointVariable(filename, schema.lon_name, ...
            schema, grid_index, 1, "");
      catch exception
         status.status(request_rows) = "invalid_source";
         status.detail(request_rows) = string(exception.identifier) + ": " ...
            + string(exception.message);
         continue
      end
      if ~isscalar(surface_type) || ~isfinite(surface_type) ...
            || ~isscalar(sampled_lat) || ~isfinite(sampled_lat) ...
            || ~isscalar(sampled_lon) || ~isfinite(sampled_lon)
         status.status(request_rows) = "invalid_source";
         status.detail(request_rows) = ...
            "SRF/LAT/LON point values are missing or nonscalar";
         continue
      end

      for k = reshape(request_rows, 1, [])
         % Match the UTC calendar date exactly; absence never selects a nearest
         % snapshot and duplicate native dates are rejected as ambiguous.
         match = find(schema.source_days == requested_days(k));
         if isempty(match)
            status.status(k) = "out_of_window";
            status.detail(k) = "requested UTC date is absent from source file";
            continue
         elseif numel(match) ~= 1
            status.status(k) = "invalid_source_date";
            status.detail(k) = "source file contains duplicate UTC dates";
            continue
         end
         status.source_datetime(k) = schema.source_datetimes(match);

         % MAR SRF=4 is permanent ice. Non-ice columns are intentionally omitted
         % because their all-zero profile arrays are not physical density data.
         if surface_type ~= 4
            status.status(k) = "non_ice";
            status.detail(k) = "MAR SRF is not permanent ice (SRF=4)";
            continue
         end

         % Read only one time/grid hyperslab while retaining the OUTLAY axis.
         try
            density = readPointVariable(filename, schema.ro1_name, schema, ...
               grid_index, match, schema.outlay_dim);
         catch exception
            status.status(k) = "invalid_source";
            status.detail(k) = string(exception.identifier) + ": " ...
               + string(exception.message);
            continue
         end
         [valid, detail] = validateDensity(density, schema.outlay);
         if ~valid
            status.status(k) = "invalid_profile";
            status.detail(k) = detail;
            continue
         end

         % Preserve one explicit identity per site/date and repeat provenance on
         % every level so a standalone table remains self-describing.
         id = makeProfileId(kwargs.site_id, grid_index, ...
            schema.source_datetimes(match));
         profile_count = profile_count + 1;
         profile_parts{profile_count} = makeProfileRows(id, ...
            schema.source_datetimes(match), density, schema, filename, ...
            grid_index, kwargs, sampled_lat, sampled_lon);
         status.status(k) = "selected";
         status.detail(k) = "exact-date RO1 snapshot selected";
         status.profile_id(k) = id;
         status.selected_rows(k) = numel(density);

         % Dynamic-layer diagnostics are optional and cannot invalidate a valid
         % public RO1 profile. Read/schema failures remain visible in the nested
         % diagnostic record rather than changing selection status.
         diagnostic = readDynamicDiagnostic(filename, schema, grid_index, ...
            match, density);
         qa_count = qa_count + 1;
         qa_parts{qa_count} = struct( ...
            'profile_id', id, ...
            'datetime', schema.source_datetimes(match), ...
            'source_file', filename, ...
            'diagnostic', diagnostic);
      end
   end

   % Concatenate only selected compact snapshots and stamp one stable public
   % profile contract after table concatenation preserves row metadata.
   if profile_count > 0
      profiles = vertcat(profile_parts{1:profile_count});
   end
   profiles = stampProfileMetadata(profiles);
   if qa_count > 0
      dynamic_qa = vertcat(qa_parts{1:qa_count});
   end
end

function dates = ensureUtc(dates)
   %ENSUREUTC Interpret timezone-free requests as UTC, then convert zoned input.
   dates.TimeZone = 'UTC';
end

function [filename, state, detail] = locateYearFile(source_dir, yyyy)
   %LOCATEYEARFILE Require one exact yearly MAR source file.
   matches = dir(fullfile(source_dir, sprintf('*-%d.nc', yyyy)));
   if isempty(matches)
      filename = "";
      state = "source_missing";
      detail = compose("no MAR *-%d.nc source file", yyyy);
   elseif numel(matches) > 1
      filename = "";
      state = "ambiguous_source";
      detail = compose("multiple MAR *-%d.nc source files", yyyy);
   else
      filename = string(fullfile(matches(1).folder, matches(1).name));
      state = "ok";
      detail = "";
   end
end

function schema = inspectSource(filename)
   %INSPECTSOURCE Resolve exact variable/dimension names and compact axes.
   info = ncinfo(filename);
   schema = struct('profile_available', false, 'detail', "");
   schema.time_name = variableName(info, "TIME");
   schema.time_component_names = [ ...
      variableName(info, "YYYY"), variableName(info, "MM"), ...
      variableName(info, "DD"), variableName(info, "HH"), ...
      variableName(info, "MIN")];
   schema.ro1_name = variableName(info, "RO1");
   schema.outlay_name = variableName(info, "OUTLAY");
   schema.lon_name = variableName(info, "LON");
   schema.lat_name = variableName(info, "LAT");
   schema.srf_name = variableName(info, "SRF");
   profile_product = [schema.ro1_name, schema.outlay_name];
   if all(profile_product == "")
      % A genuinely reduced source omits the complete optional profile product.
      schema.detail = "missing optional profile variables: RO1, OUTLAY";
      return
   elseif any(profile_product == "")
      error('icemodel:forcing:readMarDensitySnapshots:invalidPublicSchema', ...
         'RO1 and OUTLAY must be both present or both absent.')
   end
   grid_variables = [schema.lon_name, schema.lat_name, schema.srf_name];
   if any(grid_variables == "")
      labels = ["LON", "LAT", "SRF"];
      error('icemodel:forcing:readMarDensitySnapshots:invalidPublicSchema', ...
         'A RO1/OUTLAY source is missing grid variable(s): %s.', ...
         strjoin(labels(grid_variables == ""), ', '))
   end
   has_time_components = schema.time_component_names ~= "";
   if schema.time_name == ""
      error('icemodel:forcing:readMarDensitySnapshots:invalidTimeAxis', ...
         'A RO1/OUTLAY source is missing the authoritative TIME coordinate.')
   elseif any(has_time_components) && ~all(has_time_components)
      error('icemodel:forcing:readMarDensitySnapshots:partialTimeComponents', ...
         'YYYY/MM/DD/HH/MIN must be all present or all absent.')
   end

   % LON defines the native two-dimensional point-array orientation used by
   % caller grid indices and every point hyperslab. Require the exact native
   % X/Y grid set on all three static fields before any point value is read.
   lon_info = ncinfo(filename, schema.lon_name);
   lat_info = ncinfo(filename, schema.lat_name);
   srf_info = ncinfo(filename, schema.srf_name);
   schema.grid_dims = string({lon_info.Dimensions.Name});
   grid_keys = upper(schema.grid_dims);
   if numel(schema.grid_dims) ~= 2 ...
         || numel(unique(grid_keys)) ~= 2 ...
         || sum(startsWith(grid_keys, "X")) ~= 1 ...
         || sum(startsWith(grid_keys, "Y")) ~= 1 ...
         || ~sameDimensions(lon_info.Dimensions, lat_info.Dimensions) ...
         || ~sameDimensions(lon_info.Dimensions, srf_info.Dimensions)
      error('icemodel:forcing:readMarDensitySnapshots:invalidGridAxes', ...
         'LON, LAT, and SRF must share exactly one native X and one Y axis.')
   end
   schema.grid_shape = double([lon_info.Dimensions.Length]);

   % Validate the native TIME axis metadata before any coordinate/profile read.
   time_info = ncinfo(filename, schema.time_name);
   requireExactDimensions(time_info, "TIME", ...
      'icemodel:forcing:readMarDensitySnapshots:invalidTimeAxis')
   schema.time_dim = "TIME";

   % OUTLAY and its declared bounds are preserved exactly and validated only
   % for shape, positive-down ordering, and required 0-20 m support.
   outlay_info = ncinfo(filename, schema.outlay_name);
   requireExactDimensions(outlay_info, "OUTLAY", ...
      'icemodel:forcing:readMarDensitySnapshots:invalidOutlayAxis')
   requireNativeUnits(outlay_info, "m", ...
      'icemodel:forcing:readMarDensitySnapshots:invalidOutlayUnits')
   schema.outlay_dim = "OUTLAY";
   bounds_name = outlayBoundsName(info, outlay_info);
   if bounds_name == ""
      error('icemodel:forcing:readMarDensitySnapshots:invalidOutlayBounds', ...
         'A RO1/OUTLAY source is missing its declared OUTLAY bounds variable.')
   end
   bounds_info = ncinfo(filename, bounds_name);
   validateOutlayBoundsDimensions(bounds_info)

   % RO1 must contain exactly time, OUTLAY, and both native grid dimensions;
   % an unexpected sector/level axis would change the public product meaning.
   ro1_info = ncinfo(filename, schema.ro1_name);
   expected = [schema.grid_dims, schema.outlay_dim, schema.time_dim];
   requireExactDimensions(ro1_info, expected, ...
      'icemodel:forcing:readMarDensitySnapshots:invalidRo1Shape')
   requireNativeUnits(ro1_info, ["kg/m3", "kg m-3"], ...
      'icemodel:forcing:readMarDensitySnapshots:invalidRo1Units')

   % Only after every public axis/unit guard passes, decode source time and read
   % the compact fixed-depth coordinate/bounds used to stamp selected rows.
   [schema.source_datetimes, schema.time_dim, ...
      schema.source_time_variables] = sourceDatetimes(filename, schema);
   schema.source_days = dateshift(schema.source_datetimes, 'start', 'day');
   schema.outlay = double(ncread(filename, schema.outlay_name));
   schema.outlay = schema.outlay(:);
   schema.outlay_bounds = orientBounds( ...
      double(ncread(filename, bounds_name)), numel(schema.outlay));
   validateOutlay(schema.outlay, schema.outlay_bounds)

   % Optional dynamic variables are resolved without affecting public RO1
   % availability; their exact shapes are checked only when diagnostics run.
   schema.dzsn1_name = variableName(info, "DZSN1");
   schema.rosn1_name = variableName(info, "ROSN1");
   schema.shsn3_name = variableName(info, "SHSN3");
   schema.dynamic_available = all([schema.dzsn1_name, ...
      schema.rosn1_name, schema.shsn3_name] ~= "");
   schema.mar_version = marVersion(filename);
   schema.profile_available = true;
   schema.detail = "";
end

function diagnostic = readDynamicDiagnostic(filename, schema, grid_index, ...
      time_index, ro1)
   %READDYNAMICDIAGNOSTIC Run optional dynamic-layer QA without blocking RO1.
   if ~schema.dynamic_available
      diagnostic = icemodel.forcing.helpers.marDynamicProfileQa( ...
         [], [], [], ro1, schema.outlay);
      return
   end

   % Dynamic QA accepts only the documented native axes and units. This blocks
   % unknown dimensions from silently defaulting to index one in the point read.
   try
      dz_info = ncinfo(filename, schema.dzsn1_name);
      rho_info = ncinfo(filename, schema.rosn1_name);
      shsn_info = ncinfo(filename, schema.shsn3_name);
      layer_dim = "SNOLAY";
      requireExactDimensions(dz_info, ...
         [schema.grid_dims, layer_dim, schema.time_dim], ...
         'icemodel:forcing:readMarDensitySnapshots:invalidDynamicShape')
      requireExactDimensions(rho_info, ...
         [schema.grid_dims, layer_dim, schema.time_dim], ...
         'icemodel:forcing:readMarDensitySnapshots:invalidDynamicShape')
      requireExactDimensions(shsn_info, ...
         [schema.grid_dims, "SECTOR", schema.time_dim], ...
         'icemodel:forcing:readMarDensitySnapshots:invalidDynamicShape')
      dz_names = string({dz_info.Dimensions.Name});
      rho_names = string({rho_info.Dimensions.Name});
      shsn_names = string({shsn_info.Dimensions.Name});
      dz_layer_count = double(dz_info.Dimensions( ...
         strcmpi(dz_names, layer_dim)).Length);
      rho_layer_count = double(rho_info.Dimensions( ...
         strcmpi(rho_names, layer_dim)).Length);
      sector_count = double(shsn_info.Dimensions( ...
         strcmpi(shsn_names, "SECTOR")).Length);
      if dz_layer_count < 1 || rho_layer_count ~= dz_layer_count
         error('icemodel:forcing:readMarDensitySnapshots:invalidDynamicShape', ...
            'DZSN1 and ROSN1 must share a nonempty SNOLAY axis.')
      elseif sector_count < 1
         error('icemodel:forcing:readMarDensitySnapshots:invalidDynamicShape', ...
            'SHSN3 SECTOR must include permanent-ice index 1.')
      end
      requireNativeUnits(dz_info, "m", ...
         'icemodel:forcing:readMarDensitySnapshots:invalidDynamicUnits')
      requireNativeUnits(rho_info, ["kg/m3", "kg m-3"], ...
         'icemodel:forcing:readMarDensitySnapshots:invalidDynamicUnits')
      requireNativeUnits(shsn_info, "m", ...
         'icemodel:forcing:readMarDensitySnapshots:invalidDynamicUnits')
      dz = readPointVariable(filename, schema.dzsn1_name, schema, ...
         grid_index, time_index, layer_dim);
      rho = readPointVariable(filename, schema.rosn1_name, schema, ...
         grid_index, time_index, layer_dim);
      shsn3 = readPointVariable(filename, schema.shsn3_name, schema, ...
         grid_index, time_index, "", "SECTOR");
      diagnostic = icemodel.forcing.helpers.marDynamicProfileQa( ...
         dz, rho, shsn3, ro1, schema.outlay);
   catch exception
      diagnostic = icemodel.forcing.helpers.marDynamicProfileQa( ...
         [], [], [], ro1, schema.outlay);
      diagnostic.status = "source_read_error";
      diagnostic.detail = string(exception.identifier) + ": " ...
         + string(exception.message);
   end
end

function values = readPointVariable(filename, variable, schema, grid_index, ...
      time_index, retained_dim, selected_dim)
   %READPOINTVARIABLE Read one native grid/time hyperslab by dimension name.
   if nargin < 7
      % Most point variables have no additional selected scalar dimension.
      selected_dim = "";
   end
   variable_info = ncinfo(filename, variable);
   names = string({variable_info.Dimensions.Name});
   sizes = double([variable_info.Dimensions.Length]);
   start = ones(1, numel(names));
   count = ones(1, numel(names));
   for d = 1:numel(names)
      % Dimension-name routing avoids assumptions about NetCDF-versus-MATLAB
      % array order and keeps the native grid index contract explicit.
      if strcmpi(names(d), schema.grid_dims(1))
         start(d) = grid_index(1);
      elseif strcmpi(names(d), schema.grid_dims(2))
         start(d) = grid_index(2);
      elseif strcmpi(names(d), schema.time_dim)
         start(d) = time_index;
      elseif retained_dim ~= "" && strcmpi(names(d), retained_dim)
         count(d) = sizes(d);
      elseif selected_dim ~= "" && strcmpi(names(d), selected_dim)
         % SHSN3 sector 1 is the documented permanent-ice diagnostic.
         start(d) = 1;
      end
   end
   values = double(ncread(filename, variable, start, count));
   values = values(:);
end

function [valid, detail] = validateDensity(density, outlay)
   %VALIDATEDENSITY Enforce conservative physical and shape constraints.
   density = double(density(:));
   if numel(density) ~= numel(outlay)
      valid = false;
      detail = "RO1 level count does not match OUTLAY";
   elseif any(~isfinite(density)) || all(density == 0)
      valid = false;
      detail = "ice-cell RO1 is missing or all zero";
   elseif any(density < 250 | density > 1000)
      valid = false;
      detail = "RO1 density lies outside 250-1000 kg m-3";
   else
      valid = true;
      detail = "";
   end
end

function rows = makeProfileRows(id, source_datetime, density, schema, ...
      filename, grid_index, kwargs, sampled_lat, sampled_lon)
   %MAKEPROFILEROWS Build one self-describing exact-date RO1 profile table.
   n = numel(schema.outlay);
   rows = table( ...
      repmat(id, n, 1), ...
      repmat(source_datetime, n, 1), ...
      schema.outlay, ...
      schema.outlay_bounds(:, 1), ...
      schema.outlay_bounds(:, 2), ...
      double(density(:)), ...
      repmat("model_output", n, 1), ...
      repmat("mar3.11", n, 1), ...
      repmat("RO1", n, 1), ...
      repmat(filename, n, 1), ...
      repmat(schema.mar_version, n, 1), ...
      repmat(schema.source_time_variables, n, 1), ...
      repmat("fixed-depth RO1; no sector dimension", n, 1), ...
      repmat("down", n, 1), ...
      repmat(kwargs.sample_method, n, 1), ...
      repmat(grid_index(1), n, 1), ...
      repmat(grid_index(2), n, 1), ...
      repmat(kwargs.requested_location(1), n, 1), ...
      repmat(kwargs.requested_location(2), n, 1), ...
      repmat(sampled_lat, n, 1), ...
      repmat(sampled_lon, n, 1), ...
      'VariableNames', profileVariableNames());
end

function profiles = emptyProfiles()
   %EMPTYPROFILES Return a typed zero-row public profile table.
   string_column = strings(0, 1);
   numeric_column = zeros(0, 1);
   profile_time = NaT(0, 1, 'TimeZone', 'UTC');
   profiles = table(string_column, profile_time, numeric_column, ...
      numeric_column, numeric_column, numeric_column, string_column, ...
      string_column, string_column, string_column, string_column, ...
      string_column, string_column, string_column, string_column, ...
      numeric_column, ...
      numeric_column, numeric_column, numeric_column, numeric_column, ...
      numeric_column, 'VariableNames', profileVariableNames());
end

function names = profileVariableNames()
   %PROFILEVARIABLENAMES Centralize the stable public table column order.
   names = {'profile_id', 'datetime', 'depth', 'depth_lower_bound', ...
      'depth_upper_bound', 'density', 'product_role', 'source_id', ...
      'source_variable', 'source_file', 'mar_version', ...
      'source_time_variables', 'sector_semantics', 'depth_positive', ...
      'sample_method', 'grid_i', 'grid_j', ...
      'requested_lat_wgs84', 'requested_lon_wgs84', ...
      'sampled_lat_wgs84', 'sampled_lon_wgs84'};
end

function profiles = stampProfileMetadata(profiles)
   %STAMPPROFILEMETADATA Attach units, descriptions, and bundle provenance.
   names = string(profiles.Properties.VariableNames);
   units = repmat({''}, 1, width(profiles));
   depth_units = ismember(names, ["depth", "depth_lower_bound", ...
      "depth_upper_bound"]);
   units(depth_units) = repmat({'m'}, 1, nnz(depth_units));
   units(names == "density") = {'kg m-3'};
   latitude_units = ismember(names, ["requested_lat_wgs84", ...
      "sampled_lat_wgs84"]);
   longitude_units = ismember(names, ["requested_lon_wgs84", ...
      "sampled_lon_wgs84"]);
   units(latitude_units) = repmat( ...
      {'degrees_north'}, 1, nnz(latitude_units));
   units(longitude_units) = repmat( ...
      {'degrees_east'}, 1, nnz(longitude_units));
   profiles.Properties.VariableUnits = units;

   % Descriptions distinguish public source values from QA-only reconstruction.
   descriptions = repmat({''}, 1, width(profiles));
   descriptions{names == "density"} = ...
      'MAR modelled snow/firn/ice density from source RO1';
   descriptions{names == "depth"} = ...
      'source OUTLAY fixed depth, positive down';
   descriptions{names == "depth_lower_bound"} = ...
      'source OUTLAY lower bound';
   descriptions{names == "depth_upper_bound"} = ...
      'source OUTLAY upper bound';
   profiles.Properties.VariableDescriptions = descriptions;
   profiles.Properties.UserData = struct( ...
      'format', 'subsurface_profile_bundle', ...
      'product_name', 'MAR modelled snow/firn/ice density', ...
      'public_variable', 'RO1', ...
      'depth_coordinate', 'OUTLAY with source bounds', ...
      'profile_grouping', 'profile_id + datetime', ...
      'dynamic_layer_role', 'DZSN1/ROSN1 diagnostic only', ...
      'correction_applied', false, ...
      'time_selection', 'exact UTC calendar date; no extrapolation', ...
      'ignored_time_variable', ...
      'DATE float32 YYYYMMDDHH is quantized and is never decoded');
end

function status = emptyStatus(requested_days)
   %EMPTYSTATUS Return one explicit selection-status row per requested date.
   n = numel(requested_days);
   status = table(requested_days, strings(n, 1), strings(n, 1), ...
      strings(n, 1), NaT(n, 1, 'TimeZone', 'UTC'), strings(n, 1), ...
      zeros(n, 1), 'VariableNames', {'requested_date', 'status', 'detail', ...
      'source_file', 'source_datetime', 'profile_id', 'selected_rows'});
end

function dynamic_qa = emptyDynamicQa()
   %EMPTYDYNAMICQA Return a stable empty wrapper for per-profile diagnostics.
   dynamic_qa = struct('profile_id', {}, 'datetime', {}, ...
      'source_file', {}, 'diagnostic', {});
end

function name = variableName(info, requested)
   %VARIABLENAME Resolve a source variable case-insensitively.
   names = string({info.Variables.Name});
   match = find(strcmpi(names, requested), 1);
   if isempty(match)
      name = "";
   else
      name = names(match);
   end
end

function name = outlayBoundsName(info, outlay_info)
   %OUTLAYBOUNDSNAME Prefer the source-declared bounds attribute.
   name = "";
   attributes = string({outlay_info.Attributes.Name});
   match = find(strcmpi(attributes, "bounds"), 1);
   if ~isempty(match)
      candidate = string(outlay_info.Attributes(match).Value);
      name = variableName(info, candidate);
   end
   if name == ""
      name = variableName(info, "OUTLAY_bnds");
   end
end

function bounds = orientBounds(bounds, level_count)
   %ORIENTBOUNDS Return source bounds as OUTLAY-by-two without changing values.
   bounds = squeeze(double(bounds));
   if isequal(size(bounds), [2, level_count])
      bounds = bounds';
   elseif ~isequal(size(bounds), [level_count, 2])
      error('icemodel:forcing:readMarDensitySnapshots:invalidOutlayBounds', ...
         'OUTLAY bounds must be 2-by-level or level-by-2.')
   end
end

function validateOutlay(outlay, bounds)
   %VALIDATEOUTLAY Require ordered source levels with complete 0-20 m support.
   endpoint_tolerance = 100 * eps(max(20, max(abs(outlay))));
   valid = all(isfinite(outlay)) && all(diff(outlay) > 0) ...
      && abs(outlay(1)) <= endpoint_tolerance ...
      && abs(outlay(end) - 20) <= endpoint_tolerance ...
      && all(isfinite(bounds), 'all') ...
      && all(bounds(:, 1) < bounds(:, 2)) ...
      && all(bounds(:, 1) <= outlay & outlay <= bounds(:, 2));
   if ~valid
      error('icemodel:forcing:readMarDensitySnapshots:invalidOutlay', ...
         'OUTLAY/bounds must be monotonic and span exact 0-20 m support.')
   end
end

function validateOutlayBoundsDimensions(variable_info)
   %VALIDATEOUTLAYBOUNDSDIMENSIONS Require OUTLAY plus one two-value bound axis.
   names = string({variable_info.Dimensions.Name});
   sizes = double([variable_info.Dimensions.Length]);
   outlay_index = find(strcmpi(names, "OUTLAY"));
   bound_index = find(~strcmpi(names, "OUTLAY"));
   valid = numel(names) == 2 && numel(unique(lower(names))) == 2 ...
      && isscalar(outlay_index) && isscalar(bound_index) ...
      && sizes(bound_index) == 2;
   if ~valid
      error('icemodel:forcing:readMarDensitySnapshots:invalidOutlayBounds', ...
         'OUTLAY bounds must contain exactly OUTLAY and one length-two axis.')
   end
end

function [dates, time_dim, provenance] = sourceDatetimes(filename, schema)
   %SOURCEDATETIMES Decode native TIME and verify exact component timestamps.
   has_components = all(schema.time_component_names ~= "");
   % TIME values plus units retain sub-day placement without integer packing.
   time_info = ncinfo(filename, schema.time_name);
   requireExactDimensions(time_info, "TIME", ...
      'icemodel:forcing:readMarDensitySnapshots:invalidTimeAxis')
   time_dim = "TIME";
   units = attributeValue(time_info, "units");
   dates = decodeTimeValues(ncread(filename, schema.time_name), units);
   provenance = schema.time_name + " (" + units + ")";

   % Full RUH2 files carry both representations. Require exact agreement so a
   % stale component vector cannot silently relabel the native coordinate.
   if has_components
      component_dates = decodeTimeComponents(filename, ...
         schema.time_component_names, time_dim);
      if numel(component_dates) ~= numel(dates) ...
            || any(abs(seconds(component_dates - dates)) > 0.5)
         error('icemodel:forcing:readMarDensitySnapshots:timeMismatch', ...
            'TIME and YYYY/MM/DD/HH/MIN timestamps disagree.')
      end
      provenance = provenance + "; checked against " ...
         + strjoin(schema.time_component_names, "/");
   end

   % Native timestamps must be finite. Duplicate calendar dates remain visible
   % to the per-request ambiguity status rather than being silently deduplicated.
   if any(isnat(dates))
      error('icemodel:forcing:readMarDensitySnapshots:invalidTimeValues', ...
         'Native MAR timestamps contain NaT values.')
   end
   dates = dates(:);
end

function dates = decodeTimeValues(values, units)
   %DECODETIMEVALUES Decode a CF-style elapsed-time coordinate in UTC.
   match = regexp(units, ...
      '^(seconds|minutes|hours|days) since (.+)$', 'tokens', 'once');
   if isempty(match)
      error('icemodel:forcing:readMarDensitySnapshots:invalidTimeUnits', ...
         'Unsupported MAR TIME units: %s', units)
   end
   origin = decodeTimeOrigin(string(match{2}));
   values = double(values(:));
   if any(~isfinite(values))
      error('icemodel:forcing:readMarDensitySnapshots:invalidTimeValues', ...
         'TIME contains nonfinite values.')
   end

   % Apply the declared elapsed unit exactly; DATE is intentionally untouched.
   switch string(match{1})
      case "seconds"
         dates = origin + seconds(values);
      case "minutes"
         dates = origin + minutes(values);
      case "hours"
         dates = origin + hours(values);
      case "days"
         dates = origin + days(values);
   end
end

function origin = decodeTimeOrigin(text)
   %DECODETIMEORIGIN Parse the stable date/time forms used by MAR TIME units.
   formats = ["yyyy-MM-dd HH:mm:ss", "yyyy-MM-dd HH:mm", "yyyy-MM-dd"];
   origin = NaT(1, 1, 'TimeZone', 'UTC');
   for format = formats
      % Try explicit formats so locale settings cannot reinterpret the origin.
      try
         origin = datetime(text, 'InputFormat', format, 'TimeZone', 'UTC');
         break
      catch
         % The next supported format may match a shorter declared origin.
      end
   end
   if isnat(origin)
      error('icemodel:forcing:readMarDensitySnapshots:invalidTimeUnits', ...
         'Cannot parse MAR TIME origin: %s', text)
   end
end

function dates = decodeTimeComponents(filename, names, time_dim)
   %DECODETIMECOMPONENTS Build exact UTC timestamps from five native channels.
   values = cell(1, numel(names));
   for k = 1:numel(names)
      % Each exact component must be a one-dimensional integer time channel.
      info = ncinfo(filename, names(k));
      if numel(info.Dimensions) ~= 1 ...
            || ~strcmpi(string(info.Dimensions(1).Name), time_dim)
         error('icemodel:forcing:readMarDensitySnapshots:invalidTimeAxis', ...
            '%s does not share the native time dimension.', names(k))
      end
      values{k} = double(ncread(filename, names(k)));
      values{k} = values{k}(:);
      if any(~isfinite(values{k}) | values{k} ~= round(values{k}))
         error('icemodel:forcing:readMarDensitySnapshots:invalidTimeValues', ...
            '%s contains noninteger or nonfinite values.', names(k))
      end
   end
   lengths = cellfun(@numel, values);
   if any(lengths ~= lengths(1))
      error('icemodel:forcing:readMarDensitySnapshots:invalidTimeAxis', ...
         'YYYY/MM/DD/HH/MIN lengths differ.')
   end
   try
      dates = datetime(values{1}, values{2}, values{3}, values{4}, ...
         values{5}, zeros(size(values{1})), 'TimeZone', 'UTC');
   catch exception
      error('icemodel:forcing:readMarDensitySnapshots:invalidTimeValues', ...
         'Cannot decode YYYY/MM/DD/HH/MIN: %s', exception.message)
   end
end

function value = attributeValue(variable_info, requested)
   %ATTRIBUTEVALUE Read one required NetCDF variable attribute as text.
   if isempty(variable_info.Attributes)
      % Attribute-free NetCDF variables are numeric empty in ncinfo output.
      match = [];
   else
      names = string({variable_info.Attributes.Name});
      match = find(strcmpi(names, requested), 1);
   end
   if isempty(match)
      error('icemodel:forcing:readMarDensitySnapshots:missingTimeUnits', ...
         '%s has no %s attribute.', variable_info.Name, requested)
   end
   value = string(variable_info.Attributes(match).Value);
end

function requireExactDimensions(variable_info, expected, identifier)
   %REQUIREEXACTDIMENSIONS Reject missing, duplicated, unknown, or extra axes.
   actual = string({variable_info.Dimensions.Name});
   expected = string(expected);
   valid = numel(actual) == numel(expected) ...
      && numel(unique(lower(actual))) == numel(actual) ...
      && all(ismember(lower(expected), lower(actual)));
   if ~valid
      error(identifier, '%s must have exactly dimensions {%s}; found {%s}.', ...
         variable_info.Name, strjoin(expected, ', '), strjoin(actual, ', '))
   end
end

function requireNativeUnits(variable_info, accepted, identifier)
   %REQUIRENATIVEUNITS Validate the finite set of source-faithful MAR spellings.
   if isempty(variable_info.Attributes)
      % ncinfo represents an attribute-free variable as numeric empty, so do
      % not dot-index it before raising the source-specific unit diagnostic.
      match = [];
   else
      names = string({variable_info.Attributes.Name});
      match = find(strcmpi(names, "units"), 1);
   end
   if isempty(match)
      error(identifier, '%s has no native units attribute.', variable_info.Name)
   end
   units = strtrim(string(variable_info.Attributes(match).Value));
   if ~any(units == string(accepted))
      error(identifier, '%s native units are %s; expected one of {%s}.', ...
         variable_info.Name, units, strjoin(string(accepted), ', '))
   end
end

function tf = sameDimensions(left, right)
   %SAMEDIMENSIONS Compare NetCDF dimension-name/length sets without reordering.
   left_names = string({left.Name});
   right_names = string({right.Name});
   left_sizes = double([left.Length]);
   right_sizes = double([right.Length]);
   tf = numel(left_names) == numel(right_names) ...
      && numel(unique(lower(left_names))) == numel(left_names) ...
      && numel(unique(lower(right_names))) == numel(right_names);
   for k = 1:numel(left_names)
      % Match lengths by name because NetCDF/MATLAB axis order is not semantic.
      match = find(strcmpi(right_names, left_names(k)));
      tf = tf && isscalar(match) && left_sizes(k) == right_sizes(match);
   end
end

function version = marVersion(filename)
   %MARVERSION Derive the source-version token from the yearly filename.
   [~, basename] = fileparts(filename);
   token = regexp(basename, 'MARv\d+(?:\.\d+)+', 'match', 'once');
   if isempty(token)
      version = "unknown";
   else
      version = string(token);
   end
end

function id = makeProfileId(site_id, grid_index, source_datetime)
   %MAKEPROFILEID Build a stable site/grid plus UTC-date public profile key.
   if strlength(strtrim(site_id)) == 0
      prefix = compose("mar_i%d_j%d", grid_index(1), grid_index(2));
   else
      prefix = regexprep(strtrim(site_id), '[^A-Za-z0-9_-]+', '_');
   end
   id = prefix + "_mar3p11_" + string(source_datetime, 'yyyyMMdd');
end
