function [inventory, status] = gcnetInventory(source_dir, kwargs)
   %GCNETINVENTORY Index Vandecrux/GC-Net products without loading arrays.
   %
   %  inventory = icemodel.verification.setup.gcnetInventory(source_dir)
   %  [inventory, status] = icemodel.verification.setup.gcnetInventory( ...
   %     source_dir, stations=["dye2", "sum"])
   %
   % Role
   %  Lightweight discovery layer for the Vandecrux/GC-Net Arctic Data cache.
   %  It records product files, NetCDF header variables, dimensions, coordinate
   %  time/depth summaries, canonical comparison variables, station metadata,
   %  and XML provenance while avoiding large array reads.

   arguments
      source_dir (1, 1) string = ""
      kwargs.stations (1, :) string = ["DYE_2", "Summit"]
      kwargs.products (1, :) string ...
         {icemodel.verification.validators.mustBeGcnetProductSelection( ...
         kwargs.products)} = ...
         icemodel.verification.setup.gcnetProductNames()
   end

   if source_dir == ""
      source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
         'verification', 'gcnet'));
   end

   % The fetch helper supplies status rows and DOI expectations; this inventory
   % remains tolerant so downstream tooling can see partial caches.
   [~, status] = icemodel.verification.setup.fetchGcnet( ...
      cache_dir=source_dir, stations=kwargs.stations, ...
      products=kwargs.products, strict=false, silent=true);

   stations = canonicalStations(kwargs.stations);
   records = repmat(emptyRecord(), 1, 0);
   n_records = 0;
   provenance = provenanceByProduct(status, kwargs.products);

   % Build one record per station/product/file class. The fetch status owns the
   % registry-backed classification, so inventory must consume its resolved
   % station/suffix/file bindings without matching paths a second time.
   for product = reshape(kwargs.products, 1, [])
      resolved = resolvedProductEntries(status, product);
      for entry = reshape(resolved, 1, [])
         n_records = n_records + 1;
         records(n_records) = fileRecord(entry.filename, product, ...
            entry.station, entry.suffix, provenance.(char(product)));
      end
   end

   inventory = struct( ...
      'source_dir', source_dir, ...
      'stations', stations, ...
      'products', kwargs.products, ...
      'records', records);
end

function entries = resolvedProductEntries(status, product)
   %RESOLVEDPRODUCTENTRIES Reuse fetch validation's exact identity decisions.
   rows = status(string({status.product}) == product);
   if isempty(rows)
      prototype = struct('station', "", 'suffix', "", 'filename', "");
      entries = repmat(prototype, 1, 0);
      return
   end

   % An incomplete product may still contain uniquely resolved classes that are
   % useful for tolerant inventory. Ambiguous and partial names never enter the
   % status row's resolved_files and therefore cannot reach header inspection.
   entries = [rows.resolved_files];
   [~, keep] = unique([entries.filename], 'stable');
   entries = entries(keep);
end

function stations = canonicalStations(stations)
   %CANONICALSTATIONS Map RetMIP and filename aliases to Vandecrux station ids.
   stations = icemodel.forcing.helpers.gcnetVandecruxStation(stations);
   stations = unique(stations, 'stable');
end

function row = emptyRecord()
   %EMPTYRECORD Prototype for one station/product NetCDF header record.
   row = struct( ...
      'station', "", ...
      'aliases', strings(1, 0), ...
      'product', "", ...
      'kind', "", ...
      'canonical_variables', strings(1, 0), ...
      'filename', "", ...
      'variables', strings(1, 0), ...
      'units', struct(), ...
      'dimensions', struct(), ...
      'depth', struct('variable', "", 'units', "", 'levels', NaN, ...
      'sample_min_m', NaN, 'sample_max_m', NaN, 'sample_policy', ""), ...
      'period', struct('start', NaT, 'end', NaT), ...
      'cadence', "", ...
      'time_units', "", ...
      'site_location', struct(), ...
      'provenance', struct());
end

function record = fileRecord(filename, product, station, suffix, provenance)
   %FILERECORD Read one NetCDF header and convert it to inventory metadata.
   info = ncinfo(filename);
   variables = string({info.Variables.Name});

   % Variable units are attributes, so reading them does not touch payload data.
   units = struct();
   for v = reshape(variables, 1, [])
      units.(matlab.lang.makeValidName(v)) = ...
         icemodel.forcing.helpers.readNetcdfAttribute(filename, v, 'units');
   end

   dimensions = struct();
   for k = 1:numel(info.Dimensions)
      name = matlab.lang.makeValidName(info.Dimensions(k).Name);
      dimensions.(name) = info.Dimensions(k).Length;
   end

   [period, cadence, time_units] = timeSummary(filename, variables);
   record = emptyRecord();
   record.station = station;
   station_info = icemodel.forcing.helpers.gcnetVandecruxStationMetadata( ...
      station);
   record.aliases = station_info.aliases;
   record.product = product;
   record.kind = productKind(product, suffix);
   record.canonical_variables = canonicalVariables(record.kind);
   record.filename = string(filename);
   record.variables = variables(:);
   record.units = units;
   record.dimensions = dimensions;
   record.depth = depthSummary(filename, variables);
   record.period = period;
   record.cadence = cadence;
   record.time_units = time_units;
   record.site_location = station_info.site_location;
   record.provenance = provenance;
end

function [period, cadence, units] = timeSummary(filename, variables)
   %TIMESUMMARY Read only first/last time coordinates and infer cadence.
   period = struct('start', NaT, 'end', NaT);
   cadence = "";
   units = "";
   if ~any(variables == "time")
      return
   end

   info = ncinfo(filename, 'time');
   units = icemodel.forcing.helpers.readNetcdfAttribute( ...
      filename, "time", "units");
   n_time = info.Size(1);
   if n_time == 0
      return
   end

   % Only the first, second, and last coordinates are needed for the inventory.
   first = double(ncread(filename, 'time', 1, 1));
   last = double(ncread(filename, 'time', n_time, 1));
   period.start = icemodel.forcing.helpers.gcnetTime(first, units);
   period.end = icemodel.forcing.helpers.gcnetTime(last, units);
   if n_time > 1
      second = double(ncread(filename, 'time', 2, 1));
      cadence = string(hours( ...
         icemodel.forcing.helpers.gcnetTime(second, units) ...
         - period.start)) + "hr";
   end
end

function summary = depthSummary(filename, variables)
   %DEPTHSUMMARY Read small depth samples without loading full state arrays.
   summary = struct('variable', "", 'units', "", 'levels', NaN, ...
      'sample_min_m', NaN, 'sample_max_m', NaN, 'sample_policy', "");
   if ~any(variables == "Depth")
      return
   end

   info = ncinfo(filename, 'Depth');
   dims = string({info.Dimensions.Name});
   sizes = double([info.Dimensions.Length]);
   summary.variable = "Depth";
   summary.units = icemodel.forcing.helpers.readNetcdfAttribute( ...
      filename, "Depth", "units");
   if any(dims == "level")
      summary.levels = sizes(dims == "level");
   end

   samples = readDepthSamples(filename, dims, sizes);
   if isempty(samples)
      return
   end
   summary.sample_min_m = min(samples, [], 'omitnan');
   summary.sample_max_m = max(samples, [], 'omitnan');
   summary.sample_policy = ...
      "first/last time columns only; full depth array not loaded";
end

function samples = readDepthSamples(filename, dims, sizes)
   %READDEPTHSAMPLES Return first/last time depth samples for coverage hints.
   samples = [];
   if isequal(dims, ["level", "time"])
      n_level = sizes(dims == "level");
      n_time = sizes(dims == "time");
      first = double(ncread(filename, 'Depth', [1 1], [n_level 1]));
      last = double(ncread(filename, 'Depth', [1 n_time], [n_level 1]));
      samples = [first(:); last(:)];
   elseif isequal(dims, ["time", "level"])
      n_level = sizes(dims == "level");
      n_time = sizes(dims == "time");
      first = double(ncread(filename, 'Depth', [1 1], [1 n_level]));
      last = double(ncread(filename, 'Depth', [n_time 1], [1 n_level]));
      samples = [first(:); last(:)];
   elseif isequal(dims, "level")
      samples = double(ncread(filename, 'Depth'));
   end
end

function provenance = provenanceByProduct(status, products)
   %PROVENANCEBYPRODUCT Read XML selected by the shared fetch-status authority.
   provenance = struct();
   for product = reshape(products, 1, [])
      rows = status(string({status.product}) == product);
      if isempty(rows)
         filename = "";
      else
         % Fetch validation already applies the optional-provenance lookup. Reuse
         % that path so data and XML have one cache-status authority and scan.
         filename = string(rows(1).provenance_file);
      end
      provenance.(char(product)) = readXmlProvenance(filename, product);
   end
end

function provenance = readXmlProvenance(filename, product)
   %READXMLPROVENANCE Extract stable DOI/title/citation snippets from EML XML.
   provenance = struct('product', product, 'doi', "", 'title', "", ...
      'citation', "", 'xml_file', string(filename));
   if filename == "" || ~isfile(filename)
      return
   end

   text = string(fileread(filename));
   provenance.doi = icemodel.verification.setup.regexpOnce( ...
      text, 'packageId="doi:([^"]+)"');
   provenance.title = icemodel.verification.setup.regexpOnce( ...
      text, '<title>(.*?)</title>');
   provenance.citation = icemodel.verification.setup.regexpOnce(text, ...
      '<para>(.*?Vandecrux.*?2020.*?</para>)');
end

function kind = productKind(product, suffix)
   %PRODUCTKIND Assign a finer-grained kind within broad product groups.
   if product ~= "simulated_firn"
      kind = product;
      return
   end

   suffix = string(suffix);
   if contains(suffix, "T_ice")
      kind = "simulated_temperature";
   elseif contains(suffix, "rho")
      kind = "simulated_density";
   elseif contains(suffix, "slwc")
      kind = "simulated_liquid_water";
   else
      kind = "simulated_compaction";
   end
end

function variables = canonicalVariables(kind)
   %CANONICALVARIABLES Map product kinds to comparison variable names.
   switch string(kind)
      case "firn_temperature"
         variables = "subsurface_temperature";
      case "simulated_temperature"
         variables = "subsurface_temperature";
      case "simulated_density"
         variables = "density";
      case "simulated_liquid_water"
         variables = "lwc";
      case "simulated_compaction"
         variables = "compaction_rate";
      otherwise
         variables = strings(1, 0);
   end
end
