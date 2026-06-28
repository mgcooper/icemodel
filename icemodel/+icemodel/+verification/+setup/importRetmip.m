function manifest = importRetmip(source_dir, kwargs)
   %IMPORTRETMIP Build the RetMIP family manifest.
   %
   %  manifest = icemodel.verification.setup.importRetmip(source_dir)
   %  manifest = icemodel.verification.setup.importRetmip(source_dir, ...
   %     dry_run=true)
   %
   % Role
   %  First-pass RetMIP staging hook. It records RetMIP protocol cases and their
   %  source associations without pretending the 3-hour protocol variables are a
   %  normal icemodel met timetable.

   arguments
      source_dir (1, 1) string = ""
      kwargs.case_ids (1, :) string = strings(1, 0)
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.dry_run (1, 1) logical = true
      kwargs.strict (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
   end

   % Validate the cache in non-strict mode by default so dry-run manifests can be
   % built on machines without the large RetMIP archives.
   [source_dir, cache_status] = icemodel.verification.setup.fetchRetmip( ...
      cache_dir=resolveCacheDir(source_dir), strict=kwargs.strict, silent=true);
   cases = icemodel.verification.setup.retmipCaseMetadata(kwargs.case_ids);

   % Resolve the family root even for dry runs; non-dry runs use it to stage
   % protocol userdata bundles beside the manifest.
   eval_root = icemodel.verification.helpers.evaluationDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   family_root = fullfile(eval_root, "retmip");

   % Convert RetMIP metadata rows to the existing firn manifest schema.
   entries = repmat(emptyEntry(), 1, numel(cases));
   for k = 1:numel(cases)
      entries(k) = manifestEntry(cases(k), source_dir, cache_status, ...
         family_root, kwargs.dry_run);
   end

   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "retmip", "10.22008/FK2/GZ3CSN;10.22008/FK2/CVPUJL", ...
      "https://doi.org/10.22008/FK2/GZ3CSN", ...
      "retmip-protocol+outputs", string(datetime('today')), entries);

   if kwargs.dry_run
      return
   end

   % Write only when explicitly requested; non-dry runs stage real protocol
   % userdata bundles and then merge the family manifest.
   icemodel.helpers.ensureDirExists(family_root);
   manifest_file = fullfile(family_root, "manifest.json");
   manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
      manifest_file, manifest, requested_ids=[entries.case_id], ...
      overwrite_family=kwargs.overwrite_family);
end

function cache_dir = resolveCacheDir(source_dir)
   %RESOLVECACHEDIR Resolve the RetMIP cache directory.
   if strlength(source_dir) > 0
      cache_dir = source_dir;
   else
      cache_dir = string(fullfile(icemodel.getpath('data'), ...
         'verification', 'retmip'));
   end
end

function entry = emptyEntry()
   %EMPTYENTRY Prototype RetMIP manifest entry.
   values = { ...
      ''
      'firn_observational'
      ''
      ''
      'unknown'
      {'firn'}
      'none'
      struct('lat_wgs84', NaN, 'lon_wgs84', NaN, ...
         'x_epsg3413', NaN, 'y_epsg3413', NaN, 'elev_m', NaN)
      struct('start', '', 'end', '')
      ''
      {}
      {'retmip_protocol'}
      {'tsfc', 'melt', 'snowf_subl', 'density', 'subsurface_temperature', 'lwc'}
      struct()
      struct()
      '3hr'
      ''};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
end

function entry = manifestEntry(c, source_dir, cache_status, family_root, dry_run)
   %MANIFESTENTRY Build one RetMIP firn case entry.
   files = retmipFiles(source_dir, c);
   evaluation_file = "";
   if ~dry_run
      evaluation_file = writeProtocolBundle(family_root, c, files);
   end

   retmip = struct();
   retmip.kind = 'protocol_userdata';
   retmip.staged = ~dry_run;
   retmip.source_dir = char(source_dir);
   retmip.protocol_timestep = char(c.protocol_timestep);
   retmip.surface_file = char(files.surface);
   retmip.profile_files = cellstr(files.profiles);
   retmip.model_output_files = cellstr(files.outputs);
   retmip.model_output_variables = cellstr(modelOutputVariables(files.outputs));
   retmip.cache_status = cache_status;
   if string(c.case_id) == "summit"
      retmip.native_met_status = 'pending_gcnet_import';
   end

   colocation = struct();
   colocation.retmip = retmip;
   colocation.source_association = c.source_association;
   if string(c.case_id) == "summit"
      colocation.native_met = struct('staged', false, ...
         'status', 'pending_gcnet_import');
   end

   comparison_variables = comparisonVariables(files.outputs);

   values = { ...
      char(c.case_id)
      'firn_observational'
      char(c.site_id)
      char(c.site_name)
      char(c.surface_zone)
      cellstr(string(c.eval_target))
      char(c.permafrost_zone)
      c.site_location
      c.period
      char(evaluation_file)
      {'retmip'}
      {'retmip_protocol'}
      cellstr(comparison_variables)
      struct('protocol_id', c.protocol_id, ...
         'retmip_station_id', c.retmip_station_id)
      colocation
      char(c.protocol_timestep)
      sprintf(['RetMIP protocol case; standard met source label is retmip ' ...
         'only after full meteorological source coverage is confirmed.'])};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
end

function files = retmipFiles(source_dir, c)
   %RETMIPFILES Resolve protocol/profile/output files for one RetMIP case.
   forcing_dir = fullfile(source_dir, "forcing");
   outputs_dir = fullfile(source_dir, "outputs");
   token = filenameToken(c.case_id);
   profile_token = profileFilenameToken(c.case_id);

   % The surface protocol file is required for a non-dry import; profiles are
   % optional by variable because not every case ships every initial state axis.
   files = struct( ...
      'surface', firstMatch(forcing_dir, "RetMIP_surface_forcing_" + token + ".tab"), ...
      'profiles', profileMatches(forcing_dir, profile_token), ...
      'outputs', outputMatches(outputs_dir, token));
end

function relpath = writeProtocolBundle(family_root, c, files)
   %WRITEPROTOCOLBUNDLE Stage one RetMIP observations.mat protocol bundle.
   if strlength(files.surface) == 0
      error('icemodel:verification:importRetmip:missingSurfaceFile', ...
         'missing RetMIP surface protocol file for %s', c.case_id);
   end

   % Parse the official surface series and any initial state profile tables into
   % one data-only target bundle. Model outputs remain indexed as file products.
   [surface, surface_meta] = ...
      icemodel.verification.setup.readRetmipProtocolTable(files.surface);
   profiles = struct();
   profile_meta = struct();
   for k = 1:numel(files.profiles)
      [profile, meta] = ...
         icemodel.verification.setup.readRetmipProfileTable(files.profiles(k));
      name = profileName(files.profiles(k));
      profiles.(name) = profile;
      profile_meta.(name) = meta;
   end

   targets = struct( ...
      'format', 'retmip_protocol_bundle', ...
      'data', struct('surface', surface, 'profiles', profiles, ...
         'model_output_files', files.outputs), ...
      'metadata', struct('surface', surface_meta, ...
         'profiles', profile_meta, ...
         'source_kind', 'RetMIP protocol userdata and comparison products'));

   case_root = fullfile(family_root, c.case_id);
   icemodel.helpers.ensureDirExists(case_root);
   relpath = fullfile(c.case_id, "observations.mat");
   save(fullfile(family_root, relpath), 'targets');
end

function token = filenameToken(case_id)
   %FILENAMETOKEN Map filename-safe case ids to RetMIP source filename tokens.
   switch string(case_id)
      case "kanu"
         token = "KANU";
      case "dye2_long"
         token = "Dye-2_long";
      case "dye2_2016"
         token = "Dye-2_16";
      case "summit"
         token = "Summit";
      case "fa"
         token = "FA";
      otherwise
         error('icemodel:verification:importRetmip:unknownCase', ...
            'unknown RetMIP case id: %s', case_id);
   end
end

function token = profileFilenameToken(case_id)
   %PROFILEFILENAMETOKEN Map case ids to RetMIP profile filename tokens.
   token = filenameToken(case_id);
   if string(case_id) == "kanu"
      token = "KAN-U";
   end
end

function filename = firstMatch(folder, pattern)
   %FIRSTMATCH Return the first matching file or an empty string.
   hits = dir(fullfile(folder, pattern));
   if isempty(hits)
      filename = "";
   else
      filename = string(fullfile(hits(1).folder, hits(1).name));
   end
end

function files = profileMatches(folder, token)
   %PROFILEMATCHES Return initial density/temperature/LWC files for one case.
   kinds = ["density", "temperature", "lwc"];
   files = strings(1, numel(kinds));
   n_files = 0;
   for k = 1:numel(kinds)
      hit = firstMatch(folder, ...
         "RetMIP_initial_firn_" + kinds(k) + "_" + token + ".tab");
      if strlength(hit) > 0
         n_files = n_files + 1;
         files(n_files) = hit;
      end
   end
   files = files(1:n_files);
end

function files = outputMatches(folder, token)
   %OUTPUTMATCHES Return RetMIP model-output NetCDF files for one case.
   hits = [dir(fullfile(folder, "*_" + token + "_3hourly_*.nc")); ...
      dir(fullfile(folder, "*_" + upper(token) + "_3hourly_*.nc"))];
   files = strings(1, numel(hits));
   for k = 1:numel(hits)
      files(k) = string(fullfile(hits(k).folder, hits(k).name));
   end
   files = unique(files, 'stable');
end

function variables = comparisonVariables(output_files)
   %COMPARISONVARIABLES Return RetMIP userdata plus discovered output variables.
   protocol_vars = ["tsfc", "melt", "snowf_subl", "density", ...
      "subsurface_temperature", "lwc"];
   variables = unique([protocol_vars, modelOutputVariables(output_files).'], ...
      'stable');
end

function variables = modelOutputVariables(output_files)
   %MODELOUTPUTVARIABLES Read model-output headers when they are valid NetCDFs.
   if isempty(output_files)
      variables = strings(0, 1);
      return
   end

   collected = cell(1, numel(output_files));
   for k = 1:numel(output_files)
      try
         inventory = ...
            icemodel.verification.setup.retmipOutputInventory(output_files(k));
         collected{k} = inventory.variables(:);
      catch
         % Invalid or placeholder output files are still useful path evidence.
         collected{k} = strings(0, 1);
      end
   end
   nonempty = ~cellfun(@isempty, collected);
   if ~any(nonempty)
      variables = strings(0, 1);
      return
   end
   variables = vertcat(collected{nonempty});
   variables = unique(variables, 'stable');
end

function name = profileName(filename)
   %PROFILENAME Build a valid struct field for a RetMIP profile file.
   [~, stem] = fileparts(filename);
   if contains(stem, "density")
      name = 'density';
   elseif contains(stem, "temperature")
      name = 'temperature';
   elseif contains(stem, "lwc")
      name = 'lwc';
   else
      name = matlab.lang.makeValidName(stem);
   end
end
