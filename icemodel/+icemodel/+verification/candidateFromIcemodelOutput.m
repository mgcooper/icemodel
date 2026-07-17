function candidate = candidateFromIcemodelOutput(ice1, ice2, opts, case_manifest)
   %CANDIDATEFROMICEMODELOUTPUT Convert icemodel outputs for verification.
   %
   %  candidate = icemodel.verification.candidateFromIcemodelOutput( ...
   %     ice1, ice2, opts, case_manifest)
   %
   % This adapter is the stable handoff point between model execution and the
   % verification metrics. Real snow-model development should add model outputs
   % to ICE1/ICE2, then map them here only when names or units differ from the
   % verification schema.

   arguments
      ice1 (1, 1) struct
      ice2 (1, 1) struct
      opts (1, 1) struct
      case_manifest (1, 1) struct
   end

   switch case_manifest.case_type
      case "esm_site"
         candidate = timeseriesCandidateFromIce1(ice1, ice2, opts, case_manifest);
      case "firn_observational"
         candidate = firnCandidateFromIce1(ice1, ice2, opts, case_manifest);
      case "synthetic_process"
         candidate = experimentBundleCandidateFromIce1(ice1, opts);
      otherwise
         error('unsupported verification case type: %s', ...
            case_manifest.case_type)
   end
end

function candidate = firnCandidateFromIce1(ice1, ice2, opts, case_manifest)
   %FIRNCANDIDATEFROMICE1 Map ICE1/ICE2 fields to firn verification targets.
   %
   % Firn-observational cases compare against PROMICE station Data. The
   % comparison axes the staged manifests declare are the surface-energy /
   % ablation series (ablation, snow_depth, tsfc) plus the subsurface
   % thermistor profile (tice1..tice8 = T(z,t) sampled at thermistor depths)
   % and, when staged, density rho(z,t) and smb (surface mass balance). Variables
   % are mapped only as far as the staged cases require; porosity / saturation /
   % runoff are deferred with the firn physics.

   if ~isfield(ice1, "Time")
      error('icemodel output ice1 must contain Time for firn verification')
   end

   variable_names = string(case_manifest.comparison_variables);

   % RetMIP protocol targets use a protocol-bundle schema, not the SUMup-style
   % profile bundle, even though both declare density/subsurface-temperature.
   if isRetmipProtocolCase(case_manifest)
      candidate = retmipProtocolCandidateFromOutput( ...
         ice1, ice2, opts, variable_names);
      return
   end

   % SUMup firn-observational cases declare profile-bundle comparison axes
   % (density rho(z), subsurface_temperature T(z,t), smb) rather than
   % the PROMICE thermistor/ablation timeseries. Dispatch on the declared axes:
   % a case that asks for the profile-bundle variables gets a profile-bundle
   % candidate; the PROMICE thermistor timeseries stays the default firn path.
   if any(ismember(variable_names, ["density", "subsurface_temperature", "smb"]))
      candidate = firnProfileCandidateFromIce2(ice1, ice2, opts, case_manifest);
      return
   end

   data = timetable('RowTimes', ice1.Time);
   Tf = icemodel.physicalConstant('Tf');

   for name = variable_names(:)'
      tice_tok = regexp(char(name), '^tice(\d+)$', 'tokens', 'once');
      if isfield(ice1, name)
         % Direct passthrough when the model already emits the named series.
         data.(name) = ice1.(name);
      elseif name == "tsfc" && isfield(ice1, "Tsfc")
         % Surface temperature: icemodel stores Kelvin Tsfc; firn obs is
         % degrees C (PROMICE tsfc). Convert once here.
         data.tsfc = ice1.Tsfc - Tf;
      elseif name == "snow_depth" && isfield(ice1, "snow_depth_m")
         data.snow_depth = ice1.snow_depth_m;
      elseif name == "tice10m"
         % PROMICE tice10m is standardized at 10 m below the evolving surface,
         % independent of the station's depth-tagged thermistor string.
         values = ticeFromIce2AtDepth(ice2, opts, 10);
         if ~isempty(values)
            data.tice10m = values;
         end
      elseif ~isempty(tice_tok)
         % Subsurface thermistor: sample the column T(z,t) at the k-th
         % staged thermistor depth and convert to degrees C.
         values = ticeFromIce2(ice2, opts, case_manifest, ...
            str2double(tice_tok{1}));
         if ~isempty(values)
            data.(name) = values;
         end
      end
   end

   candidate = struct( ...
      "format", "timeseries", ...
      "data", data, ...
      "metadata", metadata(opts, "icemodel_output"));
end

function candidate = firnProfileCandidateFromIce2(ice1, ice2, opts, case_manifest)
   %FIRNPROFILECANDIDATEFROMICE2 Map ICE2 column state to SUMup firn profiles.
   %
   % SUMup firn cases compare against firn observation profiles: a density
   % profile rho(z) and a subsurface temperature profile T(z,t). icemodel
   % carries these as column state in ice2 (T = column temperature [K], f_liq /
   % density depending on the run). Map only the axes the staged case declares;
   % an axis whose source column is unavailable is left out so it is reported as
   % a missing-candidate diagnostic (soft gate), never fabricated. The candidate
   % format is "subsurface_profile_bundle", matching the observation target so
   % the soft firn comparison can align them by depth.

   variable_names = string(case_manifest.comparison_variables);
   Tf = icemodel.physicalConstant('Tf');

   bundle = struct('format', 'subsurface_profile_bundle');

   % Build the depth axis from the column node spacing when available.
   if isfield(opts, 'dz_thermal') && ~isempty(opts.dz_thermal)
      dz = opts.dz_thermal;
   else
      dz = NaN;
   end

   if ismember("density", variable_names)
      rho = columnProfile(ice2, "ro_sno", "density");
      if ~isempty(rho)
         bundle.density = datedProfileTable( ...
            rho, dz, ice1.Time, "density", "kg m-3");
      end
   end

   if ismember("subsurface_temperature", variable_names)
      if isfield(ice2, 'T') && ~isempty(ice2.T)
         % Column temperature T(z,t) in degrees C, depth-indexed.
         bundle.subsurface_temperature = datedProfileTable( ...
            ice2.T - Tf, dz, ice1.Time, ...
            "subsurface_temperature", "degC");
      end
   end

   if ismember("smb", variable_names) && isfield(ice1, "smb")
      bundle.smb = timetable(ice1.smb(:), 'RowTimes', ice1.Time(:), ...
         'VariableNames', {'smb'});
   end

   candidate = struct( ...
      "format", "subsurface_profile_bundle", ...
      "data", bundle, ...
      "metadata", metadata(opts, "icemodel_output"));
end

function tf = isRetmipProtocolCase(case_manifest)
   %ISRETMIPPROTOCOLCASE True when a firn case compares to RetMIP protocol data.
   tf = (isfield(case_manifest, 'dataset_family') ...
      && string(case_manifest.dataset_family) == "retmip") ...
      || (isfield(case_manifest, 'eval_sources') ...
      && any(string(case_manifest.eval_sources) == "retmip_protocol"));
end

function candidate = retmipProtocolCandidateFromOutput( ...
      ice1, ice2, opts, variable_names)
   %RETMIPPROTOCOLCANDIDATEFROMOUTPUT Build a RetMIP-shaped model candidate.
   data = struct('surface', timetable(), 'profiles', struct());
   if isfield(ice1, "Time")
      data.surface = timetable('RowTimes', ice1.Time);
   end
   for name = reshape(variable_names, 1, [])
      if isfield(ice1, name)
         data.surface.(char(name)) = ice1.(char(name));
      elseif name == "tsfc" && isfield(ice1, "Tsfc")
         data.surface.tsfc = ice1.Tsfc;
      end
   end
   if ismember("density", variable_names)
      rho = columnProfile(ice2, "ro_sno", "density");
      if ~isempty(rho)
         data.profiles.density = table( ...
            depthAxis(size(rho, 1), profileDz(opts)), rho(:, 1), ...
            'VariableNames', {'depth', 'density'});
      end
   end
   if ismember("subsurface_temperature", variable_names) ...
         && isfield(ice2, "T") && ~isempty(ice2.T)
      data.profiles.subsurface_temperature = table( ...
         depthAxis(size(ice2.T, 1), profileDz(opts)), ...
         ice2.T(:, 1), ...
         'VariableNames', {'depth', 'subsurface_temperature'});
   end
   candidate = struct( ...
      "format", "retmip_protocol_bundle", ...
      "data", data, ...
      "metadata", metadata(opts, "icemodel_output"));
end

function dz = profileDz(opts)
   %PROFILEDZ Return thermal-grid spacing when the run options carry it.
   if isfield(opts, 'dz_thermal') && ~isempty(opts.dz_thermal)
      dz = opts.dz_thermal;
   else
      dz = NaN;
   end
end

function values = columnProfile(ice2, primary, fallback)
   %COLUMNPROFILE Return a column-state field by primary or fallback name.
   values = [];
   if isfield(ice2, primary) && ~isempty(ice2.(primary))
      values = ice2.(primary);
   elseif isfield(ice2, fallback) && ~isempty(ice2.(fallback))
      values = ice2.(fallback);
   end
end

function tbl = datedProfileTable(values, dz, model_times, varname, unit)
   %DATEDPROFILETABLE Preserve every model profile column with its UTC time.
   values = double(values);
   if ~ismatrix(values) || isempty(values)
      error('icemodel:verification:candidateFromIcemodelOutput:badProfile', ...
         '%s model profiles must be a nonempty depth-by-time matrix', varname)
   end

   % Model column state is depth-by-time. A single stored column retains the
   % first corresponding output timestamp; multi-column histories must map
   % one-to-one to the model time axis so no profile can be silently discarded.
   model_times = icemodel.verification.setup.ensureUtc(model_times(:));
   n_profiles = size(values, 2);
   if n_profiles == 1
      profile_times = model_times(1);
   elseif n_profiles == numel(model_times)
      profile_times = model_times;
   else
      error('icemodel:verification:candidateFromIcemodelOutput:profileTimeSize', ...
         ['%s has %d profile columns but the model output has %d ' ...
         'timestamps'], varname, n_profiles, numel(model_times))
   end

   % Tall rows keep profile identity, timestamp, depth, value, and units
   % available to the reusable date/identity matcher.
   nz = size(values, 1);
   profile_id = repmat("icemodel", nz * n_profiles, 1);
   profile_time = reshape(repelem(profile_times(:), nz), [], 1);
   depth = repmat(depthAxis(nz, dz), n_profiles, 1);
   tbl = table(profile_id, profile_time, depth, values(:), ...
      'VariableNames', {'profile_id', 'datetime', 'depth', char(varname)});
   tbl.Properties.VariableUnits = {'', '', 'm', char(unit)};
end

function z = depthAxis(nz, dz)
   %DEPTHAXIS Node-center depth axis for an nz-node column at spacing dz.
   if isnan(dz)
      z = (0:nz - 1)';
   else
      z = (0:nz - 1)' * dz;
   end
end

function values = ticeFromIce2(ice2, opts, case_manifest, k)
   %TICEFROMICE2 Sample column T at the k-th manifested thermistor depth.
   %
   % Returns Celsius column-T at depth manifest.observation_variables
   % .thermistor_depths_m(k), or [] when the depth/grid is unavailable.
   % Mirrors soilTempFromIce2 (esm_site soil temps) for the firn thermistor
   % string. When the manifest does not record explicit depths, the function
   % returns [] so the variable is reported as missing rather than fabricated.

   values = [];
   if ~isfield(ice2, 'T') || isempty(ice2.T)
      return
   end
   obsvars = case_manifest.observation_variables;
   if ~isstruct(obsvars) || ~isfield(obsvars, 'thermistor_depths_m')
      return
   end
   depths = obsvars.thermistor_depths_m;
   if k < 1 || k > numel(depths)
      return
   end
   dz = opts.dz_thermal;
   zidx = max(1, min(size(ice2.T, 1), round(depths(k) / dz) + 1));
   Tf = icemodel.physicalConstant('Tf');
   values = ice2.T(zidx, :)' - Tf;
end

function values = ticeFromIce2AtDepth(ice2, opts, depth_m)
   %TICEFROMICE2ATDEPTH Sample Kelvin column T at a depth below the surface.

   values = [];
   if ~isfield(ice2, 'T') || isempty(ice2.T) ...
         || ~isfield(opts, 'dz_thermal') || isempty(opts.dz_thermal)
      return
   end
   zidx = max(1, min(size(ice2.T, 1), round(depth_m / opts.dz_thermal) + 1));
   values = ice2.T(zidx, :)';
end

function candidate = timeseriesCandidateFromIce1(ice1, ice2, opts, case_manifest)
   %TIMESERIESCANDIDATEFROMICE1 Map ICE1 fields to ESM verification targets.

   if ~isfield(ice1, "Time")
      error('icemodel output ice1 must contain Time for site verification')
   end

   variable_names = case_manifest.comparison_variables;
   data = timetable('RowTimes', ice1.Time);

   for name = variable_names(:)'
      if isfield(ice1, name)
         data.(name) = ice1.(name);
      elseif name == "snow_depth_m" && isfield(ice1, "snow_depth")
         data.snow_depth_m = ice1.snow_depth;
      elseif name == "swe_kg_m2" ...
            && isfield(ice1, "snow_depth_m") ...
            && isfield(ice1, "snow_density_kg_m3")
         data.swe_kg_m2 = ice1.snow_depth_m .* ice1.snow_density_kg_m3;
      elseif name == "swe_kg_m2" ...
            && isfield(ice1, "snow_depth") ...
            && isfield(ice1, "snow_density_kg_m3")
         data.swe_kg_m2 = ice1.snow_depth .* ice1.snow_density_kg_m3;
      elseif name == "surface_temp_C" && isfield(ice1, "Tsfc")
         Tf = icemodel.physicalConstant('Tf');
         data.surface_temp_C = ice1.Tsfc - Tf;
      else
         soil_tok = regexp(char(name), '^soil_temp_(\d+)_C$', 'tokens', 'once');
         if ~isempty(soil_tok)
            values = soilTempFromIce2(ice2, opts, case_manifest, ...
               str2double(soil_tok{1}));
            if ~isempty(values)
               data.(name) = values;
            end
         end
      end
   end

   candidate = struct( ...
      "format", "timeseries", ...
      "data", data, ...
      "metadata", metadata(opts, "icemodel_output"));
end

function candidate = experimentBundleCandidateFromIce1(ice1, opts)
   %EXPERIMENTBUNDLECANDIDATEFROMICE1 Map ICE1 process-benchmark payloads.

   if ~isfield(ice1, "verification_experiments")
      error(['icemodel output ice1 must contain verification_experiments ' ...
         'for synthetic process verification'])
   end

   candidate = struct( ...
      "format", "experiment_bundle", ...
      "experiments", ice1.verification_experiments, ...
      "metadata", metadata(opts, "icemodel_output"));
end

function values = soilTempFromIce2(ice2, opts, case_manifest, k)
   %SOILTEMPFROMICE2 Sample column T at the k-th manifested soil depth.
   %
   % Returns Celsius column-T at depth manifest.observation_variables
   % .soil_depths_m(k), or [] if the depth/grid is unavailable. icemodel
   % does not model a separate soil layer; the column temperature is
   % sampled at the requested depth using opts.dz_thermal spacing.

   values = [];
   if ~isfield(ice2, 'T') || isempty(ice2.T)
      return
   end
   obsvars = case_manifest.observation_variables;
   if ~isstruct(obsvars) || ~isfield(obsvars, 'soil_depths_m')
      return
   end
   depths = obsvars.soil_depths_m;
   if k < 1 || k > numel(depths)
      return
   end
   dz = opts.dz_thermal;
   zidx = max(1, min(size(ice2.T, 1), round(depths(k) / dz) + 1));
   Tf = icemodel.physicalConstant('Tf');
   values = ice2.T(zidx, :)' - Tf;
end

function info = metadata(opts, source)
   %METADATA Record enough provenance to diagnose candidate generation.

   info = struct( ...
      "source", source, ...
      "smbmodel", string(opts.smbmodel), ...
      "sitename", string(opts.sitename), ...
      "simyears", opts.simyears);

   if isfield(opts, "verification_candidate_source")
      info.verification_candidate_source = opts.verification_candidate_source;
   end
end
