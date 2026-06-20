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
   % and, when staged, density rho(z,t) and accumulation/SMB. Variables are
   % mapped only as far as the staged cases require; porosity / saturation /
   % runoff are deferred with the firn physics.

   if ~isfield(ice1, "Time")
      error('icemodel output ice1 must contain Time for firn verification')
   end

   variable_names = string(case_manifest.comparison_variables);
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
