function compatibility = comparisonCompatibility(case_ids, kwargs)
   %COMPARISONCOMPATIBILITY Derive staged verification comparison pairs.
   %
   %  compatibility = icemodel.verification.comparisonCompatibility()
   %  compatibility = icemodel.verification.comparisonCompatibility("kanm", ...
   %     dataset_family="promice")
   %
   % Role
   %  Additive discovery helper for staged verification artifacts. It reads
   %  family manifests, inspects MAT-file headers and NetCDF headers where files
   %  are available, and reports possible model-data, model-model, and data-data
   %  comparison pairs from common canonical variables. Manifest fields remain
   %  complete provenance declarations; this helper must not be used to
   %  under-populate manifests.
   %
   % See also: icemodel.verification.listcases,
   %  icemodel.verification.comparecase

   arguments
      case_ids (1, :) string = strings(1, 0)
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
      kwargs.dataset_family (1, 1) string = ""
   end

   cases = icemodel.verification.listcases( ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      input_data_root=kwargs.input_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename, ...
      dataset_family=kwargs.dataset_family);
   cases = attachInputRoot(cases, kwargs.input_data_root);
   cases = filterCases(cases, case_ids);
   if isempty(cases)
      compatibility = struct('cases', repmat(emptyCase(), 0, 1), ...
         'pairs', repmat(emptyPair(), 0, 1));
      return
   end

   reports = repmat(emptyCase(), numel(cases), 1);
   pair_groups = cell(numel(cases), 1);
   for k = 1:numel(cases)
      reports(k) = caseCompatibility(cases(k));
      pair_groups{k} = reports(k).pairs(:);
   end
   pair_groups = pair_groups(~cellfun(@isempty, pair_groups));
   if isempty(pair_groups)
      all_pairs = repmat(emptyPair(), 0, 1);
   else
      all_pairs = vertcat(pair_groups{:});
   end
   compatibility = struct('cases', reports, 'pairs', all_pairs);
end

%% Local functions
function cases = attachInputRoot(cases, input_data_root)
   %ATTACHINPUTROOT Carry an explicit input root into artifact resolution.
   if isempty(cases) || input_data_root == ""
      return
   end
   [cases.input_data_root] = deal(char(input_data_root));
end

function cases = filterCases(cases, case_ids)
   %FILTERCASES Apply an optional case-id selector to listed manifests.
   if isempty(cases) || isempty(case_ids)
      return
   end
   ids = string({cases.case_id});
   keep = ismember(ids, case_ids);
   cases = cases(keep);
end

function report = caseCompatibility(case_manifest)
   %CASECOMPATIBILITY Build artifact and pair records for one case.
   artifacts = collectArtifacts(case_manifest);
   pairs = derivePairs(case_manifest, artifacts);
   report = emptyCase();
   report.case_id = string(case_manifest.case_id);
   report.dataset_family = string(case_manifest.dataset_family);
   report.declared = declaredMetadata(case_manifest);
   report.artifacts = artifacts;
   report.pairs = pairs;
end

function report = emptyCase()
   %EMPTYCASE Prototype for one case compatibility report.
   report = struct( ...
      'case_id', "", ...
      'dataset_family', "", ...
      'declared', struct(), ...
      'artifacts', emptyArtifact(), ...
      'pairs', emptyPair());
end

function metadata = declaredMetadata(case_manifest)
   %DECLAREDMETADATA Preserve manifest declarations separately from discovery.
   metadata = struct( ...
      'comparison_variables', fieldStrings(case_manifest, ...
      'comparison_variables'), ...
      'observation_variables', icemodel.verification.helpers.fieldOr( ...
      case_manifest, ...
      'observation_variables', struct()), ...
      'forcing_sources', fieldStrings(case_manifest, 'forcing_sources'), ...
      'eval_sources', fieldStrings(case_manifest, 'eval_sources'), ...
      'colocation', icemodel.verification.helpers.fieldOr( ...
      case_manifest, 'colocation', struct()));
end

function artifacts = collectArtifacts(case_manifest)
   %COLLECTARTIFACTS Gather staged artifact records for one manifest case.
   artifact_parts = cell(maxArtifacts(case_manifest), 1);
   n_parts = 0;
   if isfield(case_manifest, 'evaluation_path') ...
         && strlength(string(case_manifest.evaluation_path)) > 0
      n_parts = n_parts + 1;
      artifact_parts{n_parts} = artifactRecord(case_manifest, ...
         primaryEvalSource(case_manifest), "data", "evaluation", ...
         string(case_manifest.evaluation_path));
   end
   reference_path = topLevelReferencePath(case_manifest);
   if reference_path ~= ""
      n_parts = n_parts + 1;
      artifact_parts{n_parts} = artifactRecord(case_manifest, ...
         "reference", "model", "reference", reference_path);
   end

   if ~isfield(case_manifest, 'colocation') ...
         || ~isstruct(case_manifest.colocation)
      artifact_parts = artifact_parts(1:n_parts);
      artifacts = combineArtifacts(artifact_parts);
      return
   end

   names = string(fieldnames(case_manifest.colocation));
   for source = reshape(names, 1, [])
      leg = case_manifest.colocation.(char(source));
      if ~isstruct(leg) || ~stagedLeg(leg)
         continue
      end
      n_parts = n_parts + 1;
      artifact_parts{n_parts} = legArtifacts(case_manifest, source, ...
         leg);
   end
   artifact_parts = artifact_parts(1:n_parts);
   artifacts = combineArtifacts(artifact_parts);
end

function n = maxArtifacts(case_manifest)
   %MAXARTIFACTS Return the maximum artifact chunks for one manifest case.
   n = 2;
   if isfield(case_manifest, 'colocation') && isstruct(case_manifest.colocation)
      n = n + numel(fieldnames(case_manifest.colocation));
   end
end

function pathname = topLevelReferencePath(case_manifest)
   %TOPLEVELREFERENCEPATH Return a case-level reference artifact path if present.
   pathname = "";
   values = fieldStrings(case_manifest, 'reference_path');
   values = values(strlength(values) > 0);
   if isempty(values)
      values = fieldStrings(case_manifest, 'reference_file');
      values = values(strlength(values) > 0);
   end
   if ~isempty(values)
      pathname = values(1);
   end
end

function tf = stagedLeg(leg)
   %STAGEDLEG True when a colocation leg represents staged data.
   tf = false;
   if isfield(leg, 'staged')
      tf = logical(leg.staged);
   end
end

function artifacts = legArtifacts(case_manifest, source, leg)
   %LEGARTIFACTS Convert one colocation source leg into artifact records.
   artifact_source = legSourceLabel(source, leg);
   specs = legArtifactSpecs(source);

   artifact_parts = cell(maxLegArtifacts(leg, specs), 1);
   n_parts = 0;
   for spec = reshape(specs, 1, [])
      values = legFileValues(leg, spec.field);
      for value = reshape(values, 1, [])
         n_parts = n_parts + 1;
         artifact_parts{n_parts} = artifactRecord(case_manifest, ...
            artifact_source, artifactRole(artifact_source, spec.kind), ...
            spec.kind, value);
      end
   end

   if n_parts == 0 && isRcmSource(source)
      artifacts = repmat(emptyArtifact(), 0, 1);
   elseif n_parts == 0
      artifacts = artifactRecord(case_manifest, artifact_source, ...
         sourceRole(artifact_source), "metadata", "");
   else
      artifact_parts = artifact_parts(1:n_parts);
      artifacts = vertcat(artifact_parts{:});
   end
end

function specs = legArtifactSpecs(source)
   %LEGARTIFACTSPECS Return eval artifact fields for a colocation leg.
   specs = [ ...
      struct('field', "obs_file", 'kind', "observation")
      struct('field', "evaluation_file", 'kind', "observation")
      struct('field', "met_files", 'kind', "met")
      struct('field', "data_files", 'kind', "userdata")
      struct('field', "model_output_files", 'kind', "model_output")
      struct('field', "reference_file", 'kind', "model_output")];

   % RCM met files are selectable forcing artifacts, not evaluation artifacts.
   % Until met-channel comparison is explicitly designed, RCM comparison
   % discovery is driven only by userdata/Data or model-output files.
   if isRcmSource(source)
      specs = specs(string({specs.field}) ~= "met_files");
   end
end

function tf = isRcmSource(source)
   %ISRCMSOURCE True for internal RCM source keys and public product ids.
   rcm_sources = [icemodel.verification.namelists.rcmsources(), ...
      icemodel.verification.namelists.rcmProductIds()];
   tf = ismember(string(source), rcm_sources);
end

function label = legSourceLabel(source, leg)
   %LEGSOURCELABEL Return the public source id for one colocation leg.
   source = string(source);
   if isstruct(leg) && isfield(leg, 'source_id') ...
         && strlength(string(leg.source_id)) > 0
      label = string(leg.source_id);
   elseif ismember(source, icemodel.verification.namelists.rcmsources())
      label = icemodel.verification.namelists.rcmProductIds(source);
   else
      label = source;
   end
end

function n = maxLegArtifacts(leg, specs)
   %MAXLEGARTIFACTS Return the maximum file records for one colocation leg.
   n = 0;
   for spec = reshape(specs, 1, [])
      n = n + numel(legFileValues(leg, spec.field));
   end
end

function artifacts = combineArtifacts(artifact_parts)
   %COMBINEARTIFACTS Preserve the empty artifact prototype for empty cases.
   if isempty(artifact_parts)
      artifacts = repmat(emptyArtifact(), 0, 1);
   else
      artifacts = vertcat(artifact_parts{:});
      artifacts = deduplicateFileArtifacts(artifacts);
   end
end

function artifacts = deduplicateFileArtifacts(artifacts)
   %DEDUPLICATEFILEARTIFACTS Keep the first record for each nonempty file/kind.
   keep = true(numel(artifacts), 1);
   seen = strings(numel(artifacts), 1);
   n_seen = 0;
   for k = 1:numel(artifacts)
      if strlength(string(artifacts(k).path)) == 0
         continue
      end
      key = dedupKind(artifacts(k).kind) + "|" + string(artifacts(k).path);
      if any(seen(1:n_seen) == key)
         keep(k) = false;
      else
         n_seen = n_seen + 1;
         seen(n_seen) = key;
      end
   end
   artifacts = artifacts(keep);
end

function kind = dedupKind(kind)
   %DEDUPKIND Treat top-level and colocation observation records as one file kind.
   kind = string(kind);
   if any(kind == ["evaluation", "observation"])
      kind = "observation";
   end
end

function values = legFileValues(leg, fieldname)
   %LEGFILEVALUES Return one colocation file field as a string row vector.
   values = strings(1, 0);
   fieldname = char(fieldname);
   if ~isfield(leg, fieldname)
      return
   end
   raw = leg.(fieldname);
   if isempty(raw)
      return
   end
   values = reshape(string(raw), 1, []);
end

function artifact = emptyArtifact()
   %EMPTYARTIFACT Prototype for one artifact compatibility record.
   artifact = struct( ...
      'source', "", ...
      'role', "", ...
      'kind', "", ...
      'path', "", ...
      'exists', false, ...
      'artifact_variables', strings(1, 0), ...
      'declared_variables', strings(1, 0), ...
      'variables', strings(1, 0), ...
      'evidence', "");
end

function artifact = artifactRecord(case_manifest, source, role, kind, pathname)
   %ARTIFACTRECORD Build one artifact record from headers and declarations.
   source = string(source);
   kind = string(kind);
   pathname = resolvePath(case_manifest, pathname, kind);
   [artifact_variables, evidence, exists] = artifactVariables(pathname);
   declared_variables = declaredVariables(case_manifest, source, kind);
   if exists
      variables = unique([artifact_variables(:); declared_variables(:)], ...
         'stable');
   else
      variables = strings(0, 1);
   end
   artifact = struct( ...
      'source', source, ...
      'role', string(role), ...
      'kind', kind, ...
      'path', pathname, ...
      'exists', exists, ...
      'artifact_variables', reshape(artifact_variables, 1, []), ...
      'declared_variables', reshape(declared_variables, 1, []), ...
      'variables', reshape(variables, 1, []), ...
      'evidence', evidence);
end

function source = primaryEvalSource(case_manifest)
   %PRIMARYEVALSOURCE Pick the most specific evaluation source declaration.
   eval_sources = fieldStrings(case_manifest, 'eval_sources');
   obs = eval_sources(endsWith(eval_sources, "_obs") ...
      | eval_sources == "retmip_protocol");
   if ~isempty(obs)
      source = obs(1);
   else
      source = string(case_manifest.dataset_family);
   end
end

function role = sourceRole(source)
   %SOURCEROLE Classify source legs for comparison-pair classes.
   rcm_sources = [icemodel.verification.namelists.rcmsources(), ...
      icemodel.verification.namelists.rcmProductIds()];
   if ismember(string(source), rcm_sources)
      role = "model";
   else
      role = "data";
   end
end

function role = artifactRole(source, kind)
   %ARTIFACTROLE Classify a concrete artifact for comparison-pair classes.
   if ismember(string(kind), ["model_output", "reference"])
      role = "model";
   else
      role = sourceRole(source);
   end
end

function pathname = resolvePath(case_manifest, pathname, kind)
   %RESOLVEPATH Resolve portable manifest paths to local files when possible.
   pathname = string(pathname);
   if pathname == "" || isfile(pathname)
      return
   end

   candidates = string.empty(1, 0);
   if isfield(case_manifest, 'input_data_root') ...
         && strlength(string(case_manifest.input_data_root)) > 0
      input_root = string(case_manifest.input_data_root);
      candidates = [candidates, ...
         fullfile(input_root, pathname), ...
         fullfile(input_root, 'met', pathname), ...
         fullfile(input_root, 'userdata', pathname)];
   end
   if isfield(case_manifest, 'family_root')
      family_root = string(case_manifest.family_root);
      eval_root = string(fileparts(family_root));
      data_root = string(fileparts(eval_root));
      candidates = [candidates, ...
         fullfile(family_root, pathname), ...
         fullfile(eval_root, pathname), ...
         fullfile(data_root, pathname), ...
         fullfile(data_root, 'input', 'met', pathname), ...
         fullfile(data_root, 'input', 'userdata', pathname)];
   end

   if kind == "met" && isfield(case_manifest, 'family_root')
      data_root = string(fileparts(fileparts(string(case_manifest.family_root))));
      candidates(end + 1) = fullfile(data_root, 'input', 'met', pathname);
   elseif kind == "userdata" && isfield(case_manifest, 'family_root')
      data_root = string(fileparts(fileparts(string(case_manifest.family_root))));
      candidates(end + 1) = fullfile(data_root, 'input', 'userdata', pathname);
   end

   hit = find(isfile(candidates), 1);
   if ~isempty(hit)
      pathname = candidates(hit);
   end
end

function [variables, evidence, exists] = artifactVariables(pathname)
   %ARTIFACTVARIABLES Read MAT/NetCDF headers without loading payload arrays.
   variables = strings(1, 0);
   evidence = "missing_artifact";
   exists = false;
   if pathname == "" || ~isfile(pathname)
      return
   end

   exists = true;
   [~, ~, ext] = fileparts(pathname);
   switch lower(string(ext))
      case ".nc"
         info = ncinfo(pathname);
         names = string({info.Variables.Name});
         variables = canonicalizeVariables(names);
         evidence = "netcdf_header";
      case ".mat"
         info = whos('-file', pathname);
         names = string({info.name});
         variables = canonicalizeVariables(names);
         payload_variables = matPayloadVariables(pathname, info);
         if ~isempty(payload_variables)
            variables = unique([payload_variables(:); variables(:)], ...
               'stable');
            evidence = "mat_payload";
         else
            evidence = "mat_header";
         end
      otherwise
         evidence = "file_exists";
   end
end

function variables = matPayloadVariables(pathname, info)
   %MATPAYLOADVARIABLES Inspect staged MAT containers for table channels.
   variables = strings(0, 1);
   names = string({info.name});
   containers = intersect(["met", "Data", "targets", "observations", ...
      "reference"], ...
      names, 'stable');
   for container = reshape(containers, 1, [])
      if matVariableBytes(info, container) > 5e6
         continue
      end
      loaded = load(pathname, char(container));
      variables = unique([variables(:); ...
         variablesFromValue(loaded.(char(container)))], 'stable');
   end
end

function bytes = matVariableBytes(info, name)
   %MATVARIABLEBYTES Return saved byte size for one MAT variable header row.
   hit = find(string({info.name}) == name, 1);
   if isempty(hit) || ~isfield(info(hit), 'bytes')
      bytes = Inf;
   else
      bytes = double(info(hit).bytes);
   end
end

function variables = variablesFromValue(value)
   %VARIABLESFROMVALUE Return canonical channel names from container metadata.
   variables = strings(0, 1);
   if istimetable(value) || istable(value)
      variables = canonicalizeVariables(string(value.Properties.VariableNames));
   elseif isstruct(value)
      fields = string(fieldnames(value));
      variables = canonicalizeVariables(fields);
      for k = 1:numel(value)
         for field = reshape(fields, 1, [])
            child = value(k).(char(field));
            if istable(child) || istimetable(child) || isstruct(child)
               variables = unique([variables(:); variablesFromValue(child)], ...
                  'stable');
            end
         end
      end
   end
end

function variables = declaredVariables(case_manifest, source, kind)
   %DECLAREDVARIABLES Return manifest-declared comparison variables for hints.
   variables = strings(0, 1);
   if ismember(kind, ["evaluation", "observation", "metadata", "userdata"])
      variables = canonicalizeVariables( ...
         fieldStrings(case_manifest, 'comparison_variables'));
   end
   if kind == "met"
      variables = unique([variables(:); "tair"; "swd"; "lwd"; "ppt"], ...
         'stable');
   end
   if source == "retmip_protocol"
      variables = unique([variables(:); "tsfc"; "melt"; "snowf_subl"], ...
         'stable');
   end
end

function variables = canonicalizeVariables(names)
   %CANONICALIZEVARIABLES Map source/header names to canonical comparison ids.
   names = reshape(string(names), 1, []);
   variables = strings(numel(names), 1);
   n_variables = 0;
   for name = names
      mapped = canonicalName(name);
      if mapped ~= ""
         n_variables = n_variables + 1;
         variables(n_variables, 1) = mapped;
      end
   end
   variables = variables(1:n_variables);
   variables = unique(variables, 'stable');
end

function mapped = canonicalName(name)
   %CANONICALNAME Return the comparison variable represented by one name.
   token = lower(regexprep(string(name), '[\s\-]', '_'));
   if ~isempty(regexp(token, "^tice\d+$", "once"))
      mapped = "subsurface_temperature";
      return
   end
   soil_tokens = regexp(token, "^soil_temp_(\d+)_c$", "tokens", "once");
   if ~isempty(soil_tokens)
      mapped = "soil_temp_" + string(soil_tokens{1}) + "_C";
      return
   end

   switch token
      case {"smb", "surface_mass_balance"}
         mapped = "smb";
      case {"density", "rho"}
         mapped = "density";
      case {"subsurface_temperature", "t_firn", "t_ice", ...
            "tice10m", "tice"}
         mapped = "subsurface_temperature";
      case {"lwc", "slwc", "liquid_water"}
         mapped = "lwc";
      case {"compaction", "compaction_rate"}
         mapped = "compaction_rate";
      case {"tsfc", "tsurf", "surface_temperature", "surface_temp_c"}
         mapped = "tsfc";
      case {"snow_depth", "snow_depth_m", "snd", "snowdepth"}
         mapped = "snow_depth";
      case {"swe", "swe_kg_m2", "snw"}
         mapped = "swe";
      case {"melt"}
         mapped = "melt";
      case {"snowf_subl"}
         mapped = "snowf_subl";
      case {"snow_liquid_water_storage_m", "bottom_outflow_mps"}
         mapped = token;
      case {"tair", "swd", "lwd", "ppt", "rainf", "snowf"}
         mapped = token;
      otherwise
         mapped = "";
   end
end

function pairs = derivePairs(case_manifest, artifacts)
   %DERIVEPAIRS Return all artifact pairs with at least one common variable.
   max_pairs = numel(artifacts) * max(numel(artifacts) - 1, 0) / 2;
   pairs = repmat(emptyPair(), max_pairs, 1);
   n_pairs = 0;
   for i = 1:numel(artifacts)
      for j = i + 1:numel(artifacts)
         common = intersect(artifacts(i).variables, artifacts(j).variables, ...
            'stable');
         if isempty(common)
            continue
         end
         n_pairs = n_pairs + 1;
         pairs(n_pairs, 1) = pairRecord(case_manifest, artifacts(i), ...
            artifacts(j), common);
      end
   end
   pairs = pairs(1:n_pairs);
end

function pair = emptyPair()
   %EMPTYPAIR Prototype for one possible comparison pair.
   pair = struct( ...
      'case_id', "", ...
      'dataset_family', "", ...
      'class', "", ...
      'source_a', "", ...
      'source_b', "", ...
      'kind_a', "", ...
      'kind_b', "", ...
      'variables', strings(1, 0), ...
      'path_a', "", ...
      'path_b', "");
end

function pair = pairRecord(case_manifest, a, b, variables)
   %PAIRRECORD Build one comparison-pair record.
   pair = emptyPair();
   pair.case_id = string(case_manifest.case_id);
   pair.dataset_family = string(case_manifest.dataset_family);
   pair.class = pairClass(a.role, b.role);
   pair.source_a = a.source;
   pair.source_b = b.source;
   pair.kind_a = a.kind;
   pair.kind_b = b.kind;
   pair.variables = reshape(variables, 1, []);
   pair.path_a = a.path;
   pair.path_b = b.path;
end

function class = pairClass(role_a, role_b)
   %PAIRCLASS Classify a pair as model-data, model-model, or data-data.
   roles = sort([string(role_a), string(role_b)]);
   if isequal(roles, ["data", "model"])
      class = "model-data";
   elseif all(roles == "model")
      class = "model-model";
   else
      class = "data-data";
   end
end

function value = fieldStrings(s, fieldname)
   %FIELDSTRINGS Return a struct field as a column string array.
   if isfield(s, fieldname)
      value = reshape(string(s.(fieldname)), [], 1);
   else
      value = strings(0, 1);
   end
end
