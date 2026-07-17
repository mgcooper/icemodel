function bundle = loadColocatedData(manifest, source, kwargs)
   %LOADCOLOCATEDDATA Assemble a timeseries bundle from staged per-source files.
   %
   %  bundle = icemodel.verification.helpers.loadColocatedData(manifest, "promice")
   %  bundle = icemodel.verification.helpers.loadColocatedData(manifest, "racmo2.3p3")
   %
   % Inputs
   %  manifest         Resolved firn case manifest. Its `colocation` field
   %                   records which source legs were staged and the per-year
   %                   Data filenames under data/input/userdata.
   %  source           Source id to load. Product ids such as "mar3.11" and
   %                   "racmo2.3p3" resolve to their colocation legs. "promice"
   %                   loads the eval target Data, used only for legacy fixtures
   %                   staged before the bundled observations.mat eval contract;
   %                   freshly staged cases bundle observations.mat instead.
   %  input_data_root  Optional base input-data root. When blank, the standard
   %                   chain (inputDataRoot) resolves it.
   %
   % Outputs
   %  bundle   Struct with the verification target/reference contract:
   %             .format   "timeseries"
   %             .data     vertically concatenated per-year Data timetable
   %             .metadata provenance struct
   %
   % Role
   %  Operational helper for the firn lane. The forcing/reference side is never
   %  bundled: colocation is recorded as metadata (available sources + per-leg
   %  windows) pointing at individual per-source userdata files. This helper
   %  reconstitutes a comparison bundle on demand from those files, so
   %  comparecase/plotcase keep their timeseries contract for the RACMO
   %  reference (and the PROMICE-obs target on legacy fixtures) without a
   %  committed data bundle.
   %
   % See also: icemodel.verification.comparecase,
   %  icemodel.verification.helpers.resolveCandidateBundle,
   %  icemodel.verification.setup.importPromiceSites

   arguments
      manifest (1, 1) struct
      source (1, 1) string
      kwargs.input_data_root (1, 1) string = ""
   end

   % Resolve userdata beside the manifest's own eval tree when no input root is
   % explicit. This keeps committed demo manifests and top-level research
   % manifests self-consistent even though the global default root is now
   % top-level data/.
   ud_dir = fullfile(resolveInputRoot(manifest, kwargs.input_data_root), ...
      'userdata');
   [source_field, source_label] = resolveColocationSource(manifest, source);

   if source_field == ""
      bundle = emptyBundle(sprintf('source "%s" not in manifest colocation', ...
         source));
      return
   end
   leg = manifest.colocation.(char(source_field));
   model_output = loadModelOutput(leg, ud_dir);
   if (~isfield(leg, 'data_files') || isempty(leg.data_files)) ...
         && isempty(model_output)
      reason = 'leg not staged';
      if isfield(leg, 'reason') && strlength(string(leg.reason)) > 0
         reason = char(string(leg.reason));
      end
      bundle = emptyBundle(reason);
      return
   end

   files = strings(0, 1);
   if isfield(leg, 'data_files') && ~isempty(leg.data_files)
      files = string(leg.data_files);
   end
   files = files(strlength(files) > 0);
   data = timetable.empty;
   for f = reshape(files, 1, [])
      pathname = fullfile(ud_dir, f);
      if ~isfile(pathname)
         continue
      end
      loaded = load(pathname, 'Data');
      if ~isfield(loaded, 'Data') || ~istimetable(loaded.Data)
         continue
      end
      if isempty(data)
         data = loaded.Data;
      else
         data = outerjoinrows(data, loaded.Data);
      end
   end

   if isempty(data) && isempty(model_output)
      bundle = emptyBundle('no per-year Data files resolved on disk');
      return
   end

   if ~isempty(data)
      data = subsetToLegWindow(data, leg);
   end
   if isempty(data) && isempty(model_output)
      bundle = emptyBundle('no samples resolved inside the staged leg window');
      return
   end

   if isempty(model_output)
      bundle = struct('format', 'timeseries', 'data', data, ...
         'metadata', struct('source', char(source_label), ...
         'reconstituted_from', 'per-year userdata files'));
   else
      % Profile model output remains the primary bundle. Exact non-profile
      % channels stay available as one-variable timetables for shared soft
      % comparisons such as SUMup SMB.
      bundle = model_output;
      if ~isempty(data) && isstruct(bundle.data)
         for name = reshape(string(data.Properties.VariableNames), 1, [])
            if ~isfield(bundle.data, name)
               bundle.data.(char(name)) = data(:, name);
            end
         end
      end
      bundle.metadata.source = char(source_label);
      bundle.metadata.reconstituted_from = ...
         'model-output sidecar plus per-year userdata files';
   end
end

%% Local helpers
function bundle = loadModelOutput(leg, userdata_dir)
   %LOADMODELOUTPUT Load the first valid optional profile-model sidecar.
   bundle = [];
   if ~isfield(leg, 'model_output_files') || isempty(leg.model_output_files)
      return
   end
   for file = reshape(string(leg.model_output_files), 1, [])
      pathname = fullfile(userdata_dir, file);
      if ~isfile(pathname)
         continue
      end
      info = whos('-file', pathname);
      if ~any(string({info.name}) == "reference")
         continue
      end
      candidate = icemodel.verification.helpers.loadArtifact(pathname, "reference");
      if isstruct(candidate) && isfield(candidate, 'format') ...
            && isfield(candidate, 'data')
         bundle = candidate;
         return
      end
   end
end

function data = subsetToLegWindow(data, leg)
   %SUBSETTOLEGWINDOW Restrict a broad cached artifact to its manifest window.
   if ~isfield(leg, 'window') || ~isstruct(leg.window)
      return
   end

   [t1, t2] = icemodel.verification.setup.periodBounds(leg.window);
   if isnat(t1) || isnat(t2)
      return
   end

   row_times = data.Properties.RowTimes;
   if ~isa(row_times, 'datetime')
      return
   end

   % Manifest timestamps are stored as UTC strings. Match the timetable's
   % row-time zone before comparing so fixtures with unzoned row times and
   % production UTC artifacts both filter through the same path.
   if isempty(row_times.TimeZone)
      t1.TimeZone = '';
      t2.TimeZone = '';
   elseif isempty(t1.TimeZone)
      t1.TimeZone = row_times.TimeZone;
      t2.TimeZone = row_times.TimeZone;
   end

   data = data(row_times >= t1 & row_times <= t2, :);
end

function [source_field, source_label] = resolveColocationSource(manifest, source)
   %RESOLVECOLOCATIONSOURCE Map public source ids to manifest colocation fields.
   source = string(source);
   source_field = "";
   source_label = source;
   if ~isfield(manifest, 'colocation') || ~isstruct(manifest.colocation)
      return
   end

   fields = string(fieldnames(manifest.colocation));
   if any(fields == source)
      source_field = source;
      source_label = manifestSourceLabel(manifest.colocation, source);
      return
   end

   for f = reshape(fields, 1, [])
      label = manifestSourceLabel(manifest.colocation, f);
      if label == source
         source_field = f;
         source_label = label;
         return
      end
   end

   obs_labels = ["promice_obs", "sumup_obs", "imau_obs", ...
      "research_site_obs", "gcnet_obs"];
   obs_fields = ["promice", "sumup", "imau", "research_site", "gcnet"];
   hit = find(obs_labels == source, 1);
   if ~isempty(hit) && any(fields == obs_fields(hit))
      source_field = obs_fields(hit);
      source_label = source;
   end
end

function label = manifestSourceLabel(colocation, source)
   %MANIFESTSOURCELABEL Return the public source id for a colocation field.
   source = string(source);
   name = char(source);
   if isfield(colocation, name) && isstruct(colocation.(name)) ...
         && isfield(colocation.(name), 'source_id') ...
         && strlength(string(colocation.(name).source_id)) > 0
      label = string(colocation.(name).source_id);
   elseif ismember(source, icemodel.verification.namelists.rcmsources())
      label = icemodel.verification.namelists.rcmProductIds(source);
   else
      label = source;
   end
end

function input_root = resolveInputRoot(manifest, explicit_root)
   %RESOLVEINPUTROOT Derive input/ from an explicit root or manifest location.
   if ~isblanktext(explicit_root)
      input_root = explicit_root;
      return
   end

   if isfield(manifest, 'input_data_root') ...
         && strlength(string(manifest.input_data_root)) > 0
      input_root = string(manifest.input_data_root);
      return
   end

   if isfield(manifest, 'family_root') && strlength(string(manifest.family_root)) > 0
      eval_root = fileparts(string(manifest.family_root));
      data_root = fileparts(eval_root);
      input_root = fullfile(data_root, 'input');
      return
   end

   if isfield(manifest, 'manifest_path') && strlength(string(manifest.manifest_path)) > 0
      family_root = fileparts(string(manifest.manifest_path));
      eval_root = fileparts(family_root);
      data_root = fileparts(eval_root);
      input_root = fullfile(data_root, 'input');
      return
   end

   % Fall back to the setup default for ad hoc structs without manifest paths.
   input_root = icemodel.verification.helpers.inputDataRoot();
end

function bundle = emptyBundle(note)
   %EMPTYBUNDLE Explicit empty timeseries bundle with a reason note.
   bundle = struct('format', 'timeseries', 'data', timetable.empty, ...
      'metadata', struct('note', char(note)));
end

function combined = outerjoinrows(a, b)
   %OUTERJOINROWS Concatenate two Data timetables on a shared variable set.
   %
   % Per-year files share the same variable schema, so a row concatenation on
   % the union of times reconstructs the multi-year record. Reorder b's columns
   % to a's order and drop columns absent from a so vertcat succeeds.
   common = intersect(a.Properties.VariableNames, ...
      b.Properties.VariableNames, 'stable');
   combined = [a(:, common); b(:, common)];
   combined = sortrows(combined);
end
