function candidate = resolveCandidateBundle(manifest, kwargs)
   %RESOLVECANDIDATEBUNDLE Resolve the comparison bundle for one case.
   %
   %  candidate = ...
   %     icemodel.verification.helpers.resolveCandidateBundle(manifest)
   %  candidate = ...
   %     icemodel.verification.helpers.resolveCandidateBundle(manifest, ...
   %     candidate=bundle)
   %
   % Inputs
   %  manifest         Resolved case manifest from loadmanifest.
   %  candidate        Optional in-memory candidate artifact.
   %  candidate_file   Optional MAT file containing `candidate` or `reference`.
   %
   % Outputs
   %  candidate   Candidate bundle used by comparecase and plotcase.
   %
   % Role
   %  Operational helper that defines candidate precedence. It never mutates
   %  staged setup artifacts.

   arguments
      manifest
      kwargs.candidate = []
      kwargs.candidate_file (1, 1) string = ""
   end

   % Caller-supplied in-memory candidates take precedence for model-development
   % workflows that have already assembled outputs.
   if ~isempty(kwargs.candidate)
      candidate = kwargs.candidate;
      return
   end

   % Candidate files are for actual model/synthetic outputs. With no supplied
   % candidate file, the staged reference is loaded below for smoke comparisons.
   if ~isblanktext(kwargs.candidate_file)
      data = load(kwargs.candidate_file, "candidate");
      if ~isfield(data, "candidate")
         error(['candidate file must contain variable "candidate" ' ...
            'with the verification artifact contract'])
      end
      candidate = data.candidate;
      return
   end

   % With no supplied model output, compare against a staged smoke reference.
   % Firn cases bundle no reference.mat, so pick a declared staged model source
   % and reconstitute its per-year userdata files. RACMO remains preferred when
   % present to preserve the old smoke-reference behavior, but MAR/MERRA-only
   % cases still resolve a usable default candidate.
   if isfield(manifest, 'reference_path') ...
         && strlength(string(manifest.reference_path)) > 0 ...
         && isfile(manifest.reference_path)
      candidate = icemodel.verification.helpers.loadArtifact( ...
         manifest.reference_path, "reference");
   else
      source = defaultCandidateSource(manifest);
      if source == ""
         candidate = emptyCandidate("no staged model candidate source");
      else
         candidate = icemodel.verification.helpers.loadColocatedData( ...
            manifest, source);
      end
   end
end

function source = defaultCandidateSource(manifest)
   %DEFAULTCANDIDATESOURCE Pick a runnable model source from manifest metadata.
   source = "";
   if isfield(manifest, 'eval_sources')
      eval_sources = reshape(string(manifest.eval_sources), [], 1);
   else
      eval_sources = strings(0, 1);
   end

   preferred = icemodel.verification.namelists.rcmProductIds( ...
      ["racmo", "mar", "merra"]);
   for candidate = reshape(preferred, 1, [])
      if any(eval_sources == candidate) && hasStagedData(manifest, candidate)
         source = candidate;
         return
      end
   end

   legacy = ["racmo", "mar", "merra"];
   for candidate = reshape(legacy, 1, [])
      if hasStagedData(manifest, candidate)
         source = candidate;
         return
      end
   end
end

function tf = hasStagedData(manifest, source)
   %HASSTAGEDDATA True when a colocation leg has userdata or model output.
   tf = false;
   if ~isfield(manifest, 'colocation') || ~isstruct(manifest.colocation)
      return
   end
   source = string(source);
   fields = string(fieldnames(manifest.colocation));
   field = "";
   if any(fields == source)
      field = source;
   else
      for f = reshape(fields, 1, [])
         if sourceLabel(manifest.colocation, f) == source
            field = f;
            break
         end
      end
   end
   if field == ""
      return
   end
   leg = manifest.colocation.(char(field));
   has_data = isfield(leg, 'data_files') && ~isempty(leg.data_files);
   has_model_output = isfield(leg, 'model_output_files') ...
      && ~isempty(leg.model_output_files);
   tf = isstruct(leg) && isfield(leg, 'staged') && logical(leg.staged) ...
      && (has_data || has_model_output);
end

function label = sourceLabel(colocation, source)
   %SOURCELABEL Return the manifest source id for a colocation field.
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

function candidate = emptyCandidate(note)
   %EMPTYCANDIDATE Return an explicit empty timeseries candidate bundle.
   candidate = struct('format', 'timeseries', 'data', timetable.empty, ...
      'metadata', struct('note', char(note)));
end
