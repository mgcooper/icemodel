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

   % With no supplied model output, compare against the staged smoke reference.
   candidate = icemodel.verification.helpers.loadArtifact( ...
      manifest.reference_path, "reference");
end
