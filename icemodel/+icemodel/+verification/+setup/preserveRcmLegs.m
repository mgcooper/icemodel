function incoming = preserveRcmLegs( ...
      prior, incoming, sources, period, kwargs)
   %PRESERVERCMLEGS Keep compatible staged RCM legs after a failed refresh.
   %
   %  incoming = icemodel.verification.setup.preserveRcmLegs( ...
   %     prior, incoming, sources, period, met_outdir=..., ...
   %     userdata_outdir=..., method="nearest", point=[lat lon])
   %
   % For each explicitly refreshed source, remove a skipped/degraded incoming
   % leg when the prior staged leg still resolves on disk and overlaps the case
   % period. mergeColocation then leaves that prior leg unchanged. An exact
   % overwrite=false reuse also keeps the prior leg so a diagnostic reuse note
   % cannot churn otherwise identical manifest bytes. A successful overwrite or
   % materially changed request replaces the prior leg, while a failed/degraded
   % overwrite preserves a compatible prior fallback. Concrete metadata conflicts
   % remain visible replacement decisions. Missing metadata alone remains
   % compatible only through an already-referenced prior leg; directory discovery
   % never attaches it.

   arguments
      prior (1, 1) struct
      incoming (1, 1) struct
      sources (1, :) string
      period (1, 1) struct
      kwargs.met_outdir (1, 1) string = ""
      kwargs.userdata_outdir (1, 1) string = ""
      kwargs.method (1, 1) string = "nearest"
      kwargs.point (1, 2) double = [NaN, NaN]
      kwargs.overwrite (1, 1) logical = false
   end

   % Preserve only sources present in both graphs; omitted sources are already
   % additive and successful incoming legs must replace their prior versions.
   for source = reshape(sources, 1, [])
      name = char(source);
      if ~isfield(prior, name) || ~isfield(incoming, name)
         continue
      end
      old = prior.(name);
      new = incoming.(name);
      [preserve, replace_prior_artifacts] = ...
         shouldPreserve(old, new, period, kwargs);
      metadata_conflict = isMetadataConflict(new);
      explicit_replacement = isfield(new, 'replace_prior_artifacts') ...
         && isscalar(new.replace_prior_artifacts) ...
         && logical(new.replace_prior_artifacts);
      if preserve && ~metadata_conflict && ~explicit_replacement
         incoming = rmfield(incoming, name);
      elseif replace_prior_artifacts || metadata_conflict ...
            || explicit_replacement
         % Tell the additive manifest merge that this requested leg is a
         % validated replacement, not an ordinary failed patch. The merge
         % consumes this transient flag before writing the manifest.
         new.replace_prior_artifacts = true;
         incoming.(name) = new;
      end
   end
end

function [tf, replace_prior_artifacts] = ...
      shouldPreserve(old, new, period, kwargs)
   %SHOULDPRESERVE Classify failed refreshes and exact no-op cache reuse.
   replace_prior_artifacts = false;
   old_staged = isstruct(old) && isfield(old, 'staged') ...
      && logical(old.staged);
   new_leg = isstruct(new) && isfield(new, 'staged');
   failed_refresh = new_leg ...
      && (~logical(new.staged) || dropsPriorData(old, new));
   exact_reuse = new_leg && logical(new.staged) && ~kwargs.overwrite ...
      && sameSuccessfulReuse(old, new);
   candidate = old_staged && (failed_refresh || exact_reuse);
   if ~candidate
      tf = false;
      return
   end
   [tf, replace_prior_artifacts] = ...
      priorLegCompatible(old, period, kwargs);
end

function tf = sameSuccessfulReuse(old, new)
   %SAMESUCCESSFULREUSE Detect one byte-logical no-op staged-leg replay.

   % Preserve provenance only when the staging identity, exact artifact
   % references, and resolved coverage are unchanged. Any wider/different path,
   % source product, method, or window remains a material incoming patch.
   identity_fields = ["kind", "source", "source_id", "sample_method"];
   tf = true;
   for field = identity_fields
      tf = tf && sameTextField(old, new, field);
   end
   tf = tf && sameReferenceField(old, new, "met_files") ...
      && sameReferenceField(old, new, "data_files") ...
      && sameCoverageWindow(old, new) ...
      && sameForcingReadiness(old, new);
end

function tf = sameForcingReadiness(old, new)
   %SAMEFORCINGREADINESS Compare exact-artifact diagnostics across JSON types.
   fields = ["forcing_ready", "forcing_ready_reason", ...
      "forcing_complete_windows"];
   for field = fields
      name = char(field);
      if isfield(old, name) ~= isfield(new, name)
         % A newly diagnosed or removed readiness field is a material patch.
         tf = false;
         return
      end
   end
   if ~isfield(old, 'forcing_ready')
      % RACMO and legacy non-met legs legitimately carry no readiness group.
      tf = true;
      return
   end

   % Logical and text scalars can change representation after JSON decoding;
   % compare their manifest meaning rather than their MATLAB storage class.
   old_ready = old.forcing_ready;
   new_ready = new.forcing_ready;
   tf = isscalar(old_ready) && isscalar(new_ready) ...
      && logical(old_ready) == logical(new_ready) ...
      && sameTextField(old, new, "forcing_ready_reason");
   if ~tf
      return
   end

   % Empty windows round-trip as numeric [] while populated windows round-trip
   % from strings to chars. Normalize both valid portable schemas explicitly.
   old_windows = old.forcing_complete_windows;
   new_windows = new.forcing_complete_windows;
   if isempty(old_windows) || isempty(new_windows)
      tf = isempty(old_windows) && isempty(new_windows);
      return
   end
   required = {'start_time', 'end_time', 'sample_count'};
   if ~isstruct(old_windows) || ~isstruct(new_windows) ...
         || ~all(isfield(old_windows, required)) ...
         || ~all(isfield(new_windows, required)) ...
         || numel(old_windows) ~= numel(new_windows)
      tf = false;
      return
   end
   tf = isequal(reshape(string({old_windows.start_time}), [], 1), ...
      reshape(string({new_windows.start_time}), [], 1)) ...
      && isequal(reshape(string({old_windows.end_time}), [], 1), ...
      reshape(string({new_windows.end_time}), [], 1)) ...
      && isequal(reshape(double([old_windows.sample_count]), [], 1), ...
      reshape(double([new_windows.sample_count]), [], 1));
end

function tf = sameTextField(old, new, field)
   %SAMETEXTFIELD Require one identical optional scalar text field.
   name = char(field);
   old_has = isfield(old, name);
   new_has = isfield(new, name);
   tf = old_has == new_has;
   if tf && old_has
      old_value = string(old.(name));
      new_value = string(new.(name));
      tf = isscalar(old_value) && isscalar(new_value) ...
         && isequal(old_value, new_value);
   end
end

function tf = sameReferenceField(old, new, field)
   %SAMEREFERENCEFIELD Require the same exact ordered artifact references.
   name = char(field);
   old_has = isfield(old, name);
   new_has = isfield(new, name);
   tf = old_has == new_has;
   if tf && old_has
      tf = isequal(reshape(string(old.(name)), [], 1), ...
         reshape(string(new.(name)), [], 1));
   end
end

function tf = sameCoverageWindow(old, new)
   %SAMECOVERAGEWINDOW Compare normalized resolved staging coverage.
   tf = isfield(old, 'window') && isstruct(old.window) ...
      && isfield(new, 'window') && isstruct(new.window);
   if ~tf
      return
   end
   try
      [old_start, old_end] = ...
         icemodel.verification.setup.periodBounds(old.window);
      [new_start, new_end] = ...
         icemodel.verification.setup.periodBounds(new.window);
      tf = isequaln(old_start, new_start) && isequaln(old_end, new_end);
   catch
      % Malformed coverage is a material patch and must remain visible.
      tf = false;
   end
end

function [tf, replace_prior_artifacts] = ...
      priorLegCompatible(leg, period, kwargs)
   %PRIORLEGCOMPATIBLE Classify prior files for preservation or replacement.
   replace_prior_artifacts = false;
   files = priorLegFiles(leg, kwargs);
   if isempty(files) || any(~isfile(files))
      % A requested refresh has proven that the prior runtime references are
      % unusable. They must not be resurrected by the additive merge.
      tf = false;
      replace_prior_artifacts = true;
      return
   end
   if ~isfield(leg, 'window') || ~isstruct(leg.window)
      tf = false;
      return
   end

   % A declared prior method must agree with this refresh. Legacy legs without
   % the field remain eligible because their manifest reference is provenance.
   if isfield(leg, 'sample_method') ...
         && strlength(string(leg.sample_method)) > 0 ...
         && string(leg.sample_method) ~= kwargs.method
      tf = false;
      replace_prior_artifacts = true;
      return
   end

   % Referenced legacy files may sit outside normal alias-window discovery.
   % Honor missing metadata as unknown, but reject concrete saved method/point
   % evidence that disagrees with this refresh.
   if ~priorFilesMatchMetadata(files, kwargs.method, kwargs.point)
      tf = false;
      replace_prior_artifacts = true;
      return
   end

   % A bounded overlap is required; unbounded/stale metadata must be rebuilt.
   [case_start, case_end] = ...
      icemodel.verification.setup.periodBounds(period);
   [leg_start, leg_end] = ...
      icemodel.verification.setup.periodBounds(leg.window);
   tf = ~isnat(case_start) && ~isnat(case_end) ...
      && ~isnat(leg_start) && ~isnat(leg_end) ...
      && leg_start <= case_end && leg_end >= case_start;
end

function tf = priorFilesMatchMetadata(files, method, point)
   %PRIORFILESMATCHMETADATA Reject concrete conflicts in referenced artifacts.
   tf = true;
   check_point = all(isfinite(point));
   for file = reshape(files, 1, [])
      metadata = ...
         icemodel.verification.setup.readRcmArtifactMetadata(file);
      if isempty(fieldnames(metadata))
         continue
      end

      % A populated saved method is concrete provenance; blanks remain unknown.
      if isfield(metadata, 'sample_method') ...
            && strlength(string(metadata.sample_method)) > 0 ...
            && string(metadata.sample_method) ~= method
         tf = false;
         return
      end
      if ~check_point
         continue
      end

      % Current metadata carries both coordinates. A partial pair is malformed,
      % while a wholly absent pair remains legacy unknown provenance.
      has_lat = isfield(metadata, 'lat_wgs84');
      has_lon = isfield(metadata, 'lon_wgs84');
      if has_lat ~= has_lon
         tf = false;
         return
      end
      if has_lat
         saved_point = [double(metadata.lat_wgs84), ...
            double(metadata.lon_wgs84)];
         if ~all(isfinite(saved_point)) ...
               || any(abs(saved_point - point) > 1e-8)
            tf = false;
            return
         end
      end
   end
end

function files = priorLegFiles(leg, kwargs)
   %PRIORLEGFILES Resolve prior manifest met/data refs to full paths.
   files = strings(0, 1);
   if isfield(leg, 'met_files') && ~isempty(leg.met_files)
      files = [files; fullfile(kwargs.met_outdir, ...
         reshape(string(leg.met_files), [], 1))];
   end
   if isfield(leg, 'data_files') && ~isempty(leg.data_files)
      files = [files; fullfile(kwargs.userdata_outdir, ...
         reshape(string(leg.data_files), [], 1))];
   end
end

function tf = dropsPriorData(old, new)
   %DROPSPRIORDATA True when a refresh degrades met+Data to met-only.
   tf = isfield(old, 'data_files') && ~isempty(old.data_files) ...
      && (~isfield(new, 'data_files') || isempty(new.data_files));
end

function tf = isMetadataConflict(leg)
   %ISMETADATACONFLICT True when skipped cache metadata is incompatible.
   tf = isfield(leg, 'reason') ...
      && contains(string(leg.reason), "artifact metadata does not match");
end
