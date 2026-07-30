function write_artifacts = prepareCaseRoot( ...
      case_root, overwrite, artifacts, requested_case)
   %PREPARECASEROOT Create one case folder and select artifacts needing writes.
   %
   %  write_artifacts = icemodel.verification.setup.prepareCaseRoot( ...
   %     case_root, overwrite, artifacts, requested_case)
   %
   % Inputs
   %  case_root   Folder where one case's setup artifacts are staged.
   %  overwrite   When false, existing requested artifacts are additive no-ops.
   %              When true, every requested artifact may be replaced and any
   %              existing replacement is reported with a warning.
   %  artifacts   Relative filenames requested by the caller. Omit only for
   %              legacy folder-level checks.
   %  requested_case  Optional struct with the requested .period and, when
   %              applicable, .site_location and .artifact_metadata. An
   %              existing fixed-name artifact is current only when the prior
   %              manifest case covers that period, matches the requested
   %              point, and has no concrete scalar provenance conflict.
   %
   % Outputs
   %  write_artifacts   Logical mask: true for missing or explicitly replaced
   %                    artifacts, false for current non-overwrite targets.
   %
   % Role
   %  Setup helper used only by import/update tooling. Repeated ordinary calls
   %  are idempotent; explicit overwrite remains the setup-refresh control.

   % Optional inputs preserve the legacy folder-level call while letting
   % artifact-aware callers provide durable currentness evidence.
   if nargin < 3
      artifacts = strings(1, 0);
   end
   if nargin < 4
      requested_case = struct();
   end
   artifacts = reshape(string(artifacts), 1, []);

   % A new site can always be added without a separate mkdir branch.
   icemodel.helpers.ensureDirExists(case_root);
   if ~isempty(artifacts)
      paths = fullfile(case_root, artifacts);
      existing = isfile(paths);
      current = existing;
      % A changed requested case contract makes a fixed-name artifact stale;
      % narrower covered windows remain reusable without replacement.
      if any(existing) && ~overwrite && ~isempty(fieldnames(requested_case))
         case_current = priorCaseCovers(case_root, requested_case);
         for k = find(existing)
            % Read only the already-staged target metadata. Importers build
            % requested metadata before this helper, and dry runs never enter
            % this fixed-artifact staging path.
            current(k) = case_current && fixedArtifactIdentityMatches( ...
               paths(k), requested_case);
         end
      end
      replace = existing & (overwrite | ~current);
      % Make requested replacement visible without warning on first writes or
      % on ordinary additive calls.
      if any(replace)
         warning('icemodel:verification:prepareCaseRoot:overwrite', ...
            'Replacing %d existing staged artifact(s) under %s.', ...
            nnz(replace), case_root);
      end
      write_artifacts = ~existing | replace;
      return
   end

   % Legacy callers without explicit filenames get one folder-level decision:
   % an empty folder needs a write, while an existing non-overwrite folder is a
   % no-op instead of the former hard error.
   files = dir(fullfile(case_root, '*'));
   names = string({files.name});
   names = names(~ismember(names, [".", ".."]));
   % Legacy folder-level replacement receives the same visible signal.
   if overwrite && ~isempty(names)
      warning('icemodel:verification:prepareCaseRoot:overwrite', ...
         'Replacing existing staged artifacts under %s.', case_root);
   end
   write_artifacts = overwrite || isempty(names);
end

function tf = fixedArtifactIdentityMatches(filename, requested_case)
   %FIXEDARTIFACTIDENTITYMATCHES Compare saved and requested scalar provenance.
   expected = requested_case;
   if isfield(requested_case, 'artifact_metadata')
      expected = requested_case.artifact_metadata;
   end
   if ~isstruct(expected) || ~isscalar(expected)
      tf = false;
      return
   end

   % Verification observation bundles save one `targets` struct. Loading that
   % staged variable may read its small payload but never reopens a raw source;
   % metadata-free legacy artifacts remain unknown-compatible, while unreadable
   % or malformed present targets require repair.
   candidate = struct();
   try
      saved = load(filename, 'targets');
   catch
      tf = false;
      return
   end
   if isfield(saved, 'targets') ...
         && (~isstruct(saved.targets) || ~isscalar(saved.targets))
      tf = false;
      return
   end
   if isfield(saved, 'targets') && isfield(saved.targets, 'metadata')
      if ~isstruct(saved.targets.metadata) ...
            || ~isscalar(saved.targets.metadata)
         tf = false;
         return
      end
      candidate = saved.targets.metadata;
   end
   tf = icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
      candidate, expected);
end

function tf = priorCaseCovers(case_root, requested_case)
   %PRIORCASECOVERS True when the prior manifest case covers this request.
   tf = false;
   family_root = fileparts(case_root);
   manifest_file = fullfile(family_root, 'manifest.json');
   % An artifact without a prior manifest case cannot be assumed current.
   if ~isfile(manifest_file)
      return
   end
   manifest = jsondecode(fileread(manifest_file));
   % Skip-only or empty manifests carry no reusable case contract.
   if ~isfield(manifest, 'cases') || isempty(manifest.cases)
      return
   end
   [~, case_id] = fileparts(case_root);
   ids = string({manifest.cases.case_id});
   hit = find(ids == string(case_id), 1);
   % A manually present artifact outside the manifest is stale for this call.
   if isempty(hit)
      return
   end
   prior_case = manifest.cases(hit);
   tf = periodCovers(prior_case, requested_case) ...
      && locationMatches(prior_case, requested_case);
end

function tf = periodCovers(prior_case, requested_case)
   %PERIODCOVERS True when a bounded prior period contains the requested one.
   tf = false;
   % Both sides need an explicit period to prove enclosing coverage.
   if ~isfield(prior_case, 'period') || ~isfield(requested_case, 'period')
      return
   end
   [prior_start, prior_end] = ...
      icemodel.verification.setup.periodBounds(prior_case.period);
   [request_start, request_end] = ...
      icemodel.verification.setup.periodBounds(requested_case.period);
   % periodBounds returns unzoned NaT for blank bounds and UTC for concrete
   % bounds; normalize all four scalars before concatenation or comparison.
   prior_start.TimeZone = 'UTC';
   prior_end.TimeZone = 'UTC';
   request_start.TimeZone = 'UTC';
   request_end.TimeZone = 'UTC';
   % Unbounded periods require a rebuild because containment is unknowable.
   if any(isnat([prior_start, prior_end, request_start, request_end]))
      return
   end
   tf = prior_start <= request_start && prior_end >= request_end;
end

function tf = locationMatches(prior_case, requested_case)
   %LOCATIONMATCHES Apply point compatibility only when one is requested.
   tf = true;
   % Fixed-site families need only period coverage and omit this field.
   if ~isfield(requested_case, 'site_location')
      return
   end
   requested = requested_case.site_location;
   % Nonfinite placeholder locations do not constrain artifact reuse.
   if ~hasFinitePoint(requested)
      return
   end
   % A finite requested point requires equivalent prior location metadata.
   if ~isfield(prior_case, 'site_location') ...
         || ~hasFinitePoint(prior_case.site_location)
      tf = false;
      return
   end
   prior = prior_case.site_location;
   tf = abs(double(prior.lat_wgs84) - double(requested.lat_wgs84)) <= 1e-8 ...
      && abs(double(prior.lon_wgs84) - double(requested.lon_wgs84)) <= 1e-8;
end

function tf = hasFinitePoint(location)
   %HASFINITEPOINT True for a location with finite WGS84 latitude/longitude.
   tf = isstruct(location) && isfield(location, 'lat_wgs84') ...
      && isfield(location, 'lon_wgs84') ...
      && all(isfinite([double(location.lat_wgs84), ...
      double(location.lon_wgs84)]));
end
