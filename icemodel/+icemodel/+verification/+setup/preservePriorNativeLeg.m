function [colocation, identity_conflict] = preservePriorNativeLeg( ...
      colocation, prior_case, metadata, build_native_forcing, kwargs)
   %PRESERVEPRIORNATIVELEG Merge prior native artifacts into a fresh case leg.
   %
   %  [colocation, identity_conflict] = ...
   %     icemodel.verification.setup.preservePriorNativeLeg(colocation, ...
   %     prior_case, metadata, build_native_forcing, family="imau", ...
   %     overwrite_family=false)
   %
   % Role
   %  Shared observation-refresh rule for native-AWS families: when a case's
   %  observations are rebuilt without rebuilding its native forcing leg, the
   %  prior manifest's forcing-owned references (met/data files, readiness,
   %  windows, artifact provenance) remain valid only while the upstream
   %  product identity is unchanged. This helper verifies that identity
   %  (product kind, DOI fields, producer family/station, and station point)
   %  and either merges the prior native references into the fresh colocation
   %  or reports a conflict so the caller can mark the leg for an explicit
   %  native rebuild.
   %
   % Name-value
   %  family : string. Family leg name inside the colocation struct (e.g.
   %     "imau", "ktransect").
   %  overwrite_family : logical. Family replacement discards prior legs, so
   %     preservation is skipped.
   %  obs_source_id_field : string. observation_variables field naming the
   %     producer station for the prior case ("hourly_site" for the IMAU
   %     hourly network; "station_site" for the K-transect annual series).
   %
   % Returns
   %  colocation : struct. Input colocation with the prior native-forcing
   %     references merged in when the identity check passed.
   %  identity_conflict : logical. True when any identity signal changed; the
   %     caller owns the family-specific conflict annotation.
   %
   % See also: icemodel.verification.setup.importImau,
   %  icemodel.verification.setup.importKtransect

   arguments
      colocation (1, 1) struct
      prior_case struct
      metadata (1, 1) struct
      build_native_forcing (1, 1) logical
      kwargs.family (1, 1) string
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.obs_source_id_field (1, 1) string = "hourly_site"
   end

   family = char(kwargs.family);
   identity_conflict = false;
   if build_native_forcing || kwargs.overwrite_family || isempty(prior_case)
      return
   end
   if ~isfield(prior_case, 'colocation') ...
         || ~isfield(prior_case.colocation, family)
      return
   end
   prior = prior_case.colocation.(family);

   % Native artifacts remain reusable only across the same stable upstream
   % product. Cache paths and source filenames may change during an observation
   % refresh, so identity rests on the product kind and DOI fields.
   identity_fields = ["kind", "doi", "bundle_doi"];
   for fieldname = identity_fields
      name = char(fieldname);
      if isfield(prior, name) && isfield(colocation.(family), name) ...
            && string(prior.(name)) ~= string(colocation.(family).(name))
         identity_conflict = true;
         return
      end
   end

   % Compare canonical and native producer names through the shared scalar
   % identity rules. Downstream source_association describes a related product
   % and is not the producer identity for these met/Data artifacts.
   fresh_identity = struct('family', metadata.source_family, ...
      'source_id', metadata.station);
   prior_metadata = struct();
   if isfield(prior, 'artifact_metadata') && ~isempty(prior.artifact_metadata)
      if ~isstruct(prior.artifact_metadata) ...
            || ~isscalar(prior.artifact_metadata)
         identity_conflict = true;
         return
      end
      prior_metadata = prior.artifact_metadata;
      if ~icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
            prior_metadata, fresh_identity)
         identity_conflict = true;
         return
      end
      native_identity = struct();
      if isfield(prior_metadata, 'source_family')
         native_identity.family = prior_metadata.source_family;
      end
      if isfield(prior_metadata, 'station')
         native_identity.source_id = prior_metadata.station;
      end
      if ~icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
            native_identity, fresh_identity)
         identity_conflict = true;
         return
      end
   end
   if isfield(prior_case, 'observation_variables') ...
         && isstruct(prior_case.observation_variables)
      prior_obs = prior_case.observation_variables;
      obs_identity = struct();
      if isfield(prior_obs, 'source_family')
         obs_identity.family = prior_obs.source_family;
      end
      if isfield(prior_obs, char(kwargs.obs_source_id_field))
         obs_identity.source_id = prior_obs.(char(kwargs.obs_source_id_field));
      end
      if ~icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
            obs_identity, fresh_identity)
         identity_conflict = true;
         return
      end
   end

   % Point semantics remain local: every known prior case/artifact point must
   % agree with the fresh producer, while missing or nonfinite legacy points are
   % compatible. Accept direct or nested artifact metadata location records.
   fresh_point = [metadata.site_location.lat_wgs84, ...
      metadata.site_location.lon_wgs84];
   prior_points = nan(3, 2);
   n_points = 0;
   if isfield(prior_case, 'site_location') ...
         && isstruct(prior_case.site_location) ...
         && all(isfield(prior_case.site_location, ...
         ["lat_wgs84", "lon_wgs84"]))
      n_points = n_points + 1;
      prior_points(n_points, :) = [ ...
         prior_case.site_location.lat_wgs84, ...
         prior_case.site_location.lon_wgs84];
   end
   if all(isfield(prior_metadata, ["lat_wgs84", "lon_wgs84"]))
      n_points = n_points + 1;
      prior_points(n_points, :) = [prior_metadata.lat_wgs84, ...
         prior_metadata.lon_wgs84];
   end
   if isfield(prior_metadata, 'site_location') ...
         && isstruct(prior_metadata.site_location) ...
         && all(isfield(prior_metadata.site_location, ...
         ["lat_wgs84", "lon_wgs84"]))
      n_points = n_points + 1;
      prior_points(n_points, :) = [ ...
         prior_metadata.site_location.lat_wgs84, ...
         prior_metadata.site_location.lon_wgs84];
   end
   for n = 1:n_points
      if all(isfinite(prior_points(n, :))) && all(isfinite(fresh_point)) ...
            && any(abs(prior_points(n, :) - fresh_point) > 1e-8)
         identity_conflict = true;
         return
      end
   end

   % Retain only forcing-owned runtime references, status, and artifact
   % provenance. Fresh source/window/evaluation fields intentionally win.
   native_fields = ["met_files", "met_file_identities", "data_files", ...
      "forcing_ready", ...
      "forcing_ready_reason", "forcing_complete_windows", ...
      "artifact_metadata"];
   for fieldname = native_fields
      name = char(fieldname);
      if isfield(prior, name)
         colocation.(family).(name) = prior.(name);
      end
   end
end
