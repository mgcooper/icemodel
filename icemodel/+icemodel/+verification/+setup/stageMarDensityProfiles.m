function state = stageMarDensityProfiles(state, alive_idx, source, kwargs)
   %STAGEMARDENSITYPROFILES Add optional MAR RO1 profiles to SUMup cases.
   %
   %  state = icemodel.verification.setup.stageMarDensityProfiles( ...
   %     state, alive_idx, source, family_root=..., userdata_outdir=..., ...)
   %
   % The helper runs after one MAR colocation leg has been merged but before the
   % family manifest is persisted. It selects only SUMup density observation
   % dates, writes one additive model-output sidecar per case, and records the
   % sidecar on the MAR leg. Absence or failed refresh never changes forcing
   % readiness and never removes an existing valid sidecar.

   arguments
      state (1, :) struct
      alive_idx (1, :) double {mustBeInteger, mustBePositive}
      source (1, 1) string
      kwargs.family_root (1, 1) string
      kwargs.userdata_outdir (1, 1) string
      kwargs.mar_dir (1, 1) string
      kwargs.overwrite (1, 1) logical = false
   end

   % Other RCM sources have no MAR RO1 product and remain untouched.
   if source ~= "mar" || isempty(alive_idx)
      return
   end

   % Resolve the same default MAR root as the primary builder; yearly source
   % files share one fixed grid.
   mar_dir = kwargs.mar_dir;
   if mar_dir == ""
      mar_dir = string(getenv("ICEMODEL_MAR_DIR"));
      if mar_dir == ""
         mar_dir = "/Volumes/S03/DATA/greenland/mar3p11/RUH2";
      end
   end
   source_files = dir(fullfile(mar_dir, '*-*.nc'));
   if isempty(source_files)
      % A cached primary MAR leg can survive without the native cache. Record
      % the optional-product outcome explicitly while preserving any earlier
      % sidecar reference and the independent forcing-ready state.
      for idx = reshape(alive_idx, 1, [])
         if ~isfield(state(idx).colocation, 'mar') ...
               || ~isstruct(state(idx).colocation.mar) ...
               || ~isfield(state(idx).colocation.mar, 'staged') ...
               || ~state(idx).colocation.mar.staged
            continue
         end
         leg = state(idx).colocation.mar;
         leg = recordModelOutputFailure(leg, 'profile_not_available', ...
            'Native MAR RO1 source files are unavailable; primary forcing data are unchanged.');
         state(idx).colocation.mar = leg;
      end
      return
   end
   grid_file = string(fullfile(source_files(1).folder, source_files(1).name));
   % Only source cell identity is needed here; avoid reading unrelated static
   % fields required by the forcing builder's full grid contract.
   grid = struct('LON', double(ncread(grid_file, 'LON')), ...
      'LAT', double(ncread(grid_file, 'LAT')));
   output_dir = fullfile(kwargs.userdata_outdir, "mar3.11");
   if ~isfolder(output_dir)
      mkdir(output_dir)
   end

   for idx = reshape(alive_idx, 1, [])
      % Optional profile staging requires a staged nearest-cell MAR leg and a
      % native SUMup density observation table with explicit dates.
      if ~isfield(state(idx).colocation, 'mar')
         continue
      end
      leg = state(idx).colocation.mar;
      if ~isstruct(leg) || ~isfield(leg, 'staged') || ~leg.staged
         continue
      end
      method = string(fieldOr(leg, 'sample_method', "nearest"));
      observation_file = fullfile( ...
         kwargs.family_root, state(idx).evaluation_file_rel);
      if method ~= "nearest" || ~isfile(observation_file)
         leg = recordModelOutputFailure(leg, 'profile_not_available', ...
            'MAR density profiles require nearest sampling and staged SUMup observations.');
         state(idx).colocation.mar = leg;
         continue
      end
      loaded_target = load(observation_file, 'targets');
      if ~isfield(loaded_target, 'targets') ...
            || ~isfield(loaded_target.targets, 'data') ...
            || ~isfield(loaded_target.targets.data, 'density') ...
            || ~istable(loaded_target.targets.data.density) ...
            || ~ismember('datetime', ...
            loaded_target.targets.data.density.Properties.VariableNames)
         leg = recordModelOutputFailure(leg, 'profile_not_available', ...
            'SUMup density observations have no profile dates.');
         state(idx).colocation.mar = leg;
         continue
      end
      requested_dates = loaded_target.targets.data.density.datetime;
      requested_dates = requested_dates(~isnat(requested_dates));
      if isempty(requested_dates)
         leg = recordModelOutputFailure(leg, 'profile_not_available', ...
            'SUMup density observations have no finite profile dates.');
         state(idx).colocation.mar = leg;
         continue
      end

      % Recover the exact MAR cell from the staged data artifact provenance.
      [grid_index, provenance_status] = artifactGridIndex( ...
         leg, kwargs.userdata_outdir, grid);
      if provenance_status ~= "ok"
         leg = recordModelOutputFailure(leg, char(provenance_status), ...
            'MAR data artifact lacks an unambiguous nearest-cell identity.');
         state(idx).colocation.mar = leg;
         continue
      end

      % Read exact-date public RO1 snapshots. Dynamic layers stay nested QA and
      % never replace or correct the public density profile.
      [profiles, selection_status, dynamic_qa] = ...
         icemodel.forcing.helpers.readMarDensitySnapshots( ...
         requested_dates, grid_index, source_dir=mar_dir, ...
         site_id=string(state(idx).site_id), sample_method=method, ...
         requested_location=double(state(idx).point));

      relative_file = fullfile("mar3.11", ...
         string(state(idx).case_id) + "_mar3.11_density_profiles.mat");
      output_file = fullfile(kwargs.userdata_outdir, relative_file);
      if isempty(profiles)
         % A failed/empty refresh is additive: keep an earlier valid sidecar and
         % record the new outcome without deleting model-output references.
         if isfile(output_file)
            status = 'preserved_after_empty_refresh';
         else
            statuses = unique(string(selection_status.status));
            status = char(strjoin(statuses, ','));
         end
         note = char(strjoin(unique(string(selection_status.detail)), '; '));
         leg = recordModelOutputFailure(leg, status, note);
         state(idx).colocation.mar = leg;
         continue
      end

      % Merge by explicit profile identity so surgical imports replace only
      % requested valid dates and preserve all unrequested historical groups.
      merged_profiles = profiles;
      if isfile(output_file)
         prior = load(output_file, 'reference');
         has_prior_table = isfield(prior, 'reference') ...
               && isstruct(prior.reference) ...
               && isfield(prior.reference, 'data') ...
               && isfield(prior.reference.data, 'density') ...
               && istable(prior.reference.data.density);
         if has_prior_table
            prior_profiles = prior.reference.data.density;
            compatible = isequal( ...
               prior_profiles.Properties.VariableNames, ...
               profiles.Properties.VariableNames) ...
               && isequal(varfun(@class, prior_profiles, ...
               'OutputFormat', 'cell'), varfun(@class, profiles, ...
               'OutputFormat', 'cell'));
            if compatible
               replace_ids = unique(string(profiles.profile_id));
               keep = ~ismember(string(prior_profiles.profile_id), replace_ids);
               merged_profiles = [prior_profiles(keep, :); profiles];
            elseif ~kwargs.overwrite
               leg = recordModelOutputFailure(leg, ...
                  'incompatible_existing_sidecar', ...
                  'Existing MAR profile sidecar was preserved because its payload is incompatible.');
               state(idx).colocation.mar = leg;
               continue
            end
         elseif ~kwargs.overwrite
            leg = recordModelOutputFailure(leg, ...
               'incompatible_existing_sidecar', ...
               'Existing MAR profile sidecar was preserved because its payload is incompatible.');
            state(idx).colocation.mar = leg;
            continue
         end
      end
      merged_profiles = sortrows(merged_profiles, ...
         {'datetime', 'profile_id', 'depth'});

      % Validate the public artifact before the atomic write.
      groups = icemodel.verification.helpers.profileGroups(merged_profiles);
      for group = reshape(groups, 1, [])
         depth = double(group.data.depth);
         density = double(group.data.density);
         assert(all(diff(depth) > 0) && min(depth) >= 0 && max(depth) <= 20, ...
            'MAR RO1 profile depth support is invalid.')
         assert(all(isfinite(density) & density >= 250 & density <= 1000), ...
            'MAR RO1 density is outside the accepted source-backed range.')
      end

      reference = struct('format', 'subsurface_profile_bundle', ...
         'data', struct('format', 'subsurface_profile_bundle', ...
         'density', merged_profiles), ...
         'metadata', struct('source_id', 'mar3.11', ...
         'product_name', 'MAR modelled snow/firn/ice density', ...
         'public_variable', 'RO1', ...
         'profile_grouping', 'profile_id + datetime', ...
         'selection_status', selection_status, 'dynamic_qa', dynamic_qa));
      artifact_metadata = reference.metadata;
      temporary_file = string(tempname(output_dir)) + ".mat";
      save(temporary_file, 'reference', 'artifact_metadata')
      movefile(temporary_file, output_file, 'f')

      % Model output is discoverable without entering forcing data_files or
      % changing the advisory forcing_ready contract.
      leg.model_output_files = char(relative_file);
      leg.model_output_format = 'subsurface_profile_bundle';
      leg.model_output_variables = {'density'};
      leg.model_output_status = 'staged';
      leg.model_output_note = sprintf( ...
         '%d MAR RO1 profile group(s), exact SUMup UTC dates.', numel(groups));
      state(idx).colocation.mar = leg;
   end
end

function leg = recordModelOutputFailure(leg, status, note)
   %RECORDMODELOUTPUTFAILURE Record failure only when no sidecar is inherited.
   % A forcing-only replay cannot rebuild every optional profile sidecar. Keep
   % an inherited reference and all of its metadata as one atomic artifact.
   has_reference = (isfield(leg, 'model_output_file') ...
      && ~isempty(leg.model_output_file)) ...
      || (isfield(leg, 'model_output_files') ...
      && ~isempty(leg.model_output_files));
   if has_reference
      return
   end

   % Without an existing artifact, expose the optional failure diagnostically.
   leg.model_output_status = status;
   leg.model_output_note = note;
end

function value = fieldOr(s, name, fallback)
   %FIELDOR Return one nonblank struct field or its fallback value.
   value = fallback;
   if isfield(s, name) && ~isempty(s.(name)) ...
         && strlength(string(s.(name))) > 0
      value = s.(name);
   end
end

function [grid_index, status] = artifactGridIndex(leg, userdata_outdir, grid)
   %ARTIFACTGRIDINDEX Recover exact MAR [X Y] identity from saved provenance.
   grid_index = [NaN NaN];
   status = "missing_grid_provenance";
   if ~isfield(leg, 'data_files') || isempty(leg.data_files)
      return
   end
   files = string(leg.data_files);
   for file = reshape(files, 1, [])
      pathname = fullfile(userdata_outdir, file);
      if ~isfile(pathname)
         continue
      end
      artifact = load(pathname, 'artifact_metadata');
      if ~isfield(artifact, 'artifact_metadata') ...
            || ~isstruct(artifact.artifact_metadata)
         continue
      end
      metadata = artifact.artifact_metadata;
      if isfield(metadata, 'grid_start') && isfield(metadata, 'grid_count') ...
            && isequal(double(metadata.grid_count(:)'), [1 1])
         grid_index = double(metadata.grid_start(:)');
         status = "ok";
         return
      end
      if ~isfield(metadata, 'source_lat_wgs84') ...
            || ~isfield(metadata, 'source_lon_wgs84')
         continue
      end
      distance = sind((grid.LAT - metadata.source_lat_wgs84) / 2).^2 ...
         + cosd(metadata.source_lat_wgs84) .* cosd(grid.LAT) ...
         .* sind((grid.LON - metadata.source_lon_wgs84) / 2).^2;
      [minimum, linear_index] = min(distance(:));
      if ~isfinite(minimum)
         continue
      end
      [grid_index(1), grid_index(2)] = ind2sub(size(grid.LON), linear_index);
      status = "ok";
      return
   end
end
