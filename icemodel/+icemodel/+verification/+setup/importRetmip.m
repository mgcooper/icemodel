function manifest = importRetmip(source_dir, kwargs)
   %IMPORTRETMIP Build the RetMIP family manifest.
   %
   %  manifest = icemodel.verification.setup.importRetmip(source_dir)
   %  manifest = icemodel.verification.setup.importRetmip(source_dir, ...
   %     dry_run=true)
   %
   %  Role
   %    RetMIP staging hook. It records RetMIP protocol cases as evaluation
   %    userdata and stages confirmed native meteorological sources separately
   %    under the standard icemodel input/met and input/userdata layout. Optional
   %    MAR/MERRA/RACMO legs are delegated to the shared dataset-family RCM
   %    staging helper after protocol/native products are safely persisted.
   %    forcing_sources selects runtime sources requested by the current call.
   %    Ordinary calls preserve omitted existing legs; overwrite_family=true
   %    deliberately replaces the whole family state.
   %    build_observations=false is a guarded non-dry fast path: requested cases
   %    must already exist in the target manifest, whose protocol/native entry is
   %    reused while selected forcing is attached.
   %
   %  Default roots
   %    source_dir="" reads <repo>/data/verification/retmip. With no output_root,
   %    evaluation artifacts go to <repo>/data/eval/retmip/<case_id>/ and native
   %    met/userdata go to <repo>/data/input/{met,userdata}/retmip/. Explicit
   %    source_dir, output_root, evaluation_data_root, and input_data_root
   %    overrides are honored as-is.
   %
   %  startdate/enddate optionally clamp each requested RetMIP protocol period.
   %  This is intended for short preview staging; omit both for full production
   %  artifacts.
   %  Model met defaults to dt_out="15m"; pass dt_out="" for native cadence.
   %  Data/userdata defaults to hourly at the shared writer boundary.
   %  Native met schema completion is fixed at the importer boundary: absent
   %  required channels are retained as NaN placeholders.
   %  Call the source-family met builder directly for strict validation.
   %
   %  Name-value
   %    case_ids : string vector  RetMIP protocol cases; blank selects all.
   %    forcing_sources : string vector  Native/RCM legs requested by this call.
   %    startdate, enddate : paired optional protocol/forcing clamp.
   %    output_root : string  Convenience root for eval and input artifacts.
   %    build_observations : logical  Build or reuse protocol/native entries.
   %    build_forcing : logical  Stage requested runtime legs (default false).
   %    overwrite : logical  Refresh requested case artifacts (default false).
   %    overwrite_family : logical  Replace the whole family (default false).
   %    skip_missing : logical  Record source-local skips instead of failing.
   %    dry_run : logical  Return metadata without source or output writes.
   %
   %  Returns
   %    manifest : struct  Final or dry-run family manifest.
   %
   %  See also: icemodel.verification.setup.importPromiceSites,
   %    icemodel.verification.setup.importImau

   arguments
      source_dir (1, 1) string = ""
      kwargs.case_ids (1, :) string = strings(1, 0)
      kwargs.forcing_sources (1, :) string ...
         {icemodel.verification.validators.mustBeRcmSourceSelection( ...
         kwargs.forcing_sources, "retmip")} = "retmip"
      kwargs.startdate = ""
      kwargs.enddate = ""
      kwargs.promice_dir (1, 1) string = ""
      kwargs.gcnet_dir (1, 1) string = ""
      kwargs.samimi_dir (1, 1) string = ""
      kwargs.imau_dir (1, 1) string = ""
      kwargs.mar_dir (1, 1) string = ""
      kwargs.merra_dir (1, 1) string = ""
      kwargs.racmo_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.output_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.input_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
      kwargs.dt_out (1, 1) string ...
         {mustBeMember(kwargs.dt_out, ["", "15m"])} = "15m"
      kwargs.overwrite (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.skip_missing (1, 1) logical = true
      kwargs.dry_run (1, 1) logical = false
      kwargs.build_forcing (1, 1) logical = false
      kwargs.build_observations (1, 1) logical = true
   end

   forcing_sources = ...
      icemodel.verification.setup.normalizeForcingSources( ...
      kwargs.forcing_sources, kwargs.build_forcing);
   kwargs.forcing_sources = forcing_sources;

   % Validate the optional clamp before any cache or staging side effect.
   [window_start, window_end, window_enabled] = ...
      icemodel.internal.pairedWindow( ...
      kwargs.startdate, kwargs.enddate);
   requested_window = struct('enabled', window_enabled, ...
      'start', window_start, 'end', window_end);

   % Resolve the family identity and requested runtime source sets once.
   dataset_family = "retmip";
   build_native_forcing = kwargs.build_forcing ...
      && ismember(dataset_family, forcing_sources);
   rcm_sources = intersect(forcing_sources, ...
      icemodel.verification.namelists.rcmsources(), "stable");
   build_rcm_forcing = kwargs.build_forcing && ~isempty(rcm_sources);

   cases = icemodel.verification.setup.retmipCaseCatalog(kwargs.case_ids);
   requested_ids = string({cases.case_id});

   % Resolve output roots and paths before raw sources. Forcing-only calls can
   % reuse the existing manifest without requiring observation/native caches.
   [evaluation_data_root, input_root] = ...
      icemodel.verification.setup.resolveStagingRoots( ...
      output_root=kwargs.output_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      input_data_root=kwargs.input_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
   [family_root, manifest_file, met_outdir, userdata_outdir] = ...
      icemodel.verification.setup.datasetFamilyStagingPaths( ...
      evaluation_data_root, input_root, dataset_family);

   prior_cases = struct([]);
   coverage = struct();
   reuse_sources = strings(1, 0);

   if ~kwargs.dry_run
      prior_cases = ...
         icemodel.verification.setup.loadPriorDatasetFamilyCases( ...
         manifest_file, overwrite_family=kwargs.overwrite_family, ...
         build_observations=kwargs.build_observations);
   end

   if ~kwargs.dry_run && kwargs.build_observations
      icemodel.helpers.ensureDirExists(family_root);
   end
   if ~kwargs.dry_run && kwargs.build_forcing
      icemodel.helpers.ensureDirExists(met_outdir);
      icemodel.helpers.ensureDirExists(userdata_outdir);
   end

   % Resolve RCM coverage only for a real requested build.
   if ~kwargs.dry_run && build_rcm_forcing
      reuse_sources = rcm_sources;
      coverage = icemodel.verification.setup.promiceSourceCoverage( ...
         rcm_sources, struct('mar', kwargs.mar_dir, ...
         'merra', kwargs.merra_dir, 'racmo', kwargs.racmo_dir));
   end

   if ~kwargs.dry_run && ~kwargs.build_observations ...
         && ~build_native_forcing
      % Reuse the staged case entry so an RCM-only attachment does not require
      % source caches.
      [state, alive, skipped] = ...
         icemodel.verification.setup.reuseDatasetFamilyCases( ...
         manifest_file, requested_ids, emptyState(), ...
         dataset_family=dataset_family, ...
         overwrite_family=kwargs.overwrite_family, ...
         forcing_sources=reuse_sources, coverage=coverage, ...
         startdate=kwargs.startdate, enddate=kwargs.enddate);
      for k = 1:numel(state)
         % The forcing-only path must retain the same collision-safe cache name
         % as a full RetMIP import without changing the manifest case id.
         state(k).storage_alias = ...
            icemodel.verification.setup.rcmStorageAlias( ...
            dataset_family, state(k).case_id);
      end
   else
      % Validate caches only when building observations or native runtime files.
      % Dry runs remain metadata-only; optional skips stay quiet while required
      % RetMIP products print their retrieval guidance before failing.
      cache_status = struct();
      if ~kwargs.dry_run && kwargs.build_observations
         strict_cache = ~kwargs.skip_missing;
         [source_dir, cache_status] = icemodel.verification.setup.fetchRetmip( ...
            cache_dir=icemodel.forcing.helpers.verificationSourceDir( ...
            source_dir, "retmip"), ...
            products=["forcing", "outputs"], strict=strict_cache, ...
            silent=kwargs.skip_missing, ...
            create_cache_dir=true);
      end
      if ~kwargs.dry_run && ~kwargs.skip_missing
         preflightCases( ...
            cases, source_dir, kwargs, requested_window, ...
            build_native_forcing);
      end

      stage_callback = @(~, n) stageRetmipCase( ...
         cases(n), source_dir, cache_status, family_root, input_root, kwargs, ...
         coverage, rcm_sources, requested_window, ...
         icemodel.verification.setup.priorCaseById( ...
         prior_cases, cases(n).case_id), build_native_forcing);

      % Stage each requested case into importer state.
      [state, alive, skipped] = ...
         icemodel.verification.setup.stageDatasetFamilyCases( ...
         1:numel(requested_ids), emptyState(), stage_callback, ...
         skip_missing=kwargs.skip_missing, ...
         warning_id="icemodel:verification:importRetmip:caseSkipped", ...
         label_callback=@(~, n) requested_ids(n));
   end

   % Record source provenance once so dry-run and persisted manifests match.
   source_doi = "10.22008/FK2/GZ3CSN;10.22008/FK2/CVPUJL";
   source_url = "https://doi.org/10.22008/FK2/GZ3CSN";
   source_version = "retmip-protocol+outputs";
   retrieval_date = string(datetime('today'));

   entry_callback = @icemodel.verification.setup.stateCaseEntry;

   % Return metadata without writing case or manifest artifacts.
   if kwargs.dry_run
      manifest = icemodel.verification.setup.runDatasetFamilyDryRun( ...
         state, alive, dataset_family=dataset_family, ...
         requested_ids=requested_ids, skipped=skipped, ...
         source_doi=source_doi, source_url=source_url, ...
         source_version=source_version, retrieval_date=retrieval_date, ...
         entry_callback=entry_callback);
      return
   end

   % Persist case entries first, then attach requested RCM source legs.
   [manifest, ~] = icemodel.verification.setup.runDatasetFamilyImport( ...
      state, alive, dataset_family=dataset_family, ...
      manifest_file=manifest_file, requested_ids=requested_ids, ...
      skipped=skipped, ...
      source_doi=source_doi, source_url=source_url, ...
      source_version=source_version, retrieval_date=retrieval_date, ...
      overwrite_family=kwargs.overwrite_family, overwrite=kwargs.overwrite, ...
      entry_callback=entry_callback, ...
      build_forcing=build_rcm_forcing, forcing_sources=rcm_sources, ...
      leg_callback=@(s, src) s.leg.(char(src)), ...
      met_outdir=met_outdir, userdata_outdir=userdata_outdir, ...
      mar_dir=kwargs.mar_dir, merra_dir=kwargs.merra_dir, ...
      racmo_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir, ...
      method="nearest", dt_out=kwargs.dt_out);
end

function entry = emptyEntry()
   %EMPTYENTRY Prototype RetMIP manifest entry.
   case_values = { ...
      ''
      'firn_observational'
      ''
      ''
      'unknown'
      {'firn'}
      'none'
      struct('lat_wgs84', NaN, 'lon_wgs84', NaN, ...
      'x_epsg3413', NaN, 'y_epsg3413', NaN, 'elev_m', NaN)
      struct('start', '', 'end', '')
      ''
      {}
      {'retmip_protocol'}
      {'tsfc', 'melt', 'snowf_subl', 'density', 'subsurface_temperature', 'lwc'}
      struct()
      struct()
      '3hr'
      ''};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(case_values);
end

function s = emptyState()
   %EMPTYSTATE Prototype RetMIP staging state.
   s = struct('case_id', "", 'storage_alias', "", ...
      'point', [NaN NaN], ...
      'evaluation_file_rel', "", ...
      'entry', emptyEntry(), 'period', struct('start', '', 'end', ''), ...
      'colocation', struct(), 'leg', struct(), 'reuse_entry', false, ...
      'dry_run', false);
end

function s = stageRetmipCase(c, source_dir, cache_status, family_root, ...
      input_root, kwargs, coverage, rcm_sources, requested_window, prior_case, ...
      build_native_forcing)
   %STAGERETMIPCASE Stage one RetMIP case and return importer state.
   c = clampCasePeriod(c, requested_window);
   entry = manifestEntry(c, source_dir, cache_status, family_root, input_root, ...
      kwargs.promice_dir, kwargs.gcnet_dir, kwargs.samimi_dir, ...
      kwargs.imau_dir, kwargs.dry_run, kwargs.build_observations, ...
      kwargs.build_forcing, build_native_forcing, kwargs.overwrite_family, ...
      kwargs.skip_missing, kwargs.overwrite, ...
      kwargs.dt_out, prior_case);
   leg = struct();
   if ~kwargs.dry_run && kwargs.build_forcing && ~isempty(rcm_sources)
      [t1, t2] = icemodel.verification.setup.periodBounds(c.period);
      leg = icemodel.verification.setup.resolveLegWindows( ...
         rcm_sources, coverage, t1, t2);
   end
   s = struct('case_id', string(c.case_id), 'storage_alias', ...
      icemodel.verification.setup.rcmStorageAlias("retmip", c.case_id), ...
      'point', [c.site_location.lat_wgs84, c.site_location.lon_wgs84], ...
      'evaluation_file_rel', string(entry.evaluation_file), ...
      'entry', entry, 'period', entry.period, ...
      'colocation', entry.colocation, 'leg', leg, 'reuse_entry', false, ...
      'dry_run', kwargs.dry_run);
end

function c = clampCasePeriod(c, requested_window)
   %CLAMPCASEPERIOD Intersect a RetMIP case period with the requested window.
   c.period_is_clamped = false;
   if ~requested_window.enabled
      return
   end

   [t1, t2] = icemodel.verification.setup.periodBounds(c.period);
   t1 = max(t1, requested_window.start);
   t2 = min(t2, requested_window.end);
   if t2 < t1
      error('icemodel:verification:importRetmip:emptyWindow', ...
         'requested RetMIP window does not overlap %s', c.case_id)
   end
   c.period = icemodel.verification.setup.manifestWindow(t1, t2);
   c.period_is_clamped = true;
end

function preflightCases(cases, source_dir, kwargs, requested_window, ...
      build_native_forcing)
   %PREFLIGHTCASES Validate every strict RetMIP case before staging begins.
   % Parse all protocol bundles first so a later malformed case cannot follow an
   % earlier case-root or manifest mutation.
   if kwargs.build_observations
      for k = 1:numel(cases)
         c = clampCasePeriod(cases(k), requested_window);
         files = retmipFiles(source_dir, c);
         prepareProtocolBundle(c, files);
      end
   end

   % Native meteorological caches are runtime forcing inputs, not protocol
   % observation prerequisites, so validate them only after every bundle parses.
   if ~build_native_forcing
      return
   end

   for k = 1:numel(cases)
      assoc = cases(k).source_association;
      family = string(assoc.family);
      if ~ismember(family, ["promice", "gcnet_vandecrux", "samimi", "imau"])
         continue
      end
      [available, resolved_dir, reason] = nativeSourceStatus( ...
         family, kwargs.promice_dir, ...
         kwargs.gcnet_dir, kwargs.samimi_dir, kwargs.imau_dir, ...
         string(assoc.source_id));
      if ~available
         throwMissingNative(reason, family);
      end
      c = clampCasePeriod(cases(k), requested_window);
      buildNativeMet(family, string(assoc.source_id), resolved_dir, c.period);
   end
end

function throwMissingNative(reason, family)
   %THROWMISSINGNATIVE Raise the RetMIP missing-native-source error.
   if family == "gcnet_vandecrux"
      error('icemodel:verification:importRetmip:missingGcnetVandecrux', ...
         '%s', reason)
   end
   error('icemodel:verification:importRetmip:missingNativeSource', ...
      '%s', reason)
end

function entry = manifestEntry(c, source_dir, cache_status, family_root, ...
      input_root, promice_dir, gcnet_dir, samimi_dir, imau_dir, dry_run, ...
      build_observations, build_forcing, build_native_forcing, ...
      overwrite_family, skip_missing, overwrite, dt_out, prior_case)
   %MANIFESTENTRY Build one RetMIP firn case entry.
   if ~dry_run && ~build_observations
      if isempty(prior_case)
         error('icemodel:verification:importRetmip:missingPriorCase', ...
            ['build_observations=false requires case "%s" in the target ' ...
            'RetMIP manifest.'], c.case_id);
      end
      entry = prior_case;
      if overwrite_family
         entry = ...
            icemodel.verification.setup.prepareReplacementCaseEntry( ...
            entry, "retmip");
      end
      [entry.colocation, ~, native_obs, native_timestep] = stageNativeSource( ...
         entry.colocation, c, input_root, promice_dir, gcnet_dir, samimi_dir, ...
         imau_dir, false, build_forcing, build_native_forcing, skip_missing, ...
         overwrite, dt_out, prior_case);
      if ~isempty(fieldnames(native_obs))
         entry.observation_variables.native_source = native_obs;
      end
      entry.native_timestep = char(native_timestep);
      return
   end

   % Dry runs use the catalog contract without discovering protocol files.
   files = struct('surface', "", 'profiles', strings(0, 1), ...
      'outputs', strings(0, 1));
   if ~dry_run
      files = retmipFiles(source_dir, c);
   end
   evaluation_file = "";
   if ~dry_run
      % Parse and stamp the complete bundle before prepareCaseRoot can create,
      % clear, or replace any case artifact.
      targets = prepareProtocolBundle(c, files);
      requested_case = struct('period', c.period, ...
         'site_location', c.site_location, ...
         'artifact_metadata', protocolArtifactIdentity(c));
      write_observation = icemodel.verification.setup.prepareCaseRoot( ...
         fullfile(family_root, c.case_id), overwrite, "observations.mat", ...
         requested_case);
      evaluation_file = fullfile(c.case_id, "observations.mat");
      if write_observation
         evaluation_file = saveProtocolBundle(family_root, c, targets);
      end
   end

   retmip = struct();
   retmip.kind = 'protocol_userdata';
   retmip.staged = ~dry_run;
   retmip.source_dir = char(source_dir);
   retmip.protocol_timestep = char(c.protocol_timestep);
   retmip.surface_file = char(files.surface);
   retmip.profile_files = cellstr(files.profiles(:));
   retmip.model_output_files = cellstr(files.outputs(:));
   retmip.model_output_variables = ...
      cellstr(reshape(modelOutputVariables(files.outputs), [], 1));
   retmip.cache_status = cache_status;

   colocation = struct();
   colocation.retmip = retmip;
   colocation.source_association = c.source_association;
   [colocation, ~, native_obs, native_timestep] = stageNativeSource( ...
      colocation, c, input_root, promice_dir, gcnet_dir, samimi_dir, ...
      imau_dir, dry_run, build_forcing, build_native_forcing, skip_missing, ...
      overwrite, dt_out, prior_case);

   comparison_variables = comparisonVariables(files);
   if dry_run
      forcing_sources = strings(0, 1);
      eval_sources = "retmip_protocol";
   else
      [forcing_sources, eval_sources] = ...
         icemodel.verification.setup.colocationSourceLists(colocation);
   end
   observation_variables = struct('protocol_id', c.protocol_id, ...
      'retmip_station_id', c.retmip_station_id);
   if ~isempty(fieldnames(native_obs))
      observation_variables.native_source = native_obs;
   end

   case_values = { ...
      char(c.case_id)
      'firn_observational'
      char(c.site_id)
      char(c.site_name)
      char(c.surface_zone)
      cellstr(string(c.eval_target))
      char(c.permafrost_zone)
      c.site_location
      c.period
      char(evaluation_file)
      cellstr(forcing_sources)
      cellstr(eval_sources)
      cellstr(comparison_variables)
      observation_variables
      colocation
      char(native_timestep)
      caseNote(c, colocation)};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(case_values);
end

function [colocation, native_vars, native_obs, native_timestep] = stageNativeSource( ...
      colocation, c, input_root, promice_dir, gcnet_dir, samimi_dir, ...
      imau_dir, dry_run, build_forcing, build_native_forcing, skip_missing, ...
      overwrite, dt_out, prior_case)
   %STAGENATIVESOURCE Stage the confirmed native source for one RetMIP case.
   native_vars = strings(0, 1);
   native_obs = struct();
   native_timestep = string(c.protocol_timestep);
   assoc = c.source_association;
   family = string(assoc.family);
   if ~ismember(family, ["promice", "gcnet_vandecrux", "samimi", "imau"])
      return
   end

   % RetMIP uses the runtime source label `retmip`, while this provenance record
   % identifies the underlying native source that supplies it.
   station = string(assoc.source_id);
   source = struct('family', family, ...
      'source_id', station, 'relationship', assoc.relationship);
   colocation.retmip.native_source = source;

   % Dry runs publish selector intent without reading native source caches.
   if dry_run
      status = "not_requested";
      reason = "retmip was not requested in forcing_sources";
      if ~build_forcing
         status = "forcing_disabled";
         reason = "native forcing disabled because build_forcing=false";
      elseif build_native_forcing
         status = "not_checked";
         reason = "native source cache not read during dry run";
      end
      colocation = recordNativeStatus( ...
         colocation, false, status, reason, family);
      return
   end

   % A protocol-only import neither probes nor writes native runtime artifacts.
   % Ordinary refreshes retain a compatible previously staged native leg;
   % identity changes leave a fresh unstaged leg that requires a forcing rebuild.
   if ~build_native_forcing
      [colocation, native_obs, native_timestep, preserved, identity_conflict] = ...
         preservePriorNativeSource(colocation, prior_case, native_timestep);
      if preserved
         return
      end
      status = "not_requested";
      reason = "retmip was not requested in forcing_sources";
      if ~build_forcing
         status = "forcing_disabled";
         reason = "native forcing disabled because build_forcing=false";
      end
      if identity_conflict
         status = "identity_changed_requires_rebuild";
         reason = "native source identity changed; rerun with " + ...
            "build_forcing=true and forcing_sources including retmip";
      end
      colocation = recordNativeStatus( ...
         colocation, false, status, reason, family);
      return
   end

   [available, resolved_dir, reason] = nativeSourceStatus( ...
      family, promice_dir, gcnet_dir, samimi_dir, imau_dir, station, dry_run);
   colocation.retmip.native_source.cache_dir = char(resolved_dir);
   if ~available
      colocation = handleNativeMissing(colocation, reason, skip_missing, family);
      return
   end

   try
      % Build both products from one source read so the met and userdata files
      % share identical time windows and provenance.
      [met, metadata, Data] = buildNativeMet( ...
         family, station, resolved_dir, c.period);
      % Keep the protocol evaluation period intact while recording the exact
      % native forcing support used by runtime window resolution.
      native_window = icemodel.verification.setup.manifestWindow( ...
         met.Time(1), met.Time(end));
      met_outdir = fullfile(input_root, 'met');
      userdata_outdir = fullfile(input_root, 'userdata');
      met_files = icemodel.forcing.helpers.writemet(met, c.case_id, ...
         "retmip", outdir=met_outdir, naming="window", dt_out=dt_out, ...
         overwrite=overwrite);
      % Diagnose the exact returned path because no-overwrite staging may
      % select an existing exact or broader enclosing artifact.
      [forcing_ready, forcing_ready_reason, forcing_complete_windows] = ...
         icemodel.verification.setup.metArtifactReadiness(met_files);
      data_files = icemodel.forcing.helpers.writeuserdata(Data, c.case_id, ...
         "retmip", outdir=userdata_outdir, naming="window", ...
         overwrite=overwrite);

      colocation.retmip.kind = 'protocol_userdata_and_native_met';
      colocation.retmip.native_met_status = 'staged';
      colocation.retmip.window = native_window;
      colocation.retmip.met_files = ...
         icemodel.verification.setup.relpaths(met_files, met_outdir);
      colocation.retmip.data_files = ...
         icemodel.verification.setup.relpaths(data_files, userdata_outdir);
      colocation.retmip.forcing_ready = logical(forcing_ready);
      colocation.retmip.forcing_ready_reason = char(forcing_ready_reason);
      colocation.retmip.forcing_complete_windows = forcing_complete_windows;
      colocation.retmip.native_source.station = metadata.station;
      colocation.retmip.native_source.source_file = metadata.source_file;
      colocation.retmip.native_source.precip_policy = metadata.precip_policy;
      if isfield(metadata, 'lwd_policy')
         colocation.retmip.native_source.lwd_policy = metadata.lwd_policy;
      end
      colocation.native_met = struct('staged', true, 'status', "staged", ...
         'source_family', family, 'source_id', station, ...
         'forcing_ready', logical(forcing_ready), ...
         'forcing_ready_reason', string(forcing_ready_reason));
      colocation.native_met.forcing_complete_windows = ...
         forcing_complete_windows;

      native_vars = string(Data.Properties.VariableNames);
      native_obs = nativeObservationMetadata(metadata, met, Data);
      native_timestep = nativeTimestepLabel(met);

   catch err
      if ~skip_missing || ~isSkippableNativeBuildError(err)
         rethrow(err)
      end
      colocation = handleNativeMissing(colocation, err.message, true, family);
   end
end

function [colocation, native_obs, native_timestep, preserved, identity_conflict] = ...
      preservePriorNativeSource(colocation, prior_case, native_timestep)
   %PRESERVEPRIORNATIVESOURCE Retain compatible native state on obs-only refresh.
   native_obs = struct();
   identity_conflict = false;
   preserved = hasPriorNativeSource(prior_case);
   if ~preserved
      return
   end

   % Missing legacy identity remains compatible, but every concrete native
   % family/source/relationship fact must match the fresh case association.
   prior_source = struct();
   if isfield(prior_case.colocation.retmip, 'native_source') ...
         && isstruct(prior_case.colocation.retmip.native_source)
      prior_source = prior_case.colocation.retmip.native_source;
   end
   fresh_source = colocation.retmip.native_source;
   preserved = icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
      prior_source, fresh_source);
   identity_conflict = ~preserved;
   if identity_conflict
      return
   end

   % Refresh protocol-source fields while retaining every native status, file,
   % policy, and provenance field owned by the earlier forcing build.
   fresh = colocation.retmip;
   retained = prior_case.colocation.retmip;
   protocol_fields = ["staged", "source_dir", "protocol_timestep", ...
      "surface_file", "profile_files", "model_output_files", ...
      "model_output_variables", "cache_status"];
   for fieldname = protocol_fields
      name = char(fieldname);
      if isfield(fresh, name)
         retained.(name) = fresh.(name);
      end
   end
   colocation.retmip = retained;
   if isfield(prior_case.colocation, 'native_met')
      colocation.native_met = prior_case.colocation.native_met;
   end
   if isfield(prior_case, 'observation_variables') ...
         && isfield(prior_case.observation_variables, 'native_source')
      native_obs = prior_case.observation_variables.native_source;
   end
   if isfield(prior_case, 'native_timestep') ...
         && strlength(string(prior_case.native_timestep)) > 0
      native_timestep = string(prior_case.native_timestep);
   end
end

function tf = hasPriorNativeSource(prior_case)
   %HASPRIORNATIVESOURCE True when a prior RetMIP leg owns runtime file refs.
   tf = ~isempty(fieldnames(prior_case)) ...
      && isfield(prior_case, 'colocation') ...
      && isfield(prior_case.colocation, 'retmip') ...
      && isstruct(prior_case.colocation.retmip) ...
      && isfield(prior_case.colocation.retmip, 'met_files') ...
      && ~isempty(prior_case.colocation.retmip.met_files) ...
      && isfield(prior_case.colocation.retmip, 'data_files') ...
      && ~isempty(prior_case.colocation.retmip.data_files);
end

function label = nativeTimestepLabel(met)
   %NATIVETIMESTEPLABEL Return the runtime cadence label for staged native met.
   label = "hourly";
   if height(met) < 2
      return
   end
   dt = seconds(met.Time(2) - met.Time(1));
   if abs(dt - 900) < 1
      label = "15m";
   elseif abs(dt - 3600) < 1
      label = "1hr";
   elseif abs(dt - 10800) < 1
      label = "3hr";
   else
      label = string(duration(0, 0, dt, 'Format', 'hh:mm:ss'));
   end
end

function tf = isSkippableNativeBuildError(err)
   %ISSKIPPABLENATIVEBUILDERROR True only for absent native data/window errors.
   id = string(err.identifier);
   tf = endsWith(id, ":fileNotFound") || endsWith(id, ":emptyWindow");
end

function [available, source_dir, reason] = nativeSourceStatus( ...
      family, promice_dir, gcnet_dir, samimi_dir, imau_dir, station, dry_run)
   %NATIVESOURCESTATUS Return cheap cache presence for one native source.
   if nargin < 7
      dry_run = false;
   end
   switch family
      case "promice"
         [source_dir, status] = icemodel.verification.setup.fetchPromice( ...
            cache_dir=promice_dir, stations=station, products="hour", ...
            strict=false, silent=true, create_cache_dir=~dry_run);
         available = all([status.present]);
         reason = "";
         if ~available
            missing_products = string({status(~[status.present]).product});
            reason = sprintf('missing PROMICE source product for %s: %s', ...
               station, strjoin(missing_products, ', '));
         end
      case "gcnet_vandecrux"
         [source_dir, status] = icemodel.verification.setup.fetchGcnet( ...
            cache_dir=gcnet_dir, stations=station, products="surface", ...
            strict=false, silent=true, create_cache_dir=~dry_run);
         row = status(string({status.product}) == "surface");
         available = ~isempty(row) && logical(row.present);
         reason = "";
         if ~available
            missing = strings(1, 0);
            if ~isempty(row)
               missing = string(row.missing_files);
            end
            reason = sprintf( ...
               'missing GC-Net/Vandecrux surface product for %s: %s', ...
               station, strjoin(missing, ', '));
         end
      case "samimi"
         source_dir = icemodel.forcing.helpers.verificationSourceDir( ...
            samimi_dir, ["retmip", "samimi"]);
         try
            icemodel.forcing.helpers.locateSamimiDye2Workbook(source_dir);
            available = true;
         catch err
            if ~strcmp(err.identifier, ...
                  'icemodel:forcing:locateSamimiDye2Workbook:fileNotFound')
               rethrow(err)
            end
            available = false;
         end
         reason = "";
         if ~available
            reason = sprintf('missing Samimi Dye-2 AWS workbook under %s', ...
               source_dir);
         end
      case "imau"
         source_dir = icemodel.forcing.helpers.verificationSourceDir( ...
            imau_dir, "imau");
         available = imauHourlyFilePresent(source_dir, station);
         reason = "";
         if ~available
            reason = sprintf('missing IMAU hourly source product for %s under %s', ...
               station, source_dir);
         end
   end
end

function [met, metadata, Data] = buildNativeMet(family, station, source_dir, ...
      period)
   %BUILDNATIVEMET Build a native RetMIP met/Data pair by source family.
   switch family
      case "gcnet_vandecrux"
         [met, metadata, Data] = icemodel.forcing.buildGcnetVandecruxMet( ...
            station, source_dir=source_dir, ...
            startdate=period.start, enddate=period.end, ...
            fillwithmissing=true);
      case "promice"
         [met, metadata] = icemodel.forcing.buildPromiceMet( ...
            station, source_dir=source_dir, ...
            startdate=period.start, enddate=period.end, ...
            fillwithmissing=true);
         [Data, data_metadata] = icemodel.forcing.buildPromiceData( ...
            station, source_dir=source_dir, ...
            startdate=period.start, enddate=period.end);
         metadata.data_metadata = data_metadata;
         metadata.source_family = "promice";
         metadata.station = station;
         metadata.source_file = data_metadata.source_file;
      case "samimi"
         [met, metadata, Data] = icemodel.forcing.buildSamimiDye2Met( ...
            source_dir=source_dir, startdate=period.start, ...
            enddate=period.end, fillwithmissing=true);
      case "imau"
         [met, metadata, Data] = icemodel.forcing.buildImauHourlyMet( ...
            station, source_dir=source_dir, startdate=period.start, ...
            enddate=period.end, fillwithmissing=true);
   end
end

function tf = imauHourlyFilePresent(source_dir, station)
   %IMAUHOURLYFILEPRESENT Check that the requested hourly station exists.
   try
      icemodel.forcing.helpers.locateImauHourlyFile(source_dir, station);
      tf = true;
   catch err
      if strcmp(err.identifier, ...
            'icemodel:forcing:locateImauHourlyFile:fileNotFound')
         tf = false;
      else
         rethrow(err)
      end
   end
end

function colocation = handleNativeMissing(colocation, reason, skip_missing, ...
      family)
   %HANDLENATIVEMISSING Either record or raise a missing-native-source failure.
   if ~skip_missing
      if family == "gcnet_vandecrux"
         error('icemodel:verification:importRetmip:missingGcnetVandecrux', ...
            '%s', reason)
      else
         error('icemodel:verification:importRetmip:missingNativeSource', ...
            '%s', reason)
      end
   end
   colocation = recordNativeStatus(colocation, false, ...
      "missing_" + family, reason, family);
end

function colocation = recordNativeStatus(colocation, staged, status, reason, ...
      family)
   %RECORDNATIVESTATUS Store native-met status in the RetMIP colocation record.
   colocation.retmip.native_met_status = char(status);
   if strlength(string(reason)) > 0
      colocation.retmip.native_met_skipped_reason = char(reason);
   end
   source_id = "";
   if isfield(colocation.retmip, 'native_source')
      source_id = colocation.retmip.native_source.source_id;
   end
   colocation.native_met = struct('staged', logical(staged), ...
      'status', status, 'source_family', family, ...
      'source_id', source_id, 'reason', string(reason));
end

function obs = nativeObservationMetadata(metadata, met, Data)
   %NATIVEOBSERVATIONMETADATA Keep compact native-source variable provenance.
   obs = struct( ...
      'family', metadata.source_family, ...
      'source_id', metadata.station, ...
      'source_file', metadata.source_file, ...
      'met_variables', {cellstr(string(met.Properties.VariableNames))}, ...
      'data_variables', {cellstr(string(Data.Properties.VariableNames))}, ...
      'precip_policy', metadata.precip_policy);
   if isfield(metadata, 'lwd_policy')
      obs.lwd_policy = metadata.lwd_policy;
   end
end

function note = caseNote(~, colocation)
   %CASENOTE Explain the RetMIP protocol/native-source split for the manifest.
   note = "RetMIP protocol case; protocol userdata remains retmip_protocol.";
   if ~isfield(colocation.retmip, 'native_source')
      note = note + ...
         " Standard met source label is retmip only after native source staging.";
      note = char(note);
      return
   end

   status = "";
   if isfield(colocation.retmip, 'native_met_status')
      status = string(colocation.retmip.native_met_status);
   end
   family = string(colocation.retmip.native_source.family);
   if status == "staged"
      note = note + " Native met/userdata staged from " + family ...
         + " under source label retmip.";
   else
      note = note + " Native " + family + " source not staged.";
   end
   note = char(note);
end

function files = retmipFiles(source_dir, c)
   %RETMIPFILES Resolve protocol/profile/output files for one RetMIP case.
   forcing_dirs = [fullfile(source_dir, "forcing"); source_dir];
   outputs_dirs = [fullfile(source_dir, "outputs"); source_dir];
   token = filenameToken(c.case_id);
   profile_token = profileFilenameToken(c.case_id);

   % The surface protocol file is required for a non-dry import; profiles are
   % optional by variable because not every case ships every initial state axis.
   files = struct( ...
      'surface', firstMatch(forcing_dirs, ...
      "RetMIP_surface_forcing_" + token + ".tab"), ...
      'profiles', profileMatches(forcing_dirs, profile_token), ...
      'outputs', outputMatches(outputs_dirs, token));
end

function targets = prepareProtocolBundle(c, files)
   %PREPAREPROTOCOLBUNDLE Parse one RetMIP bundle without filesystem mutation.
   assertProtocolBundleAvailable(files, c);

   % Parse the official surface series and any initial state profile tables into
   % one data-only target bundle. Model outputs remain indexed as file products.
   [surface, surface_meta] = ...
      icemodel.verification.setup.readRetmipProtocolTable(files.surface);
   if isfield(c, 'period_is_clamped') && c.period_is_clamped
      surface = filterProtocolSurface(surface, surface_meta, c.period);
   end
   profiles = struct();
   profile_meta = struct();
   for k = 1:numel(files.profiles)
      [profile, meta] = ...
         icemodel.verification.setup.readRetmipProfileTable(files.profiles(k));
      name = profileName(files.profiles(k));
      profiles.(name) = profile;
      profile_meta.(name) = meta;
   end

   metadata = protocolArtifactIdentity(c);
   metadata.surface = surface_meta;
   metadata.profiles = profile_meta;
   metadata.source_kind = ...
      'RetMIP protocol userdata and comparison products';
   targets = struct( ...
      'format', 'retmip_protocol_bundle', ...
      'data', struct('surface', surface, 'profiles', profiles, ...
      'model_output_files', files.outputs(:)), ...
      'metadata', metadata);
   targets = icemodel.verification.setup.stampArtifactMetadata(targets);
end

function relpath = saveProtocolBundle(family_root, c, targets)
   %SAVEPROTOCOLBUNDLE Persist a prepared bundle after root preparation.
   % prepareCaseRoot has already authorized the write and created the case root.
   relpath = fullfile(c.case_id, "observations.mat");
   save(fullfile(family_root, relpath), 'targets');
end

function metadata = protocolArtifactIdentity(c)
   %PROTOCOLARTIFACTIDENTITY Return stable RetMIP bundle provenance identity.
   metadata = struct( ...
      'source_family', 'retmip', ...
      'station', char(string(c.retmip_station_id)), ...
      'doi', '10.22008/FK2/GZ3CSN', ...
      'product', char(string(c.protocol_id)), ...
      'schema', 'retmip_protocol_bundle');
end

function surface = filterProtocolSurface(surface, metadata, period)
   %FILTERPROTOCOLSURFACE Restrict RetMIP protocol rows to the staged period.
   [t1, t2] = icemodel.verification.setup.periodBounds(period);
   time_name = char(metadata.time_variable);
   keep = surface.(time_name) >= t1 & surface.(time_name) <= t2;
   surface = surface(keep, :);
   if isempty(surface)
      error('icemodel:verification:importRetmip:emptyProtocolWindow', ...
         'RetMIP protocol surface has no rows in requested window')
   end
end

function assertProtocolBundleAvailable(files, c)
   %ASSERTPROTOCOLBUNDLEAVAILABLE Fail before mutating a skipped case root.
   if strlength(files.surface) == 0
      error('icemodel:verification:importRetmip:missingSurfaceFile', ...
         'missing RetMIP surface protocol file for %s', c.case_id);
   end
end

function token = filenameToken(case_id)
   %FILENAMETOKEN Map filename-safe case ids to RetMIP source filename tokens.
   switch string(case_id)
      case "kanu"
         token = "KANU";
      case "dye2_long"
         token = "Dye-2_long";
      case "dye2_2016"
         token = "Dye-2_16";
      case "summit"
         token = "Summit";
      case "fa"
         token = "FA";
      otherwise
         error('icemodel:verification:importRetmip:unknownCase', ...
            'unknown RetMIP case id: %s', case_id);
   end
end

function token = profileFilenameToken(case_id)
   %PROFILEFILENAMETOKEN Map case ids to RetMIP profile filename tokens.
   token = filenameToken(case_id);
   if string(case_id) == "kanu"
      token = "KAN-U";
   end
end

function filename = firstMatch(folders, pattern)
   %FIRSTMATCH Return the first matching file or an empty string.
   filename = "";
   for folder = reshape(folders, 1, [])
      hits = dir(fullfile(folder, pattern));
      if ~isempty(hits)
         filename = string(fullfile(hits(1).folder, hits(1).name));
         return
      end
   end
end

function files = profileMatches(folders, token)
   %PROFILEMATCHES Return initial density/temperature/LWC files for one case.
   kinds = ["density", "temperature", "lwc"];
   files = strings(1, numel(kinds));
   n_files = 0;
   for k = 1:numel(kinds)
      hit = firstMatch(folders, ...
         "RetMIP_initial_firn_" + kinds(k) + "_" + token + ".tab");
      if strlength(hit) > 0
         n_files = n_files + 1;
         files(n_files) = hit;
      end
   end
   files = files(1:n_files);
end

function files = outputMatches(folders, token)
   %OUTPUTMATCHES Return RetMIP model-output NetCDF files for one case.
   hit_groups = cell(numel(folders), 1);
   n_groups = 0;
   for folder = reshape(folders, 1, [])
      found = [dir(fullfile(folder, "*_" + token + "_3hourly_*.nc")); ...
         dir(fullfile(folder, "*_" + upper(token) + "_3hourly_*.nc"))];
      if isempty(found)
         continue
      end
      n_groups = n_groups + 1;
      hit_groups{n_groups} = found(:);
   end
   if n_groups == 0
      files = strings(0, 1);
      return
   end
   hits = vertcat(hit_groups{1:n_groups});
   files = strings(1, numel(hits));
   for k = 1:numel(hits)
      files(k) = string(fullfile(hits(k).folder, hits(k).name));
   end
   files = unique(files, 'stable');
end

function variables = comparisonVariables(files)
   %COMPARISONVARIABLES Return variables present in the protocol target bundle.
   variables = ["tsfc", "melt", "snowf_subl"];
   for filename = reshape(string(files.profiles), 1, [])
      variables(end + 1) = string(profileName(filename)); %#ok<AGROW>
   end
   variables = unique(variables, 'stable');
end

function variables = modelOutputVariables(output_files)
   %MODELOUTPUTVARIABLES Read model-output headers when they are valid NetCDFs.
   if isempty(output_files)
      variables = strings(0, 1);
      return
   end

   collected = cell(1, numel(output_files));
   for k = 1:numel(output_files)
      try
         inventory = ...
            icemodel.verification.setup.retmipOutputInventory(output_files(k));
         collected{k} = inventory.variables(:);
      catch
         % Invalid or placeholder output files are still useful path evidence.
         collected{k} = strings(0, 1);
      end
   end
   nonempty = ~cellfun(@isempty, collected);
   if ~any(nonempty)
      variables = strings(0, 1);
      return
   end
   variables = vertcat(collected{nonempty});
   variables = unique(variables, 'stable');
end

function name = profileName(filename)
   %PROFILENAME Build a valid struct field for a RetMIP profile file.
   [~, stem] = fileparts(filename);
   if contains(stem, "density")
      name = 'density';
   elseif contains(stem, "temperature")
      name = 'subsurface_temperature';
   elseif contains(stem, "lwc")
      name = 'lwc';
   else
      name = matlab.lang.makeValidName(stem);
   end
end
