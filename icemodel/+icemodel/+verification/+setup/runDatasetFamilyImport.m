function [manifest, state] = runDatasetFamilyImport(state, alive, kwargs)
   %RUNDATASETFAMILYIMPORT Persist native state and optional requested RCMs.
   %
   %  [manifest, state] = icemodel.verification.setup.runDatasetFamilyImport( ...
   %     state, alive, dataset_family=..., manifest_file=..., ...
   %     entry_callback=@makeEntry)
   %
   % Neutral orchestration shared by firn-evaluation importers. Family code
   % builds native observations/state; this helper persists that state first,
   % delegates explicitly requested RCM work one source at a time, and persists
   % after each source. Cache discovery, derivation, and artifact writes belong
   % only to stageDatasetRcmForcing -> stageRcmForcing. Consequently,
   % build_forcing=false performs no RCM discovery, attachment, or input write.
   % Ordinary manifest writes remain additive through writeFamilyManifestMerge.
   % Model met defaults to dt_out="15m"; shared userdata writers default to 1h.
   % RCM selectors and nearest/natural interpolation are validated before the
   % first native-manifest persistence, keeping invalid calls non-mutating.
   %
   % See also: icemodel.verification.setup.buildDatasetFamilyManifest,
   %  icemodel.verification.setup.stageDatasetRcmForcing,
   %  icemodel.verification.setup.writeFamilyManifestMerge

   arguments
      state (1, :) struct
      alive (1, :) logical
      kwargs.dataset_family (1, 1) string
      kwargs.manifest_file (1, 1) string
      kwargs.requested_ids (1, :) string = strings(1, 0)
      kwargs.skipped (1, :) struct = struct('site', {}, 'reason', {})
      kwargs.source_doi (1, 1) string = ""
      kwargs.source_url (1, 1) string = ""
      kwargs.source_version (1, 1) string = ""
      kwargs.retrieval_date (1, 1) string = string(datetime('today'))
      kwargs.overwrite_family (1, 1) logical = false
      kwargs.overwrite (1, 1) logical = false
      kwargs.entry_callback = []
      kwargs.build_forcing (1, 1) logical = false
      kwargs.forcing_sources (1, :) string ...
         {icemodel.verification.validators.mustBeRcmSourceSelection( ...
         kwargs.forcing_sources)} = strings(1, 0)
      kwargs.leg_callback = []
      kwargs.after_source_callback = []
      kwargs.met_outdir (1, 1) string = ""
      kwargs.userdata_outdir (1, 1) string = ""
      kwargs.mar_dir (1, 1) string = ""
      kwargs.merra_dir (1, 1) string = ""
      kwargs.racmo_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.method (1, 1) string ...
         {mustBeMember(kwargs.method, ["nearest", "natural"])} = "nearest"
      kwargs.dt_out (1, 1) string ...
         {mustBeMember(kwargs.dt_out, ["", "15m"])} = "15m"
   end

   % Persist native observation state before optional RCM work so interruption
   % cannot lose a completed family-specific import.
   persist = @(st) icemodel.verification.setup.buildDatasetFamilyManifest(st, ...
      alive, dataset_family=kwargs.dataset_family, ...
      manifest_file=kwargs.manifest_file, ...
      requested_ids=kwargs.requested_ids, skipped=kwargs.skipped, ...
      source_doi=kwargs.source_doi, source_url=kwargs.source_url, ...
      source_version=kwargs.source_version, ...
      retrieval_date=kwargs.retrieval_date, ...
      overwrite_family=kwargs.overwrite_family, ...
      entry_callback=kwargs.entry_callback);
   manifest = persist(state);

   % A disabled/empty forcing request ends here without reading input caches or
   % mutating importer state with omitted source legs.
   if ~kwargs.build_forcing || isempty(kwargs.forcing_sources)
      return
   end

   % Carry additive manifest legs into the in-memory state before a requested
   % source refresh. This is preservation only; stageRcmForcing remains the one
   % owner of cache discovery, derivation, and writes.
   state = carryPersistedColocation(state, alive, manifest);
   state = icemodel.verification.setup.stageDatasetRcmForcing(state, alive, ...
      dataset_family=kwargs.dataset_family, ...
      forcing_sources=kwargs.forcing_sources, ...
      leg_callback=kwargs.leg_callback, ...
      after_source_callback=kwargs.after_source_callback, ...
      persist_callback=persist, ...
      met_outdir=kwargs.met_outdir, ...
      userdata_outdir=kwargs.userdata_outdir, ...
      mar_dir=kwargs.mar_dir, merra_dir=kwargs.merra_dir, ...
      racmo_dir=kwargs.racmo_dir, modis_dir=kwargs.modis_dir, ...
      method=kwargs.method, dt_out=kwargs.dt_out, ...
      overwrite=kwargs.overwrite);
   manifest = persist(state);
end

function state = carryPersistedColocation(state, alive, manifest)
   %CARRYPERSISTEDCOLOCATION Seed requested refreshes with additive prior legs.
   if isempty(state) || ~isfield(state, 'colocation') ...
         || ~isfield(manifest, 'cases') || isempty(manifest.cases)
      return
   end

   % Match each live importer record by its stable manifest case identity.
   cases = manifest.cases;
   ids = string({cases.case_id});
   for k = find(alive)
      id = stateCaseId(state(k));
      hit = find(ids == id, 1);
      if isempty(hit) || ~isfield(cases(hit), 'colocation') ...
            || ~isstruct(cases(hit).colocation)
         continue
      end
      state(k).colocation = ...
         icemodel.verification.setup.mergeColocation( ...
         cases(hit).colocation, state(k).colocation);
   end
end

function id = stateCaseId(state)
   %STATECASEID Return the stable case identity used by the family manifest.
   if isfield(state, 'case_id') && strlength(string(state.case_id)) > 0
      id = string(state.case_id);
   else
      id = "";
   end
end
