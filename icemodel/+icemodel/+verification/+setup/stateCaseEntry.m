function entry = stateCaseEntry(state)
   %STATECASEENTRY Refresh one staged state record into a manifest case entry.
   %
   %  entry = icemodel.verification.setup.stateCaseEntry(state)
   %
   % Atomic states carry their complete entry. Firn states either construct the
   % entry from the common staging schema or refresh a reused entry after RCM
   % legs have been attached.

   % Construct a fresh firn entry only after every requested leg is attached so
   % all analogous importers share one canonical state-to-manifest conversion.
   entry = state.entry;
   if isempty(fieldnames(entry))
      [forcing_sources, eval_sources] = ...
         icemodel.verification.setup.colocationSourceLists(state.colocation);
      if state.dry_run && isempty(eval_sources)
         eval_sources = state.eval_source;
      end
      case_values = { ...
         char(state.case_id)
         'firn_observational'
         char(state.site_id)
         char(state.site_name)
         char(state.surface_zone)
         cellstr(state.eval_target)
         char(state.permafrost_zone)
         state.site_location
         state.period
         state.evaluation_file_rel
         cellstr(forcing_sources(:))
         cellstr(eval_sources(:))
         cellstr(state.comparison_variables(:))
         state.observation_variables
         state.colocation
         char(state.native_timestep)
         char(state.notes)};
      entry = ...
         icemodel.verification.setup.makeFirnCaseManifestEntry(case_values);
      return
   end

   % Preserve atomic or dry-run entries exactly because they have no final staged
   % colocation graph to recompute.
   if ~isfield(state, 'colocation') ...
         || (isfield(state, 'dry_run') && state.dry_run)
      return
   end

   % Persisted entries derive source inventories from the final staged legs so
   % repeated imports cannot leave stale forcing or evaluation source lists.
   entry.colocation = state.colocation;
   [forcing_sources, eval_sources] = ...
      icemodel.verification.setup.colocationSourceLists(entry.colocation);
   entry.forcing_sources = cellstr(forcing_sources);
   entry.eval_sources = cellstr(eval_sources);
end
