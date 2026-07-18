function entry = stateCaseEntry(state)
   %STATECASEENTRY Refresh one staged state record into a manifest case entry.
   %
   %  entry = icemodel.verification.setup.stateCaseEntry(state)
   %
   % Atomic states and dry-run states already carry their complete entry.
   % Persisted firn states replace colocation and derive the public source lists
   % from those staged legs.

   % Preserve atomic or dry-run entries exactly because they have no final staged
   % colocation graph to recompute.
   entry = state.entry;
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
