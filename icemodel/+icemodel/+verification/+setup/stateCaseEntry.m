function entry = stateCaseEntry(state)
   %STATECASEENTRY Refresh one staged state record into a manifest case entry.
   %
   %  entry = icemodel.verification.setup.stateCaseEntry(state)
   %
   % Dry-run states already carry their complete metadata-only entry. Persisted
   % states replace colocation and derive the public source lists from those
   % staged legs.

   % Preserve source-light dry-run entries exactly as the family adapter built
   % them because no staged colocation facts are available to recompute.
   entry = state.entry;
   if state.dry_run
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
