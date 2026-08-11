function simyears = outputYears(opts)
   %OUTPUTYEARS Return the simulation years retained in saved output.
   %
   %  simyears = icemodel.outputYears(opts)
   %
   % OPTS.SIMYEARS lists the forcing years in run order. The model runs the
   % first OPTS.N_SPINUP_YEARS years only for spinup, and the saved and
   % postprocessed output excludes them.

   simyears = opts.simyears(opts.n_spinup_years+1:end);
end
