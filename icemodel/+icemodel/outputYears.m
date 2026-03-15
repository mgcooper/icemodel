function simyears = outputYears(opts)
   %OUTPUTYEARS Return the simulation years retained in saved output.
   %
   %  simyears = icemodel.outputYears(opts)
   %
   % OPTS.SIMYEARS lists the forcing years in run order. The first
   % OPTS.N_SPINUP_YEARS years are run only for spinup and are excluded from
   % saved/postprocessed output.

   simyears = opts.simyears(opts.n_spinup_years+1:end);
end
