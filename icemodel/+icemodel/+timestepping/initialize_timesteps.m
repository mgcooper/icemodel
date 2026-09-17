function [metstep, substep, numsteps, dt_new, numyears, numspinup, ...
      force_advance_streak_dt] = initialize_timesteps(opts, Time)
   %INITIALIZE_TIMESTEPS Initialize the model timestep counters.
   %
   %  [metstep, substep, numsteps, dt_new, numyears, numspinup, ...
   %     force_advance_streak_dt] = ...
   %     icemodel.timestepping.initialize_timesteps(opts, Time)
   %
   % NUMSTEPS is the number of forcing steps per simulated year. DT_NEW
   % starts at the full forcing-step length OPTS.DT. NUMSPINUP is the number
   % of leading simulated years excluded from saved output.
   % FORCE_ADVANCE_STREAK_DT is the elapsed time [s] of consecutive forced
   % advances; icemodel.timestepping.checksubstep updates it, and it starts
   % at zero.
   %
   % See also: icemodel, skinmodel, icemodel.timestepping.nexttimestep
   %
   %#codegen

   narginchk(0, 2)

   % Start at the first forcing step and substep with no forced advances.
   metstep = 1;
   substep = 1;
   force_advance_streak_dt = 0.0;

   if nargin == 0
      return
   end

   % Compute the number of forcing steps per simulated year.
   assert(mod(numel(Time), opts.numyears) == 0)
   numsteps = numel(Time) / opts.numyears;
   dt_new = opts.dt;

   % Compute the number of leading spinup years. The model runs the forcing
   % years in opts.simyears in order. The saved and postprocessed output
   % excludes the first numspinup years.
   numspinup = opts.n_spinup_years;
   assert(numspinup < opts.numyears)
   numyears = numel(opts.simyears);
end
