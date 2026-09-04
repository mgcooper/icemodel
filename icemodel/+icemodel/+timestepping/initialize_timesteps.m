function [metstep, substep, numsteps, dt_new, numyears, numspinup] = ...
      initialize_timesteps(opts, Time)
   %INITIALIZE_TIMESTEPS Initialize the model timestep counters.
   %
   %  [metstep, substep, numsteps, dt_new, numyears, numspinup] = ...
   %     icemodel.timestepping.initialize_timesteps(opts, Time)
   %
   % NUMSTEPS is the number of forcing steps per simulated year. DT_NEW
   % starts at the full forcing-step length OPTS.DT. NUMSPINUP is the number
   % of leading simulated years excluded from saved output.
   %
   % See also: icemodel, skinmodel, icemodel.timestepping.nexttimestep
   %
   %#codegen

   narginchk(0, 2)

   metstep = 1;
   substep = 1;

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
