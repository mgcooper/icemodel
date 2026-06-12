function [metstep, substep, numsteps, maxsubstep, dt_new, dt_FULL_STEP, ...
      numyears, numspinup, simyears] = initialize_timesteps(opts, Time)
   % Initialize timestep counters and the full-step integration contract.
   %
   %#codegen

   narginchk(0, 2)

   metstep = 1;
   substep = 1;

   if nargin == 0
      return
   end

   % Compute the timestep and the total number of model timesteps
   dt_FULL_STEP = opts.dt;
   assert(mod(numel(Time), opts.numyears) == 0)
   numsteps = numel(Time) / opts.numyears;
   maxsubstep = opts.dt; % allow 1 sec dt min
   minsubstep = 1;

   % dt_min = dt_FULL_STEP / maxsubstep; % keep for reference
   dt_new = dt_FULL_STEP / minsubstep;

   % Compute the number of leading spinup years. The forcing years in
   % opts.simyears are run in order; the first numspinup years are excluded
   % from saved/postprocessed output.
   numspinup = opts.n_spinup_years;
   assert(numspinup < opts.numyears)
   simyears = opts.simyears(:);
   numyears = numel(simyears);
end
