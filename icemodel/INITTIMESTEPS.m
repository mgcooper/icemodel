function [metstep, substep, numsteps, maxsubstep, dt_new, dt_FULL_STEP, ...
      numyears, numspinup, simyears] = INITTIMESTEPS(opts, Time)
   %INITTIMESTEP initialize timestepping
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
   numsteps = numel(Time) / opts.numyears;
   maxsubstep = opts.dt; % allow 1 sec dt min
   minsubstep = 1;

   dt_min = dt_FULL_STEP / maxsubstep;
   dt_new = dt_FULL_STEP / minsubstep;

   % Compute the number of spinup years and a vector of simulation years
   numspinup = opts.spinup_loops;
   simyears = opts.simyears(:);
   if numspinup > 1
      simyears = [repmat(simyears(1), numspinup,1); simyears(2:end)];
   end
   numyears = numel(simyears);
end
