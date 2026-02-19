function [metiter, subiter, maxiter, maxsubiter, dt_new, dt_FULL_STEP, ...
      numyears, numspinup, simyears] = INITTIMESTEPS(opts, Time)
   %INITTIMESTEP initialize timestepping
   %
   %#codegen

   narginchk(0, 2)

   metiter = 1;
   subiter = 1;

   if nargin == 0
      return
   end

   % Compute the timestep and the total number of model timesteps
   dt_FULL_STEP = opts.dt;
   maxiter = numel(Time) / opts.numyears;
   maxsubiter = opts.dt; % allow 1 sec dt min
   minsubiter = 1;

   dt_min = dt_FULL_STEP / maxsubiter;
   dt_new = dt_FULL_STEP / minsubiter;

   % Compute the number of spinup years and a vector of simulation years
   numspinup = opts.spinup_loops;
   simyears = opts.simyears(:);
   if numspinup > 1
      simyears = [repmat(simyears(1), numspinup,1); simyears(2:end)];
   end
   numyears = numel(simyears);
end
