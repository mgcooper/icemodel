function [iter, subiter, maxiter, maxsubiter, dt, dt_min, dt_max, dt_new, ...
      numyears, simyears, numspinup] = INITTIMESTEPS(opts, Time)
   %INITTIMESTEP initialize timestepping

   narginchk(0, 2)

   iter = 1;
   subiter = 1;

   if nargin == 0
      return
   end

   % compute the timestep and the total number of model timesteps
   dt = seconds(Time(2)-Time(1));
   maxiter = numel(Time)/opts.numyears;
   maxsubiter = 200;
   minsubiter = 1;

   dt_min = dt/maxsubiter;  % or: dt_min = 5; maxsubiter  = dt/dt_min;
   dt_max = dt/minsubiter;  % or: dt_min = 5; maxsubiter  = dt/dt_min;
   dt_new = dt_max;

   % this is
   numspinup = opts.spinup_loops;
   simyears = opts.simyears;
   if numspinup > 1
      simyears = [repmat(simyears(1),numspinup,1) simyears(2:end)];
   end
   numyears = numel(simyears);
end
