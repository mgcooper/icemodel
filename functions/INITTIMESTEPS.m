function [iter, subiter, itime, maxiter, maxsubiter, dt, dt_min, ...
   dt_max, dt_new] = INITTIMESTEPS(opts, Time)
%INITTIMESTEP initialize timestepping

dt          =  opts.dt;
maxiter     =  opts.maxiter;
maxsubiter  =  200;
minsubiter  =  1;
dt_min      =  dt/maxsubiter;  % or: dt_min = 5; maxsubiter  = dt/dt_min;
dt_max      =  dt/minsubiter;  % or: dt_min = 5; maxsubiter  = dt/dt_min;
dt_new      =  dt_max;
iter        =  1;
subiter     =  1;
itime       =  Time(1);   % initialize model time