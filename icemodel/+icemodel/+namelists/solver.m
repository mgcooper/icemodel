function list = solver()
   %SOLVER Return the supported icemodel solver ids.
   %
   %  list = icemodel.namelists.solver()
   %
   % Solver 0 is the Dirichlet single sweep, 1 the coupled Dirichlet
   % iterations, 2 the Robin single sweep, and 3 the coupled Robin iterations
   % (icemodel.setopts).

   list = [0 1 2 3];
end
