function numcells = maxcells(memoryGB, numlayers, numtimesteps, numvars, dtype)
   %MAXCELLS Calculate the maximum number of gridcells for icemodel nc file.
   %
   % NUMCELLS = MAXCELLS(MEMORYGB, NUMLAYERS, NUMTIMESTEPS, NUMVARS, DTYPE)
   % Calculates the maximum number of grid cells that can be processed given a
   % memory constraint.
   %
   % Args:
   % - memoryGB: The memory limit in GB.
   % - numLayers: Number of layers in each array.
   % - numTimesteps: Number of time steps in each array.
   % - numVars: Number of variables, with one array per variable.
   % - dataType: The data type of the array elements ('single' or 'double').
   %
   % Returns:
   % - The maximum number of grid cells that fit within the memory limit.

   % Convert memory limit from GB to bytes
   memoryBytes = memoryGB * 1024^3; % 1GB = 1024^3 bytes

   % Calculate the bytes per element based on data type
   switch dtype
      case 'single'
         bytesPerElement = 4; % single precision (float) uses 4 bytes
      case 'double'
         bytesPerElement = 8; % double precision uses 8 bytes
      otherwise
         error('Unsupported data type. Choose ''single'' or ''double''.');
   end

   % Calculate bytes for one variable
   bytesPerVar = numlayers * numtimesteps * bytesPerElement;

   % Calculate total bytes for all variables
   totalBytesForAllVars = bytesPerVar * numvars;

   % Calculate the maximum number of grid cells that can fit into memory
   numcells = floor(memoryBytes / totalBytesForAllVars);
end
