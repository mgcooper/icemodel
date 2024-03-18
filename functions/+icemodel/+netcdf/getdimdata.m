function dims = getdimdata(Z, dz, T, dt, opts)
   %GETDIMDATA Get dimensions of icemodel simulation data
   %
   %  DIMS = GETDIMDATA(Z, DZ, T, DT) Returns the spatial and temporal
   %  dimensions given icemodel vertical domain depth Z, grid spacing DZ,
   %  number of timesteps per year T, and timestep DT.
   %
   %  Note: If Z = 0, dz = 0, then dims.depth is an empty 1x0 vector and
   %  numel(dims.depth) = 0. See icemodel.netcdf.getdimsize where a max(1, ...)
   %  statement is used to set dimsize.depth = 1 in this case.
   %
   % See also:

   arguments
      Z  (1, 1) double {mustBeNumeric} = 0
      dz (1, 1) double {mustBeNumeric} = 0
      T  (1, 1) double {mustBeNumeric} = 8760
      dt (1, 1) double {mustBeNumeric} = 3600
      opts.whichmask = "icemask"
   end

   % -----------------------------  Below here sets the actual dimension values

   % Load the grid coordinates and elevation.
   [dims.x_easting, dims.y_northing, dims.elevation, ...
      dims.latitude, dims.longitude] = loadIceMask(opts.whichmask, ...
      varnames=["X", "Y", "Elev", "Lat", "Lon"]);

   % Define the grid cell index
   dims.gridcell = 1:numel(dims.elevation);

   % Define the layer depth dimension
   dims.depth = dz/2:dz:Z;

   % Define the time dimension
   dims.time = 0:dt:dt*(T-1);

   % -----------------------------  Below here returns their sizes

   % dims.numcells = numel(dims.x_easting);    % number of grid cells
   % dims.numlayers = numel(dims.depth);       % number of vertical layers
   % dims.numtimesteps = numel(dims.time);     % number of timesteps (hours) per year

   % Use this to set the dims directly from the data
   % [dims.numlayers, dims.numtimesteps] = size(data);

end
