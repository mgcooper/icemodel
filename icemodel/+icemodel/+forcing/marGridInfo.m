function grid = marGridInfo(filename)
   %MARGRIDINFO Read the MAR grid coordinates and static fields.
   %
   %  grid = icemodel.forcing.marGridInfo(filename)
   %
   % Reads the MAR NetCDF coordinate grids and static surface fields and
   % returns them in both the native MAR projection and EPSG:3413:
   %
   %    grid.LON, grid.LAT - geographic cell centers [degrees], nx x ny
   %    grid.Xnat, grid.Ynat - MAR native polar-stereo cell centers [m],
   %                         a REGULAR grid (the 15 km MAR projection axes)
   %    grid.X, grid.Y     - cell centers reprojected to EPSG:3413 [m]
   %                         (for output metadata only)
   %    grid.elev          - surface height SH [m]
   %    grid.slope         - surface slope SLO [m/m]
   %    grid.srf           - MAR surface type (4 = ice sheet)
   %
   % All fields keep the native NetCDF orientation (x down rows, y across
   % columns), so linear indices align with hyperslab reads from
   % icemodel.forcing.readMar3p11.
   %
   % grid.Xnat/grid.Ynat come straight from the MAR projection axes
   % (variables named for the LON/LAT dimensions, e.g. X10_105/Y21_199, in
   % km), so they form an exactly regular grid suitable for conservative
   % area-weighted remap. grid.X/grid.Y are the SAME cells reprojected to
   % EPSG:3413, which is non-uniform (a curvilinear stereo-to-stereo map)
   % and therefore only used for reporting site coordinates.
   %
   % See also: icemodel.forcing.readMar3p11,
   %  icemodel.forcing.helpers.psnProjection

   arguments
      filename (1, 1) string
   end

   grid.LON = double(ncread(filename, 'LON'));
   grid.LAT = double(ncread(filename, 'LAT'));
   grid.elev = double(ncread(filename, 'SH'));
   grid.slope = double(ncread(filename, 'SLO'));
   grid.srf = double(ncread(filename, 'SRF'));

   % Native MAR projection axes: the LON dimensions name the 1-D coordinate
   % variables (km). ncinfo returns them in MATLAB (column-major) order, so
   % the first is the row (x) axis and the second the column (y) axis,
   % matching grid.LON's [nx ny] layout.
   loninfo = ncinfo(filename, 'LON');
   dimnames = {loninfo.Dimensions.Name};
   xkm = double(ncread(filename, dimnames{1}));
   ykm = double(ncread(filename, dimnames{2}));
   [grid.Xnat, grid.Ynat] = ndgrid(xkm * 1000, ykm * 1000);   % km -> m

   proj = icemodel.forcing.helpers.psnProjection();
   [grid.X, grid.Y] = projfwd(proj, grid.LAT, grid.LON);
end
