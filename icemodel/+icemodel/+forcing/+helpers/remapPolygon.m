function series = remapPolygon(X, Y, block, P, kwargs)
   %REMAPPOLYGON Conservative area-weighted remap of a grid block to a polygon.
   %
   %  series = icemodel.forcing.helpers.remapPolygon(X, Y, block, P)
   %  series = ... remapPolygon(_, cellareas=..., validmask=..., ...
   %     normalizationdomain="valid", infillmasked=true, missingvalue=...)
   %
   % Thin wrapper over the exactremap toolbox that area-weight-averages a
   % gridded data block onto a polygon (the conservative alternative to the
   % equal-weight cell mean in icemodel.forcing.helpers.gridLocation). The
   % weight is the exact overlap area of the polygon with each grid cell.
   %
   % Inputs
   %  X, Y  : grid cell-center coordinates (N x M), in a metric CRS
   %          consistent with P (the builders use EPSG:3413 metres). The
   %          icemodel builders deliver these in ndgrid layout (X varies
   %          down dim 1); exactremap expects meshgrid layout (X along
   %          dim 2), so this wrapper transposes when needed (see below).
   %  block : data, cells-by-time, cells flattened column-major to match
   %          X(:)/Y(:) (as produced by reshape(ncread_slab, prod(count), []))
   %  P     : polyshape (or vertex array) in the same CRS as X, Y
   %
   % Name-value
   %  cellareas : N x M true cell areas [m^2]. When omitted, exactremap
   %      derives areas from the coordinates. RACMO ships gridarea (km^2);
   %      pass it as gridarea*1e6. Areas matter for area-sums; for the
   %      area-weighted MEAN returned here they are second-order (the
   %      weights are a ratio), but passing true areas is exact.
   %  validmask : N x M logical, true = valid cell. Invalid cells (e.g.
   %      off-ice) are infilled from valid neighbours (infillmasked=true,
   %      the default) before averaging; cells that remain invalid follow
   %      the normalizationdomain / missingvalue policy.
   %  infillmasked : extend valid data across masked cells (default true).
   %  normalizationdomain : "polygon" (default; result represents the whole
   %      polygon) or "valid" (renormalise over valid-overlap area only).
   %  missingvalue : value assigned to masked area when not infilled and
   %      normalizationdomain="polygon" (e.g. 0). Default unset.
   %
   % Output
   %  series : time-by-1 area-weighted average over the polygon.
   %
   % exactremap contract (verified 2026-06-17):
   %  - input grids must be in meshgrid / north-up orientation (dim 1 = Y
   %    rows, dim 2 = X columns); exactremap normalises row/column FLIPS but
   %    not the ndgrid/meshgrid TRANSPOSE, so this wrapper transposes when
   %    the incoming grid is ndgrid (X increasing down dim 1).
   %  - data must be shaped like the grid: a 2-D (single time) or 3-D
   %    N-by-M-by-T array, OR fully-flattened coordinate lists
   %    (V(:), X(:), Y(:) all consistent). A V(:) column paired with 2-D
   %    X, Y is a mixed form that mis-aligns, so the caller's cells-by-time
   %    block is reshaped to a 3-D meshgrid-oriented array here.
   %
   % Dependency: the exactremap toolbox must be on the MATLAB path
   % (icemodel.internal.installRequiredFiles handles this in a later pass;
   % until then add toolboxes/libraries/spatial/exactremap/toolbox).
   %
   % See also: exactremap, icemodel.forcing.helpers.gridLocation

   arguments
      X (:, :) double
      Y (:, :) double
      block (:, :) double
      P
      kwargs.cellareas double = []
      kwargs.validmask = []
      kwargs.infillmasked (1, 1) logical = true
      kwargs.normalizationdomain (1, 1) string ...
         {mustBeMember(kwargs.normalizationdomain, ["polygon", "valid"])} ...
         = "polygon"
      kwargs.missingvalue double = []
      kwargs.usegeocoords (1, 1) logical = false
   end

   if isempty(which('exactremap'))
      error('icemodel:forcing:remapPolygon:exactremapMissing', ...
         ['the exactremap toolbox is not on the MATLAB path; add ' ...
         'toolboxes/libraries/spatial/exactremap/toolbox'])
   end

   % Reorient ndgrid -> meshgrid if X increases predominantly down dim 1
   % (the icemodel builder convention). exactremap normalises flips but not
   % the transpose, so align here and reshuffle the data/masks identically.
   [ni, nm] = size(X);
   if gridIsNdgrid(X)
      X = X.';
      Y = Y.';
      if ~isempty(kwargs.cellareas), kwargs.cellareas = kwargs.cellareas.'; end
      if ~isempty(kwargs.validmask), kwargs.validmask = kwargs.validmask.'; end
      % block rows are column-major over the original (ni x nm) slab; map to
      % the transposed (nm x ni) layout, then reshape to N-by-M-by-T below.
      V = permute(reshape(block, ni, nm, []), [2 1 3]);
   else
      V = reshape(block, ni, nm, []);
   end

   % Assemble the exactremap name-value options. UseGeoCoords=true treats
   % X, Y as longitude/latitude and computes ellipsoidal overlap areas (the
   % correct conservative weighting for a native geographic grid, e.g.
   % MERRA-2); false treats them as planar metric coordinates.
   args = {'UseGeoCoords', kwargs.usegeocoords};
   if ~isempty(kwargs.cellareas)
      args = [args, {'CellAreas', kwargs.cellareas}];
   end
   if ~isempty(kwargs.validmask)
      args = [args, {'ValidCellsMask', logical(kwargs.validmask), ...
         'NormalizationDomain', kwargs.normalizationdomain}];
      if kwargs.infillmasked
         args = [args, {'InfillMasked', true}];
      end
      if ~isempty(kwargs.missingvalue)
         args = [args, {'MissingValue', kwargs.missingvalue}];
      end
   end

   % exactremap returns the per-polygon area-weighted average for each of
   % the T layers of the N-by-M-by-T block.
   series = exactremap(V, X, Y, P, 'areaavg', args{:});
   series = series(:);
end

%% Local functions
function tf = gridIsNdgrid(X)
   %GRIDISNDGRID True if X increases predominantly down dim 1 (ndgrid layout).
   % For a meshgrid grid X is (near) constant down rows; for ndgrid it is
   % (near) constant across columns. Robust for the regular polar-stereo
   % grids (MAR, MERRA) and the gently rotated FGRN11/RACMO curvilinear grid.
   if size(X, 1) < 2 && size(X, 2) < 2
      tf = false;   % degenerate single-cell slab: nothing to transpose
      return
   end
   d1 = mean(abs(diff(X, 1, 1)), 'all');   % variation down rows   (dim 1)
   d2 = mean(abs(diff(X, 1, 2)), 'all');   % variation across cols (dim 2)
   tf = d1 > d2;
end
