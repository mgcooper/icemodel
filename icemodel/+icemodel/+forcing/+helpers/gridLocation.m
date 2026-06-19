function [start, count, collapse, inslab, loctype] = gridLocation( ...
      X, Y, location, method, kwargs)
   %GRIDLOCATION Map a point or polygon onto a grid hyperslab + collapse rule.
   %
   %  [start, count, collapse, inslab, loctype] = ...
   %     icemodel.forcing.helpers.gridLocation(X, Y, location)
   %  [...] = gridLocation(X, Y, location, method)
   %
   % Shared cell-selection logic for the gridded-source builders (MAR,
   % RACMO, MERRA). Given the projected cell-center grids X, Y [m,
   % EPSG:3413] and a target location, returns the bounding hyperslab over
   % the two grid dimensions and a COLLAPSE function handle that reduces a
   % hyperslab data block (cells x time, cells flattened column-major to
   % match an ncread/reshape of start..start+count-1) to the target series
   % (time x 1):
   %
   %  - location = [x y] (projected point, 1x2 metres):
   %      method "nearest" (default) -> the single nearest cell.
   %      method "natural"           -> a padded neighbourhood collapsed by
   %                                    natural-neighbour interpolation at the
   %                                    point (icemodel.forcing.helpers.interpRcm).
   %  - location = polyshape (vertices in EPSG:3413 metres):
   %      remap="conservative" (default) -> exact overlap-area-weighted remap
   %      via the exactremap toolbox (helpers.remapPolygon), honouring optional
   %      true cell areas and a valid (ice) mask; remap="equal" -> plain mean
   %      of in-polygon cell centres (no external dependency, the interim
   %      fallback when exactremap is unavailable).
   %
   % The "natural" point option exists because near steep terrain (e.g. the
   % Greenland western slope) a station sits between cells of differing model
   % elevation, and nearest-cell extraction yields a step discontinuity at
   % cell boundaries; natural-neighbour blends the surrounding cells
   % continuously. (Note: this is horizontal blending only; it does not
   % lapse-correct to the station elevation.)
   %
   % Outputs
   %  start, count - hyperslab start [i j] (1-based) and extent [ni nj]
   %  collapse     - function handle: (cells x time block) -> (time x 1)
   %  inslab       - slab-relative linear index(es) used for site metadata
   %                 (the nearest cell for a point; the in-polygon cells for
   %                 a polygon), to be averaged with helpers.slabMean
   %  loctype      - "point" or "polygon"
   %
   % See also: icemodel.forcing.helpers.slabMean,
   %  icemodel.forcing.helpers.interpRcm, icemodel.forcing.buildMarData

   arguments
      X (:, :) double
      Y (:, :) double
      location
      method (1, 1) string {mustBeMember(method, ...
         ["nearest", "natural"])} = "nearest"
      kwargs.remap (1, 1) string {mustBeMember(kwargs.remap, ...
         ["equal", "conservative"])} = "conservative"
      kwargs.cellareas double = []
      kwargs.validmask = []
      kwargs.normalizationdomain (1, 1) string ...
         {mustBeMember(kwargs.normalizationdomain, ["polygon", "valid"])} ...
         = "polygon"
      kwargs.missingvalue double = []
      kwargs.usegeocoords (1, 1) logical = false
   end

   if isnumeric(location)
      assert(isequal(size(location), [1 2]), ...
         'point location must be a 1x2 projected [x y] in metres')
      [~, idx] = min(hypot(X(:) - location(1), Y(:) - location(2)));
      [i, j] = ind2sub(size(X), idx);
      loctype = "point";

      switch method
         case "nearest"
            start = [i j];
            count = [1 1];
            inslab = 1;
            collapse = @(block) reshape(block, 1, []).';

         case "natural"
            % Padded neighbourhood around the nearest cell (5x5 interior;
            % clamped at domain edges). Natural-neighbour interpolation to
            % the point reuses the polygon-side interpolation machinery.
            pad = 2;
            [start, count, inslab] = paddedSlab(size(X), i, j, pad);
            [xn, yn] = slabCoords(X, Y, start, count);
            xq = location(1);
            yq = location(2);
            collapse = @(block) icemodel.forcing.helpers.interpRcm( ...
               xn, yn, block, xq, yq, method="natural");
      end

   elseif isa(location, 'polyshape')
      loctype = "polygon";

      % Slab covering the polygon bounding box plus a 2-cell pad, so every
      % polygon-overlapping cell (including boundary cells whose centres
      % fall outside the polygon) and the edge-reconstruction neighbours
      % are present. Using only centre-in-polygon cells underfills small
      % catchments and breaks the conservative overlap geometry.
      [bx, by] = boundingbox(location);
      within = X >= bx(1) & X <= bx(2) & Y >= by(1) & Y <= by(2);
      if any(within(:))
         [ii, jj] = find(within);
      else   % sub-cell polygon: seed from the cell nearest the centroid
         [cx, cy] = centroid(location);
         [~, idx] = min(hypot(X(:) - cx, Y(:) - cy));
         [ii, jj] = ind2sub(size(X), idx);
      end
      pad = 2;
      start = [max(1, min(ii) - pad), max(1, min(jj) - pad)];
      iend = [min(size(X, 1), max(ii) + pad), min(size(X, 2), max(jj) + pad)];
      count = [iend(1) - start(1) + 1, iend(2) - start(2) + 1];

      rows = start(1):start(1) + count(1) - 1;
      cols = start(2):start(2) + count(2) - 1;
      Xs = X(rows, cols);
      Ys = Y(rows, cols);

      % In-polygon cells within the slab (site metadata + equal weights).
      inslab = find(isinterior(location, Xs(:), Ys(:)));
      if isempty(inslab)   % sub-cell polygon: nearest cell to the centroid
         [cx, cy] = centroid(location);
         [~, inslab] = min(hypot(Xs(:) - cx, Ys(:) - cy));
      end

      switch kwargs.remap
         case "equal"
            % Interim equal-weight mean of in-polygon cell centres.
            weights = zeros(numel(Xs), 1);
            weights(inslab) = 1 / numel(inslab);
            collapse = @(block) (weights.' * block).';

         case "conservative"
            % Exact overlap-area-weighted remap via exactremap over the
            % slab (cells outside the polygon get zero weight). Optional
            % true cell areas and a valid (ice) mask slice to the slab;
            % masked cells are inpainted from valid neighbours (see
            % helpers.remapPolygon).
            ca = sliceOpt(kwargs.cellareas, rows, cols);
            vm = sliceOpt(kwargs.validmask, rows, cols);
            P = location;
            nd = kwargs.normalizationdomain;
            mv = kwargs.missingvalue;
            collapse = @(block) icemodel.forcing.helpers.remapPolygon( ...
               Xs, Ys, block, P, cellareas=ca, validmask=vm, ...
               normalizationdomain=nd, missingvalue=mv, ...
               usegeocoords=kwargs.usegeocoords);
      end

   else
      error('icemodel:forcing:gridLocation:badLocation', ...
         'location must be a projected [x y] point or a polyshape')
   end
end

%% Local functions
function [start, count, inslab] = paddedSlab(gridsize, i, j, pad)
   %PADDEDSLAB Bounding hyperslab of cell (i,j) padded by PAD, clamped.
   i0 = max(1, i - pad);
   i1 = min(gridsize(1), i + pad);
   j0 = max(1, j - pad);
   j1 = min(gridsize(2), j + pad);
   start = [i0 j0];
   count = [i1 - i0 + 1, j1 - j0 + 1];
   inslab = sub2ind(count, i - i0 + 1, j - j0 + 1);   % nearest cell in slab
end

function [xn, yn] = slabCoords(X, Y, start, count)
   %SLABCOORDS Column-major cell coordinates over the hyperslab.
   rows = start(1):start(1) + count(1) - 1;
   cols = start(2):start(2) + count(2) - 1;
   xs = X(rows, cols);
   ys = Y(rows, cols);
   xn = xs(:);
   yn = ys(:);
end

function S = sliceOpt(M, rows, cols)
   %SLICEOPT Slice an optional full-grid field to the slab, or pass [] through.
   if isempty(M)
      S = [];
   else
      S = M(rows, cols);
   end
end
