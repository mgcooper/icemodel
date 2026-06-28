function [is_colocated, anchor, distance_km] = anchorColocation(x3413, y3413, kwargs)
   %ANCHORCOLOCATION Flag a point as co-located with the nearest anchor.
   %
   %  [tf, anchor, d_km] = ...
   %     icemodel.verification.setup.anchorColocation(x3413, y3413)
   %  [tf, anchor, d_km] = ...
   %     icemodel.verification.setup.anchorColocation(x3413, y3413, ...
   %        anchors=anchors, threshold_km=7.5)
   %
   % Inputs
   %  x3413, y3413   Query point coordinates in EPSG:3413.
   %
   % Name-value
   %  threshold_km   Co-location distance threshold in km (default 7.5).
   %  anchors        Optional struct array with x_epsg3413/y_epsg3413 fields.
   %                 When omitted, staged mixed anchors are read from manifests.
   %
   % Outputs
   %  is_colocated   true when the nearest anchor falls within threshold_km.
   %  anchor         Nearest anchor with distance_km attached, or struct([]).
   %  distance_km    Distance to nearest anchor in km, or Inf when unavailable.
   %
   % Role
   %  Generic setup helper for nearest-anchor provenance across PROMICE, RetMIP,
   %  IMAU, SUMup, and research_site staging. Pure geometry; mutates nothing.
   %
   % See also: icemodel.verification.setup.mixedAnchorCatalog,
   %  icemodel.verification.helpers.sumupColocation

   arguments
      x3413 (1, 1) double
      y3413 (1, 1) double
      kwargs.threshold_km (1, 1) double {mustBePositive} = 7.5
      kwargs.anchors = []
      kwargs.output_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
   end

   anchors = kwargs.anchors;
   if isempty(anchors)
      anchors = icemodel.verification.setup.mixedAnchorCatalog( ...
         output_root=kwargs.output_root, ...
         evaluation_data_root=kwargs.evaluation_data_root, ...
         icemodel_config_casename=kwargs.icemodel_config_casename);
   end

   % Missing or incomplete staged catalogs mean there is no usable distance
   % anchor; return a non-colocated result instead of fabricating provenance.
   [anchors, ax, ay] = finiteCoordinateAnchors(anchors);
   if isempty(anchors)
      is_colocated = false;
      anchor = struct([]);
      distance_km = Inf;
      return
   end

   % EPSG:3413 planar distance is appropriate for the few-km Greenland anchor
   % threshold used by the firn-evaluation staging workflow.
   d_km = hypot(ax - x3413, ay - y3413) / 1000;
   [distance_km, idx] = min(d_km);
   is_colocated = distance_km <= kwargs.threshold_km;
   anchor = anchors(idx);
   anchor.distance_km = distance_km;
end

function [anchors, ax, ay] = finiteCoordinateAnchors(anchors)
   %FINITECOORDINATEANCHORS Keep only anchors with projected coordinates.
   if isempty(anchors) || ~isfield(anchors, 'x_epsg3413') ...
         || ~isfield(anchors, 'y_epsg3413')
      anchors = struct([]);
      ax = [];
      ay = [];
      return
   end

   ax = [anchors.x_epsg3413];
   ay = [anchors.y_epsg3413];
   keep = isfinite(ax) & isfinite(ay);
   anchors = anchors(keep);
   ax = ax(keep);
   ay = ay(keep);
end
