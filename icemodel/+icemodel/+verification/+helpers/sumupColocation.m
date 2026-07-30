function [is_colocated, anchor, distance_km] = sumupColocation(x3413, y3413, kwargs)
   %SUMUPCOLOCATION Flag a SUMup point as co-located with a mixed anchor.
   %
   %  [tf, anchor, d_km] = ...
   %     icemodel.verification.helpers.sumupColocation(x3413, y3413)
   %  [tf, anchor, d_km] = ...
   %     icemodel.verification.helpers.sumupColocation(x3413, y3413, ...
   %        threshold_km=7.5)
   %
   % Inputs
   %  x3413, y3413   SUMup point coordinates in EPSG:3413 (NSIDC Sea Ice
   %                 Polar Stereographic North, the projected CRS the firn
   %                 staging driver records for every site).
   %
   % Name-value
   %  threshold_km   Co-location distance threshold in km (default 7.5).
   %  anchors        Optional struct array of anchor points with fields site,
   %                 x_epsg3413, y_epsg3413, and optionally family/source_id.
   %                 When omitted, anchors are resolved from the top-level
   %                 staged firn/research manifests.
   %
   % Outputs
   %  is_colocated   true when the SUMup point is within threshold_km of any
   %                 anchor.
   %  anchor         The nearest anchor struct (site + coords + distance), or
   %                 an empty struct when no anchor coordinates are available.
   %  distance_km    Distance to the nearest anchor in km (Inf when no anchor
   %                 coordinates are available).
   %
   % Co-location threshold rationale
   %  7.5 km is chosen as the default. PROMICE AWS sites sit on the GrIS
   %  ablation-to-percolation transect where surface-mass-balance gradients
   %  with elevation are steep, so the co-location radius must be small enough
   %  that a SUMup firn core and the PROMICE anchor share broadly the same
   %  accumulation regime, yet large enough to capture the sparse SUMup
   %  coverage near each anchor. 7.5 km is also comfortably inside one MAR
   %  (~7.5-15 km) / RACMO (~5.5-11 km) grid cell, so the co-located MAR/RACMO
   %  forcing extracted at the SUMup point is effectively the same model cell
   %  as the anchor's - the co-location is meaningful at the RCM resolution the
   %  reference is built from. EPSG:3413 is a conformal projection at GrIS
   %  latitudes, so planar distance is an accurate proxy for great-circle
   %  distance over a few km.
   %
   % Role
   %  Setup/staging helper used by importSumup to record which SUMup points
   %  fall near PROMICE, RetMIP, IMAU, or research_site anchors.
   %  Pure geometry; reads staged manifests but mutates nothing.
   %
   % See also: icemodel.verification.setup.importSumup,
   %  icemodel.verification.setup.promiceSiteCatalog

   arguments
      x3413 (1, 1) double
      y3413 (1, 1) double
      kwargs.threshold_km (1, 1) double {mustBePositive} = 7.5
      kwargs.anchors = []
      kwargs.output_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = ""
   end

   % Keep the legacy helper name for importSumup callers while sharing the
   % generic mixed-anchor geometry with research_site and future families.
   [is_colocated, anchor, distance_km] = ...
      icemodel.verification.setup.anchorColocation(x3413, y3413, ...
      threshold_km=kwargs.threshold_km, anchors=kwargs.anchors, ...
      output_root=kwargs.output_root, ...
      evaluation_data_root=kwargs.evaluation_data_root, ...
      icemodel_config_casename=kwargs.icemodel_config_casename);
end
