function [is_colocated, anchor, distance_km] = sumupColocation(x3413, y3413, kwargs)
   %SUMUPCOLOCATION Flag a SUMup point as co-located with a PROMICE anchor.
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
   %  anchors        Optional Nx3 struct array of anchor points with fields
   %                 site, x_epsg3413, y_epsg3413. When omitted, the anchors
   %                 are resolved from the committed firn/promice manifest
   %                 (the staged KAN transect), falling back to the curated
   %                 promicesiteinfo catalog ids with no coordinates (in which
   %                 case co-location cannot be evaluated and returns false).
   %
   % Outputs
   %  is_colocated   true when the SUMup point is within threshold_km of any
   %                 anchor.
   %  anchor         The nearest anchor struct (site + coords + distance), or
   %                 an empty struct when no anchor is within range.
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
   %  fall near a PROMICE anchor (and thus share a curated co-located bundle).
   %  Pure geometry; reads the committed manifest but mutates nothing.
   %
   % See also: icemodel.verification.setup.importSumup,
   %  icemodel.verification.helpers.promicesiteinfo

   arguments
      x3413 (1, 1) double
      y3413 (1, 1) double
      kwargs.threshold_km (1, 1) double {mustBePositive} = 7.5
      kwargs.anchors = []
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
   end

   anchors = kwargs.anchors;
   if isempty(anchors)
      anchors = anchorsFromPromiceManifest( ...
         kwargs.evaluation_data_root, kwargs.icemodel_config_casename);
   end

   % With no anchor coordinates available, co-location cannot be evaluated.
   if isempty(anchors)
      is_colocated = false;
      anchor = struct([]);
      distance_km = Inf;
      return
   end

   ax = [anchors.x_epsg3413];
   ay = [anchors.y_epsg3413];
   d_km = hypot(ax - x3413, ay - y3413) / 1000;
   [distance_km, idx] = min(d_km);

   is_colocated = distance_km <= kwargs.threshold_km;
   if is_colocated
      anchor = anchors(idx);
      anchor.distance_km = distance_km;
   else
      anchor = struct([]);
   end
end

%% Local helpers
function anchors = anchorsFromPromiceManifest(evaluation_data_root, casename)
   %ANCHORSFROMPROMICEMANIFEST Read anchor coords from the promice manifest.
   %
   % The PROMICE anchor coordinates are not stored in a static list; the
   % staging driver records them (x_epsg3413/y_epsg3413) in the committed
   % promice manifest. Read them from there so the co-location reference
   % never drifts from the staged anchors. Returns [] when no manifest exists.

   anchors = struct('site', {}, 'x_epsg3413', {}, 'y_epsg3413', {});

   eval_root = icemodel.verification.helpers.evaluationDataRoot( ...
      "evaluation_data_root", evaluation_data_root, ...
      "icemodel_config_casename", casename);
   manifest_file = fullfile(eval_root, "promice", "manifest.json");
   if exist(manifest_file, 'file') ~= 2
      return
   end

   manifest = jsondecode(fileread(manifest_file));
   cases = manifest.cases;
   for n = 1:numel(cases)
      loc = cases(n).site_location;
      anchors(end + 1) = struct( ...
         'site', string(cases(n).site_id), ...
         'x_epsg3413', loc.x_epsg3413, ...
         'y_epsg3413', loc.y_epsg3413); %#ok<AGROW>
   end
end
