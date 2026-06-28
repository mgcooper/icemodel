function anchors = mixedAnchorCatalog(kwargs)
   %MIXEDANCHORCATALOG Read staged firn/research anchors from manifests.
   %
   %  anchors = icemodel.verification.setup.mixedAnchorCatalog()
   %  anchors = icemodel.verification.setup.mixedAnchorCatalog( ...
   %     output_root=output_root)
   %  anchors = icemodel.verification.setup.mixedAnchorCatalog( ...
   %     evaluation_data_root=eval_root)
   %
   % Returns a flat struct array with family, site/source ids, WGS84/EPSG:3413
   % coordinates, period, and available met/userdata/eval source labels. Missing
   % family manifests are ignored so callers can use it during partial staging.
   %
   % Role
   %  Setup catalog helper for mixed-anchor provenance. It reads staged manifests
   %  only; it does not fetch, stage, or validate source data.
   %
   % See also: icemodel.verification.setup.anchorColocation

   arguments
      kwargs.output_root (1, 1) string = ""
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
   end

   families = ["promice", "retmip", "imau", "research_site"];
   eval_root = resolveEvalRoot(kwargs);
   rows = cell(1, numel(families));
   n_rows = 0;

   % Each family can be staged independently; skip absent manifests rather than
   % turning a partial eval tree into a hard failure for colocation metadata.
   for family = families
      manifest_file = fullfile(eval_root, family, "manifest.json");
      if ~isfile(manifest_file)
         continue
      end
      add = familyAnchors(manifest_file, family);
      for k = 1:numel(add)
         n_rows = n_rows + 1;
         rows{n_rows} = add{k};
      end
   end

   if n_rows == 0
      anchors = emptyAnchor();
   else
      anchors = [rows{1:n_rows}];
   end
end

function eval_root = resolveEvalRoot(kwargs)
   %RESOLVEEVALROOT Honor output_root=<root>/eval before normal eval roots.
   if kwargs.output_root ~= ""
      eval_root = fullfile(kwargs.output_root, "eval");
   else
      eval_root = icemodel.verification.helpers.evaluationDataRoot( ...
         evaluation_data_root=kwargs.evaluation_data_root, ...
         icemodel_config_casename=kwargs.icemodel_config_casename);
   end
end

function rows = familyAnchors(manifest_file, family)
   %FAMILYANCHORS Convert one family manifest into anchor rows.
   manifest = jsondecode(fileread(manifest_file));
   if ~isfield(manifest, 'cases') || isempty(manifest.cases)
      rows = {};
      return
   end

   cases = manifest.cases;
   rows = cell(1, numel(cases));
   for k = 1:numel(cases)
      rows{k} = caseAnchor(cases(k), family);
   end
end

function a = caseAnchor(c, family)
   %CASEANCHOR Pull the common anchor fields from one manifest case.
   loc = fieldOr(c, 'site_location', struct());
   colocation = fieldOr(c, 'colocation', struct());

   a = anchorPrototype();
   a.family = string(family);
   a.site = string(fieldOr(c, 'site_id', fieldOr(c, 'case_id', "")));
   a.source_id = sourceId(c, colocation);
   a.case_id = string(fieldOr(c, 'case_id', a.source_id));
   a.lat_wgs84 = scalarField(loc, 'lat_wgs84');
   a.lon_wgs84 = scalarField(loc, 'lon_wgs84');
   a.x_epsg3413 = scalarField(loc, 'x_epsg3413');
   a.y_epsg3413 = scalarField(loc, 'y_epsg3413');
   a.elev_m = scalarField(loc, 'elev_m');
   a.period = fieldOr(c, 'period', struct('start', '', 'end', ''));
   a.met_sources = string(fieldOr(c, 'forcing_sources', strings(0, 1)));
   a.userdata_sources = userdataSources(colocation);
   a.eval_sources = string(fieldOr(c, 'eval_sources', strings(0, 1)));
end

function a = emptyAnchor()
   %EMPTYANCHOR Prototype catalog row with stable field order.
   a = anchorPrototype();
   a = a([]);
end

function a = anchorPrototype()
   %ANCHORPROTOTYPE Scalar catalog row used before field assignment.
   a = struct('family', "", 'site', "", 'source_id', "", 'case_id', "", ...
      'lat_wgs84', NaN, 'lon_wgs84', NaN, ...
      'x_epsg3413', NaN, 'y_epsg3413', NaN, 'elev_m', NaN, ...
      'period', struct('start', '', 'end', ''), ...
      'met_sources', strings(0, 1), ...
      'userdata_sources', strings(0, 1), ...
      'eval_sources', strings(0, 1));
end

function id = sourceId(c, colocation)
   %SOURCEID Prefer explicit source-association ids, then case ids.
   id = "";
   if isfield(colocation, 'source_association') ...
         && isfield(colocation.source_association, 'source_id')
      id = string(colocation.source_association.source_id);
   elseif isfield(c, 'source_id')
      id = string(c.source_id);
   elseif isfield(c, 'case_id')
      id = string(c.case_id);
   end
end

function sources = userdataSources(colocation)
   %USERDATASOURCES List staged colocation legs that carry userdata/Data files.
   fields = string(fieldnames(colocation));
   keep = false(numel(fields), 1);
   for k = 1:numel(fields)
      leg = colocation.(fields(k));
      keep(k) = isstruct(leg) && isfield(leg, 'staged') ...
         && logical(leg.staged) && isfield(leg, 'data_files') ...
         && ~isempty(leg.data_files);
   end
   sources = fields(keep);
end

function v = fieldOr(s, name, default)
   %FIELDOR Read a struct field with a default for mixed manifest schemas.
   if isstruct(s) && isfield(s, name)
      v = s.(name);
   else
      v = default;
   end
end

function v = scalarField(s, name)
   %SCALARFIELD Read one numeric field, returning NaN when absent.
   if isstruct(s) && isfield(s, name)
      v = double(s.(name));
   else
      v = NaN;
   end
end
