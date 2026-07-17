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
      kwargs.icemodel_config_casename (1, 1) string = ""
   end

   families = setdiff(icemodel.verification.namelists.firndatasetfamily(), ...
      "sumup", 'stable');
   manifest_files = anchorManifestFiles(kwargs, families);
   rows = cell(1, numel(families));
   n_rows = 0;

   % Each family can be staged independently; skip absent manifests rather than
   % turning a partial eval tree into a hard failure for colocation metadata.
   for manifest_file = reshape(manifest_files, 1, [])
      [~, family] = fileparts(fileparts(manifest_file));
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

function files = anchorManifestFiles(kwargs, families)
   %ANCHORMANIFESTFILES List available anchor family manifests.
   if kwargs.output_root ~= ""
      files = icemodel.verification.helpers.familyManifestFiles( ...
         evaluation_data_root=fullfile(kwargs.output_root, "eval"));
   else
      files = icemodel.verification.helpers.familyManifestFiles( ...
         evaluation_data_root=kwargs.evaluation_data_root, ...
         icemodel_config_casename=kwargs.icemodel_config_casename);
   end
   file_families = strings(numel(files), 1);
   for k = 1:numel(files)
      [~, file_families(k)] = fileparts(fileparts(files(k)));
   end
   files = files(ismember(file_families, families));
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
   loc = icemodel.verification.helpers.fieldOr(c, 'site_location', struct());
   colocation = icemodel.verification.helpers.fieldOr(c, ...
      'colocation', struct());

   a = anchorPrototype();
   a.family = string(family);
   a.site = string(icemodel.verification.helpers.fieldOr(c, 'site_id', ...
      icemodel.verification.helpers.fieldOr(c, 'case_id', "")));
   a.source_id = sourceId(c, colocation);
   a.case_id = string(icemodel.verification.helpers.fieldOr(c, ...
      'case_id', a.source_id));
   a.lat_wgs84 = scalarField(loc, 'lat_wgs84');
   a.lon_wgs84 = scalarField(loc, 'lon_wgs84');
   a.x_epsg3413 = scalarField(loc, 'x_epsg3413');
   a.y_epsg3413 = scalarField(loc, 'y_epsg3413');
   [a.x_epsg3413, a.y_epsg3413] = fillProjectedCoordinates( ...
      a.lat_wgs84, a.lon_wgs84, a.x_epsg3413, a.y_epsg3413);
   a.elev_m = scalarField(loc, 'elev_m');
   a.period = icemodel.verification.helpers.fieldOr(c, 'period', ...
      struct('start', '', 'end', ''));
   a.surface_zone = string(icemodel.verification.helpers.fieldOr(c, ...
      'surface_zone', ""));
   a.eval_target = string(icemodel.verification.helpers.fieldOr(c, ...
      'eval_target', strings(0, 1)));
   a.permafrost_zone = string(icemodel.verification.helpers.fieldOr(c, ...
      'permafrost_zone', ""));
   a.met_sources = string(icemodel.verification.helpers.fieldOr(c, ...
      'forcing_sources', strings(0, 1)));
   a.userdata_sources = userdataSources(colocation);
   a.eval_sources = string(icemodel.verification.helpers.fieldOr(c, ...
      'eval_sources', strings(0, 1)));
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
      'surface_zone', "", 'eval_target', strings(0, 1), ...
      'permafrost_zone', "", ...
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

function v = scalarField(s, name)
   %SCALARFIELD Read one numeric field, returning NaN when absent.
   if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
      v = double(s.(name));
      if ~isscalar(v)
         v = NaN;
      end
   else
      v = NaN;
   end
end

function [x, y] = fillProjectedCoordinates(lat, lon, x, y)
   %FILLPROJECTEDCOORDINATES Derive EPSG:3413 when manifests omit it.
   if all(isfinite([lat lon])) && any(~isfinite([x y]))
      proj = icemodel.forcing.helpers.psnProjection();
      [x, y] = projfwd(proj, lat, lon);
   end
end
