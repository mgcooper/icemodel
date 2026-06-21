function info = promicesiteinfo(site)
   %PROMICESITEINFO Catalog of PROMICE anchor sites for firn evaluation.
   %
   %  info = icemodel.verification.helpers.promicesiteinfo()
   %  info = icemodel.verification.helpers.promicesiteinfo("KAN_M")
   %
   %  Returns the curated catalog of PROMICE automatic-weather-station anchor
   %  sites used to stage firn/snow-evaluation cases. With no arguments, returns
   %  a struct array; with one site argument (canonical id "KAN_M" or compact
   %  alias "kanm") returns the matching scalar struct.
   %
   %  The catalog carries TWO independent descriptors per site:
   %    surface_zone   the glaciological substrate REGIME the site sits in
   %                   (icemodel.verification.namelists.surfacezone vocabulary:
   %                   ablation/percolation/wet_snow/dry_snow/accumulation/
   %                   land/tundra/unknown). This is "where is the site".
   %    eval_target    a string ARRAY of which model CAPABILITIES the site
   %                   exercises (icemodel.verification.namelists.evaltarget
   %                   vocabulary: seasonal_snow/bare_ice/firn/ablation). This is
   %                   "what does a case here test". A site can have several.
   %
   %  Each entry has fields:
   %    site            canonical PROMICE station id ("KAN_M")
   %    alias           compact lowercase alias used in filenames ("kanm")
   %    long_name       display-friendly name
   %    surface_zone    glaciological zone (see above)
   %    eval_target     capability descriptor string array (see above)
   %    classification  "authoritative" for the curated KAN anchors, or
   %                    "first_pass" for the elevation/location_type-derived
   %                    non-KAN classifications (see provenance below)
   %    has_recipe      true when buildPromiceData carries a curated
   %                    service-window ablation recipe for the site
   %    models          co-located model set available at the site, in the
   %                    canonical order [promice, mar, merra, racmo]
   %    note            short provenance note
   %
   %  Site coordinates are NOT stored here: they are read live from the L3
   %  NetCDF metadata by readPromiceAws (latitude / longitude variables) and
   %  converted to EPSG:3413 by the staging driver, so the committed catalog
   %  never drifts from the source files.
   %
   %  ===========================================================================
   %  CLASSIFICATION PROVENANCE
   %  ---------------------------------------------------------------------------
   %  KAN transect (AUTHORITATIVE, reviewed with the user):
   %    KAN_L  surface_zone="ablation"     eval_target=["seasonal_snow","bare_ice"]
   %    KAN_M  surface_zone="ablation"     eval_target=["seasonal_snow","bare_ice"]
   %    KAN_U  surface_zone="percolation"  eval_target=["seasonal_snow","firn"]
   %
   %  All OTHER PROMICE stations: FIRST PASS ONLY, derived programmatically from
   %  AWS_sites_metadata.csv location_type + altitude_installation by
   %  firstPassZone() below. These are NOT authoritative and are flagged
   %  classification="first_pass" for USER REVIEW. The mapping is coarse:
   %    location_type "tundra"               -> surface_zone "tundra"
   %    location_type "bedrock"/"not Greenland" -> surface_zone "land"
   %    location_type "local glacier"        -> surface_zone "ablation"
   %      (marginal glaciers; eval_target ["seasonal_snow","bare_ice"])
   %    location_type "ice sheet", by elevation band (SW-GrIS firn-zone proxy):
   %       elev <  1500 m -> "ablation"     eval ["seasonal_snow","bare_ice"]
   %       1500-2000 m    -> "percolation"  eval ["seasonal_snow","firn"]
   %       elev >= 2000 m -> "accumulation" eval ["seasonal_snow","firn"]
   %    no elevation in the CSV               -> "unknown" (eval empty)
   %  The elevation bands are a single-transect heuristic; the firn-line
   %  elevation varies by region (NE/N GrIS differs from SW), so these need an
   %  expert pass before any non-KAN case is treated as authoritative.
   %  ===========================================================================
   %
   %  Role
   %    Catalog metadata used by the firn staging driver and report tooling.
   %    Not an argument validator - the driver accepts any PROMICE station id.
   %
   % See also: icemodel.verification.setup.importPromiceSites,
   %  icemodel.verification.namelists.surfacezone,
   %  icemodel.verification.namelists.evaltarget,
   %  icemodel.forcing.readPromiceAws, icemodel.forcing.buildPromiceData,
   %  icemodel.verification.helpers.snowmipinfo

   arguments
      site (1, 1) string = ""
   end

   % Canonical co-located model order. PROMICE is the station anchor; MAR and
   % MERRA are point met sources whose albedo is swapped downstream; RACMO is a
   % Data-only SMB/eval source (never met).
   models = ["promice", "mar", "merra", "racmo"];

   % Curated AUTHORITATIVE KAN anchors. Order is south-to-north up the KAN
   % transect so the staged sites read as an elevation transect.
   anchors = struct( ...
      'site',        {"KAN_L", "KAN_M", "KAN_U"}, ...
      'alias',       {"kanl",  "kanm",  "kanu"}, ...
      'long_name',   {"KAN-L", "KAN-M", "KAN-U"}, ...
      'surface_zone', {"ablation", "ablation", "percolation"}, ...
      'eval_target', { ...
         ["seasonal_snow"; "bare_ice"], ...
         ["seasonal_snow"; "bare_ice"], ...
         ["seasonal_snow"; "firn"]}, ...
      'classification', {"authoritative", "authoritative", "authoritative"}, ...
      'has_recipe',  {true, true, false}, ...
      'models',      {models, models, models}, ...
      'note', { ...
         "Lower ablation zone (~670 m); curated KAN_L recipe. KAN anchor.", ...
         "Upper ablation zone / bare ice (~1270 m); curated KAN_M recipe.", ...
         "Lower percolation zone (~1840 m); generic recipe. KAN anchor."});

   % First-pass non-KAN catalog, derived from AWS_sites_metadata.csv. These are
   % advisory and flagged classification="first_pass" for user review.
   firstpass = firstPassCatalog(models);

   catalog = [reshape(anchors, 1, []), reshape(firstpass, 1, [])];

   if site == ""
      info = catalog;
      return
   end

   wanted = lower(erase(site, "_"));
   match = string({catalog.alias}) == wanted;
   if ~any(match)
      valid = strjoin(string({catalog.site}), ', ');
      error('icemodel:verification:promicesiteinfo:unknown', ...
         ['unknown PROMICE anchor site "%s". Cataloged sites: %s. ' ...
         '(Any L3 station id can still be staged via importPromiceSites.)'], ...
         site, valid);
   end
   info = catalog(match);
end

%% Local helpers
function catalog = firstPassCatalog(models)
   %FIRSTPASSCATALOG Build the first-pass non-KAN catalog from the AWS CSV.
   %
   % Reads AWS_sites_metadata.csv once and classifies every non-KAN station via
   % firstPassZone(). Returns an empty struct array if the CSV is unavailable.
   catalog = repmat(emptyEntry(), 1, 0);

   csv = locateAwsCsv();
   if csv == "" || ~isfile(csv)
      return
   end
   T = readtable(csv, 'TextType', 'string');
   if ~all(ismember({'site_id', 'location_type'}, T.Properties.VariableNames))
      return
   end

   kan_anchors = ["KAN_L"; "KAN_M"; "KAN_U"];
   for n = 1:height(T)
      site = string(T.site_id(n));
      if ismember(site, kan_anchors) || strlength(site) == 0
         continue
      end
      loctype = lower(strtrim(string(T.location_type(n))));
      elev = NaN;
      if ismember('altitude_installation', T.Properties.VariableNames)
         elev = double(T.altitude_installation(n));
      end
      [zone, target] = firstPassZone(loctype, elev);

      entry = emptyEntry();
      entry.site = site;
      entry.alias = lower(erase(site, "_"));
      entry.long_name = site;
      entry.surface_zone = zone;
      entry.eval_target = target;
      entry.classification = "first_pass";
      entry.has_recipe = false;
      entry.models = models;
      entry.note = sprintf( ...
         "FIRST PASS (review): location_type=%s, elev=%s m.", ...
         loctype, num2str(elev));
      catalog(end + 1) = entry; %#ok<AGROW>
   end
end

function [zone, target] = firstPassZone(loctype, elev)
   %FIRSTPASSZONE Coarse first-pass zone+target from location_type + elevation.
   %
   % See the CLASSIFICATION PROVENANCE block in the function header for the full
   % mapping rationale. This is a heuristic, NOT authoritative.
   seasonal_ice = ["seasonal_snow"; "bare_ice"];
   seasonal_firn = ["seasonal_snow"; "firn"];

   switch loctype
      case "tundra"
         zone = "tundra";
         target = strings(0, 1);
         return
      case {"bedrock", "not greenland"}
         zone = "land";
         target = strings(0, 1);
         return
      case "local glacier"
         zone = "ablation";
         target = seasonal_ice;
         return
   end

   % "ice sheet" (and any other on-ice type): band by elevation.
   if ~isfinite(elev)
      zone = "unknown";
      target = strings(0, 1);
   elseif elev < 1500
      zone = "ablation";
      target = seasonal_ice;
   elseif elev < 2000
      zone = "percolation";
      target = seasonal_firn;
   else
      zone = "accumulation";
      target = seasonal_firn;
   end
end

function entry = emptyEntry()
   %EMPTYENTRY One catalog entry with the canonical field order.
   entry = struct('site', "", 'alias', "", 'long_name', "", ...
      'surface_zone', "", 'eval_target', strings(0, 1), ...
      'classification', "", 'has_recipe', false, 'models', strings(1, 0), ...
      'note', "");
end

function csv = locateAwsCsv()
   %LOCATEAWSCSV Resolve AWS_sites_metadata.csv under the staged product.
   csv = string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'promice', 'AWS_sites_metadata.csv'));
end
