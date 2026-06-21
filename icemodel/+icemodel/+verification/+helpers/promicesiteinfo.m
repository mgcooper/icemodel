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
   %  The catalog carries THREE independent descriptors per site:
   %    surface_zone   the glaciological substrate REGIME the site sits in
   %                   (icemodel.verification.namelists.surfacezone vocabulary:
   %                   ablation/percolation/wet_snow/dry_snow/accumulation/
   %                   land/tundra/unknown). This is "where is the site".
   %    eval_target    a string ARRAY of which model CAPABILITIES the site
   %                   exercises (icemodel.verification.namelists.evaltarget
   %                   vocabulary: seasonal_snow/bare_ice/firn/ablation). This is
   %                   "what does a case here test". A site can have several.
   %    permafrost_zone the permafrost EXTENT class of the GROUND, ORTHOGONAL to
   %                   surface_zone (icemodel.verification.namelists.permafrostzone
   %                   vocabulary: continuous/discontinuous/sporadic/isolated/
   %                   none/unknown). Ice-sheet/glacier sites sit on ice, not
   %                   permafrost ground, so they carry "none"; off-ice land/
   %                   tundra sites carry the Brown et al. (1997) extent class.
   %
   %  Each entry has fields:
   %    site            canonical PROMICE station id ("KAN_M")
   %    alias           compact lowercase alias used in filenames ("kanm")
   %    long_name       display-friendly name
   %    surface_zone    glaciological zone (see above)
   %    eval_target     capability descriptor string array (see above)
   %    permafrost_zone permafrost extent class (see above)
   %    classification  "authoritative" for the data-derived classifications and
   %                    the curated KAN anchors, or "unknown" where no dataset
   %                    could resolve the site (e.g. missing coordinates).
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
   %  CLASSIFICATION PROVENANCE (AUTHORITATIVE - data-derived, hard-coded)
   %  ---------------------------------------------------------------------------
   %  The surface_zone and permafrost_zone values below are HARD-CODED results of
   %  spatially sampling three reference datasets at each site's installation
   %  lon/lat. The analysis tool is test/interactive/classify_site_facies.m (which
   %  requires /Volumes/S03); its results are baked in here so the committed
   %  catalog has NO S03 runtime dependency. Re-run that tool to refresh.
   %
   %  surface_zone (PRIMARY signal: MODIS MOD10A1 albedo facies, ISMIP6 1 km grid,
   %    /Volumes/S03/DATA/greenland/snowlines/MOD10A1_albedo_classify_<YR>.nc):
   %    Per-site multi-year bare-ice fraction f_bare = N(class=bare-ice) /
   %    N(class in {bare-ice, snow}) over the 92-slice melt-season composites.
   %      f_bare >= 0.40           -> ablation     (frequently snow-free)
   %      0.10 <= f_bare < 0.40    -> percolation  (mixed bare-ice / snow margin)
   %      f_bare <  0.10, elev<2500-> accumulation (snow all season; MODIS surface
   %                                                 albedo CANNOT see subsurface
   %                                                 percolation, so percolation is
   %                                                 NOT separable from accum here)
   %      f_bare <  0.10, elev>=2500 -> dry_snow   (cold high interior)
   %    This method REPRODUCES the KAN anchors: KAN_L f_bare=0.91 -> ablation,
   %    KAN_M 0.67 -> ablation. KAN_U f_bare=0.02 reads "snow" at the surface and
   %    the method alone would call it "accumulation"; it is the lower PERCOLATION
   %    zone by firn-core ground truth, so KAN_U is PINNED to the user-authoritative
   %    "percolation" (a known limit of surface remote sensing, not a method error).
   %    location_type tundra -> tundra ; bedrock/not-Greenland -> land ; local
   %    glacier -> ablation (marginal). MERRA-2 SMB (51x112 ~50 km grid) is too
   %    coarse to resolve the SW ablation margin (positive at KAN_L/M) so it is a
   %    diagnostic ONLY and never used. The 510x1400 ELAmask.mat carries no
   %    georeferencing and is unused.
   %
   %  permafrost_zone (Brown et al. 1997 circum-arctic permafrost, EXTENT code,
   %    /Volumes/S03/DATA/interface/GIS/Brown_Permafrost_Map/*.shp, Lambert
   %    Azimuthal Equal Area). Off-ice sites only; ice-sheet/glacier -> "none".
   %      C->continuous  D->discontinuous  S->sporadic  I->isolated.
   %
   %  Remaining "unknown": ZAC_A/L/U carry no coordinates in
   %  AWS_sites_metadata.csv, so no dataset could be sampled.
   %  ===========================================================================
   %
   %  Role
   %    Catalog metadata used by the firn staging driver and report tooling.
   %    Not an argument validator - the driver accepts any PROMICE station id.
   %    Uncataloged station ids fall back to firstPassZone() (a coarse
   %    elevation/location_type heuristic) flagged classification="first_pass".
   %
   % See also: icemodel.verification.setup.importPromiceSites,
   %  icemodel.verification.namelists.surfacezone,
   %  icemodel.verification.namelists.evaltarget,
   %  icemodel.verification.namelists.permafrostzone,
   %  icemodel.forcing.readPromiceAws, icemodel.forcing.buildPromiceData,
   %  icemodel.verification.helpers.snowmipinfo,
   %  test/interactive/classify_site_facies.m

   arguments
      site (1, 1) string = ""
   end

   % Canonical co-located model order. PROMICE is the station anchor; MAR and
   % MERRA are point met sources whose albedo is swapped downstream; RACMO is a
   % Data-only SMB/eval source (never met).
   models = ["promice", "mar", "merra", "racmo"];

   catalog = buildCatalog(models);

   if site == ""
      info = catalog;
      return
   end

   wanted = lower(erase(site, "_"));
   match = string({catalog.alias}) == wanted;
   if ~any(match)
      % Uncataloged station: fall back to the coarse first-pass heuristic so the
      % staging driver still gets a (flagged) zone for any L3 station id.
      info = firstPassEntry(site, models);
      return
   end
   info = catalog(match);
end

%% Local helpers
function catalog = buildCatalog(models)
   %BUILDCATALOG Hard-coded authoritative PROMICE catalog (see header provenance).
   %
   % Columns: site, alias, long_name, surface_zone, eval_target, permafrost_zone,
   % has_recipe, note. eval_target uses sentinel codes expanded by evalTarget():
   %   "si" = ["seasonal_snow";"bare_ice"]   (seasonal snow + bare ice)
   %   "sf" = ["seasonal_snow";"firn"]       (seasonal snow + firn)
   %   ""   = empty (off-ice land/tundra/unknown surfaces)
   % Order: KAN transect first (anchors), then the rest alphabetically.
   rows = {
   % site      surf_zone        et    pfz                recipe note
     "KAN_L",  "ablation",      "si", "none",            true,  "Lower ablation zone (~679 m); curated KAN_L recipe. KAN anchor. MODIS f_bare=0.91."
     "KAN_M",  "ablation",      "si", "none",            true,  "Upper ablation / bare ice (~1272 m); curated KAN_M recipe. KAN anchor. MODIS f_bare=0.67."
     "KAN_U",  "percolation",   "sf", "none",            false, "Lower percolation zone (~1845 m); KAN anchor. MODIS reads snow (f_bare=0.02); percolation by firn-core ground truth."
     "CEN",    "accumulation",  "sf", "none",            false, "Accumulation (MODIS f_bare=0.00, snow all season)."
     "CP1",    "accumulation",  "sf", "none",            false, "Accumulation (MODIS f_bare=0.00)."
     "DY2",    "accumulation",  "sf", "none",            false, "Accumulation (MODIS f_bare=0.02)."
     "EGP",    "dry_snow",      "sf", "none",            false, "Dry-snow interior (MODIS f_bare=0.00, elev 2663 m)."
     "FRE",    "ablation",      "si", "none",            false, "Marginal local glacier -> ablation (MODIS f_bare=0.70)."
     "HUM",    "accumulation",  "sf", "none",            false, "Accumulation (MODIS f_bare=0.00)."
     "JAR",    "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.70)."
     "KAN_B",  "tundra",        "",   "continuous",      false, "Off-ice tundra; Brown EXTENT continuous."
     "KAN_T",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.92)."
     "KPC_L",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.68)."
     "KPC_U",  "accumulation",  "sf", "none",            false, "Accumulation (MODIS f_bare=0.09)."
     "LYN_L",  "ablation",      "si", "none",            false, "Marginal local glacier -> ablation (MODIS f_bare=0.58)."
     "LYN_T",  "ablation",      "si", "none",            false, "Marginal local glacier -> ablation (MODIS f_bare=0.75)."
     "MIT",    "ablation",      "si", "none",            false, "Marginal local glacier -> ablation (MODIS f_bare=0.70)."
     "NAE",    "dry_snow",      "sf", "none",            false, "Dry-snow interior (MODIS f_bare=0.00, elev 2624 m)."
     "NAU",    "accumulation",  "sf", "none",            false, "Accumulation (MODIS f_bare=0.00)."
     "NEM",    "accumulation",  "sf", "none",            false, "Accumulation (MODIS f_bare=0.00)."
     "NSE",    "accumulation",  "sf", "none",            false, "Accumulation (MODIS f_bare=0.00)."
     "NUK_B",  "tundra",        "",   "isolated",        false, "Off-ice tundra; Brown EXTENT isolated."
     "NUK_K",  "ablation",      "si", "none",            false, "Marginal local glacier -> ablation (MODIS f_bare=0.81)."
     "NUK_L",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.90)."
     "NUK_N",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.83)."
     "NUK_P",  "tundra",        "",   "isolated",        false, "Off-ice tundra; Brown EXTENT isolated."
     "NUK_U",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.75)."
     "ORO",    "land",          "",   "sporadic",        false, "Off-ice (not Greenland); Brown EXTENT sporadic (nearest, low conf)."
     "QAS_A",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.48)."
     "QAS_L",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.78)."
     "QAS_M",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.73)."
     "QAS_U",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.55)."
     "RED_L",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.53)."
     "SCO_L",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.94)."
     "SCO_U",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.76)."
     "SDL",    "accumulation",  "sf", "none",            false, "Accumulation (MODIS f_bare=0.00)."
     "SDM",    "dry_snow",      "sf", "none",            false, "Dry-snow interior (MODIS f_bare=0.00, elev 2879 m)."
     "SER_B",  "land",          "",   "discontinuous",   false, "Off-ice bedrock; Brown EXTENT discontinuous."
     "SWC",    "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.43)."
     "TAS_A",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.51)."
     "TAS_L",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.83)."
     "TAS_U",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.70)."
     "THU_L",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.54)."
     "THU_L2", "percolation",   "sf", "none",            false, "Percolation margin (MODIS f_bare=0.38, mixed)."
     "THU_U",  "percolation",   "sf", "none",            false, "Percolation margin (MODIS f_bare=0.24, mixed)."
     "TUN",    "accumulation",  "sf", "none",            false, "Accumulation (MODIS f_bare=0.00)."
     "UPE_L",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.75)."
     "UPE_U",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.60)."
     "UWN",    "land",          "",   "sporadic",        false, "Off-ice (not Greenland); Brown EXTENT sporadic."
     "WEG_B",  "tundra",        "",   "continuous",      false, "Off-ice tundra; Brown EXTENT continuous."
     "WEG_L",  "ablation",      "si", "none",            false, "Ablation (MODIS f_bare=0.55)."
     "ZAC_A",  "unknown",       "",   "unknown",         false, "No coordinates in metadata; unresolved."
     "ZAC_L",  "unknown",       "",   "unknown",         false, "No coordinates in metadata; unresolved."
     "ZAC_U",  "unknown",       "",   "unknown",         false, "No coordinates in metadata; unresolved."
   };

   n = size(rows, 1);
   catalog = repmat(emptyEntry(), 1, n);
   for k = 1:n
      site = rows{k, 1};
      classification = "authoritative";
      if rows{k, 2} == "unknown"
         classification = "unknown";
      end
      catalog(k) = struct( ...
         'site',           site, ...
         'alias',          lower(erase(site, "_")), ...
         'long_name',      replace(site, "_", "-"), ...
         'surface_zone',   rows{k, 2}, ...
         'eval_target',    evalTarget(rows{k, 3}), ...
         'permafrost_zone',rows{k, 4}, ...
         'classification', classification, ...
         'has_recipe',     rows{k, 5}, ...
         'models',         models, ...
         'note',           rows{k, 6});
   end
end

function target = evalTarget(code)
   %EVALTARGET Expand a catalog eval_target sentinel code to a string array.
   switch code
      case "si"
         target = ["seasonal_snow"; "bare_ice"];
      case "sf"
         target = ["seasonal_snow"; "firn"];
      otherwise
         target = strings(0, 1);
   end
end

function entry = firstPassEntry(site, models)
   %FIRSTPASSENTRY Coarse first-pass entry for an UNCATALOGED station.
   %
   % Used only when a station id is not in the hard-coded authoritative catalog.
   % Reads location_type + elevation from AWS_sites_metadata.csv (if present) and
   % classifies with the legacy elevation-band heuristic, flagged "first_pass".
   loctype = "";
   elev = NaN;
   csv = locateAwsCsv();
   if csv ~= "" && isfile(csv)
      T = readtable(csv, 'TextType', 'string');
      row = string(T.site_id) == site;
      if any(row)
         if ismember('location_type', T.Properties.VariableNames)
            loctype = lower(strtrim(string(T.location_type(find(row, 1)))));
         end
         if ismember('altitude_installation', T.Properties.VariableNames)
            elev = double(T.altitude_installation(find(row, 1)));
         end
      end
   end
   [zone, target, pfz] = firstPassZone(loctype, elev);

   entry = emptyEntry();
   entry.site = site;
   entry.alias = lower(erase(site, "_"));
   entry.long_name = replace(site, "_", "-");
   entry.surface_zone = zone;
   entry.eval_target = target;
   entry.permafrost_zone = pfz;
   entry.classification = "first_pass";
   entry.has_recipe = false;
   entry.models = models;
   entry.note = sprintf( ...
      "FIRST PASS (uncataloged): location_type=%s, elev=%s m.", ...
      loctype, num2str(elev));
end

function [zone, target, pfz] = firstPassZone(loctype, elev)
   %FIRSTPASSZONE Coarse first-pass zone+target from location_type + elevation.
   %
   % Legacy elevation-band heuristic retained ONLY as a fallback for uncataloged
   % stations. NOT authoritative. permafrost_zone is left "unknown" (off-ice) or
   % "none" (on-ice) since the heuristic cannot sample Brown.
   seasonal_ice = ["seasonal_snow"; "bare_ice"];
   seasonal_firn = ["seasonal_snow"; "firn"];

   switch loctype
      case "tundra"
         zone = "tundra"; target = strings(0, 1); pfz = "unknown"; return
      case {"bedrock", "not greenland"}
         zone = "land"; target = strings(0, 1); pfz = "unknown"; return
      case "local glacier"
         zone = "ablation"; target = seasonal_ice; pfz = "none"; return
   end

   pfz = "none"; % "ice sheet" (and any other on-ice type): band by elevation.
   if ~isfinite(elev)
      zone = "unknown"; target = strings(0, 1); pfz = "unknown";
   elseif elev < 1500
      zone = "ablation"; target = seasonal_ice;
   elseif elev < 2000
      zone = "percolation"; target = seasonal_firn;
   else
      zone = "accumulation"; target = seasonal_firn;
   end
end

function entry = emptyEntry()
   %EMPTYENTRY One catalog entry with the canonical field order.
   entry = struct('site', "", 'alias', "", 'long_name', "", ...
      'surface_zone', "", 'eval_target', strings(0, 1), ...
      'permafrost_zone', "", 'classification', "", 'has_recipe', false, ...
      'models', strings(1, 0), 'note', "");
end

function csv = locateAwsCsv()
   %LOCATEAWSCSV Resolve AWS_sites_metadata.csv under the staged product.
   csv = string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'promice', 'AWS_sites_metadata.csv'));
end
