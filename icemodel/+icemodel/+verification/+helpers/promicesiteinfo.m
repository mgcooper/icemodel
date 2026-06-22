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
   %    site_type       PROMICE/GC-Net AWS readme Table 1 "Site type"
   %                    (Ablation/Accumulation/Bedrock). This is the
   %                    AUTHORITATIVE surface-height-channel discriminator:
   %                    Ablation sites ship z_ice_surf + snow_height; the
   %                    others ship only z_surf_combined. buildPromiceData
   %                    branches on z_ice_surf presence, which agrees with
   %                    this field. NOTE site_type (data-product class) and
   %                    surface_zone (glaciological facies) are distinct:
   %                    KAN_U is site_type=Accumulation (no z_ice_surf) but
   %                    surface_zone=percolation (firn-core truth).
   %    stations        string array of the composing AWS for the site (readme
   %                    Table 1 / AWS_sites_metadata.csv `stations`). A site
   %                    merges multiple AWS over time, so a station transition
   %                    can produce an expected step-shift in the surface or
   %                    subsurface series (e.g. CEN ~2021, MIT early).
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
   %  surface_zone (PRIMARY signal: MODIS end-of-summer BARE-ICE EXTENT 2000-2018,
   %    /Volumes/S03/DATA/greenland/racmo2p3/annual/
   %      MODIS_Bare_Ice_Extent_2000-2018.nc, binary 0/1 annual mask on the native
   %    RACMO2.3p2 FGRN055 ~1 km grid). The end-of-summer snowline ~ the ELA, so
   %    this is the proper ablation-vs-accumulation discriminator. Per site we take
   %    the bare-ice FREQUENCY f_bare = N(year bare)/N(valid years), max over a 3x3
   %    native-cell window (so a margin site whose nearest cell lands off-ice does
   %    not falsely read 0):
   %      f_bare >= 0.50           -> ablation     (bare ice in a majority of years)
   %    In the accumulation area (f_bare < 0.50) the firn facies is refined with
   %    SUMup density co-location:
   %      SUMup_2025 density profile <= 15 km -> percolation (firn observed)
   %      otherwise                           -> accumulation (facies unresolved)
   %    (A former elev >= 2500 m & f_bare==0 -> dry_snow branch was removed: the
   %    elevation cutoff did not generalize, so the three former dry_snow sites
   %    EGP/NAE/SDM collapse to accumulation. "dry_snow" stays in the surfacezone
   %    vocabulary but is currently unused.)
   %    This method REPRODUCES the KAN anchors: KAN_L f_bare=1.00 -> ablation,
   %    KAN_M 1.00 -> ablation. KAN_U f_bare=0.00 reads snow-covered at the surface
   %    every year and the surface signal alone would call it accumulation; it is
   %    the lower PERCOLATION zone by firn-core ground truth (SUMup density 0.2 km),
   %    so KAN_U is PINNED to the user-authoritative "percolation" (consistent with
   %    the SUMup firn evidence; surface remote sensing cannot see subsurface melt).
   %    location_type tundra -> tundra ; bedrock/not-Greenland -> land ; local
   %    glacier -> ablation (marginal).
   %
   %  permafrost_zone (Obu et al. 2019 ESA GlobPermafrost / UiO PEX permafrost
   %    zones, EXTENT field, WGS84 polygons, /Volumes/S03/DATA/interface/permafrost/
   %    obu/UiO_PEX_PERZONES_5/wgs/UiO_PEX_PERZONES_5.0_20181128_2000_2016_NH.shp).
   %    Off-ice sites only; ice-sheet/glacier -> "none". Point-in-polygon on
   %    lon/lat; EXTENT Cont/Discon/Spora/Isol ->
   %    continuous/discontinuous/sporadic/isolated (matching
   %    activelayer.readobuzones' parsing); an off-ice site outside all permafrost
   %    polygons -> "none" (permafrost-free ground).
   %
   %    NOTE: the staged activelayer.readobuzones reader could not be invoked on
   %    this machine (missing helpers parseFileName/dealout/activate); the SAME Obu
   %    shapefile is read directly with shaperead applying that reader's exact
   %    EXTENT->zone mapping. Replaces the v1 Brown et al. (1997) source.
   %
   %  ZAC_A/L/U (A.P. Olsen / Zackenberg GlacioBasis transect, NE Greenland)
   %  carry no installation coordinates in AWS_sites_metadata.csv; their lon/lat
   %  are sourced from the per-station L3 NetCDF global latitude/longitude
   %  attributes (data/verification/promice/hour/ZAC_*_hour.nc, the same
   %  coordinates readPromiceAws reads). They are local-glacier AWS -> ablation;
   %  on glacier ice so permafrost_zone=none (the surrounding Zackenberg tundra
   %  is continuous permafrost, but the sites sit on ice). No site now remains
   %  "unknown".
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
     "KAN_L",  "ablation",      "si", "none",            true,  "Lower ablation zone (~679 m); curated KAN_L recipe. KAN anchor. MODIS bare-ice freq=1.00."
     "KAN_M",  "ablation",      "si", "none",            true,  "Upper ablation / bare ice (~1272 m); curated KAN_M recipe. KAN anchor. MODIS bare-ice freq=1.00."
     "KAN_U",  "percolation",   "sf", "none",            false, "Lower percolation zone (~1845 m); KAN anchor. MODIS bare-ice freq=0.00 (snow-covered surface); percolation by firn-core truth, consistent with SUMup density 0.2 km."
     "CEN",    "percolation",   "sf", "none",            false, "Percolation (MODIS bare-ice freq=0.00; SUMup density 0.0 km, firn observed; elev 1872 m)."
     "CP1",    "percolation",   "sf", "none",            false, "Percolation (MODIS bare-ice freq=0.00; SUMup density 0.3 km, firn observed; elev 1951 m)."
     "DY2",    "percolation",   "sf", "none",            false, "Percolation (MODIS bare-ice freq=0.00; SUMup density 0.2 km, firn observed; elev 2113 m)."
     "EGP",    "accumulation",  "sf", "none",            false, "Accumulation interior (MODIS bare-ice freq=0.00, elev 2663 m; high-interior site, no SUMup firn co-location)."
     "FRE",    "ablation",      "si", "none",            false, "Marginal local glacier -> ablation."
     "HUM",    "percolation",   "sf", "none",            false, "Percolation (MODIS bare-ice freq=0.00; SUMup density 0.1 km, firn observed; elev 1967 m)."
     "JAR",    "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "KAN_B",  "tundra",        "",   "continuous",      false, "Off-ice tundra; Obu EXTENT continuous."
     "KAN_T",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00, window-max recovered margin signal)."
     "KPC_L",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "KPC_U",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=0.68)."
     "LYN_L",  "ablation",      "si", "none",            false, "Marginal local glacier -> ablation."
     "LYN_T",  "ablation",      "si", "none",            false, "Marginal local glacier -> ablation."
     "MIT",    "ablation",      "si", "none",            false, "Marginal local glacier -> ablation."
     "NAE",    "accumulation",  "sf", "none",            false, "Accumulation interior (MODIS bare-ice freq=0.00, elev 2624 m; high-interior site, no SUMup firn co-location)."
     "NAU",    "percolation",   "sf", "none",            false, "Percolation (MODIS bare-ice freq=0.00; SUMup density 0.2 km, firn observed; elev 2335 m < 2500)."
     "NEM",    "percolation",   "sf", "none",            false, "Percolation (MODIS bare-ice freq=0.00; SUMup density 0.1 km, firn observed; elev 2451 m < 2500)."
     "NSE",    "percolation",   "sf", "none",            false, "Percolation (MODIS bare-ice freq=0.00; SUMup density 0.1 km, firn observed; elev 2375 m < 2500)."
     "NUK_B",  "tundra",        "",   "discontinuous",   false, "Off-ice tundra; Obu EXTENT discontinuous."
     "NUK_K",  "ablation",      "si", "none",            false, "Marginal local glacier -> ablation."
     "NUK_L",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "NUK_N",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "NUK_P",  "tundra",        "",   "sporadic",        false, "Off-ice tundra; Obu EXTENT sporadic."
     "NUK_U",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "ORO",    "land",          "",   "none",            false, "Off-ice (not Greenland); outside Obu permafrost polygons -> none."
     "QAS_A",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=0.68)."
     "QAS_L",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "QAS_M",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=0.95)."
     "QAS_U",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=0.68)."
     "RED_L",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00, window-max recovered margin signal)."
     "SCO_L",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "SCO_U",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "SDL",    "percolation",   "sf", "none",            false, "Percolation (MODIS bare-ice freq=0.00; SUMup density 0.0 km, firn observed; elev 2459 m < 2500)."
     "SDM",    "accumulation",  "sf", "none",            false, "Accumulation interior (MODIS bare-ice freq=0.00, elev 2879 m; high-interior site, no SUMup firn co-location)."
     "SER_B",  "land",          "",   "discontinuous",   false, "Off-ice bedrock; Obu EXTENT discontinuous."
     "SWC",    "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=0.84)."
     "TAS_A",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=0.95)."
     "TAS_L",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "TAS_U",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "THU_L",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "THU_L2", "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "THU_U",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "TUN",    "percolation",   "sf", "none",            false, "Percolation (MODIS bare-ice freq=0.00; SUMup density 0.1 km, firn observed; elev 2076 m)."
     "UPE_L",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=1.00)."
     "UPE_U",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=0.95)."
     "UWN",    "land",          "",   "isolated",        false, "Off-ice (not Greenland); Obu EXTENT isolated."
     "WEG_B",  "tundra",        "",   "continuous",      false, "Off-ice tundra; Obu EXTENT continuous."
     "WEG_L",  "ablation",      "si", "none",            false, "Ablation (MODIS bare-ice freq=0.95)."
     "ZAC_A",  "ablation",      "si", "none",            false, "A.P. Olsen / Zackenberg local glacier (NE Greenland, 74.65N -21.65E, 1481 m); GlacioBasis AWS. coords from L3 nc latitude/longitude attr (CSV installation blank). Marginal local glacier -> ablation. On glacier ice -> permafrost_zone none (surrounding tundra is continuous permafrost)."
     "ZAC_L",  "ablation",      "si", "none",            false, "A.P. Olsen / Zackenberg local glacier (NE Greenland, 74.62N -21.37E, 629 m); GlacioBasis AWS. coords from L3 nc latitude/longitude attr. Marginal local glacier -> ablation. On glacier ice -> permafrost_zone none (Obu EXTENT at point continuous)."
     "ZAC_U",  "ablation",      "si", "none",            false, "A.P. Olsen / Zackenberg local glacier (NE Greenland, 74.64N -21.46E, 862 m); GlacioBasis AWS. coords from L3 nc latitude/longitude attr. Marginal local glacier -> ablation. On glacier ice -> permafrost_zone none (Obu EXTENT at point continuous)."
   };

   typemap = siteTypeMap();
   stationmap = stationsMap();

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
         'site_type',      lookupSiteType(typemap, site), ...
         'stations',       lookupStations(stationmap, site), ...
         'surface_zone',   rows{k, 2}, ...
         'eval_target',    evalTarget(rows{k, 3}), ...
         'permafrost_zone',rows{k, 4}, ...
         'classification', classification, ...
         'has_recipe',     rows{k, 5}, ...
         'models',         models, ...
         'note',           rows{k, 6});
   end
end

function st = lookupSiteType(typemap, site)
   %LOOKUPSITETYPE Readme Table 1 "Site type" for a site, "unknown" if absent.
   if isKey(typemap, site)
      st = typemap(site);
   else
      st = "unknown";
   end
end

function s = lookupStations(stationmap, site)
   %LOOKUPSTATIONS Composing AWS for a site (from AWS_sites_metadata.csv).
   if isKey(stationmap, site)
      s = stationmap{site};
   else
      s = site;   % single-station fallback (site == station)
   end
end

function typemap = siteTypeMap()
   %SITETYPEMAP Readme Table 1 "Site type" (Ablation/Accumulation/Bedrock).
   %
   % Transcribed verbatim from the PROMICE/GC-Net AWS readme Table 1
   % (data/verification/promice/AWS_data_readme.pdf). This is the data-product
   % site class and the authoritative surface-height-channel discriminator
   % (Ablation -> z_ice_surf + snow_height; Accumulation/Bedrock ->
   % z_surf_combined only).
   pairs = {
     "CEN","Accumulation"; "CP1","Accumulation"; "DY2","Accumulation"
     "EGP","Accumulation"; "FRE","Ablation";     "HUM","Accumulation"
     "JAR","Ablation";     "KAN_B","Bedrock";    "KAN_L","Ablation"
     "KAN_M","Ablation";   "KAN_T","Ablation";   "KAN_U","Accumulation"
     "KPC_L","Ablation";   "KPC_U","Ablation";   "LYN_L","Ablation"
     "LYN_T","Ablation";   "MIT","Ablation";     "NAE","Accumulation"
     "NAU","Accumulation"; "NEM","Accumulation"; "NSE","Accumulation"
     "NUK_B","Bedrock";    "NUK_K","Ablation";   "NUK_L","Ablation"
     "NUK_N","Ablation";   "NUK_U","Ablation";   "QAS_A","Ablation"
     "QAS_L","Ablation";   "QAS_M","Ablation";   "QAS_U","Ablation"
     "RED_L","Ablation";   "SCO_L","Ablation";   "SCO_U","Ablation"
     "SDL","Accumulation"; "SDM","Accumulation"; "SER_B","Bedrock"
     "SWC","Ablation";     "TAS_A","Ablation";   "TAS_L","Ablation"
     "TAS_U","Ablation";   "THU_L","Ablation";   "THU_L2","Ablation"
     "THU_U","Ablation";   "TUN","Accumulation"; "UPE_L","Ablation"
     "UPE_U","Ablation";   "WEG_B","Bedrock";    "WEG_L","Ablation"
     "ZAC_A","Ablation";   "ZAC_L","Ablation";   "ZAC_U","Ablation"
   };
   typemap = dictionary(string(pairs(:, 1)), string(pairs(:, 2)));
end

function stationmap = stationsMap()
   %STATIONSMAP site_id -> composing-AWS string array, from the sites CSV.
   stationmap = configureDictionary("string", "cell");
   csv = locateAwsCsv();
   if csv == "" || ~isfile(csv)
      return
   end
   T = readtable(csv, 'TextType', 'string');
   if ~all(ismember({'site_id', 'stations'}, T.Properties.VariableNames))
      return
   end
   for r = 1:height(T)
      sid = string(T.site_id(r));
      stns = strsplit(strtrim(string(T.stations(r))));
      stationmap(sid) = {string(stns)};
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
   entry.site_type = lookupSiteType(siteTypeMap(), site);
   entry.stations = lookupStations(stationsMap(), site);
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
      'site_type', "", 'stations', strings(1, 0), ...
      'surface_zone', "", 'eval_target', strings(0, 1), ...
      'permafrost_zone', "", 'classification', "", 'has_recipe', false, ...
      'models', strings(1, 0), 'note', "");
end

function csv = locateAwsCsv()
   %LOCATEAWSCSV Resolve AWS_sites_metadata.csv under the staged product.
   csv = string(fullfile(icemodel.internal.fullpath('data'), ...
      'verification', 'promice', 'AWS_sites_metadata.csv'));
end
