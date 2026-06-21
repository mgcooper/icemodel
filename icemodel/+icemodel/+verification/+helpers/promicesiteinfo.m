function info = promicesiteinfo(site)
   %PROMICESITEINFO Catalog of PROMICE anchor sites for firn co-location.
   %
   %  info = icemodel.verification.helpers.promicesiteinfo()
   %  info = icemodel.verification.helpers.promicesiteinfo("KAN_M")
   %
   %  Returns the curated catalog of PROMICE automatic-weather-station
   %  anchor sites used to stage co-located multi-model firn-evaluation
   %  bundles. With no arguments, returns a struct array; with one site
   %  argument (canonical id "KAN_M" or compact alias "kanm") returns the
   %  matching scalar struct.
   %
   %  The catalog only enumerates the curated anchor sites - those with a
   %  station ablation recipe and/or a firn role established in the legacy
   %  GREENLAND processing. The staging driver
   %  icemodel.verification.setup.importPromiceSites generalizes to any
   %  PROMICE station in the v3 bundle (it reads lat/lon from the NetCDF
   %  metadata and falls back to the generic ablation recipe), so this
   %  catalog is advisory metadata, not an argument-validation gate.
   %
   %  Each entry has fields:
   %    site            canonical PROMICE station id ("KAN_M")
   %    alias           compact lowercase alias used in filenames ("kanm")
   %    long_name       display-friendly name
   %    zone            firn/ablation zone classification
   %    has_recipe      true when buildPromiceData carries a curated
   %                    service-window ablation recipe for the site
   %    models          co-located model set staged at the site, in the
   %                    canonical order [promice, mar, merra, racmo]
   %    note            short provenance note
   %
   %  Site coordinates are NOT stored here: they are read live from the v3
   %  NetCDF metadata by readPromiceAws (latitude / longitude variables)
   %  and converted to EPSG:3413 by the staging driver, so the committed
   %  catalog never drifts from the source files.
   %
   %  Anchor sites (the KAN transect on the southwest GrIS):
   %    KAN_L - lower ablation zone (~670 m), curated recipe
   %    KAN_M - upper ablation zone / bare ice (~1270 m), curated recipe
   %    KAN_U - lower percolation zone (~1840 m), generic recipe
   %
   %  Role
   %    Catalog metadata used by the firn staging driver and report tooling.
   %    Not an argument validator - the driver accepts any PROMICE station id.
   %
   % See also: icemodel.verification.setup.importPromiceSites,
   %  icemodel.forcing.readPromiceAws, icemodel.forcing.buildPromiceData,
   %  icemodel.verification.helpers.snowmipinfo

   arguments
      site (1, 1) string = ""
   end

   % Canonical co-located model order for every firn bundle. PROMICE is the
   % station anchor; MAR and MERRA are point met sources whose albedo is
   % swapped downstream; RACMO is a Data-only SMB/eval source (never met).
   models = ["promice", "mar", "merra", "racmo"];

   % Curated anchor catalog. Order is south-to-north up the KAN transect so
   % the staged sites read as an elevation transect.
   catalog = struct( ...
      'site',      {"KAN_L",          "KAN_M",                 "KAN_U"}, ...
      'alias',     {"kanl",           "kanm",                  "kanu"}, ...
      'long_name', {"KAN-L",          "KAN-M",                 "KAN-U"}, ...
      'zone',      {"lower_ablation", "upper_ablation",        "lower_percolation"}, ...
      'has_recipe', {true,            true,                    false}, ...
      'models',    {models,           models,                  models}, ...
      'note', { ...
      "Lower ablation zone (~670 m); curated KAN_L service-window recipe.", ...
      "Upper ablation zone / bare ice (~1270 m); curated KAN_M service-window recipe.", ...
      "Lower percolation zone (~1840 m); generic ablation recipe."});

   if site == ""
      info = catalog;
      return
   end

   wanted = lower(erase(site, "_"));
   match = string({catalog.alias}) == wanted;
   if ~any(match)
      valid = strjoin(string({catalog.site}), ', ');
      error('icemodel:verification:promicesiteinfo:unknown', ...
         ['unknown PROMICE anchor site "%s". Curated anchors: %s. ' ...
         '(Any v3 station id can still be staged via importPromiceSites.)'], ...
         site, valid);
   end
   info = catalog(match);
end
