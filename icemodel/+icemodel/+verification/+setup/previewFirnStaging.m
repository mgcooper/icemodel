function manifests = previewFirnStaging(families, kwargs)
   %PREVIEWFIRNSTAGING Stage short build_forcing=true firn QA previews.
   %
   %  manifests = icemodel.verification.setup.previewFirnStaging()
   %  manifests = icemodel.verification.setup.previewFirnStaging("retmip")
   %  manifests = icemodel.verification.setup.previewFirnStaging( ...
   %     "all", mar_dir=mar_dir, merra_dir=merra_dir, racmo_dir=racmo_dir)
   %
   % Role
   %  User-facing setup helper for preflight QA. It stages representative,
   %  one-year-or-shorter previews with build_forcing=true into a separate
   %  output root so users can visualize observations, native met/userdata, and
   %  RCM legs before running the full final staging workflow into data/eval and
   %  data/input. It calls the production importers; there is no separate scratch
   %  conversion path. It is an optional preflight/troubleshooting workflow, not
   %  a repair, promotion, or substitute for a canonical full stage. Once a full
   %  tree has passed QA, use it again only to isolate a source/importer change.
   %
   % Defaults
   %  output_root defaults to <repo>/data/preview/firn_staging. The default
   %  representative cases are PROMICE KAN_M, RetMIP dye2_2016, IMAU S21,
   %  research_site Humphrey, and SUMup at staged mixed anchors (or Humphrey
   %  when SUMup is run alone). forcing_sources selects the RCM runtime sources
   %  requested from every importer. Pass source/cache kwargs to override them.
   %  Model met defaults to dt_out="15m"; pass dt_out="" to retain each
   %  model-met source's native cadence. Userdata defaults to hourly.
   %  overwrite=true also clears prior PNGs for the selected cases before plotting;
   %  all writes remain beneath output_root.
   %
   % See also: icemodel.verification.plotVerificationArtifacts,
   %  icemodel.verification.setup.importPromiceSites,
   %  icemodel.verification.setup.importRetmip,
   %  icemodel.verification.setup.importImau,
   %  icemodel.verification.setup.importResearchSites,
   %  icemodel.verification.setup.importSumup

   arguments
      families (1, :) string = "all"
      kwargs.output_root (1, 1) string = defaultPreviewRoot()
      kwargs.promice_dir (1, 1) string = ""
      kwargs.retmip_dir (1, 1) string = ""
      kwargs.imau_dir (1, 1) string = ""
      kwargs.sumup_dir (1, 1) string = ""
      kwargs.gcnet_dir (1, 1) string = ""
      kwargs.samimi_dir (1, 1) string = ""
      kwargs.mar_dir (1, 1) string = ""
      kwargs.merra_dir (1, 1) string = ""
      kwargs.racmo_dir (1, 1) string = ""
      kwargs.modis_dir (1, 1) string = ""
      kwargs.promice_sites (1, :) string = "KAN_M"
      kwargs.retmip_case_ids (1, :) string = "dye2_16"
      kwargs.imau_site_ids (1, :) string = "S21"
      kwargs.research_site_ids (1, :) string = "humphrey"
      kwargs.forcing_sources (1, :) string ...
         {icemodel.verification.validators.mustBeRcmSourceSelection( ...
         kwargs.forcing_sources)} = ...
         icemodel.verification.namelists.rcmsources()
      kwargs.dt_out (1, 1) string ...
         {mustBeMember(kwargs.dt_out, ["", "15m"])} = "15m"
      kwargs.overwrite (1, 1) logical = true
      kwargs.overwrite_family (1, 1) logical = true
      kwargs.skip_missing (1, 1) logical = true
      kwargs.make_plots (1, 1) logical = true
      kwargs.figure_root (1, 1) string = ""
      kwargs.visible (1, 1) logical = false
   end

   validateFamilies(families);
   families = expandFamilies(families);
   icemodel.helpers.ensureDirExists(kwargs.output_root);
   manifests = struct();

   % Stage source-owning families before SUMup so the SUMup mixed-anchor
   % preview can discover the anchors already written into the preview eval
   % tree.
   ordered = orderedPreviewFamilies(families);
   for family = ordered
      manifests.(char(family)) = stageFamily(family, kwargs);
   end

   % Plot after staging so the figures always reflect the actual manifest state,
   % including skipped RCM legs and partial source caches.
   if kwargs.make_plots
      figure_root = kwargs.figure_root;
      if figure_root == ""
         figure_root = fullfile(kwargs.output_root, "figures");
      end
      icemodel.verification.plotVerificationArtifacts( ...
         dataset_family=families, output_root=kwargs.output_root, ...
         figure_root=figure_root, save_figs=true, ...
         overwrite=kwargs.overwrite, visible=kwargs.visible);
   end
end

function validateFamilies(families)
   %VALIDATEFAMILIES Accept all or a known firn dataset family selector.
   mustBeMember(families, ["all", ...
      icemodel.verification.namelists.firndatasetfamily()]);
end

function root = defaultPreviewRoot()
   %DEFAULTPREVIEWROOT Gitignored preview tree separate from final data/input.
   root = string(fullfile(icemodel.internal.fullpath('data'), ...
      'preview', 'firn_staging'));
end

function ordered = orderedPreviewFamilies(families)
   %ORDEREDPREVIEWFAMILIES Keep SUMup after source-owning anchor families.
   source_families = setdiff(icemodel.verification.namelists.firndatasetfamily(), ...
      "sumup", 'stable');
   ordered = [intersect(source_families, families, 'stable'), ...
      intersect("sumup", families, 'stable')];
end

function families = expandFamilies(families)
   %EXPANDFAMILIES Replace "all" with the firn-development family list.
   if any(families == "all")
      families = icemodel.verification.namelists.firndatasetfamily();
   end
   families = unique(families, 'stable');
end

function manifest = stageFamily(family, kwargs)
   %STAGEFAMILY Dispatch one short-window family preview.
   rcm = rcmArgs(kwargs);
   switch family
      case "promice"
         window = oneYearWindow("2012-01-01 00:00:00", ...
            "2012-12-31 23:00:00");
         manifest = icemodel.verification.setup.importPromiceSites( ...
            kwargs.promice_dir, "sites", kwargs.promice_sites, ...
            "output_root", kwargs.output_root, "build_forcing", true, ...
            "forcing_sources", ["promice", kwargs.forcing_sources], ...
            "startdate", window.start, "enddate", window.end, ...
            "overwrite", kwargs.overwrite, ...
            "overwrite_family", kwargs.overwrite_family, ...
            "skip_missing", kwargs.skip_missing, rcm{:});

      case "retmip"
         manifest = stageRetmipPreviewCases(kwargs, rcm);

      case "imau"
         window = oneYearWindow("2014-04-12 00:00:00", ...
            "2014-12-02 21:00:00");
         manifest = icemodel.verification.setup.importImau( ...
            kwargs.imau_dir, "site_ids", kwargs.imau_site_ids, ...
            "output_root", kwargs.output_root, "dry_run", false, ...
            "build_forcing", true, ...
            "forcing_sources", kwargs.forcing_sources, ...
            "startdate", window.start, ...
            "enddate", window.end, ...
            "overwrite", kwargs.overwrite, ...
            "overwrite_family", kwargs.overwrite_family, ...
            "skip_missing", kwargs.skip_missing, rcm{:});

      case "research_site"
         window = oneYearWindow("2011-07-01 00:00:00", ...
            "2012-04-30 23:00:00");
         manifest = icemodel.verification.setup.importResearchSites( ...
            kwargs.sumup_dir, "site_ids", kwargs.research_site_ids, ...
            "output_root", kwargs.output_root, "build_forcing", true, ...
            "forcing_sources", kwargs.forcing_sources, ...
            "startdate", window.start, "enddate", window.end, ...
            "overwrite", kwargs.overwrite, ...
            "overwrite_family", kwargs.overwrite_family, ...
            "skip_missing", kwargs.skip_missing, rcm{:});

      case "sumup"
         [points, anchors, case_ids] = sumupPreviewPoints(kwargs);
         window = oneYearWindow("2011-07-01 00:00:00", ...
            "2012-04-30 23:00:00");
         manifest = icemodel.verification.setup.importSumup( ...
            kwargs.sumup_dir, "points", points, "anchors", anchors, ...
            "case_ids", case_ids, "output_root", kwargs.output_root, ...
            "build_forcing", true, "startdate", window.start, ...
            "forcing_sources", kwargs.forcing_sources, ...
            "enddate", window.end, ...
            "overwrite", kwargs.overwrite, ...
            "overwrite_family", kwargs.overwrite_family, ...
            "skip_missing", kwargs.skip_missing, rcm{:});
   end
end

function args = rcmArgs(kwargs)
   %RCMARGS Common RCM source directory kwargs for every importer.
   args = {"mar_dir", kwargs.mar_dir, ...
      "merra_dir", kwargs.merra_dir, ...
      "racmo_dir", kwargs.racmo_dir, ...
      "modis_dir", kwargs.modis_dir, ...
      "dt_out", kwargs.dt_out};
end

function manifest = stageRetmipPreviewCases(kwargs, rcm)
   %STAGERETMIPPREVIEWCASES Stage each RetMIP case with its own preview window.
   cases = icemodel.verification.setup.retmipCaseCatalog( ...
      kwargs.retmip_case_ids);
   manifest = struct();
   for k = 1:numel(cases)
      window = oneYearWindow(cases(k).period.start, cases(k).period.end);
      manifest = icemodel.verification.setup.importRetmip( ...
         kwargs.retmip_dir, "case_ids", string(cases(k).case_id), ...
         "output_root", kwargs.output_root, "dry_run", false, ...
         "build_forcing", true, "startdate", window.start, ...
         "forcing_sources", kwargs.forcing_sources, ...
         "enddate", window.end, ...
         "promice_dir", kwargs.promice_dir, ...
         "gcnet_dir", kwargs.gcnet_dir, "samimi_dir", kwargs.samimi_dir, ...
         "imau_dir", kwargs.imau_dir, ...
         "overwrite", kwargs.overwrite, ...
         "overwrite_family", kwargs.overwrite_family && k == 1, ...
         "skip_missing", kwargs.skip_missing, rcm{:});
   end
   manifest_file = fullfile(kwargs.output_root, 'eval', 'retmip', ...
      'manifest.json');
   if isfile(manifest_file)
      manifest = icemodel.verification.helpers.readFamilyManifest( ...
         manifest_file);
   end
end

function window = oneYearWindow(start_text, end_text)
   %ONEYEARWINDOW Clamp a case window to at most one year.
   t1 = icemodel.verification.setup.ensureUtc(start_text);
   t2 = icemodel.verification.setup.ensureUtc(end_text);
   tmax = t1 + calyears(1) - hours(1);
   if t2 > tmax
      t2 = tmax;
   end
   window = struct('start', ...
      icemodel.verification.setup.formatManifestTime(t1), ...
      'end', icemodel.verification.setup.formatManifestTime(t2));
end

function [points, anchors, case_ids] = sumupPreviewPoints(kwargs)
   %SUMUPPREVIEWPOINTS Prefer staged mixed anchors, otherwise Humphrey alone.
   anchors = icemodel.verification.setup.mixedAnchorCatalog( ...
      output_root=kwargs.output_root);
   points = zeros(0, 2);
   case_ids = strings(1, 0);
   if ~isempty(anchors)
      return
   end

   site = icemodel.verification.setup.researchSiteCatalog("humphrey");
   points = [site.lat_wgs84, site.lon_wgs84];
   anchors = researchSiteAnchor(site);
   case_ids = "humphrey";
end

function anchor = researchSiteAnchor(site)
   %RESEARCHSITEANCHOR Minimal Humphrey anchor for standalone SUMup previews.
   proj = icemodel.forcing.helpers.psnProjection();
   [x3413, y3413] = projfwd(proj, site.lat_wgs84, site.lon_wgs84);
   anchor = struct('family', "research_site", ...
      'site', string(site.site_id), ...
      'source_id', string(site.source_id), ...
      'case_id', string(site.case_id), ...
      'lat_wgs84', site.lat_wgs84, ...
      'lon_wgs84', site.lon_wgs84, ...
      'x_epsg3413', x3413, ...
      'y_epsg3413', y3413, ...
      'elev_m', site.elev_m);
end
