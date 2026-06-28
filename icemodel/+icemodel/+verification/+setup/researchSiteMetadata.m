function sites = researchSiteMetadata(site_ids)
   %RESEARCHSITEMETADATA Return catchall research-site metadata.
   %
   %  sites = icemodel.verification.setup.researchSiteMetadata()
   %  sites = icemodel.verification.setup.researchSiteMetadata("humphrey")
   %
   % Role
   %  Catchall site catalog for high-value research targets that do not belong
   %  cleanly to a network family such as PROMICE, IMAU, or RetMIP. These sites
   %  are source-family anchors for colocation metadata, not a special SUMup
   %  subcategory.

   arguments
      site_ids (1, :) string = strings(1, 0)
   end

   % Start with Humphrey because it is the current standalone Meyer-Hewitt target
   % and should be represented by the generic research_site family.
   sites = one("humphrey", "Humphrey percolation thermistor network");

   % A caller can request a subset; fail before writing a mislabeled manifest.
   if ~isempty(site_ids)
      keep = ismember([sites.site_id], site_ids);
      missing = setdiff(site_ids, [sites.site_id], 'stable');
      if ~isempty(missing)
         error('icemodel:verification:researchSiteMetadata:unknownSite', ...
            'unknown research_site id(s): %s', strjoin(missing, ', '));
      end
      sites = sites(keep);
   end
end

function s = one(site_id, site_name)
   %ONE Build one research-site metadata row.
   % Humphrey is a multi-site thermistor network. Use the SUMup-weighted
   % centroid for nearest-anchor bookkeeping while retaining the network
   % provenance in the note.
   s = struct( ...
      'site_id', site_id, ...
      'case_id', site_id, ...
      'site_name', site_name, ...
      'family', "research_site", ...
      'source_id', site_id, ...
      'surface_zone', "unknown", ...
      'eval_target', "firn", ...
      'permafrost_zone', "unknown", ...
      'lat_wgs84', 69.725714, ...
      'lon_wgs84', -48.190512, ...
      'x_epsg3413', NaN, ...
      'y_epsg3413', NaN, ...
      'elev_m', 1645.1, ...
      'note', "SUMup Humphrey2012/Harper2011 thermistor-network centroid; use RCM forcing first.");
end
