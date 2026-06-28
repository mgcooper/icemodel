function sites = imauSiteMetadata(site_ids)
   %IMAUSITEMETADATA Return first-pass IMAU hourly AWS site metadata.
   %
   %  sites = icemodel.verification.setup.imauSiteMetadata()
   %  sites = icemodel.verification.setup.imauSiteMetadata("S21")
   %
   % Role
   %  The first IMAU verification inventory is the hourly PANGAEA S21/S22/S23
   %  collection. Coordinates are left unset until parsed from source files so
   %  this metadata helper does not pretend to be a geodetic source of truth.

   arguments
      site_ids (1, :) string = strings(1, 0)
   end

   % Keep the hourly network separate from the daily 19-station SEB product; the
   % daily product is QA/provenance input, not a first-pass case list.
   sites = [ ...
      one("S21", "RetMIP FA meteorological source")
      one("S22", "")
      one("S23", "")];

   % A caller can request a subset; unknown ids should fail before writing a
   % manifest or staging a mislabeled case folder.
   if ~isempty(site_ids)
      keep = ismember([sites.site_id], site_ids);
      missing = setdiff(site_ids, [sites.site_id], 'stable');
      if ~isempty(missing)
         error('icemodel:verification:imauSiteMetadata:unknownSite', ...
            'unknown IMAU site id(s): %s', strjoin(missing, ', '));
      end
      sites = sites(keep);
   end
end

function s = one(site_id, note)
   %ONE Build one metadata row with a stable association field.
   assoc = struct('family', "", 'source_id', "", 'relationship', "");
   if site_id == "S21"
      assoc = struct('family', "retmip", 'source_id', "fa", ...
         'relationship', "meteorological_source_for");
   end

   s = struct( ...
      'site_id', site_id, ...
      'case_id', lower(site_id), ...
      'site_name', "IMAU " + site_id, ...
      'surface_zone', "unknown", ...
      'eval_target', "firn", ...
      'permafrost_zone', "none", ...
      'lat_wgs84', NaN, ...
      'lon_wgs84', NaN, ...
      'x_epsg3413', NaN, ...
      'y_epsg3413', NaN, ...
      'elev_m', NaN, ...
      'source_association', assoc, ...
      'note', note);
end
