function sites = selectSiteCatalogEntries(sites, site_ids, error_id, item_label)
   %SELECTSITECATALOGENTRIES Select known ids from a source-site catalog.
   %
   %  sites = icemodel.verification.setup.selectSiteCatalogEntries( ...
   %     sites, site_ids, error_id, item_label)
   %
   % Inputs
   %  sites       Struct array with a string-valued site_id field.
   %  site_ids    Requested ids. An empty array selects the full catalog.
   %  error_id    Caller-owned identifier used for unknown ids.
   %  item_label  Family-specific label used in the error message.
   %
   % Outputs
   %  sites       Selected entries in canonical catalog order.

   arguments
      sites (1, :) struct
      site_ids (1, :) string
      error_id (1, 1) string
      item_label (1, 1) string
   end

   % Empty selection means the complete catalog, matching every family API.
   if isempty(site_ids)
      return
   end

   % Fail before any caller can stage a mislabeled case, then preserve the
   % catalog's canonical order rather than the request order.
   available_ids = string({sites.site_id});
   missing = setdiff(site_ids, available_ids, 'stable');
   if ~isempty(missing)
      error(char(error_id), 'unknown %s id(s): %s', item_label, ...
         strjoin(missing, ', '));
   end
   sites = sites(ismember(available_ids, site_ids));
end
