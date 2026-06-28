function manifest = importImau(source_dir, kwargs)
   %IMPORTIMAU Build the IMAU hourly AWS family manifest.
   %
   %  manifest = icemodel.verification.setup.importImau(source_dir)
   %  manifest = icemodel.verification.setup.importImau(source_dir, dry_run=true)
   %
   % Role
   %  First-pass IMAU staging hook for the hourly PANGAEA S21/S22/S23 network.
   %  The daily 19-station product remains QA/provenance input, not a first-pass
   %  case inventory.

   arguments
      source_dir (1, 1) string = ""
      kwargs.site_ids (1, :) string = strings(1, 0)
      kwargs.evaluation_data_root (1, 1) string = ""
      kwargs.icemodel_config_casename (1, 1) string = "test"
      kwargs.dry_run (1, 1) logical = true
      kwargs.strict (1, 1) logical = false
      kwargs.overwrite_family (1, 1) logical = false
   end

   % Validate caches but allow dry-run metadata on machines without PANGAEA files.
   [source_dir, cache_status] = icemodel.verification.setup.fetchImau( ...
      cache_dir=resolveCacheDir(source_dir), strict=kwargs.strict, silent=true);
   sites = icemodel.verification.setup.imauSiteMetadata(kwargs.site_ids);

   % Convert IMAU site rows to the existing firn manifest schema.
   entries = repmat(emptyEntry(), 1, numel(sites));
   for k = 1:numel(sites)
      entries(k) = manifestEntry(sites(k), source_dir, cache_status);
   end

   manifest = icemodel.verification.setup.makeFamilyManifest( ...
      "imau", "10.1594/PANGAEA.971647;10.1594/PANGAEA.970127", ...
      "https://doi.org/10.1594/PANGAEA.971647", ...
      "hourly-s21-s22-s23+daily-qa", string(datetime('today')), entries);

   if kwargs.dry_run
      return
   end

   % Write only on explicit non-dry runs; no empty observations.mat is fabricated.
   eval_root = icemodel.verification.helpers.evaluationDataRoot( ...
      "evaluation_data_root", kwargs.evaluation_data_root, ...
      "icemodel_config_casename", kwargs.icemodel_config_casename);
   family_root = fullfile(eval_root, "imau");
   icemodel.helpers.ensureDirExists(family_root);
   manifest_file = fullfile(family_root, "manifest.json");
   manifest = icemodel.verification.setup.writeFamilyManifestMerge( ...
      manifest_file, manifest, requested_ids=[entries.case_id], ...
      overwrite_family=kwargs.overwrite_family);
end

function cache_dir = resolveCacheDir(source_dir)
   %RESOLVECACHEDIR Resolve the IMAU cache directory.
   if strlength(source_dir) > 0
      cache_dir = source_dir;
   else
      cache_dir = string(fullfile(icemodel.getpath('data'), ...
         'verification', 'imau'));
   end
end

function entry = emptyEntry()
   %EMPTYENTRY Prototype IMAU manifest entry.
   values = { ...
      ''
      'firn_observational'
      ''
      ''
      'unknown'
      {'firn'}
      'none'
      struct('lat_wgs84', NaN, 'lon_wgs84', NaN, ...
         'x_epsg3413', NaN, 'y_epsg3413', NaN, 'elev_m', NaN)
      struct('start', '', 'end', '')
      ''
      {'imau'}
      {'imau_obs'}
      {'tair', 'rh', 'wspd', 'swd', 'lwd', 'precip'}
      struct()
      struct()
      '1hr'
      ''};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
end

function entry = manifestEntry(s, source_dir, cache_status)
   %MANIFESTENTRY Build one IMAU firn case entry.
   imau = struct();
   imau.kind = 'hourly_aws';
   imau.staged = false;
   imau.source_dir = char(source_dir);
   imau.cache_status = cache_status;

   colocation = struct();
   colocation.imau = imau;
   colocation.source_association = s.source_association;

   values = { ...
      char(s.case_id)
      'firn_observational'
      char(s.site_id)
      char(s.site_name)
      char(s.surface_zone)
      cellstr(string(s.eval_target))
      char(s.permafrost_zone)
      struct('lat_wgs84', s.lat_wgs84, 'lon_wgs84', s.lon_wgs84, ...
         'x_epsg3413', s.x_epsg3413, 'y_epsg3413', s.y_epsg3413, ...
         'elev_m', s.elev_m)
      struct('start', '', 'end', '')
      ''
      {'imau'}
      {'imau_obs'}
      {'tair', 'rh', 'wspd', 'swd', 'lwd', 'precip'}
      struct('hourly_site', s.site_id)
      colocation
      '1hr'
      sprintf('IMAU hourly AWS site %s. %s', s.site_id, s.note)};
   entry = icemodel.verification.setup.makeFirnCaseManifestEntry(values);
end
