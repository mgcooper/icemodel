function filenames = locateKtransectFiles(source_dir, station)
   %LOCATEKTRANSECTFILES Find every cached K-transect annual file for a station.
   %
   %  filenames = icemodel.forcing.helpers.locateKtransectFiles( ...
   %     source_dir, "AWS9")
   %
   % Role
   %  Shared cache-layout resolver for the PANGAEA.947483 annual files, used by
   %  the Data builder and the importer preflight so both agree on what
   %  "present in the cache" means. The extracted archive layout (datasets/)
   %  is searched first, then flat and recursive fallbacks, mirroring the
   %  fetch-status glob order.
   %
   % See also: icemodel.forcing.buildKtransectData,
   %  icemodel.verification.setup.fetchKtransect

   arguments
      source_dir (1, 1) string
      station (1, 1) string
   end

   patterns = [ ...
      fullfile(source_dir, 'datasets', "K-transect_" + station + "_*.tab")
      fullfile(source_dir, "K-transect_" + station + "_*.tab")
      fullfile(source_dir, '**', "K-transect_" + station + "_*.tab")];
   for pattern = reshape(patterns, 1, [])
      hits = dir(pattern);
      hits = hits(~[hits.isdir]);
      if ~isempty(hits)
         filenames = string(fullfile({hits.folder}, {hits.name}));
         return
      end
   end
   error('icemodel:forcing:locateKtransectFiles:filesNotFound', ...
      'no K-transect annual files found for %s under %s', ...
      station, source_dir)
end
