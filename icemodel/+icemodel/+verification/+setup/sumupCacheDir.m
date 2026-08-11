function cache_dir = sumupCacheDir(source_dir)
   %SUMUPCACHEDIR Resolve the canonical SUMup verification source cache.
   %
   %  cache_dir = icemodel.verification.setup.sumupCacheDir("")
   %
   % A blank source root resolves to the repo top-level data/verification
   % tree. Source reads then do not depend on ICEMODEL_DATA_PATH.
   arguments
      source_dir (1, 1) string = ""
   end

   cache_dir = icemodel.forcing.helpers.verificationSourceDir(source_dir, ...
      "sumup");
end
