function cache_dir = resolveFetchCacheDir(cache_dir, default_cache_dir)
   %RESOLVEFETCHCACHEDIR Apply the shared fetch-cache defaulting rule.
   %
   %  cache_dir = icemodel.verification.setup.resolveFetchCacheDir( ...
   %     cache_dir, default_cache_dir)
   %
   % Fetch helpers accept cache_dir="" from higher-level option plumbing. Treat a
   % blank value the same as an omitted value so callers can forward optional roots
   % without accidentally validating the current working directory.

   cache_dir = string(cache_dir);
   if cache_dir == ""
      cache_dir = string(default_cache_dir);
   end
end
