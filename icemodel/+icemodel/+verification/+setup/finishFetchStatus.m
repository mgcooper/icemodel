function source_dir = finishFetchStatus(cache_dir, status, kwargs)
   %FINISHFETCHSTATUS Apply the shared fetch status/strict/silent contract.
   %
   % Fetchers build dataset-specific status rows and retrieval banners, then
   % delegate the common completion, strict-error, and partial-cache return flow.

   arguments
      cache_dir (1, 1) string
      status (1, :) struct
      kwargs.strict (1, 1) logical = true
      kwargs.silent (1, 1) logical = false
      kwargs.error_id (1, 1) string = ...
         "icemodel:verification:finishFetchStatus:missingSources"
      kwargs.error_label (1, 1) string = "Source"
      kwargs.banner_callback = []
      kwargs.error_callback = []
   end

   missing = string({status(~[status.present]).product});
   if isempty(missing)
      source_dir = cache_dir;
      return
   end

   % Strict failures always print retrieval guidance; silent suppresses only
   % optional, non-strict cache probes.
   if (~kwargs.silent || kwargs.strict) && ~isempty(kwargs.banner_callback)
      kwargs.banner_callback(cache_dir, status, missing);
   end

   if kwargs.strict
      if ~isempty(kwargs.error_callback)
         kwargs.error_callback(cache_dir, status, missing);
      end
      error(kwargs.error_id, '%s source cache incomplete in %s. Missing: %s', ...
         kwargs.error_label, cache_dir, strjoin(missing, ', '));
   end

   source_dir = cache_dir;
end
