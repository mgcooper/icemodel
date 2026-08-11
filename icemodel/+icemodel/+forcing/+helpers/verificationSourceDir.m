function source_dir = verificationSourceDir(source_dir, parts)
   %VERIFICATIONSOURCEDIR Resolve repo-local verification source data roots.
   %
   %  source_dir = icemodel.forcing.helpers.verificationSourceDir("", "imau")
   %  source_dir = icemodel.forcing.helpers.verificationSourceDir( ...
   %     "", ["retmip", "samimi"])
   %
   % Forcing builders read manually staged verification sources from the
   % top-level repo data tree. This helper holds that policy, so a change to
   % the staging root convention updates every builder default.

   arguments
      source_dir (1, 1) string
      parts (1, :) string
   end

   % An explicit caller root takes priority. A blank root resolves under
   % <repo>/data/verification.
   if source_dir == ""
      pieces = [{'verification'}, cellstr(parts)];
      source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
         pieces{:}));
   end
end
