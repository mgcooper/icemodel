function source_dir = verificationSourceDir(source_dir, parts)
   %VERIFICATIONSOURCEDIR Resolve repo-local verification source data roots.
   %
   %  source_dir = icemodel.forcing.helpers.verificationSourceDir("", "imau")
   %  source_dir = icemodel.forcing.helpers.verificationSourceDir( ...
   %     "", ["retmip", "samimi"])
   %
   % Forcing builders read manually staged verification sources from the top-level
   % repo data tree. Keeping this policy in one helper prevents builder-specific
   % defaults from drifting when the staging root convention changes.

   arguments
      source_dir (1, 1) string
      parts (1, :) string
   end

   % Explicit caller roots win. Blank roots resolve under <repo>/data/verification.
   if source_dir == ""
      pieces = [{'verification'}, cellstr(parts)];
      source_dir = string(fullfile(icemodel.internal.fullpath('data'), ...
         pieces{:}));
   end
end
