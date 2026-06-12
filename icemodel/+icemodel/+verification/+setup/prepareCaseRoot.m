function prepareCaseRoot(case_root, overwrite)
   %PREPARECASEROOT Create one setup case folder and guard existing artifacts.
   %
   %  icemodel.verification.setup.prepareCaseRoot(case_root, overwrite)
   %
   % Inputs
   %  case_root   Folder where one case's setup artifacts will be staged.
   %  overwrite   When false, existing files trigger an error. When true, the
   %              importer may replace the staged artifacts.
   %
   % Role
   %  Setup helper used only by import/update tooling. The overwrite flag is a
   %  legitimate setup-refresh control, not a normal verification workflow
   %  behavior.

   % Ensure the folder exists before checking contents so a first setup run can
   % proceed without a separate mkdir path.
   icemodel.helpers.ensureDirExists(case_root);
   if overwrite
      return
   end

   % Treat any existing non-dot entry as an accidental overwrite risk.
   files = dir(fullfile(case_root, '*'));
   names = string({files.name});
   names = names(~ismember(names, [".", ".."]));
   if ~isempty(names)
      error(['staged verification folder already exists for %s. ' ...
         'Rerun with overwrite=true to replace it.'], case_root)
   end
end
