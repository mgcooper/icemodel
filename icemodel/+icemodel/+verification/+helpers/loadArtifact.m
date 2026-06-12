function artifact = loadArtifact(pathname, varname)
   %LOADARTIFACT Load one named staged verification artifact from a MAT file.
   %
   %  targets = icemodel.verification.helpers.loadArtifact(path, "targets")
   %
   % Inputs
   %  pathname   MAT file containing a staged verification artifact.
   %  varname    Variable name to load from the MAT file.
   %
   % Outputs
   %  artifact   The loaded artifact value.
   %
   % Role
   %  Operational helper shared by comparecase, plotcase, and candidate
   %  resolution. It keeps MAT-file loading centralized and explicit.

   arguments
      pathname (1, :) string
      varname (1, :) string
   end

   % Load only the requested variable so accidental extra MAT-file contents do
   % not become part of the public verification contract.
   data = load(pathname, varname);
   artifact = data.(varname);
end
