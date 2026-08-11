function dirs = sourceSearchDirs(base, source)
   %SOURCESEARCHDIRS Candidate dirs for a staged file: per-source subfolder first.
   %
   %  dirs = icemodel.forcing.helpers.sourceSearchDirs(BASE, SOURCE)
   %
   % Returns the ordered candidate directories the runtime searches when
   % resolving a staged met/userdata file:
   %
   %   1. BASE/<SOURCE>  - the per-source subfolder the staging writers
   %      (icemodel.forcing.helpers.writemet / writeuserdata) create, so the
   %      flat BASE folder does not fill with files as staging grows.
   %   2. BASE           - the flat layout, kept as a fallback so committed
   %      flat fixtures and workspaces without subfolders still resolve.
   %
   % The rule is subfolder first: a staged subfolder file takes priority over a
   % flat file with the same name. This function defines that order for the
   % three runtime resolvers (icemodel.configureRun.resolveMetPaths,
   % icemodel.createMetFileNames.findEnclosingMetFile, and
   % icemodel.loadmet.resolveUserdataFile).
   %
   % Inputs
   %  base   - the flat input directory (e.g. input/met or input/userdata)
   %  source - the per-source subfolder key (opts.forcings or opts.userdata)
   %
   % Output
   %  dirs   - 1x2 cellstr {subfolder, flat} in search order
   %
   % See also: icemodel.forcing.helpers.findEnclosingWindowFile,
   %  icemodel.forcing.helpers.writemet, icemodel.forcing.helpers.writeuserdata

   base = char(base);
   dirs = {fullfile(base, char(source)), base};
end
