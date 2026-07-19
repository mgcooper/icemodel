function added = dependencies(kwargs)
   %DEPENDENCIES Add external project dependencies to the MATLAB path.
   %
   %  added = icemodel.dependencies()
   %  added = icemodel.dependencies(require=true)
   %
   %% Description
   %
   %  icemodel.dependencies() resolves the external dev-repo toolboxes that a
   %  subset of icemodel workflows need and places them on the MATLAB path.
   %  These dependencies are intentionally kept OUT of the icemodel repo; they
   %  live in sibling dev repos that the user clones separately:
   %
   %     exactremap  - conservative (overlap-area-weighted) polygon remap used
   %                   by icemodel.forcing.helpers.remapPolygon. The toolbox/
   %                   subfolder is placed on the path. Repo:
   %                   https://github.com/mgcooper/exactremap
   %     activelayer - reads the Obu (UiO PEX) permafrost-zones shapefile via
   %                   activelayer.readobuzones, used by
   %                   icemodel.verification site classification. The toolbox/
   %                   subfolder is placed on the path. Repo:
   %                   https://github.com/mgcooper/activelayer
   %     matfunclib  - shared MATLAB helpers (parseFileName, dealout,
   %                   mustBePolygonOrFile, ...) that activelayer depends on.
   %                   Placed on the path before activelayer so its helper
   %                   dependencies resolve. Repo:
   %                   https://github.com/mgcooper/matfunclib
   %
   %  This is the single, central place that wires external dependencies. The
   %  test bootstrap (icemodel.test.helpers.bootstrapTestEnvironment) and any
   %  interactive setup should CALL this function rather than embedding their
   %  own dependency-path logic, so matfunclib and friends are never
   %  bootstrapped from arbitrary locations.
   %
   %% Path resolution
   %
   %  Each dependency root is resolved in this order:
   %
   %    1. A dependency-specific environment variable (highest priority):
   %          ICEMODEL_EXACTREMAP   -> the exactremap repo root
   %          ICEMODEL_ACTIVELAYER  -> the activelayer repo root
   %          ICEMODEL_MATFUNCLIB   -> the matfunclib repo root
   %
   %    2. A shared projects-root environment variable, joined with the repo
   %       name (e.g. ICEMODEL_PROJECTS_ROOT/exactremap):
   %          ICEMODEL_PROJECTS_ROOT
   %
   %    3. The sibling projects/ layout: the parent of the icemodel repo
   %       (projects/icemodel/.. -> projects/), joined with the repo name.
   %
   %  For exactremap and activelayer the toolbox/ subfolder is added (the repos
   %  ship their code under toolbox/); matfunclib is added at its repo root.
   %  Each dependency is a clean no-op when it is already on the path (e.g.
   %  activated by the user's startup.m) or when its repo is absent.
   %
   %% Input arguments
   %
   %  require - logical (default false). When true, error if any dependency
   %            cannot be resolved instead of skipping it silently. Use this
   %            when a workflow strictly requires the dependency.
   %
   %% Output arguments
   %
   %  added - string array of the toolbox folders that were newly added to the
   %          path by this call (empty when every dependency was already on the
   %          path or absent).
   %
   % See also: icemodel.config, icemodel.test.helpers.bootstrapTestEnvironment

   arguments
      kwargs.require (1, 1) logical = false
   end

   % The icemodel repo root is the anchor for the sibling projects/ fallback.
   rootdir = icemodel.internal.fullpath();
   projects_root = fileparts(rootdir);

   added = string.empty(1, 0);

   % matfunclib must go on first so activelayer's helper dependencies resolve
   % when activelayer is added below. Its sentinel is a representative helper.
   added = [added, ensureDependency('matfunclib', ...
      'ICEMODEL_MATFUNCLIB', '', projects_root, 'parseFileName', ...
      kwargs.require)];

   % exactremap supplies the conservative polygon remap. Its sentinel is the
   % top-level exactremap function.
   added = [added, ensureDependency('exactremap', ...
      'ICEMODEL_EXACTREMAP', 'toolbox', projects_root, 'exactremap', ...
      kwargs.require)];

   % activelayer supplies the Obu permafrost-zone reader. Its sentinel is the
   % packaged reader function.
   added = [added, ensureDependency('activelayer', ...
      'ICEMODEL_ACTIVELAYER', 'toolbox', projects_root, ...
      'activelayer.readobuzones', kwargs.require)];
end

function added = ensureDependency( ...
      reponame, envvar, subfolder, projects_root, sentinel, require)
   %ENSUREDEPENDENCY Resolve one dependency root and add it to the path.
   %
   % Adds the resolved toolbox folder to the path unless the dependency is
   % already resolvable (SENTINEL is on the path) or its folder is absent.
   % Returns the folder string when newly added, or an empty string array.

   added = string.empty(1, 0);

   % A clean no-op when the dependency is already on the path.
   if ~isempty(which(sentinel))
      return
   end

   % Resolve the repo root, then append the code subfolder (toolbox/ for the
   % toolbox-style repos, empty for matfunclib which ships code at its root).
   reporoot = resolveRepoRoot(reponame, envvar, projects_root);
   if isempty(subfolder)
      tbdir = reporoot;
   else
      tbdir = fullfile(reporoot, subfolder);
   end

   % Add the folder when present; otherwise skip (or error when required).
   if isfolder(tbdir)
      addpath(genpath(tbdir))
      added = string(tbdir);
   elseif require
      error('icemodel:dependencies:missing', ...
         ['the %s dependency could not be resolved at "%s"; clone the ' ...
         '%s dev repo or set the %s environment variable'], ...
         reponame, tbdir, reponame, envvar)
   end
end

function reporoot = resolveRepoRoot(reponame, envvar, projects_root)
   %RESOLVEREPOROOT Resolve a dependency repo root from env vars or siblings.
   %
   % Priority: dependency-specific env var, then ICEMODEL_PROJECTS_ROOT joined
   % with the repo name, then the sibling projects/ layout.

   % 1. Dependency-specific environment variable (highest priority).
   reporoot = getenv(envvar);
   if ~isempty(reporoot)
      return
   end

   % 2. Shared projects-root environment variable joined with the repo name.
   shared_root = getenv('ICEMODEL_PROJECTS_ROOT');
   if ~isempty(shared_root)
      reporoot = fullfile(shared_root, reponame);
      return
   end

   % 3. Sibling projects/ layout: the parent of the icemodel repo.
   reporoot = fullfile(projects_root, reponame);
end
