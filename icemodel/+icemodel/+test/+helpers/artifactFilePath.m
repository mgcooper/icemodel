function pathname = artifactFilePath(kind, kwargs)
   %ARTIFACTFILEPATH Return the canonical artifact file path.
   %
   %  pathname = icemodel.test.helpers.artifactFilePath("perf")
   %  pathname = icemodel.test.helpers.artifactFilePath("perf", ...
   %     run_name="20260324-080229")
   %  pathname = icemodel.test.helpers.artifactFilePath("regression", ...
   %     tier="smoke", smbmodel="icemodel", solver=2)
   %
   % If RUN_NAME is empty, the most recent artifact run containing a
   % matching file is resolved by scanning test/artifacts/.

   arguments
      kind (1, :) string {mustBeMember( ...
         kind, ["perf", "regression"])}

      kwargs.tier (1, :) string {mustBeMember( ...
         kwargs.tier, ["smoke", "full"])} ...
         = "smoke"

      kwargs.smbmodel (1, :) string ...
         = "icemodel"

      kwargs.solver double ...
         = []

      kwargs.baseline_type (1, :) string {mustBeMember( ...
         kwargs.baseline_type, ["rolling", "release"])} ...
         = "rolling"

      kwargs.baseline_tag string {mustBeTextScalarOrEmpty} ...
         = string.empty()

      kwargs.run_name string {mustBeTextScalarOrEmpty} ...
         = string.empty()
   end

   [tier, smbmodel, solver, baseline_type, baseline_tag, run_name] = deal( ...
      kwargs.tier, kwargs.smbmodel, kwargs.solver, ...
      kwargs.baseline_type, kwargs.baseline_tag, kwargs.run_name);

   % Resolve baseline selector: "rolling" or blank → rolling; else release.
   if ~isblanktext(baseline_tag)
      [baseline_type, baseline_tag] = ...
         icemodel.test.helpers.resolveBaselineSelector(baseline_tag);
   end

   testdir = icemodel.getpath('test');
   artifacts_dir = fullfile(testdir, 'artifacts');

   % Build the artifact filename stem from the label components.
   filename = buildArtifactFilename( ...
      kind, tier, smbmodel, solver, baseline_type, baseline_tag);

   if ~isblanktext(run_name)
      % Use the requested run directory. If solver was not specified, try
      % a glob match to find the actual file in that directory.
      if isempty(solver)
         pathname = findInDir( ...
            fullfile(artifacts_dir, char(run_name)), filename);
      else
         pathname = fullfile(artifacts_dir, char(run_name), filename);
      end
   else
      % Scan for the most recent run that contains this artifact. When
      % solver is unspecified, use a glob pattern for the solver portion.
      pattern = buildArtifactPattern( ...
         kind, tier, smbmodel, solver, baseline_type, baseline_tag);
      pathname = findLatestArtifact(artifacts_dir, pattern);
   end
end

function filename = buildArtifactFilename( ...
      kind, tier, smbmodel, solver, baseline_type, baseline_tag)
   %BUILDARTIFACTFILENAME Construct the artifact filename from components.

   model_label = artifactSmbmodelLabel(kind, smbmodel);
   solver_label = artifactSolverLabel(solver);

   switch kind
      case "perf"
         if baseline_type == "rolling"
            baseline_label = "vs_rolling";
         else
            baseline_label = "vs_" + ...
               icemodel.test.helpers.sanitizeTag(baseline_tag);
         end
         filename = sprintf('perf_results_%s%s_%s.mat', ...
            char(tier), char(model_label + solver_label), ...
            char(baseline_label));

      case "regression"
         if baseline_type == "rolling"
            baseline_label = "rolling";
         elseif isblanktext(baseline_tag)
            baseline_label = "nobaseline";
         else
            baseline_label = char( ...
               icemodel.test.helpers.sanitizeTag(baseline_tag));
         end
         filename = sprintf('regression_report_%s%s_%s.mat', ...
            char(tier), char(model_label + solver_label), ...
            baseline_label);
   end
end

function label = artifactSmbmodelLabel(kind, smbmodel)
   %ARTIFACTSMBMODELLABEL Format the smbmodel selector for artifact filenames.
   smbmodel = string(smbmodel);
   switch kind
      case "perf"
         label = "_" + icemodel.test.helpers.smbmodelTag(smbmodel);
      case "regression"
         if any(strcmpi(smbmodel, "all"))
            label = "";
         else
            label = "_" + icemodel.test.helpers.smbmodelTag(smbmodel);
         end
   end
end

function label = artifactSolverLabel(solver)
   %ARTIFACTSOLVERLABEL Format the solver filter for artifact filenames.
   if isempty(solver)
      label = "";
   else
      label = "_s" + join(string(solver), '-');
   end
end

function pattern = buildArtifactPattern( ...
      kind, tier, smbmodel, solver, baseline_type, baseline_tag)
   %BUILDARTIFACTPATTERN Build a glob pattern for artifact filename matching.
   %
   % When solver is empty, uses a wildcard (*) for the solver portion so
   % that artifacts with any solver label are matched.

   model_label = artifactSmbmodelLabel(kind, smbmodel);

   if isempty(solver)
      % Match any solver label (including no solver label).
      solver_part = "*";
   else
      solver_part = char(artifactSolverLabel(solver));
   end

   switch kind
      case "perf"
         if baseline_type == "rolling"
            baseline_label = "vs_rolling";
         else
            baseline_label = "vs_" + ...
               icemodel.test.helpers.sanitizeTag(baseline_tag);
         end
         pattern = sprintf('perf_results_%s%s%s_%s.mat', ...
            char(tier), char(model_label), solver_part, ...
            char(baseline_label));

      case "regression"
         if isblanktext(baseline_tag) && baseline_type == "rolling"
            baseline_label = "*";
         elseif isblanktext(baseline_tag)
            baseline_label = "nobaseline";
         else
            baseline_label = char( ...
               icemodel.test.helpers.sanitizeTag(baseline_tag));
         end
         pattern = sprintf('regression_report_%s%s%s_%s.mat', ...
            char(tier), char(model_label), solver_part, ...
            baseline_label);
   end
end

function pathname = findInDir(dirpath, filename)
   %FINDINDIR Find a file in a directory, trying exact match then glob.

   exact = fullfile(dirpath, filename);
   if exist(exact, 'file') == 2
      pathname = exact;
      return
   end

   % Try glob pattern: replace the filename with a pattern that wildcards
   % the solver portion.
   [~, stem, ext] = fileparts(filename);
   files = dir(fullfile(dirpath, [char(stem) '*' char(ext)]));
   if ~isempty(files)
      pathname = fullfile(files(1).folder, files(1).name);
      return
   end

   pathname = exact;
end

function pathname = findLatestArtifact(artifacts_dir, pattern)
   %FINDLATESTARTIFACT Find the most recent run containing a matching artifact.

   folders = dir(artifacts_dir);
   folders = folders([folders.isdir]);
   names = string({folders.name});

   % Keep only timestamped run directories (yyyymmdd-HHMMSS).
   valid = matches(names, regexpPattern('^\d{8}-\d{6}$'));
   names = sort(names(valid), 'descend');

   if isempty(names)
      error('icemodel:test:noArtifactRuns', ...
         'No artifact runs found in %s', artifacts_dir)
   end

   for i = 1:numel(names)
      rundir = fullfile(artifacts_dir, char(names(i)));
      files = dir(fullfile(rundir, pattern));
      if ~isempty(files)
         pathname = fullfile(files(1).folder, files(1).name);
         return
      end
   end

   error('icemodel:test:noMatchingArtifact', ...
      'No artifact run contains files matching: %s', pattern)
end
