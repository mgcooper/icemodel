function [data, meta] = loadArtifact(kind, kwargs)
   %LOADARTIFACT Load a saved test artifact.
   %
   %  data = icemodel.test.helpers.loadArtifact("perf")
   %  [data, meta] = icemodel.test.helpers.loadArtifact("perf", ...
   %     smbmodel="skinmodel", run_name="20260324-080229")
   %  data = icemodel.test.helpers.loadArtifact("regression")
   %  data = icemodel.test.helpers.loadArtifact("perf", ...
   %     filename="/path/to/file.mat")
   %
   % Without RUN_NAME, the most recent artifact run is loaded. The loaded
   % file path is printed to the console.

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

      kwargs.filename string {mustBeTextScalarOrEmpty} ...
         = string.empty()
   end

   if ~isblanktext(kwargs.filename)
      pathname = char(kwargs.filename);
   else
      pathname = char(icemodel.test.helpers.artifactFilePath(kind, ...
         tier=kwargs.tier, smbmodel=kwargs.smbmodel, ...
         solver=kwargs.solver, baseline_type=kwargs.baseline_type, ...
         baseline_tag=kwargs.baseline_tag, run_name=kwargs.run_name));
   end

   if exist(pathname, 'file') ~= 2
      error('icemodel:test:artifactNotFound', ...
         'Artifact file not found: %s', pathname)
   end

   S = load(pathname);
   meta = struct();
   if isfield(S, 'meta')
      meta = S.meta;
   end

   switch kind
      case "perf"
         data = extractPerfData(S);
      case "regression"
         data = extractRegressionData(S);
   end

   % Print the loaded filename to the console.
   icemodel.test.helpers.printFilePath(pathname, "load");
end

function data = extractPerfData(S)
   %EXTRACTPERFDATA Extract perf artifact tables into a struct.

   data = struct();

   if isfield(S, 'case_summary')
      data.case_summary = S.case_summary;
   else
      data.case_summary = table();
   end

   if isfield(S, 'sample_detail')
      data.sample_detail = S.sample_detail;
   end

   if isfield(S, 'activity_detail')
      data.activity_detail = S.activity_detail;
   end

   if isfield(S, 'case_opts')
      data.case_opts = S.case_opts;
   end

   data.benchmark = struct('summary', table(), 'comparison', table(), ...
      'meta', struct());
   if isfield(S, 'benchmark')
      data.benchmark = S.benchmark;
   elseif isfield(S, 'benchmark_summary')
      data.benchmark.summary = S.benchmark_summary;
      if isfield(S, 'benchmark_comparison')
         data.benchmark.comparison = S.benchmark_comparison;
      end
      if isfield(S, 'benchmark_meta')
         data.benchmark.meta = S.benchmark_meta;
      end
   end
end

function data = extractRegressionData(S)
   %EXTRACTREGRESSIONDATA Extract regression artifact tables into a struct.

   data = struct();

   if isfield(S, 'report')
      data.report = S.report;
   elseif isfield(S, 'RegressionBaseline')
      data.report = S.RegressionBaseline;
   else
      data.report = table();
   end

   if isfield(S, 'case_opts')
      data.case_opts = S.case_opts;
   end
end
