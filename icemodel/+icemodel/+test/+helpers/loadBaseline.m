function [baseline, meta] = loadBaseline(kind, kwargs)
   %LOADBASELINE Load a rolling or release baseline table.
   %
   %  baseline = icemodel.test.helpers.loadBaseline("perf")
   %  baseline = icemodel.test.helpers.loadBaseline("regression")
   %  [baseline, meta] = icemodel.test.helpers.loadBaseline("perf", ...
   %     smbmodel="skinmodel", baseline_tag="v1.1")
   %  [baseline, meta] = icemodel.test.helpers.loadBaseline("perf", ...
   %     smbmodel="all")
   %  baseline = icemodel.test.helpers.loadBaseline("perf", ...
   %     filename="/path/to/file.mat")

   arguments
      kind (1, :) string {mustBeMember( ...
         kind, ["perf", "regression"])}

      kwargs.smbmodel (1, :) string ...
         = "all"

      kwargs.baseline_tag (1, :) string ...
         = "rolling"

      kwargs.simyear double ...
         = 2016

      kwargs.filename string {mustBeTextScalarOrEmpty} ...
         = string.empty()
   end

   [smbmodel, baseline_tag, simyear, filename] = deal( ...
      kwargs.smbmodel, kwargs.baseline_tag, kwargs.simyear, kwargs.filename);

   meta = struct();

   % Resolve the baseline selector.
   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector(baseline_tag);

   % If the caller provided an explicit file, load it directly.
   if ~isblanktext(filename)
      [baseline, meta] = loadOneFile(kind, char(filename), simyear);
      baseline = backfillMetadata( ...
         baseline, baseline_type, baseline_tag, smbmodel);
      icemodel.test.helpers.printFilePath(filename, "load");
      return
   end

   % smbmodel="all" is virtual: load and concatenate the per-model files.
   if smbmodel == "all"
      baseline = loadAllModels(kind, baseline_type, baseline_tag, simyear);
      baseline = backfillMetadata( ...
         baseline, baseline_type, baseline_tag, smbmodel);
      return
   end

   % Resolve the canonical baseline file path.
   pathname = icemodel.test.helpers.baselineFilePath(kind, ...
      smbmodel=smbmodel, baseline_type=baseline_type, ...
      baseline_tag=baseline_tag, simyear=simyear);

   [baseline, meta] = loadOneFile(kind, char(pathname), simyear);
   baseline = backfillMetadata( ...
      baseline, baseline_type, baseline_tag, smbmodel);

   % Print the loaded filename to the console.
   icemodel.test.helpers.printFilePath(pathname, "load");
end

function [baseline, meta] = loadOneFile(kind, pathname, simyear)
   %LOADONEFILE Load and normalize a single baseline MAT file.

   meta = struct();

   if exist(pathname, 'file') ~= 2
      baseline = table();
      return
   end

   % Load the meta struct if present (perf baselines store it).
   if any(string({whos('-file', pathname).name}) == "meta")
      S = load(pathname, 'meta');
      meta = normalizeMeta(S.meta);
   end

   % Map kind to the expected saved variable names.
   switch kind
      case "perf"
         varnames = ["PerfBaseline", "baseline"];
      case "regression"
         varnames = ["RegressionBaseline", "baseline"];
   end

   baseline = icemodel.test.helpers.loadSavedTable(pathname, varnames);

   % Perf baselines are filtered to the requested simulation year.
   if kind == "perf" && ~isempty(baseline) ...
         && ismember('simyear', baseline.Properties.VariableNames)
      baseline = baseline(baseline.simyear == simyear, :);
   end

   % Normalize shared columns.
   if ~isempty(baseline)
      baseline = normalizeColumns(baseline);
   end
end

function baseline = normalizeColumns(baseline)
   %NORMALIZECOLUMNS Normalize column types shared across baseline kinds.

   if ismember('case_id', baseline.Properties.VariableNames)
      baseline.case_id = ...
         icemodel.test.helpers.normalizeFormalCaseId(baseline.case_id);
   end

   if ismember('last_updated_utc', baseline.Properties.VariableNames) ...
         && ~isdatetime(baseline.last_updated_utc)
      baseline.last_updated_utc = datetime( ...
         baseline.last_updated_utc, 'ConvertFrom', 'datenum', ...
         'TimeZone', 'UTC');
   end
end

function baseline = backfillMetadata( ...
      baseline, baseline_type, baseline_tag, smbmodel)
   %BACKFILLMETADATA Backfill selector metadata missing from older baselines.

   if isempty(baseline)
      return
   end

   if ~ismember('baseline_type', baseline.Properties.VariableNames)
      baseline.baseline_type = repmat( ...
         baseline_type, height(baseline), 1);
   end
   baseline.baseline_type = string(baseline.baseline_type);

   if ~ismember('baseline_tag', baseline.Properties.VariableNames)
      baseline.baseline_tag = repmat( ...
         baseline_tag, height(baseline), 1);
   end
   baseline.baseline_tag = string(baseline.baseline_tag);

   if ~ismember('smbmodel_filter', baseline.Properties.VariableNames)
      baseline.smbmodel_filter = repmat( ...
         string(smbmodel), height(baseline), 1);
   end
   baseline.smbmodel_filter = string(baseline.smbmodel_filter);
end

function baseline = loadAllModels(kind, baseline_type, baseline_tag, simyear)
   %LOADALLMODELS Concatenate per-model baselines for smbmodel="all".

   models = icemodel.namelists.smbmodel("test");

   pathnames = arrayfun(@(mdl) icemodel.test.helpers.baselineFilePath( ...
      kind, smbmodel=mdl, baseline_type=baseline_type, ...
      baseline_tag=baseline_tag, simyear=simyear), ...
      models, 'UniformOutput', false);

   exists = cellfun(@(p) exist(char(p), 'file') == 2, pathnames);
   if ~any(exists)
      baseline = table();
      return
   end

   tables = cell(1, nnz(exists));
   idx = find(exists);
   for i = 1:numel(idx)

      % Load one table and backfill the metadata table.
      [tbl, ~] = loadOneFile(kind, char(pathnames{idx(i)}), simyear);
      tbl = backfillMetadata(tbl, baseline_type, baseline_tag, models(idx(i)));
      tables{i} = tbl;

      % Print the loaded filename to the console.
      icemodel.test.helpers.printFilePath(pathnames{idx(i)}, "load");
   end
   baseline = vertcat(tables{:});
end

function meta = normalizeMeta(meta)
   %NORMALIZEMETA Normalize saved baseline metadata into current string types.

   if ~isstruct(meta)
      meta = struct();
      return
   end

   fields = ["baseline_type", "baseline_tag", "smbmodel_filter", ...
      "matlab_version", "host"];
   for i = 1:numel(fields)
      if isfield(meta, fields(i))
         meta.(fields(i)) = string(meta.(fields(i)));
      end
   end
end
