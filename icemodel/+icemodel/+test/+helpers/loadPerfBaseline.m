function [baseline, meta] = loadPerfBaseline(simyear, baseline_tag, smbmodel, pathname)
   %LOADPERFBASELINE Load rolling or release performance baseline table.
   %
   %  baseline = icemodel.test.helpers.loadPerfBaseline(2016)
   %  [baseline, meta] = icemodel.test.helpers.loadPerfBaseline(2016, ...
   %     "rolling", "all")
   %  [baseline, meta] = icemodel.test.helpers.loadPerfBaseline(2016, ...
   %     "v1.1", "skinmodel")
   %  [baseline, meta] = icemodel.test.helpers.loadPerfBaseline(2016, ...
   %     "v1.1", "skinmodel", pathname)
   arguments
      simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      baseline_tag (1, :) string = "rolling"
      smbmodel (1, :) string = "all"
      pathname (1, :) string = ""
   end

   meta = struct();

   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector(baseline_tag);

   % smbmodel="all" is virtual: load and concatenate the per-model files.
   if smbmodel == "all" && isblanktext(pathname)
      baseline = loadAllModels(simyear, baseline_type, baseline_tag);
      baseline = addMetadata(baseline, baseline_type, baseline_tag, smbmodel);
      return
   end

   % Resolve the canonical baseline file unless the caller already passed
   % an explicit pathname.
   if isblanktext(pathname)
      pathname = icemodel.test.helpers.defaultBaselinePath( ...
         "perf", baseline_type, baseline_tag, smbmodel, simyear);
   end

   pathname = char(pathname);
   if exist(pathname, 'file') ~= 2
      baseline = table();
      return
   end

   if any(string({whos('-file', pathname).name}) == "meta")
      S = load(pathname, 'meta');
      meta = normalizeMeta(S.meta);
   end

   % Normalize saved MAT content back into a standard table schema.
   baseline = ...
      icemodel.test.helpers.loadSavedTable(pathname, ["PerfBaseline", "baseline"]);

   if ~isempty(baseline) ...
         && ismember('simyear', baseline.Properties.VariableNames)
      baseline = baseline(baseline.simyear == simyear, :);
   end

   if ismember('case_id', baseline.Properties.VariableNames)
      baseline.case_id = icemodel.test.helpers.normalizeFormalCaseId(baseline.case_id);
   end

   if ismember('last_updated_utc', baseline.Properties.VariableNames) ...
         && ~isdatetime(baseline.last_updated_utc)
      baseline.last_updated_utc = datetime( ...
         baseline.last_updated_utc, 'ConvertFrom', 'datenum', ...
         'TimeZone', 'UTC');
   end

   baseline = addMetadata(baseline, baseline_type, baseline_tag, smbmodel);
end

function baseline = loadAllModels(simyear, baseline_type, baseline_tag)
   %LOADALLMODELS Concatenate per-model perf baselines for smbmodel="all".
   models = icemodel.namelists.smbmodel("test");
   pathnames = arrayfun(@(mdl) icemodel.test.helpers.defaultBaselinePath( ...
      "perf", baseline_type, baseline_tag, mdl, simyear), ...
      models, 'UniformOutput', false);
   exists = cellfun(@(p) exist(char(p), 'file') == 2, pathnames);
   if ~any(exists)
      baseline = table();
   else
      tables = cellfun(@(mdl, p) icemodel.test.helpers.loadPerfBaseline( ...
         simyear, baseline_tag, string(mdl), string(p)), ...
         num2cell(models(exists)), pathnames(exists), 'UniformOutput', false);
      baseline = vertcat(tables{:});
   end
end

function baseline = addMetadata(baseline, baseline_type, baseline_tag, smbmodel)
   %ADDMETADATA Backfill selector metadata missing from older saved baselines.
   if isempty(baseline)
      return
   end

   if ~ismember('baseline_type', baseline.Properties.VariableNames)
      baseline.baseline_type = repmat(baseline_type, height(baseline), 1);
   end
   if ~ismember('baseline_tag', baseline.Properties.VariableNames)
      baseline.baseline_tag = repmat(baseline_tag, height(baseline), 1);
   end
   baseline.baseline_type = string(baseline.baseline_type);
   baseline.baseline_tag = string(baseline.baseline_tag);
   if ~ismember('smbmodel_filter', baseline.Properties.VariableNames)
      baseline.smbmodel_filter = repmat(string(smbmodel), height(baseline), 1);
   end
   baseline.smbmodel_filter = string(baseline.smbmodel_filter);
end

function meta = normalizeMeta(meta)
   %NORMALIZEMETA Normalize saved baseline metadata into current string types.

   if ~isstruct(meta)
      meta = struct();
      return
   end

   if isfield(meta, 'baseline_type')
      meta.baseline_type = string(meta.baseline_type);
   end
   if isfield(meta, 'baseline_tag')
      meta.baseline_tag = string(meta.baseline_tag);
   end
   if isfield(meta, 'smbmodel_filter')
      meta.smbmodel_filter = string(meta.smbmodel_filter);
   end
   if isfield(meta, 'matlab_version')
      meta.matlab_version = string(meta.matlab_version);
   end
   if isfield(meta, 'host')
      meta.host = string(meta.host);
   end
end
