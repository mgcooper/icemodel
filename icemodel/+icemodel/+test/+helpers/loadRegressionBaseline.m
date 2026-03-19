function baseline = loadRegressionBaseline(baseline_tag, smbmodel, pathname)
   %LOADREGRESSIONBASELINE Load rolling or versioned icemodel regression baseline table.
   %
   %  baseline = icemodel.test.helpers.loadRegressionBaseline()
   %  baseline = icemodel.test.helpers.loadRegressionBaseline("rolling", "all")
   %  baseline = icemodel.test.helpers.loadRegressionBaseline("v1.01", "skinmodel")
   %  baseline = icemodel.test.helpers.loadRegressionBaseline("v1.1", "skinmodel", pathname)
   arguments
      baseline_tag (1, :) string = ""
      smbmodel string = "all"
      pathname (1, :) string = ""
   end

   % An empty selector means there is no requested regression baseline only
   % when the caller has not already provided an explicit file to load.
   if isblanktext(pathname) && isblanktext(baseline_tag)
      baseline = table();
      return
   end

   % smbmodel="all" is virtual: load and concatenate the per-model files.
   if smbmodel == "all" && isblanktext(pathname)
      baseline = loadAllModels(baseline_tag);
      return
   end

   % Resolve the canonical baseline file unless the caller already passed
   % an explicit pathname.
   if isblanktext(pathname)
      [baseline_type, baseline_tag] = ...
         icemodel.test.helpers.resolveBaselineSelector(baseline_tag);
      pathname = icemodel.test.helpers.defaultBaselinePath( ...
         "regression", baseline_type, baseline_tag, smbmodel);
   end
   pathname = char(pathname);

   if exist(pathname, 'file') ~= 2
      error('regression baseline file not found: %s', pathname)
   end

   % Normalize saved MAT content back into a standard table schema.
   baseline = icemodel.test.helpers.loadSavedTable(...
      pathname, ["RegressionBaseline", "baseline"]);

   if ismember('case_id', baseline.Properties.VariableNames)
      baseline.case_id = icemodel.test.helpers.normalizeFormalCaseId(baseline.case_id);
   end
   if ismember('baseline_tag', baseline.Properties.VariableNames)
      baseline.baseline_tag = string(baseline.baseline_tag);
   elseif lower(baseline_tag) == "rolling"
      baseline.baseline_tag = repmat("", height(baseline), 1);
   else
      baseline.baseline_tag = repmat(string(baseline_tag), height(baseline), 1);
   end
   if ismember('baseline_type', baseline.Properties.VariableNames)
      baseline.baseline_type = string(baseline.baseline_type);
   elseif lower(baseline_tag) == "rolling"
      baseline.baseline_type = repmat("rolling", height(baseline), 1);
   else
      baseline.baseline_type = repmat("release", height(baseline), 1);
   end
   if ~ismember('smbmodel_filter', baseline.Properties.VariableNames)
      baseline.smbmodel_filter = repmat(string(smbmodel), height(baseline), 1);
   end
   baseline.smbmodel_filter = string(baseline.smbmodel_filter);
   if ismember('last_updated_utc', baseline.Properties.VariableNames) ...
         && ~isdatetime(baseline.last_updated_utc)
      baseline.last_updated_utc = datetime( ...
         baseline.last_updated_utc, 'ConvertFrom', 'datenum', ...
         'TimeZone', 'UTC');
   end
end

function baseline = loadAllModels(baseline_tag)
   %LOADALLMODELS Concatenate per-model regression baselines for smbmodel="all".
   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector(baseline_tag);
   models = icemodel.namelists.smbmodel("test");
   pathnames = arrayfun(@(mdl) icemodel.test.helpers.defaultBaselinePath( ...
      "regression", baseline_type, baseline_tag, mdl), ...
      models, 'UniformOutput', false);
   exists = cellfun(@(p) exist(char(p), 'file') == 2, pathnames);
   if ~any(exists)
      baseline = table();
   else
      tables = cellfun(@(mdl, p) icemodel.test.helpers.loadRegressionBaseline( ...
         baseline_tag, string(mdl), string(p)), ...
         num2cell(models(exists)), pathnames(exists), 'UniformOutput', false);
      baseline = vertcat(tables{:});
   end
end
