function baseline = loadRegressionBaseline(baseline_tag, smbmodel, pathname)
   %LOADREGRESSIONBASELINE Load rolling or versioned icemodel regression baseline table.
   %
   %  baseline = icemodel.test.helpers.loadRegressionBaseline()
   %  baseline = icemodel.test.helpers.loadRegressionBaseline("rolling", "all")
   %  baseline = icemodel.test.helpers.loadRegressionBaseline("v1.01", "skinmodel")
   %  baseline = icemodel.test.helpers.loadRegressionBaseline("v1.1", "skinmodel", pathname)
   arguments
      baseline_tag string = string.empty()
      smbmodel string = "all"
      pathname = string.empty()
   end

   if isempty(baseline_tag) || ...
         (isstring(baseline_tag) && all(strlength(baseline_tag) == 0))
      baseline = table();
      return
   end

   % smbmodel="all" is virtual: load and concatenate the per-model files.
   if string(smbmodel) == "all" ...
         && (isempty(pathname) || (isstring(pathname) ...
         && all(strlength(pathname) == 0)))
      baseline = loadAllModels(baseline_tag);
      return
   end

   if isempty(pathname) || (isstring(pathname) && all(strlength(pathname) == 0))
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
   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector(baseline_tag);
   models = icemodel.test.helpers.formalSmbmodels();
   tables = cell(numel(models), 1);
   k = 0;
   for i = 1:numel(models)
      pathname = icemodel.test.helpers.defaultBaselinePath( ...
         "regression", baseline_type, baseline_tag, models(i));
      if exist(char(pathname), 'file') ~= 2
         continue
      end
      k = k + 1;
      tables{k} = icemodel.test.helpers.loadRegressionBaseline( ...
         baseline_tag, models(i), pathname);
   end
   if k == 0
      baseline = table();
   else
      baseline = vertcat(tables{1:k});
   end
end
