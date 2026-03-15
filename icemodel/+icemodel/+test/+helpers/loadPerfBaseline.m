function baseline = loadPerfBaseline(simyear, baseline_tag, smbmodel, pathname)
   %LOADPERFBASELINE Load rolling or release performance baseline table.
   %
   %  baseline = icemodel.test.helpers.loadPerfBaseline(2016)
   %  baseline = icemodel.test.helpers.loadPerfBaseline(2016, "rolling", "all")
   %  baseline = icemodel.test.helpers.loadPerfBaseline(2016, "v1.1", "skinmodel")
   %  baseline = icemodel.test.helpers.loadPerfBaseline(2016, "v1.1", "skinmodel", pathname)
   arguments
      simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      baseline_tag string = "rolling"
      smbmodel string = "all"
      pathname = string.empty()
   end

   [baseline_type, baseline_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector(baseline_tag);

   % smbmodel="all" is virtual: load and concatenate the per-model files.
   if string(smbmodel) == "all" ...
         && (isempty(pathname) || (isstring(pathname) ...
         && all(strlength(pathname) == 0)))
      baseline = loadAllModels(simyear, baseline_type, baseline_tag);
      baseline = addMetadata(baseline, baseline_type, baseline_tag, smbmodel);
      return
   end

   if isempty(pathname) || (isstring(pathname) && all(strlength(pathname) == 0))
      pathname = icemodel.test.helpers.defaultBaselinePath( ...
         "perf", baseline_type, baseline_tag, smbmodel, simyear);
   end

   pathname = char(pathname);
   if exist(pathname, 'file') ~= 2
      baseline = table();
      return
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
   models = icemodel.test.helpers.formalSmbmodels();
   tables = cell(numel(models), 1);
   k = 0;
   for i = 1:numel(models)
      pathname = icemodel.test.helpers.defaultBaselinePath( ...
         "perf", baseline_type, baseline_tag, models(i), simyear);
      if exist(char(pathname), 'file') ~= 2
         continue
      end
      k = k + 1;
      tables{k} = icemodel.test.helpers.loadPerfBaseline( ...
         simyear, baseline_tag, models(i), pathname);
   end
   if k == 0
      baseline = table();
   else
      baseline = vertcat(tables{1:k});
   end
end

function baseline = addMetadata(baseline, baseline_type, baseline_tag, smbmodel)
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
