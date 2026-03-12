function baseline = loadPerfBaseline(simyear, baseline_tag, smbmodel, pathname)
   %LOADPERFBASELINE Load rolling or release performance baseline table.
   %
   %  baseline = test.helpers.loadPerfBaseline(2016)
   %  baseline = test.helpers.loadPerfBaseline(2016, "rolling", "all")
   %  baseline = test.helpers.loadPerfBaseline(2016, "v1.1", "skinmodel")
   %  baseline = test.helpers.loadPerfBaseline(2016, "v1.1", "skinmodel", pathname)
   arguments
      simyear (1, 1) double {mustBeInteger, mustBePositive} = 2016
      baseline_tag string = "rolling"
      smbmodel string = "all"
      pathname = string.empty()
   end

   [baseline_type, baseline_tag] = ...
      test.helpers.resolveBaselineSelector(baseline_tag);

   % smbmodel="all" is virtual: load and concatenate the per-model files.
   if string(smbmodel) == "all" ...
         && (isempty(pathname) || (isstring(pathname) ...
         && all(strlength(pathname) == 0)))
      baseline = loadAllModels(simyear, baseline_type, baseline_tag);
      baseline = addMetadata(baseline, baseline_type, baseline_tag, smbmodel);
      return
   end

   if isempty(pathname) || (isstring(pathname) && all(strlength(pathname) == 0))
      rootdir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
      pathname = ...
         defaultBaselinePath(rootdir, simyear, baseline_type, baseline_tag, ...
         smbmodel);
   end

   pathname = char(pathname);
   if exist(pathname, 'file') ~= 2
      baseline = table();
      return
   end

   % Normalize saved MAT content back into a standard table schema.
   baseline = ...
      test.helpers.loadSavedTable(pathname, ["PerfBaseline", "baseline"]);

   if ~isempty(baseline) ...
         && ismember('simyear', baseline.Properties.VariableNames)
      baseline = baseline(baseline.simyear == simyear, :);
   end

   if ismember('case_id', baseline.Properties.VariableNames)
      baseline.case_id = test.helpers.normalizeFormalCaseId(baseline.case_id);
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
   rootdir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
   models = test.helpers.formalSmbmodels();
   tables = cell(numel(models), 1);
   k = 0;
   for i = 1:numel(models)
      pathname = defaultBaselinePath(rootdir, simyear, baseline_type, ...
         baseline_tag, models(i));
      if exist(char(pathname), 'file') ~= 2
         continue
      end
      k = k + 1;
      tables{k} = test.helpers.loadPerfBaseline( ...
         simyear, baseline_tag, models(i), pathname);
   end
   if k == 0
      baseline = table();
   else
      baseline = vertcat(tables{1:k});
   end
end

function pathname = defaultBaselinePath(...
      rootdir, simyear, baseline_type, baseline_tag, smbmodel)

   model_tag = test.helpers.smbmodelTag(smbmodel);

   if baseline_type == "rolling"
      pathname = fullfile(rootdir, 'baselines', ...
         sprintf('perf_baseline_%d_rolling_%s.mat', simyear, model_tag));
   else
      pathname = fullfile(rootdir, 'baselines', ...
         sprintf('perf_baseline_%d_%s_%s.mat', ...
         simyear, test.helpers.sanitizeTag(baseline_tag), model_tag));
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
