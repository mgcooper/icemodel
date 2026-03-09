function baseline = snapshotBaseline(kind, baseline_tag, smbmodel, overwrite, output_file, simyear)
%SNAPSHOTBASELINE Save a release snapshot from the rolling test baseline.
%
%  baseline = test.helpers.snapshotBaseline("perf", "v1.1", "skinmodel", true, string.empty(), 2016)
%  baseline = test.helpers.snapshotBaseline("regression", "v1.1", "icemodel", true)

   arguments
      kind (1, :) string {mustBeMember(kind, ["perf", "regression"])}
      baseline_tag (1, :) string
      smbmodel (1, :) string {mustBeMember(smbmodel, ["all", "icemodel", "skinmodel"])} = "all"
      overwrite (1, 1) logical = false
      output_file string = string.empty()
      simyear double = NaN
   end

   if smbmodel == "all"
      if ~isblanktext(output_file)
         error('output_file is only supported for a single smbmodel')
      end
      models = test.helpers.formalSmbmodels();
      baselines = cell(numel(models), 1);
      for i = 1:numel(models)
         baselines{i} = test.helpers.snapshotBaseline( ...
            kind, baseline_tag, models(i), overwrite, string.empty(), simyear);
      end
      baseline = vertcat(baselines{:});
      return
   end

   if isblanktext(output_file)
      output_file = test.helpers.defaultBaselinePath( ...
         kind, "release", baseline_tag, smbmodel, simyear);
   end

   if exist(char(output_file), 'file') == 2 && ~overwrite
      error('release %s baseline already exists: %s', kind, char(output_file))
   end

   switch kind
      case "perf"
         baseline = test.helpers.loadPerfBaseline(simyear, 'rolling', smbmodel);
         if isempty(baseline)
            error('rolling perf baseline is empty or missing for simyear %d', simyear)
         end
         % Save using the conventional perf baseline variable name.
         PerfBaseline = baseline; %#ok<NASGU>
         varname = 'PerfBaseline';

      case "regression"
         baseline = test.helpers.loadRegressionBaseline('rolling', smbmodel);
         if isempty(baseline)
            error('rolling regression baseline is empty or missing')
         end
         % Save using the conventional regression baseline variable name.
         RegressionBaseline = baseline; %#ok<NASGU>
         varname = 'RegressionBaseline';
   end

   baseline.baseline_type = repmat("release", height(baseline), 1);
   baseline.baseline_tag = repmat(baseline_tag, height(baseline), 1);
   if kind == "perf"
      PerfBaseline = baseline; %#ok<NASGU>
   else
      RegressionBaseline = baseline; %#ok<NASGU>
   end
   save(char(output_file), varname);
end
