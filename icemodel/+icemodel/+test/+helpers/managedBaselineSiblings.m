function managed_files = managedBaselineSiblings(kind, baseline_selector, ...
      output_file, kwargs)
   %MANAGEDBASELINESIBLINGS Return the baseline files one build writes.
   %
   %  managed_files = icemodel.test.helpers.managedBaselineSiblings( ...
   %     "regression", baseline_selector, output_file)
   %  managed_files = icemodel.test.helpers.managedBaselineSiblings( ...
   %     "perf", baseline_selector, output_file, simyear=2016)
   %
   % A builder writes one file for each formal model, so worktreeRevision must
   % leave those files out of the source identity it records. Return every
   % managed sibling of the requested build, or the explicit override when the
   % caller named one.
   %
   % Name-value
   %  simyear  Benchmark year for perf baseline paths. Blank leaves the
   %           default to baselineFilePath. Regression paths ignore it.
   %
   % See also: icemodel.test.helpers.worktreeRevision,
   %  icemodel.test.helpers.baselineFilePath

   arguments
      kind (1, 1) string {mustBeMember(kind, ["perf", "regression"])}
      baseline_selector (1, :) string
      output_file string
      kwargs.simyear double = []
   end

   % An explicit output file names the only managed file this build writes.
   if ~isblanktext(output_file)
      managed_files = string(output_file);
      return
   end

   % Otherwise the build owns one file for each formal model.
   [baseline_type, resolved_tag] = ...
      icemodel.test.helpers.resolveBaselineSelector(baseline_selector);
   managed_models = icemodel.test.helpers.resolveRequestedSmbmodels("all");

   % Forward SIMYEAR only when the caller set it, so baselineFilePath keeps
   % ownership of the benchmark-year default.
   year_args = {};
   if ~isempty(kwargs.simyear)
      year_args = {'simyear', kwargs.simyear};
   end

   % baselineFilePath returns char for perf and string for regression, so
   % convert each result before it enters the string array.
   managed_files = strings(size(managed_models));
   for k = 1:numel(managed_models)
      managed_files(k) = string( ...
         icemodel.test.helpers.baselineFilePath(kind, ...
         'smbmodel', managed_models(k), ...
         'baseline_type', baseline_type, ...
         'baseline_tag', resolved_tag, year_args{:}));
   end
end
